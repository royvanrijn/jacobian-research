#!/usr/bin/env python3
"""Exact controls for the uniform critical-value quadratic trace chart."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import platform

import sympy as sp

HERE = Path(__file__).resolve().parent


def require(value, message):
    if not value:
        raise AssertionError(message)


def scalar(x, mod):
    return int(x) % mod if mod else Fraction(x)


def solve(columns, vector, mod=None):
    n = len(columns)
    require(all(len(c) == n for c in columns) and len(vector) == n, "Nonsquare basis")
    a = [[scalar(columns[j][i], mod) for j in range(n)] + [scalar(vector[i], mod)]
         for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j,n)
                      if (gcd(int(a[i][j]),mod) == 1 if mod else a[i][j] != 0)), None)
        require(pivot is not None, "Noninvertible basis or nonunit pivot")
        a[j],a[pivot] = a[pivot],a[j]
        inv = pow(int(a[j][j]),-1,mod) if mod else 1/a[j][j]
        a[j] = [scalar(x*inv,mod) for x in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [scalar(x-c*y,mod) for x,y in zip(a[i],a[j])]
    return [a[i][-1] for i in range(n)]


class Algebra:
    def __init__(self, polynomial, mod=None):
        self.mod = mod
        self.n = len(polynomial)-1
        self.q = [scalar(x,mod) for x in polynomial]
        require(self.q[-1] == 1, "Need monic quotient polynomial")
        self.zero = [scalar(0,mod)]*self.n
        self.one = self.zero[:]
        self.one[0] = scalar(1,mod)

    def reduce(self, value):
        value = [scalar(x,self.mod) for x in value]
        while len(value) > self.n:
            c = value.pop()
            off = len(value)-self.n
            for j in range(self.n):
                value[off+j] = scalar(value[off+j]-c*self.q[j],self.mod)
        return value + self.zero[:self.n-len(value)]

    def add(self,a,b,sign=1):
        return [scalar(x+sign*y,self.mod) for x,y in zip(a,b)]

    def mul(self,a,b):
        c = [scalar(0,self.mod)]*(2*self.n-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                c[i+j] = scalar(c[i+j]+x*y,self.mod)
        return self.reduce(c)

    def scale(self,a,c):
        return [scalar(x*c,self.mod) for x in a]

    def power(self,a,k):
        out = self.one[:]
        while k:
            if k%2:
                out = self.mul(out,a)
            a = self.mul(a,a)
            k //= 2
        return out

    def evaluate(self,coefficients,a):
        out = self.zero[:]
        for c in reversed(coefficients):
            out = self.add(self.mul(out,a),self.scale(self.one,c))
        return out

    def columns(self,a):
        return [self.mul(a,self.reduce([0]*i+[1])) for i in range(self.n)]

    def inverse(self,a):
        out = solve(self.columns(a),self.one,self.mod)
        require(self.mul(a,out) == self.one, "Failed inverse")
        return out


def charpoly(columns,mod=None):
    n = len(columns)
    matrix = sp.Matrix(n,n,lambda i,j:sp.Rational(columns[j][i]))
    coeffs = matrix.charpoly().all_coeffs()[::-1]
    return [scalar(x,mod) for x in coeffs]


def poly_mul(a,b,mod=None):
    out = [scalar(0,mod)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] = scalar(out[i+j]+x*y,mod)
    return out


def critical_polynomial(N,epsilon,mod=None):
    L = 2*(len(N)-1)-1
    R = [scalar((2*j-1)*c,mod) for j,c in enumerate(N)]
    if epsilon:
        s = [scalar(0,mod)]*(len(R)-1)
        partial = scalar(0,mod)
        for j in range(len(R)-1,0,-1):
            partial = scalar(partial+R[j],mod)
            s[j-1] = partial
        require(scalar(R[0]+s[0],mod) == 0, "No distinguished critical root at one")
    else:
        s = R
    invL = pow(L,-1,mod) if mod else Fraction(1,L)
    return [scalar(c*invL,mod) for c in s]


def quadratic_product(B,a,b,trace,norm):
    # Elements are a0+a1*u, with u²=trace*u-norm.
    aa = B.mul(a[0],b[0])
    bb = B.mul(a[1],b[1])
    const = B.add(aa,B.mul(norm,bb),-1)
    linear = B.add(B.add(B.mul(a[0],b[1]),B.mul(a[1],b[0])),B.mul(trace,bb))
    return const,linear


def quadratic_evaluate(B,coeffs,trace,norm):
    out = (B.zero[:],B.zero[:])
    U = (B.zero[:],B.one[:])
    for c in reversed(coeffs):
        out = quadratic_product(B,out,U,trace,norm)
        out = (B.add(out[0],B.scale(B.one,c)),out[1])
    return out


def recover_and_reconstruct(N,epsilon,mod=None,ambient=None,U=None):
    N = [scalar(c,mod) for c in N]
    h = ((len(N)-1)-epsilon)//2
    n = 2*h
    require(len(N)-1 == 2*h+epsilon and N[-1] == 1, "Wrong numerator degree")
    S = critical_polynomial(N,epsilon,mod)
    A = Algebra(ambient or S,mod)
    require(A.n == n, "Wrong residual critical algebra rank")
    U = A.reduce(U if U is not None else [0,1])
    require(A.evaluate(S,U) == A.zero, "Input U is not a residual critical root")
    require(charpoly(A.columns(U),mod) == S, "Input U is not primitive")
    Nu = A.evaluate(N,U)
    M = A.mul(A.mul(Nu,Nu),A.inverse(U))
    A.inverse(Nu)
    powers = [A.one]
    for _ in range(h):
        powers.append(A.mul(powers[-1],M))
    basis = powers[:h] + [A.mul(U,x) for x in powers[:h]]
    coords = solve(basis,Nu,mod)
    b,a = coords[:h],coords[h:]
    mcoords = solve(basis,powers[h],mod)
    require(not any(mcoords[h:]), "Critical values do not generate a degree-h algebra")
    q = [scalar(-x,mod) for x in mcoords[:h]] + [scalar(1,mod)]
    B = Algebra(q,mod)
    w = B.reduce([0,1])
    B.inverse(w)
    inva = B.inverse(a)
    rho = B.scale(B.mul(b,inva),-1)
    B.inverse(rho)
    norm = B.mul(rho,rho)
    inva2 = B.mul(inva,inva)
    trace = B.add(B.scale(rho,2),B.mul(w,inva2))
    delta = B.mul(w,B.add(w,B.scale(B.mul(rho,B.mul(a,a)),4)))
    B.inverse(delta)

    # Characteristic polynomial of u in the quadratic algebra A/B.
    norm_matrix = B.columns(norm)
    trace_matrix = B.columns(trace)
    columns = [[scalar(0,mod)]*h + B.reduce([0]*j+[1]) for j in range(h)]
    columns += [B.scale(norm_matrix[j],-1) + trace_matrix[j] for j in range(h)]
    recovered_S = charpoly(columns,mod)
    require(recovered_S == S, "Quadratic trace chart did not recover the actual critical polynomial")

    target = (B.scale(B.mul(a,rho),-1),a)
    require(quadratic_evaluate(B,N,trace,norm) == target, "Euler numerator failed the critical-value interpolation equation")
    nu_squared = quadratic_product(B,target,target,trace,norm)
    require(nu_squared == (B.zero,w), "nu²=u*w")

    norm_nu = B.add(B.add(B.mul(b,b),B.mul(B.mul(a,b),trace)),B.mul(B.mul(a,a),norm))
    require(norm_nu == B.scale(B.mul(rho,w),-1), "Canonical signed norm formula")

    L = 4*h+2*epsilon-1
    R = [scalar(L*x,mod) for x in (poly_mul([-1,1],S,mod) if epsilon else S)]
    require(R == [scalar((2*j-1)*c,mod) for j,c in enumerate(N)], "Euler coefficient reconstruction")
    if not mod:
        recovered_N = [Fraction(c,2*j-1) for j,c in enumerate(R)]
        require(recovered_N == N, "Characteristic-zero Euler inverse")
    if epsilon:
        require(scalar(sum(N),mod) == 0 and scalar(sum(j*c for j,c in enumerate(N)),mod) == 0,
                "Normalized double zero")
        require(scalar(sum(S),mod) != 0, "Non-simple distinguished critical root")
    norm_rho_matrix = B.columns(rho)
    norm_rho = scalar(sp.Matrix(h,h,lambda i,j:sp.Rational(norm_rho_matrix[j][i])).det(),mod)
    require(scalar(S[0]-norm_rho*norm_rho,mod) == 0, "Residual norm is a scalar square")
    require(scalar(N[0]-(-1)**(epsilon+1)*L*norm_rho*norm_rho,mod) == 0,
            "Uniform constant-term squareclass")
    t = sp.symbols("t")
    qp = sp.Poly(sum(int(x)*t**j for j,x in enumerate(q)),t,modulus=3) if mod else sp.Poly(
        sum(sp.Rational(x)*t**j for j,x in enumerate(q)),t)
    require(sp.gcd(qp,qp.diff()).degree() == 0, "Critical value algebra is not etale")
    real_test = None
    if not mod:
        second = B.add(w,B.scale(B.mul(rho,B.mul(a,a)),4))
        def trace_element(value):
            columns = B.columns(value)
            return sum(columns[i][i] for i in range(h))
        def hermite(value):
            return sp.Matrix(h,h,lambda i,j:sp.Rational(
                trace_element(B.mul(value,B.power(w,i+j)))))
        forms = [hermite(w),hermite(second)]
        minors = [[sp.factor(form[:j,:j].det()) for j in range(1,h+1)] for form in forms]
        real_test = dict(leading_principal_minors=[[str(x) for x in row] for row in minors],
                         all_critical_points_real=all(x>0 for row in minors for x in row))
    return dict(h=h,epsilon=epsilon,degree=4*h+2*epsilon,
                exact_Hermite_real_test=real_test,
                q_coefficients_ascending=list(map(str,q)),
                a_coefficients_ascending=list(map(str,a)),
                rho_coefficients_ascending=list(map(str,rho)),
                paired_product_coefficients_ascending=list(map(str,norm)),
                paired_product_is_base_constant=not any(norm[1:]),
                norm_rho=str(norm_rho), constant_term_squareclass_checked=True,
                entire_quadratic_chart_roundtrip=True,
                Euler_inverse_checked_in_characteristic_zero=not bool(mod),
                coefficient_ring=("Q" if not mod else "Z/"+str(mod)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args = parser.parse_args()
    q8 = [Fraction(-1,5488),Fraction(233,5292),Fraction(9,28),1,1]
    q10 = [Fraction(16,5625),Fraction(1804,5625),Fraction(-14006,5625),Fraction(3129,625),Fraction(-96,25),1]
    controls = [recover_and_reconstruct(q8,0),recover_and_reconstruct(q10,1)]
    require(controls[0]["paired_product_is_base_constant"] is False,"Known rational degree8 test should be nonreciprocal")
    require(controls[1]["rho_coefficients_ascending"] == ["2/15","0"],"Known rational reciprocal degree10 scale")
    require(not controls[0]["exact_Hermite_real_test"]["all_critical_points_real"],
            "Known c=t=1 degree8 negative real control")
    require(controls[1]["exact_Hermite_real_test"]["all_critical_points_real"],
            "Known degree10 positive real control")
    q8_real = [q8[0],-q8[1],q8[2],-q8[3],q8[4]]
    controls.append(recover_and_reconstruct(q8_real,0))
    require(controls[-1]["exact_Hermite_real_test"]["all_critical_points_real"],
            "Known c=t=-1 degree8 positive real control")

    hensel_path = HERE/"degree_twenty_six_three_adic_verification.json"
    hensel = json.loads(hensel_path.read_text())
    T = [1,0,0,1,1,1,2,1,1,1,0,0,1]
    local = recover_and_reconstruct(hensel["mod9_numerator_coefficients"],1,9,T,
                                   hensel["mod9_critical_root_in_A"])
    require(local["paired_product_is_base_constant"] is False,"Mod9 degree26 lift should be genuinely nonreciprocal")
    controls.append(local)

    report = dict(schema_version=1,result="PASS",python=platform.python_version(),sympy=sp.__version__,
                  checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),controls=controls,
                  local_input_sha256=sha256(hensel_path.read_bytes()).hexdigest(),
                  scope=("Written uniform equivalence is required. Exact Q degree8/10 roundtrips and "
                         "the degree26 modulo9 Hensel shadow are controlled. The modulo9 check does not "
                         "divide nonunit Euler constants and gives no rational descent or new Q3 model."),
                  formal_verification=False,independent_replay=False,external_review=False)
    text = json.dumps(report,indent=2)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(text,encoding="utf-8")
    print(text,end="")


if __name__ == "__main__":
    main()

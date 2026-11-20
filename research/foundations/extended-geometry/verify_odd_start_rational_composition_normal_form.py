#!/usr/bin/env python3
"""Exact controls for the all-odd rational composition normal form; stdlib only."""
import argparse,json,sys,hashlib
from fractions import Fraction as Q
from math import comb
from pathlib import Path
if not __debug__:raise SystemExit("Run without -O: assertion controls are required.")

def trim(p):
    p = tuple(Q(x) for x in p)
    while p and p[-1] == 0:
        p = p[:-1]
    return p


Z = ()
ONE = (Q(1),)
XX = (Q(0), Q(1))


def add(p, q):
    return trim(tuple((p[i] if i < len(p) else 0) +
                      (q[i] if i < len(q) else 0)
                      for i in range(max(len(p), len(q)))))


def scale(p, c):
    return trim(tuple(c*x for x in p))


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    if not p or not q:
        return Z
    out = [Q(0)]*(len(p)+len(q)-1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i+j] += pi*qj
    return trim(out)


def power(p, n):
    assert n >= 0
    out = ONE
    while n:
        if n & 1:
            out = mul(out, p)
        p = mul(p, p)
        n //= 2
    return out


def deriv(p):
    return trim(tuple(i*p[i] for i in range(1, len(p))))


def divmodp(p, q):
    assert q
    rem = list(p)
    quotient = [Q(0)]*max(0, len(p)-len(q)+1)
    while rem and len(rem) >= len(q):
        shift = len(rem)-len(q)
        c = rem[-1]/q[-1]
        quotient[shift] = c
        for j, qj in enumerate(q):
            rem[shift+j] -= c*qj
        rem = list(trim(rem))
    return trim(quotient), trim(rem)


def exactdiv(p, q):
    quotient, rem = divmodp(p, q)
    assert not rem
    return quotient


def monic(p):
    assert p
    return scale(p, 1/p[-1])


def gcdp(p, q):
    while q:
        p, q = q, divmodp(p, q)[1]
    return monic(p) if p else Z


def degree(p):
    return len(p)-1


def radical(p):
    return monic(exactdiv(p, gcdp(p, deriv(p))))


def max_multiplicity(p):
    count = 0
    while degree(p) > 0:
        count += 1
        p = gcdp(p, deriv(p))
    return count


class RF:
    def __init__(self, num=ONE, den=ONE):
        num, den = trim(num), trim(den)
        assert den
        g = gcdp(num, den)
        num, den = exactdiv(num, g), exactdiv(den, g)
        self.num, self.den = scale(num, 1/den[-1]), scale(den, 1/den[-1])

    def __mul__(self, other):
        return RF(mul(self.num, other.num), mul(self.den, other.den))

    def __truediv__(self, other):
        assert other.num
        return RF(mul(self.num, other.den), mul(self.den, other.num))

    def __add__(self, other):
        return RF(add(mul(self.num, other.den), mul(other.num, self.den)),
                  mul(self.den, other.den))

    def __pow__(self, n):
        if n < 0:
            return RF(power(self.den, -n), power(self.num, -n))
        return RF(power(self.num, n), power(self.den, n))

    def derivative(self):
        return RF(sub(mul(deriv(self.num), self.den),
                      mul(self.num, deriv(self.den))), power(self.den, 2))

    def at_zero(self):
        assert self.den[0] != 0
        return self.num[0]/self.den[0] if self.num else Q(0)

    def polynomial(self):
        return exactdiv(self.num, self.den)

    def __eq__(self, other):
        return self.num == other.num and self.den == other.den



def neg(f):return RF(scale(f.num,-1),f.den)
def minus(f,g):return f+neg(g)
def compose_rf(f,g):
    def evaluate(p):
        out=RF(Z)
        for c in reversed(p):out=out*g+RF((c,))
        return out
    return evaluate(f.num)/evaluate(f.den)
def regular_zero(f):return bool(f.den)and f.den[0]!=0
def deg_rf(f):return max(degree(f.num),degree(f.den))
def odd(f):return compose_rf(f,RF(scale(XX,-1)))==neg(f)
def finite_tests(a,u,T,C,D):
    A,B,M,N=u.num,u.den,T.num,T.den
    Delta=sub(mul(deriv(A),B),mul(A,deriv(B)))
    assert Delta and M
    cn=mul(mul(power(A,a-1),power(B,a+1)),power(N,a))
    cd=mul(mul(power(XX,a-1),power(M,a)),power(Delta,a))
    dn=scale(mul(mul(power(A,a+1),power(B,a+3)),power(N,a+2)),-1)
    dd=mul(mul(power(XX,a+1),power(M,a+2)),power(Delta,a+2))
    assert RF(cn,cd)==C and RF(dn,dd)==D
    return not divmodp(cn,cd)[1],not divmodp(dn,dd)[1]
def build(a,w,U,V):
    assert a>=1 and a%2==1
    assert regular_zero(w)and regular_zero(U)and regular_zero(V)
    assert w.at_zero()==U.at_zero()==V.at_zero()==0
    assert w.derivative().at_zero()==U.derivative().at_zero()==1
    assert odd(U)and odd(V)
    u=compose_rf(U,w);T=RF(ONE)+compose_rf(V,w)
    x=RF(XX);X=u/(x*T*u.derivative())
    C=x*(X**a)/u;D=neg(C*(X**2))
    assert X.at_zero()==C.at_zero()==1 and D.at_zero()==-1
    assert x*(X**a)/C==u
    assert (X**(a-1))/(C*u.derivative())==T
    ct,dt=finite_tests(a,u,T,C,D)
    assert ct and dt
    return u,T,X,C.polynomial(),D.polynomial()
def truncated_mul(p,q,n):
    return trim(mul(p,q)[:n+1])
def truncated_power(p,k,n):
    out=ONE
    while k:
        if k&1:out=truncated_mul(out,p,n)
        p=truncated_mul(p,p,n);k//=2
    return out
def raw(C,D,a,m,j=0):
    if a*m%2 or m<j:return Q(0)
    n=m-j;out=Q(0)
    for i in range(m+1):
        p=truncated_mul(truncated_power(C,m-i,n),truncated_power(D,i,n),n)
        out+=Q(comb(m,i),a*m+2*i+1)*(p[n]if len(p)>n else 0)
    return out
def hyperbolic_target(a,n,lam):
    H=(Q(1),lam)
    S=tuple(Q(comb(2*n,j+1))*lam**j for j in range(2*n))
    C=scale(mul(power(H,n+a),power(S,a-1)),Q(1,2*n)**(a-1))
    D=scale(mul(power(H,n+a+2),power(S,a+1)),-Q(1,2*n)**(a+1))
    return C,D
def hyperbolic_data(n,lam,source=None):
    x=RF(XX);H=RF((Q(1),lam));beta=n*lam/2
    w=RF((Q(2)/(n*lam),))*minus(H**n,RF(ONE))/(H**n+RF(ONE))
    U=x/minus(RF(ONE),RF((beta*beta,))*(x**2))
    V=neg(RF((2*beta,))*x)/(RF(ONE)+RF((beta*beta,))*(x**2))
    recover=RF((2,))*U*V/minus(V,RF((2*beta,))*U)
    assert recover==x
    if source is not None:w=compose_rf(w,source)
    return w,U,V
def record_fixture(name,a,w,U,V,maxm,expected=None):
    u,T,X,C,D=build(a,w,U,V)
    if expected is not None:assert (C,D)==expected
    for m in range(1,maxm+1):assert raw(C,D,a,m)==0
    degree_phase=max(degree(C)+a-1,degree(D)+a+1)
    return {"name":name,"a":a,"phase_degree":degree_phase,
            "w_degree":deg_rf(w),"U_degree":deg_rf(U),"V_degree":deg_rf(V),
            "w_polynomial":degree(w.den)==0,
            "C_degree":degree(C),"D_degree":degree(D),
            "u_degree":deg_rf(u),"T_degree":deg_rf(T),
            "raw_purity_controls":maxm,
            "C_coefficients":list(map(str,C)),"D_coefficients":list(map(str,D))}
def expect_reject(thunk):
    try:thunk()
    except AssertionError:return True
    raise AssertionError("Malformed data accepted")
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo-root",type=Path,default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("odd_start_rational_composition_verification.json"))
    args=parser.parse_args();here=Path(__file__).resolve().parent
    ip=here/"odd_start_rational_composition_inputs.json"
    lp=here/"odd_start_rational_composition_source_ledger.json"
    data=json.loads(ip.read_text());ledger=json.loads(lp.read_text())
    pins=[]
    for entry in ledger["sources"]:
        p=args.repo_root/entry["path"];assert digest(p)==entry["sha256"]
        pins.append({"path":entry["path"],"sha256":digest(p),"role":entry["role"]})
    x=RF(XX);one=RF(ONE);zero=RF(Z);records=[]
    for text in data["involution_beta_values"]:
        beta=Q(text);s=neg(x)/(one+RF((beta,))*x);h=x/(one+RF((beta/2,))*x)
        assert compose_rf(s,s)==x and compose_rf(h,s)==neg(h)
        assert h.at_zero()==0 and h.derivative().at_zero()==1
    for a in data["starting_heights"]:
        records.append(record_fixture("trivial",a,x,x,zero,data["maximum_moment"],(ONE,scale(ONE,-1))))
        for lt in data["lambda_values"]:
            lam=Q(lt);H=(Q(1),lam);w=x/RF(H)
            records.append(record_fixture("sharp_Phi1",a,w,x,zero,data["maximum_moment"],
                (power(H,a+1),scale(power(H,a+3),-1))))
            V=neg(RF((lam,))*x)
            records.append(record_fixture("affine_Phi_rational_outer",a,w,x,V,data["maximum_moment"],
                (power(H,2*a+1),scale(power(H,2*a+5),-1))))
            # The canonical outer phase is genuinely rational, not polynomial.
            outerX=one/minus(one,RF((lam,))*x)
            assert degree((outerX**a).den)>0 and degree((outerX**(a+2)).den)>0
        for n in data["hyperbolic_indices"]:
            w,U,V=hyperbolic_data(n,Q(1))
            records.append(record_fixture("hyperbolic_primitive",a,w,U,V,data["maximum_moment"],
                hyperbolic_target(a,n,Q(1))))
    for a in [3,5]:
        A=(Q(1),Q(1));B=(Q(1),Q(2,3))
        w=RF(mul(XX,power(B,2)),power(A,3))
        expected=(mul(power(B,a-2),power(A,a+3)),scale(mul(power(B,a),power(A,a+5)),-1))
        records.append(record_fixture("two_root_rational_w",a,w,x,zero,data["maximum_moment"],expected))
    source=RF(mul(XX,power((Q(1),Q(1)),2)))
    w,U,V=hyperbolic_data(2,Q(27,4),source)
    prime=record_fixture("prime59_pullback",3,w,U,V,data["maximum_moment"])
    assert prime["phase_degree"]==59 and (prime["C_degree"],prime["D_degree"])==(31,55)
    records.append(prime)
    # Same polynomial phase, redundant data: a genuine odd right factor of degree3.
    w_small=x/RF((Q(1),Q(1,2)))
    right=x*(one+RF((Q(1,12),))*(x**2))/(one+RF((Q(3,4),))*(x**2))
    _,Uo,Vo=hyperbolic_data(3,Q(1))
    Ur,Vr=compose_rf(Uo,right),compose_rf(Vo,right)
    redundant=record_fixture("nonprimitive_redundant_data",1,w_small,Ur,Vr,data["maximum_moment"],
                              hyperbolic_target(1,3,Q(1)))
    assert deg_rf(right)==3
    assert RF((2,))*Ur*Vr/minus(Vr,RF((Q(3),))*Ur)==right
    canonical_w,_,_=hyperbolic_data(3,Q(1))
    assert compose_rf(right,w_small)==canonical_w and deg_rf(canonical_w)==3
    redundant["canonical_w_degree"]=3;redundant["presented_w_degree"]=1
    records.append(redundant)
    # Defect controls enforce each independent hypothesis.
    expect_reject(lambda:build(1,x**2,x,zero))
    expect_reject(lambda:build(1,x,x**2,zero))
    badV=minus(one/RF((Q(1),Q(1))),one)
    expect_reject(lambda:build(1,x,x,badV))
    assert raw((Q(1),Q(1)),scale(power((Q(1),Q(1)),3),-1),1,2)!=0
    expect_reject(lambda:build(1,x,x,x)) # correct parity, nonpolynomial C,D
    badw=RF(mul(XX,power((Q(1),Q(2,3)),2)),power((Q(1),Q(1)),3))
    expect_reject(lambda:build(1,badw,x,zero)) # height3 example fails at height1
    names=["ODD_START_RATIONAL_COMPOSITION_NORMAL_FORM.md",Path(__file__).name,ip.name,lp.name]
    out={"status":"PASS","scope":"Exact finite controls supplement the exhaustive all-field/all-order rational-composition iff. No polynomial outer phase or polynomial w assumption.",
         "involution_straightening_controls":len(data["involution_beta_values"]),
         "fixtures":records,"raw_purity_controls":sum(r["raw_purity_controls"]for r in records),
         "malformed_controls":["ramified marked base","bad U derivative","nonodd endpoint part","parity-correct nonpolynomial source","height1 rejected counterexample"],
         "primitive_controls":"Hyperbolic identity z=2UV/(V-2betaU); explicit redundant degree3 common right factor accepted but not canonical.",
         "source_pins":pins,"packet_sha256":{name:digest(here/name)for name in names},
         "runtime":{"python":sys.version.split()[0],"standard_library_only":True},
         "assurance":"Written unrestricted proof and bounded exact controls; no external or formal verification."}
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","fixtures":len(records),"raw_purity_controls":out["raw_purity_controls"],"output":str(args.output)},sort_keys=True))
if __name__=="__main__":main()


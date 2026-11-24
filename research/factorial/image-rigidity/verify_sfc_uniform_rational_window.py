#!/usr/bin/env python3
"""Focused exact controls for the rational-observable uniform-window proof.

The all-parameter statement is proved in the note. This standard-library
replay authenticates its resultants, weighted reductions, connection,
minimal-order cofactor relation and degree estimates on the stated cases.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from math import factorial, prod
from pathlib import Path

if not __debug__:
    raise SystemExit("REJECT: optimized execution is not an approved replay")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "research/factorial/image-rigidity"
NOTE = HERE / "SFC_ALGEBRAIC_UNIFORM_WINDOWS.md"
NOTE_PIN = "5066e9fef20d387d80bf865cb07f19563bb5907ace2c8e5c0032b29028f3b2f4"
HELPER = SOURCE / "verify_sfc_monomial_linear_contact.py"
HELPER_PIN = "c80ea0aa4ade828a89d18e028072c65fc0525c6700afc661183b5cb14c2e1055"


def need(ok, message):
    if not ok:
        raise ValueError(message)


def tr(a):
    a = list(map(F, a)) or [F(0)]
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def deg(a):
    return -1 if a == [0] else len(a)-1


def add(a, b):
    return tr([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
               for i in range(max(len(a),len(b)))])


def scale(a, value):
    return tr([value*x for x in a])


def sub(a, b):
    return add(a,scale(b,-1))


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return tr(out)


def pw(a, n):
    result = [F(1)]
    for _ in range(n):
        result = mul(result,a)
    return result


def der(a):
    return tr([i*a[i] for i in range(1,len(a))])


def compose(a, b):
    result = [F(0)]
    for x in reversed(a):
        result = add(mul(result,b),[x])
    return result


def evaluate(a, value):
    result = F(0)
    for coefficient in reversed(a):
        result = result*value+coefficient
    return result


def division(a, b):
    a,b=tr(a),tr(b)
    need(b != [0], "division by zero polynomial")
    quotient=[F(0)]*max(1,len(a)-len(b)+1)
    while a != [0] and len(a)>=len(b):
        offset=len(a)-len(b)
        value=a[-1]/b[-1]
        quotient[offset]+=value
        for j,c in enumerate(b):
            a[j+offset]-=value*c
        a=tr(a)
    return tr(quotient),a


def exact(a, b):
    quotient,remainder=division(a,b)
    need(remainder == [0], "non-exact polynomial division")
    return quotient


def gcd(a, b):
    a,b=tr(a),tr(b)
    while b != [0]:
        a,b=b,division(a,b)[1]
    return scale(a,1/a[-1])


def radical(a):
    out=exact(a,gcd(a,der(a)))
    return scale(out,1/out[-1])


def determinant(matrix):
    """Fraction-free determinant in Q[t], with exact division checked."""
    n=len(matrix)
    if n==0:
        return [F(1)]
    work=[[tr(value) for value in row] for row in matrix]
    previous=[F(1)]
    sign=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if work[i][k] != [0]),None)
        if pivot is None:
            return [F(0)]
        if pivot!=k:
            work[k],work[pivot]=work[pivot],work[k]
            sign=-sign
        value=work[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                work[i][j]=exact(sub(mul(value,work[i][j]),mul(work[i][k],work[k][j])),previous)
        for i in range(k+1,n):
            work[i][k]=[F(0)]
        previous=value
    return scale(work[-1][-1],sign)


def resultant_at(a, b, value):
    f=a[:]
    f[0]-=value
    m,n=len(f)-1,len(b)-1
    matrix=[]
    for shift in range(n):
        matrix.append([[F(0)]]*shift+[[x] for x in reversed(f)]+[[F(0)]]*(n-1-shift))
    for shift in range(m):
        matrix.append([[F(0)]]*shift+[[x] for x in reversed(b)]+[[F(0)]]*(m-1-shift))
    result=determinant(matrix)
    need(deg(result)<=0,"numeric resultant remained nonconstant")
    return result[0]


def interpolate(values):
    """Lagrange interpolation on the entire integer grid 0,...,degree."""
    result=[F(0)]
    for i,value in enumerate(values):
        basis=[F(1)]
        denominator=F(1)
        for j in range(len(values)):
            if j!=i:
                basis=mul(basis,[-j,1])
                denominator*=i-j
        result=add(result,scale(basis,value/denominator))
    return result


def resultant(a,b):
    # Res_z(a(z)-t,b(z)) has degree at most deg(b) in t.
    bound=deg(b)
    values=[resultant_at(a,b,F(j)) for j in range(bound+1)]
    result=interpolate(values)
    need(deg(result)==bound,"resultant t-degree authentication")
    for j,value in enumerate(values):
        need(evaluate(result,F(j))==value,"complete resultant interpolation")
    return result


def reduce_basis(coefficients, a):
    """Reduce Q[t][z] by a(z)-t, storing polynomials in t at each z-degree."""
    d=deg(a)
    work=[tr(x) for x in coefficients]
    while len(work)>d:
        shift=len(work)-1-d
        leading=scale(work[-1],1/a[-1])
        work.pop()
        for j in range(d):
            work[j+shift]=sub(work[j+shift],scale(leading,a[j]))
        work[shift]=add(work[shift],mul([0,1],leading))
    return work+[[F(0)]]*(d-len(work))


def reconstruct(vector,a):
    result=[F(0)]
    for i,value in enumerate(vector):
        result=add(result,[F(0)]*i+compose(value,a))
    return result


def matvec(matrix,vector):
    result=[]
    for row in matrix:
        value=[F(0)]
        for entry,coefficient in zip(row,vector):
            value=add(value,mul(entry,coefficient))
        result.append(value)
    return result


def minimal_relation(vectors,d):
    """Authenticate a first-independent-block minor and full cofactor relation."""
    q=len(vectors)-1
    chosen=None
    for rows in combinations(range(d),q):
        square=[[vectors[j][i] for j in range(q)] for i in rows]
        if determinant(square) != [0]:
            chosen=list(rows)
            break
    need(chosen is not None,"first q derivative vectors must be independent")
    minors=[]
    for omitted in range(q+1):
        square=[[vectors[j][i] for j in range(q+1) if j!=omitted] for i in chosen]
        minors.append(scale(determinant(square),(-1)**omitted))
    need(minors[-1]!=[0],"leading cofactor must be nonzero")
    residual=[]
    for row in range(d):
        value=[F(0)]
        for j in range(q+1):
            value=add(value,mul(minors[j],vectors[j][row]))
        residual.append(value)
    return chosen,minors,all(value==[0] for value in residual)


def check_case(exponents,coefficients,h):
    P,A,a,_,_,D=h.build(exponents,coefficients)
    P,A,a,D=map(tr,(P,A,a,D))
    s,d=len(coefficients),deg(a)
    need(deg(D)==s and D[0]==1 and a[0]==0 and a[1]==1,"actual frame degrees")
    delta=resultant(a,D)
    quotient=exact(compose(delta,a),D)
    need(deg(quotient)==s*(d-1),"resultant quotient degree")
    U=reduce_basis([[x] for x in quotient],a)
    need(reconstruct(U,a)==quotient,"exact observable basis reconstruction")
    e=s-(s+d-1)//d
    u=max(deg(x) for x in U)
    for i,coordinate in enumerate(U):
        need(deg(coordinate)<0 or deg(coordinate)<=(s*(d-1)-i)//d,"weighted observable degree")
    need(0<=u<=e,"observable numerator maximum degree")
    if d==1:
        critical=[F(1)]
        sigma=[F(1)]
        connection=[[[F(0)]]]
    else:
        critical=resultant(a,der(a))
        sigma=radical(critical)
        quotient_connection=exact(compose(sigma,a),der(a))
        c=deg(sigma)
        need(deg(quotient_connection)==(c-1)*d+1,"connection quotient degree")
        columns=[]
        for j in range(d):
            numerator=[F(0)] if j==0 else [F(0)]*(j-1)+scale(quotient_connection,j)
            column=reduce_basis([[x] for x in numerator],a)
            need(reconstruct(column,a)==numerator,"connection basis reconstruction")
            need(all(deg(x)<=c-1 for x in column),"connection coefficient degree")
            columns.append(column)
        connection=[[columns[j][i] for j in range(d)] for i in range(d)]
    need(division(mul([0,1],delta),sigma)[1]==[0],"actual critical-value alignment")
    omega=radical(mul(delta,sigma))
    C=deg(omega)
    need(C<=s+1 and (omega[0]==0 or C<=s),"support bound for pole locations")
    logdelta=exact(mul(omega,der(delta)),delta)
    connection_scale=exact(omega,sigma)
    need(deg(logdelta)<=C-1,"logarithmic denominator degree")
    scaled_connection=[[mul(connection_scale,x) for x in row] for row in connection]
    need(all(deg(x)<=C-1 for row in scaled_connection for x in row),"scaled connection degree")
    vectors=[U]
    selected=minors=None
    for k in range(d):
        vector=vectors[-1]
        conn=matvec(scaled_connection,vector)
        multiplier=add(logdelta,scale(der(omega),k))
        next_vector=[add(mul(omega,der(x)),sub(y,mul(multiplier,x))) for x,y in zip(vector,conn)]
        need(all(deg(x)<=u+(k+1)*(C-1) for x in next_vector),"derivative vector degree")
        vectors.append(next_vector)
        rows,m,complete=minimal_relation(vectors,d)
        if complete:
            selected,minors=rows,m
            break
    need(minors is not None,"d-dimensional dependence must occur")
    q=len(vectors)-1
    M=q*u+(C-1)*q*(q+1)//2
    operator=[mul(m,pw(omega,k)) for k,m in enumerate(minors)]
    for k,p in enumerate(operator):
        need(deg(p)-k<=M,"operator coefficient degree")
        need(deg(minors[k])<=q*u+(C-1)*(q*(q+1)//2-k),"cofactor degree")
        if omega[0]==0:
            need(all(value==0 for value in p[:k]),"restored zero valuation")
    actual_window=M+(q if omega[0]!=0 else 0)
    uniform_window=d*e+s*d*(d+1)//2
    need(1<=actual_window<=q*e+s*q*(q+1)//2<=uniform_window,"uniform bound conversion")
    if d==1:
        need(q==1 and u==0 and actual_window==len(set(coefficients)),"linear distinct-coefficient window")
    # Independent actual factorial expansion authenticates the distinguished scalar branch.
    maximum=10
    inverse=h.inverse(a,maximum)
    rho=h.compose(h.reciprocal(D,maximum),inverse,maximum)
    direct=[]
    for m in range(9):
        value=h.actual_moment(exponents,coefficients,m)
        normalization=factorial(m)*prod(factorial(m*r) for r in exponents)
        normalized=F(value,normalization)
        need(normalized==h.coeff(rho,m),"actual factorial observable")
        direct.append(normalized)
    return {"exponents":exponents,"coefficients":coefficients,"s":s,"d":d,
            "P":P,"A":A,"a":a,"D_phase":D,"Delta":delta,
            "critical_resultant":critical,"sigma":sigma,"Omega":omega,"C":C,
            "U":U,"u":u,"e":e,"connection":connection,
            "derivative_vectors":vectors,"minimal_order_q":q,"independent_rows":selected,
            "cofactors":minors,"operator":operator,"M":M,
            "actual_window_bound":actual_window,"coarse_uniform_bound":uniform_window,
            "decay_exponent":{"numerator":-s,"denominator":d},
            "direct_normalized_moments":direct}


def qjson(value):
    if isinstance(value,F):
        return {"numerator":value.numerator,"denominator":value.denominator}
    raise TypeError(type(value).__name__)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note",type=Path,default=NOTE)
    parser.add_argument("--output",type=Path,default=HERE/"sfc_uniform_rational_window_verification.json")
    args=parser.parse_args()
    need(sha256(args.note.read_bytes()).hexdigest()==NOTE_PIN,"note pin mismatch")
    need(sha256(HELPER.read_bytes()).hexdigest()==HELPER_PIN,"arithmetic helper pin mismatch")
    spec=importlib.util.spec_from_file_location("retained_actual_contact_arithmetic",HELPER)
    h=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    cases=[((0,0),(1,-1)),((0,0,0),(1,1,-1)),
           ((1,0),(1,2)),((0,1),(1,2)),
           ((0,2),(1,1)),((0,0,2),(1,1,-2))]
    records=[check_case(r,c,h) for r,c in cases]
    receipt={"status":"PASS","note_sha256":NOTE_PIN,"helper_sha256":HELPER_PIN,
             "scope":"Finite exact rational-observable controls; the written proof supplies all-parameter quantifiers.",
             "case_count":len(records),"actual_factorial_moment_checks":9*len(records),
             "controls":["complete degree-bounded resultant interpolation",
                         "weighted polynomial basis reduction",
                         "actual critical-value and pole alignment",
                         "exact connection and derivative recurrence",
                         "first-independent-block minor and full cofactor relation",
                         "operator degree and zero-valuation bounds",
                         "degree-one distinct-coefficient window"],
             "cases":records,"sharp_SFC_claimed":False}
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True,default=qjson)+"\n")
    print(f"PASS: {len(records)} rational-observable cases; {9*len(records)} actual factorial moments; "
          "exact resultants, connections, minimal-order relations and degree bounds.")
    print("Receipt:",args.output)


if __name__=="__main__":
    try:
        main()
    except (ValueError,OSError) as error:
        raise SystemExit("REJECT: "+str(error)) from error

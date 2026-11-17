#!/usr/bin/env python3
"""Exact controls for QFL2; no degree search or noninvertible lift is constructed.

The general valuation, finite-type and approximation results are written
proofs. These controls distinguish bounded nonterminating scalar series,
unbounded spatial growth of an identity lift, and arithmetic center changes.
"""
from collections import defaultdict
from math import comb, prod
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
x,h=sp.symbols("x h")
u=sp.Function("u")(x)

def check_first_order(a,b,X):
    def operator(f):
        return h*a*sp.diff(f,x)+b*f
    return sp.simplify(operator(X*u)-X*operator(u)-h*u)==0

def pbw_mul(p,q,ell):
    out=defaultdict(int)
    for (a,b,c),v in p.items():
        for (d,e,f),w in q.items():
            for k in range(min(b,d)+1):
                factor=comb(b,k)*prod(range(d-k+1,d+1))
                out[a+d-k,b+e-k,c+f+k]+=v*w*factor
    return {m:c%ell for m,c in out.items() if c%ell}

def pbw_power(p,n,ell):
    out={(0,0,0):1}
    for _ in range(n):
        out=pbw_mul(out,p,ell)
    return out

def main():
    # Exact bounded-degree canonical pair over K((h)); h expansion infinite.
    g=1/(1-h)
    assert check_first_order(sp.Integer(1),g*x**2,x)
    assert sp.simplify((1-h)*g-1)==0
    assert sp.denom(g)==1-h
    assert sp.simplify(g.subs(h,2)+1)==0

    # Symmetric symbol a*p has normal PBW expression a*p+(h/2)*a'.
    a=(1-h**2*x)**2
    b=sp.expand(h*sp.diff(a,x)/2)
    X=x/(1-h**2*x)
    assert check_first_order(a,b,X)
    assert sp.simplify(a*sp.diff(X,x)-1)==0
    assert sp.expand(b+h**3*(1-h**2*x))==0
    assert X.subs(h,0)==x and a.subs(h,0)==1 and b.subs(h,0)==0
    # Only implementation controls; the geometric-series formula proves all r.
    series=sp.series(X,h,0,10).removeO().expand()
    for r in range(5):
        assert series.coeff(h,2*r)==x**(r+1)

    # Same classical reduction, different finite nonzero-h assignment.
    assert check_first_order(sp.Integer(1),h**2*x**2,x)
    center_checks=[]
    for ell in [5,7,11]:
        A={(0,1,0):1,(2,0,2):1}
        expected={(0,ell,0):1,(2*ell,0,2*ell):1}
        assert pbw_power(A,ell,ell)==expected
        center_checks.append(ell)

    result={
        "claim":"QFL2","status":"PASS",
        "bounded_nonterminating_pair":"Exact [P,X]=h on a generic test function; spatial degree2, scalar coefficient1/(1-h), regular at h=2",
        "unbounded_identity_lift":"Exact rational differential-operator identity; classical reduction identity. Written geometric-series formula gives unbounded coefficient degrees.",
        "arithmetic_shear":"Exact canonical pair with classical identity reduction; independent integer PBW powers modulo primes5,7,11",
        "prime_controls":center_checks,
        "coverage_boundary":"No bounded lift of a noninvertible classical map is asserted. Valuation/automatic nonsurjectivity, closed automorphism locus, Artin approximation and all-prime identity use the written proof.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                  [HERE/"BOUNDED_REALIZATION.md",HERE/"check_realization_controls.py"]}
    }
    (HERE/"realization_controls.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

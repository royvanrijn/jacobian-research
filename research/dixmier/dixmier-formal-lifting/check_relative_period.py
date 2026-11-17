#!/usr/bin/env python3
"""QRP1 exact controls; all-degree exhaustion is a written proof."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
x,y,b,r0,q,p,z,h,s = sp.symbols("x y b r0 q p z h s")
r=x*(1+x*y)
u=1+x*y
alpha=(4*y**2,2*x*y-1)

def theta(F):
    return sp.expand(x**2*sp.diff(F,x)-(1+2*x*y)*sp.diff(F,y))

def V(F):
    return sp.expand((1-2*x*y)*sp.diff(F,x)+4*y**2*sp.diff(F,y))

def normal(A,B=0):
    # Representatives modulo exact forms, then the complete weight recurrence.
    terms={}
    def add(a,j,c):
        if j:
            terms[a,j]=terms.get((a,j),sp.S.Zero)+c
    for (a,j),c in sp.Poly(sp.expand(A),x,y).terms():
        add(a,j,c)
    for (a,j),c in sp.Poly(sp.expand(B),x,y).terms():
        if a:
            add(a-1,j+1,-c*sp.Rational(a,j+1))
    result={}
    for (a,j),c in terms.items():
        d=a-j
        if d>=1 and j<=d:
            continue
        start=d+1 if d>=1 else max(1,-d)
        factor=sp.S.One
        for ell in range(start,j):
            factor*=-sp.Rational(ell+1,ell-d)
        result[d]=result.get(d,sp.S.Zero)+c*factor
    return {d:sp.cancel(c) for d,c in result.items() if sp.cancel(c)!=0}

def connection(A,B=0):
    vx=1-2*x*y
    vy=4*y**2
    return (sp.expand(V(A)+A*sp.diff(vx,x)+B*sp.diff(vy,x)),
            sp.expand(V(B)+A*sp.diff(vx,y)+B*sp.diff(vy,y)))

def zero(F):
    assert sp.cancel(F)==0,sp.cancel(F)

def main():
    zero(V(r)-1)
    zero(theta(r))
    zero(alpha[0]*theta(x)+alpha[1]*theta(y)-1)
    zero(theta(u)+r)
    R=x*(2-3*x*q)
    B=3*x**2*p+2*(1-3*x*q)*z
    # Normal-order commutators with positional polynomials are first-order exact.
    eta_x=sp.diff(B,p)
    eta_q=sp.diff(B,z)
    zero(R.subs(q,-2*y/3)-2*r)
    zero(eta_x/3-theta(x))
    zero((-sp.Rational(3,2)*eta_q/3).subs(q,-2*y/3)-theta(y))

    assert normal(*alpha)==normal(3*y**2)!={}
    assert normal(r*alpha[0],r*alpha[1])=={}
    zero(theta(-u)-r)
    zero(theta(-1/x)-1)
    relation_controls=0
    for a in range(6):
        for j in range(6):
            G=x**a*y**j
            assert normal(sp.diff(G,x),sp.diff(G,y))=={}
            assert normal(G*sp.diff(r,x),G*sp.diff(r,y))=={}
            assert normal(alpha[0]*theta(G),alpha[1]*theta(G))=={}
            relation_controls+=3

    eta=(y,sp.S.Zero)
    assert normal(*connection(*eta))==normal(3*y**2)
    assert normal(r*3*y**2)=={}
    for m in range(1,9):
        assert normal(*connection(y**m))==normal(sp.Rational(2*m*(2*m+1),m+1)*y**(m+1))
    for k in range(1,9):
        assert normal(r**k*y)==normal(-x**(2*k-1)*y**k)

    controls=[]
    for n in (1,2,3,5,8):
        En=sp.integrate(s**n*(s-1)**n,(s,0,sp.Symbol("U")))
        endpoint=sp.cancel(En.subs("U",1)-En.subs("U",0))
        assert endpoint==sp.Rational((-1)**n*factorial(n)**2,factorial(2*n+1))
        E=sp.expand(En.subs("U",u))
        omega=(y**(n+1),x*y**n)
        zero(r**n*omega[0]-sp.diff(E,x))
        zero(r**n*omega[1]-sp.diff(E,y))
        assert normal(r**n*omega[0],r**n*omega[1])=={}
        assert normal(r**(n-1)*omega[0],r**(n-1)*omega[1])!={}
        F=r*y**n
        assert normal(F*alpha[0],F*alpha[1])!={}
        zero(theta(E/r**n)+F)
        # Full polynomial canonical extension after the exact clearing power.
        G=-E
        shift_p=sp.diff(G,x)
        shift_w=sp.diff(G,y)
        zero(sp.diff(shift_w,x)-sp.diff(shift_p,y))
        zero(x**2*shift_p-(1+2*x*y)*shift_w-r**n*F)
        # Exact normalized Ore inner action on b: [H,b]/h=-theta(H).
        assert normal(r**n*F*alpha[0],r**n*F*alpha[1])=={}
        controls.append({"n":n,"endpoint_difference":str(endpoint),
                         "normal_form_before_last_clearing":
                             {str(k):str(v) for k,v in normal(r**(n-1)*omega[0],r**(n-1)*omega[1]).items()}})
    generic_eta=sp.cancel(y.subs(y,(r0-x)/x**2))
    assert sp.residue(generic_eta,x,0)==-1
    receipt={
        "claim":"QRP1","status":"PASS",
        "coordinate_translation_and_unimodularity":True,
        "exact_form_and_relative_relation_controls":relation_controls,
        "connection_and_cyclic_relation_controls":True,
        "torsion_samples":controls,
        "generic_eta_residue":-1,
        "canonical_extension_clearing_controls":True,
        "nonadmissible_Long_triples":"Written full-frame H1 argument and nonzero Delta_1 class; no explicit formal Long series replay.",
        "quantum_scope":"Exact Ore relation; normalized commutator with a base function is -theta(H). Innerness exhaustion is a written highest-b-order proof.",
        "coverage_boundary":"Finite symbolic controls only. Complete weight normal form, all-degree cohomology exhaustion, D-module presentation, arbitrary torsion order and global quantum vertical non-innerness and the fixed-position canonical extension criterion use the written proof. Nonadmission of the shifted all-orders Long triples uses the written infinitesimal full-frame argument, not a finite-series replay. No admitted bounded lift or DC2 conclusion follows.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("RELATIVE_PERIOD_TORSION.md","check_relative_period.py")}
    }
    (HERE/"relative_period_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()

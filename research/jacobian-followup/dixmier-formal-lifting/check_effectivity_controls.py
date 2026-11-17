#!/usr/bin/env python3
"""Exact finite controls for QEF1; the Ore-module and all-N results are written proofs."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
x,p,t=sp.symbols("x p t")

def normal(f,g):
    if f==0 or g==0:
        return sp.Integer(0)
    degree=sp.Poly(f,p).degree()
    return sp.expand(sum(t**j*sp.diff(f,p,j)*sp.diff(g,x,j)/sp.factorial(j)
                         for j in range(degree+1)))

def comm(f,g):
    return sp.expand(normal(f,g)-normal(g,f))

def power(f,n):
    out=sp.Integer(1)
    for _ in range(n):
        out=normal(out,f)
    return out

def main():
    X=x+t**2*p
    P=p
    assert comm(P,X)==t
    assert sp.expand(X-t**2*P-x)==0
    controls=[sp.Integer(1),x**2*p**3+x, x*p**5+t*p]
    for m in controls:
        G=sp.expand(normal(m,X)-normal(x,m))
        assert sp.expand(G-t*sp.diff(m,p)-t**2*p*m)==0
    Gpower=sp.Integer(1)
    for n in range(1,6):
        Gpower=sp.expand(t*sp.diff(Gpower,p)+t**2*p*Gpower)
        hermite=sum(sp.Rational(factorial(n),2**j*factorial(j)*factorial(n-2*j))*
                    t**(2*n-j)*p**(n-2*j) for j in range(n//2+1))
        assert sp.expand(Gpower-hermite)==0
        assert sp.Poly(Gpower,p).degree()==n
        assert sp.Poly(Gpower,p).LC()==t**(2*n)

    # Identities used by the generic simplicity proof of the boundary model.
    q=1-t**2*x
    imageX=x/q
    imageP=q**2*p-t**3*q
    assert sp.simplify(normal(q**-2,imageP)+t**3/q-p)==0
    for degree in range(1,5):
        lower=sp.cancel(normal(q,normal(p**degree,imageX))-normal(x,p**degree))
        lower_poly=sp.Poly(lower,p)
        assert lower_poly.degree()==degree-1
        assert sp.simplify(lower_poly.LC()-degree*t/q)==0

    seedP=p**2+x
    seedQ=-p
    assert comm(seedP,seedQ)==t
    profiles=[]
    for n in range(1,5):
        Qn=seedQ+t**2*power(seedP,n)
        assert comm(seedP,Qn)==t
        assert sp.Poly(Qn,p).degree()==2*n
        profiles.append([2,2*n])

    receipt={
        "claim":"QEF1","status":"PASS",
        "bounded_degree_one_canonical_automorphism":True,
        "PBW_graph_constraint_action":True,
        "boundary_image_contains_p":True,
        "boundary_order_lowering_controls":4,
        "graph_action_controls":len(controls),
        "Hermite_controls":5,
        "exact_order_jump_profiles":profiles,
        "coverage_boundary":"The rank-one dual criterion, algebraic dual vanishing, Ext obstruction, completion comparison, non-local-nilpotence for every n and arbitrary-order shear statement use written proofs. This checker verifies only exact symbolic identities and the listed finite controls.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in
                   [HERE/"ALGEBRAIC_DUAL_AND_POLARIZATION.md",HERE/"check_effectivity_controls.py"]}
    }
    (HERE/"effectivity_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()

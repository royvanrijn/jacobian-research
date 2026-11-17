#!/usr/bin/env python3
"""Exact symbolic identities for QSO1 and QGM1, without cover/degree searches."""
from math import comb
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent

def clean(a):
    return {i:c for i,v in a.items() if (c:=sp.simplify(v))!=0}

def add(*operators):
    out={}
    for a in operators:
        for i,c in a.items():
            out[i]=out.get(i,0)+c
    return clean(out)

def scale(a,c):
    return clean({i:c*v for i,v in a.items()})

def mul(a,b,var,weight=1):
    out={}
    for i,c in a.items():
        for j,d in b.items():
            for k in range(i+1):
                exponent=i+j-k
                coefficient=c*comb(i,k)*weight**k*sp.diff(d,var,k)
                out[exponent]=out.get(exponent,0)+coefficient
    return clean(out)

def comm(a,b,var,weight=1):
    return add(mul(a,b,var,weight),scale(mul(b,a,var,weight),-1))

def equal(a,b):
    return not add(a,scale(b,-1))

def main():
    z=sp.symbols("z")
    a,b,u,ell1,ell0=(sp.Function(n)(z) for n in ["a","b","u","ell1","ell0"])
    c0,c1=sp.symbols("c0 c1")
    top=comm({2:a},{3:b},z)
    assert sp.simplify(top[4]-(2*a*sp.diff(b,z)-3*b*sp.diff(a,z)))==0
    shifted={1:sp.Integer(1),0:-ell1/2}
    normalized=add(mul(shifted,shifted,z),scale(shifted,ell1),{0:ell0})
    assert equal(normalized,{2:1,0:ell0-ell1**2/4-sp.diff(ell1,z)/2})
    L={2:sp.Integer(1),0:u}
    Q={3:sp.Integer(1),1:sp.Rational(3,2)*u+c1,0:sp.Rational(3,4)*sp.diff(u,z)+c0}
    result=comm(L,Q,z)
    scalar=-sp.diff(u,z,3)/4-(sp.Rational(3,2)*u+c1)*sp.diff(u,z)
    assert equal(result,{0:scalar})
    w=u+sp.Rational(2,3)*c1
    assert sp.simplify(scalar+sp.diff(sp.diff(w,z,2)+3*w**2,z)/4)==0
    mu=sp.symbols("mu")
    assert sp.factor(mu*(mu-1)-12)==(mu-4)*(mu+3)
    assert sp.solve(mu*(mu-1)-12,mu)==[-3,4]
    assert sp.solve(mu*(mu-1),mu)==[0,1]
    # Sharpness: the (2,1) and (2,4) profiles really occur.
    sharpP={2:sp.Integer(1),0:z}
    sharpQ={1:sp.Integer(-1)}
    sharpQ4=add(sharpQ,mul(sharpP,sharpP,z))
    assert equal(comm(sharpP,sharpQ,z),{0:1})
    assert equal(comm(sharpP,sharpQ4,z),{0:1})
    assert max(sharpQ4)==4

    x,t=sp.symbols("x t")
    q=1-t**2*x
    imageX={0:x/q}
    imageP={1:q**2,0:-t**3*q}
    assert equal(comm(imageP,imageX,x,t),{0:t})
    assert sp.simplify(1+t**2*x/q-1/q)==0
    G={1:x**2,0:t*x}
    Gp=comm(G,{1:sp.Integer(1)},x,t)
    G2p=comm(G,Gp,x,t)
    assert equal(Gp,{1:-2*t*x,0:-t**2})
    assert equal(G2p,{1:2*t**2*x**2,0:2*t**3*x})
    assert not comm(G,G2p,x,t)
    conjugateP=add({1:sp.Integer(1)},scale(Gp,t),scale(G2p,t**2/2))
    assert equal(conjugateP,imageP)
    k=sp.symbols("k",integer=True,positive=True)
    assert equal(comm(G,{0:x**k},x,t),{0:k*t*x**(k+1)})

    files=[HERE/"ALGEBRAIC_CURVE_STRING_OBSTRUCTION.md",
           HERE/"FORMAL_GRAPH_BOUNDARY_MODEL.md",HERE/"check_curve_and_module.py"]
    receipt={
        "claims":["QSO1","QGM1"],"status":"PASS",
        "Ore_leading_coefficient_normalization":True,
        "generic_string_commutator_and_Painleve_reduction":True,
        "pole_indicial_roots":[-3,4],"regular_indicial_roots":[0,1],
        "sharp_allowed_order_profiles":[[2,1],[2,4]],
        "rational_quantum_graph_commutator":True,
        "cyclic_generator_identity":True,
        "formal_inner_conjugation_identities":True,
        "coverage_boundary":"All algebraic curves and coefficient degrees, the no-algebraic-Painleve lemma, strict infinite module chain, completion and generic holonomicity are written proofs. This checker verifies symbolic identities, not those global implications.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    }
    (HERE/"curve_module_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()

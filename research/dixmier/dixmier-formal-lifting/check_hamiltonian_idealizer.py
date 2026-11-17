#!/usr/bin/env python3
"""Exact Ore identities for QHI1; universal symbol and exclusion results are written."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
u,v,kappa,a=sp.symbols("u v kappa a", nonzero=True)

def mul(f,g):
    f,g=sp.expand(f),sp.expand(g)
    if f==0 or g==0:
        return sp.Integer(0)
    return sp.expand(sum(sp.diff(f,v,j)*sp.diff(g,u,j)/sp.factorial(j)
                         for j in range(sp.Poly(f,v).degree()+1)))

def comm(f,g):
    return sp.cancel(mul(f,g)-mul(g,f))

def zero(f):
    assert sp.cancel(f)==0,sp.cancel(f)

def degree_symbol(f):
    # A common u power clears Laurent denominators; division-ring Bernstein degree.
    f=sp.cancel(f)
    numerator,denominator=sp.fraction(f)
    polynomial=sp.Poly(numerator,u,v)
    denominator_poly=sp.Poly(denominator,u,v)
    top=polynomial.total_degree()
    bottom=denominator_poly.total_degree()
    lead_num=sum(c*u**i*v**j for (i,j),c in polynomial.terms() if i+j==top)
    lead_den=sum(c*u**i*v**j for (i,j),c in denominator_poly.terms() if i+j==bottom)
    return int(top-bottom),sp.cancel(lead_num/lead_den)

def main():
    m1,m2=u,v+1/u
    H,L=u*v,v**2-2/u**2
    multipliers=[
        (u**2,(u**2,0),(mul(v,u)+1,0)),
        (u**3,(u**3,0),(mul(v,u**2)+u,0)),
        (H,(-1,u**2),(v**2,0)),
        (L,(-v**2,2*u*v),(0,v**2))
    ]
    for r,left1,left2 in multipliers:
        for m,coeff in ((m1,left1),(m2,left2)):
            zero(mul(m,r)-mul(coeff[0],m1)-mul(coeff[1],m2))
    e,f0,h=u**2/2,-L/2,H+sp.Rational(1,2)
    zero(comm(h,e)-2*e)
    zero(comm(h,f0)+2*f0)
    zero(comm(e,f0)-h)

    ladder=[]
    current=u**3
    for index in range(4):
        deg,symbol=degree_symbol(current)
        assert deg==3
        expected=sp.prod([2*(3-j) for j in range(index)])*u**(3-index)*v**index
        zero(symbol-expected)
        ladder.append(str(symbol))
        current=comm(L,current)
    zero(current)

    # For a binary monomial, Hamiltonian uv measures the exponent difference.
    for i,j in ((0,3),(1,1),(3,0),(4,2)):
        zero(comm(H,u**i*v**j)-(i-j)*u**i*v**j)
    Q=-kappa*v/(2*a*u)
    zero(comm(a*u**2,Q)-kappa)
    assert degree_symbol(Q)[0]==0

    receipt={
        "claim":"QHI1","status":"PASS",
        "right_multiplier_membership_identities":8,
        "sl2_relations":True,
        "complete_cubic_symbol_ladder":ladder,
        "quadratic_symbol_weight_controls":4,
        "scalar_commutator_formula":True,
        "coverage_boundary":"Exact rational differential-operator identities only. Hamiltonian preservation of normalized symbol ideals, transversality, the all-degree idealizer equality and exclusion with an unrestricted partner use written proofs. No right Weyl pair or Long model is constructed.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("HAMILTONIAN_IDEALIZER_CRITERION.md","check_hamiltonian_idealizer.py")}
    }
    (HERE/"hamiltonian_idealizer_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()

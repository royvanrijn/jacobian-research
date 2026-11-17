#!/usr/bin/env python3
"""QCR1 exact controls; centralizer exhaustion and normalization optimality are written."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
AUDIT=HERE.parent/"dixmier-audit"
sys.path.insert(0,str(AUDIT))
import verify_long_pbw as poly

def projection(F,R):
    result=dict(F)
    for exponent,coefficient in F.items():
        i,j,k,l=exponent
        if k==l==0 and i==2*j:
            result=poly.add(result,poly.scale(poly.power(R,j),-coefficient/Fraction((-3)**j)))
    return result

def degree(F):
    return max(map(sum,F),default=-1)

def main():
    R,T,D,S=poly.long_images()
    assert len(D)==139
    assert degree(R)==3 and degree(D)==23
    assert max(e[2]+e[3] for e in D)==5
    assert D[(12,6,0,5)]==Fraction(144,5)
    assert D[(6,0,0,5)]==Fraction(32,405)
    assert projection(D,R)==D
    controls=[]
    for k in (1,3,8):
        central=poly.add(poly.power(R,k),poly.scale(R,Fraction(2,3)),poly.const(7))
        shifted=poly.add(D,central)
        assert projection(central,R)=={}
        assert projection(shifted,R)==D
        assert degree(shifted)==max(23,3*k)
        assert degree(projection(shifted,R))<=degree(shifted)
        assert poly.poisson(central,R)=={}
        assert poly.poisson(central,S)=={}
        assert poly.poisson(central,T)=={}
        controls.append({"R_degree":k,"completion_degree":degree(shifted)})
    trial=poly.add(D,poly.power(R,4),poly.mul(poly.var(0),poly.var(1)),poly.const(3))
    normalized=projection(trial,R)
    assert projection(normalized,R)==normalized
    assert degree(normalized)<=degree(trial)

    x,p,h,Z=sp.symbols("x p h Z")
    X=x+h**2*x**2
    inverse_factor=1/(1+2*h**2*x)
    # P has p order one, so its commutator with positional X is exact.
    assert sp.cancel(inverse_factor*sp.diff(X,x)-1)==0
    catalan=sum((-1)**k*sp.catalan(k)*h**(2*k)*Z**(k+1) for k in range(5))
    error=sp.Poly(sp.expand(catalan.subs(Z,X)-x),h)
    assert all(exponent[0]>=10 for exponent,coefficient in error.terms() if coefficient)
    assert sp.Poly(X,x).degree()==2
    assert max(e[2]+e[3] for e in normalized)==5

    receipt={
        "claim":"QCR1","status":"PASS",
        "Long_D_terms":len(D),
        "Long_D_spatial_degree":degree(D),
        "Long_D_momentum_order":5,
        "Long_D_already_normalized":True,
        "degree23_witness":[[12,6,0,5],"144/5"],
        "momentum_order5_witness":[[6,0,0,5],"32/405"],
        "central_shift_controls":controls,
        "projection_idempotence":True,
        "deformed_coordinate_exact_commutator":True,
        "Catalan_inverse_control_mod_h":10,
        "fixed_X_no_bounded_mate":"written p-order and polynomial-unit proof, including spectators",
        "coverage_boundary":"Exact finite controls only. The all-rank classical centralizer theorem uses published polynomial-dependence theory; formal centralizer exhaustion, all-degree ambiguity, projection optimality and the fixed-triple boundedness equivalence are written proofs. No admission or boundedness for an arbitrary quantum triple, arbitrary gauge invariance or bounded Long lift is established.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("CENTRALIZER_NORMAL_FORM.md","check_centralizer_normal_form.py")},
        "retained_input":{
            "path":"jacobian-followup/dixmier-audit/verify_long_pbw.py",
            "sha256":hashlib.sha256(Path(poly.__file__).read_bytes()).hexdigest(),
            "R_sha256":poly.digest(R),"D_sha256":poly.digest(D)}
    }
    (HERE/"centralizer_normal_form_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()

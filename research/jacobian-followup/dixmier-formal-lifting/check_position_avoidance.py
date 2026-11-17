#!/usr/bin/env python3
"""Exact finite controls for QPA1; universal intersections use the written proof."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s

HERE=Path(__file__).resolve().parent

def main():
    r,v,a,z,X,ss,tt=s.symbols("r v a z X ss tt")
    # r is inverted in the coefficient field; no specialization at r=0 is used.
    S=2*v+4*a-3*r*a*a
    T=(r*a**3-a*a-v*a)/2
    f=r*X**3-2*X**2+ss*X+4*tt
    actual=s.expand(f.subs({ss:S,tt:T}))
    assert s.expand(actual-(X-a)*(r*X**2+(r*a-2)*X+2*(v+a-r*a*a)))==0
    assert s.expand(actual.subs(X,a))==0
    assert s.expand(s.diff(actual,X).subs(X,a)-2*v)==0
    jac=s.det(s.Matrix([[s.diff(S,v),s.diff(S,a)],[s.diff(T,v),s.diff(T,a)]]))
    assert s.expand(jac+v)==0
    assert s.cancel((s.diff(S,a)*s.diff(T,v)-s.diff(S,v)*s.diff(T,a))/v)==1
    disc=s.discriminant(f,X)
    expected=4*ss**2-4*r*ss**3+128*tt-144*r*ss*tt-432*r*r*tt*tt
    assert s.expand(disc-expected)==0
    assert s.expand(disc.subs({ss:S,tt:T})-4*v*v*((3*r*a-2)**2-8*r*v))==0

    S0=4*z-3*r*z*z
    T0=(r*z**3-z*z)/2
    b=2/r-2*z
    v1=(2-3*r*z)**2/(2*r)
    assert s.expand(f.subs({ss:S0,tt:T0})-r*(X-z)**2*(X-b))==0
    for image,wanted in [(S,S0),(T,T0)]:
        assert s.expand(image.subs({a:z,v:0})-wanted)==0
        assert s.cancel(image.subs({a:b,v:v1})-wanted)==0
    assert s.cancel(8*r*v1-(3*r*b-2)**2)==0
    assert s.diff(v1,z)!=0
    # Generic discriminant simple/double orders, without a point census.
    assert s.rem(s.Poly(s.diff(expected,tt).subs({ss:S,tt:T}),v),s.Poly(v,v)).as_expr()!=0

    spec=importlib.util.spec_from_file_location("qpa1_ore",HERE/"check_curve_and_module.py")
    ore=importlib.util.module_from_spec(spec);spec.loader.exec_module(ore)
    x=s.symbols("x")
    # Sharpness: without the branch condition the usual pair contains a position.
    assert ore.equal(ore.comm({1:s.Integer(1)},{0:x},x),{0:1})
    # First-order mate reduction is an exact finite identity here.
    P={1:(1+x)**2}
    g=x/(1+x)
    Q=ore.add({0:g},ore.mul(P,P,x))
    assert ore.equal(ore.comm(P,Q,x),{0:1})

    # The positive-order commutator with a position is nonzero.
    u=s.Function("u")(x)
    for m in [1,2,3,5]:
        result=ore.comm({m:u},{0:x},x)
        assert ore.equal(result,{m-1:m*u})
    # A local-nilpotence negative control after inverting x.
    assert s.diff(1/x,x,5)==-120/x**6

    old=json.loads((HERE/"quadratic_mate_controls.json").read_text())
    for name,digest in old["hashes"].items():
        p=HERE/name if name!="verify_long_pbw.py" else HERE.parent/"dixmier-audit"/name
        assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
    files=[HERE/"RAMIFICATION_AND_POSITION_AVOIDANCE.md",Path(__file__),
           HERE/"QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md",HERE/"quadratic_mate_controls.json",
           HERE.parent/"dixmier-audit"/"verify_long_pbw.py",HERE/"check_curve_and_module.py"]
    receipt={
        "claim":"QPA1","status":"PASS",
        "cubic_equation_and_derivative":True,"generic_Jacobian":"-v",
        "generic_residue_cover_degree":{"value":3,"evidence":"written monic-cubic function-field argument"},"both_boundary_parameterizations":True,
        "discriminant_double_and_simple_factorization":True,
        "boundary_v1":"(2-3*r*z)^2/(2*r)","v1_is_nonconstant":True,
        "first_order_normal_form_control":True,
        "positive_operator_order_commutator_controls":[1,2,3,5],
        "localization_can_destroy_local_nilpotence_control":True,
        "retained_Long_localization_receipt_hashes":"PASS",
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "coverage_boundary":"Finite symbolic identities only. The two-branch field lemma, completed quantum intersection, full positional-image theorem, all-source low-order exclusion and ambient ad-nilpotence classification use written proofs, not finite degree tests.",
        "independent_replay":False,"formal_verification":False,"external_review":False
    }
    (HERE/"position_avoidance_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: cubic cover, discriminant, two boundary traces, first-order normal form, local-nilpotence controls and retained input hashes")
    print("Universal field and quantum-polarization statements use the written proof.")

if __name__=="__main__":
    main()

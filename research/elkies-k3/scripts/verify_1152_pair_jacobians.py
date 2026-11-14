#!/usr/bin/env python3
"""Independent rational-polynomial replay of15 paired-base Jacobian models.

The only imported arithmetic is standard-library Fraction polynomial operations.
The classical binary-quartic Jacobian theorem is an explicit external premise.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PILOT=ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2"
DATA=ROOT/"elkies-k3/data/external/hill2026/data"
DEFAULT=ROOT/"artifacts/generated-results/x1152-pair-jacobians-20261114-v1.json"
spec=importlib.util.spec_from_file_location("paired_polynomials",Path(__file__).with_name("verify_r17_correlated_genus_one_pencils.py"))
V=importlib.util.module_from_spec(spec);spec.loader.exec_module(V)
require=V.require


def invariants(cs):
    f=list(cs)+[Q(0)]*(5-len(cs));e,d,c,b,a=f
    return 12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c


def conic_pull(q,N,D):
    q=list(q)+[Q(0)]*(3-len(q))
    return V.add(V.add(V.scale(V.power(N,2),q[2]),V.scale(V.mul(N,D),q[1])),V.scale(V.power(D,2),q[0]))


def reverse_mobius(f,c,weight):
    # U^weight f(c+1/U), fixed weight >= deg f.
    out=[]
    for i,a in enumerate(f):
        term=V.mul([Q(0)]*(weight-i)+[Q(1)],V.power([Q(1),c],i))
        out=V.add(out,V.scale(term,a))
    return out


def verify(result):
    require(result["schema"]=="hill1152-paired-jacobians-2isogenies-v1","schema")
    for name,key in [("result.json","pilot_result_sha256"),("four-generic-covers.json","cover_certificate_sha256")]:
        require(sha256((PILOT/name).read_bytes()).hexdigest()==result[key],"input hash")
    pilot=json.loads((PILOT/"result.json").read_text());covercert=json.loads((PILOT/"four-generic-covers.json").read_text())
    c=Q(pilot["common_chart"]["c"])
    expectedmaps={f["index"]:[V.poly(f[k]) for k in ["N","D","W"]] for f in covercert["families"]}
    for index,published in [(13,4339),(14,36773)]:
        path=DATA/f"cover-{published}.json"
        require(sha256(path.read_bytes()).hexdigest()==result["source_hashes"][path.name],"published map hash")
        d=json.loads(path.read_text())["quadratic_basechange"]
        oldN,oldD,oldW=[V.poly(d[k]) for k in ["parameter_numerator","parameter_denominator","branch_root_numerator"]]
        expectedmaps[index]=[oldD,V.sub(oldN,V.scale(oldD,c)),oldW]
    require([r["pair"] for r in result["pairs"]]==pilot["locally_surviving_disjoint_pairs"]
            and len(result["pairs"])==15,"pair roster")
    for row in result["pairs"]:
        i,j=row["pair"]
        q,r=[V.poly(pilot["quadratic_covers"][n]["q"]) for n in [i,j]]
        require(len(q)==len(r)==3,"quadratic degrees")
        require(len(V.gcd_poly(q,r))==1,"shared branch")
        for f in [q,r]:
            require(len(V.gcd_poly(f,[k*x for k,x in enumerate(f)][1:]))==1,"multiple branch root")
        a=q[2]*r[0]+q[0]*r[2]-q[1]*r[1]/2
        discq=q[1]**2-4*q[2]*q[0];discr=r[1]**2-4*r[2]*r[0]
        b=discq*discr/16;res=a*a-4*b
        require(b and res,"singular Jacobian")
        require([str(x) for x in [0,a,0,b,0]]==row["jacobian_a_invariants"],"Jacobian coefficients")
        require(list(map(Q,row["quadratic_discriminants"]))==[discq,discr],"quadratic discriminants")
        require(Q(row["radicand_resultant"])==res,"resultant relation")
        N,D,W=[V.poly(row["conic_map"][k]) for k in ["N","D","W"]]
        require([N,D,W]==expectedmaps[i],"unbound conic map")
        require(len(V.gcd_poly(N,D))==1,"conic map common factor")
        require(conic_pull(q,N,D)==V.power(W,2),"conic identity")
        F=conic_pull(r,N,D)
        require(F==V.poly(row["paired_base_quartic"]) and len(F) in [4,5],"paired quartic")
        require(len(V.gcd_poly(F,[k*x for k,x in enumerate(F)][1:]))==1,"singular paired quartic")
        I,J=invariants(F)
        require([Q(row["quartic_I"]),Q(row["quartic_J"])]==[I,J],"quartic invariants")
        k=Q(row["jacobian_to_quartic_scale"])
        require(k and k**4*(b-a*a/3)==-27*I and k**6*(2*a**3/27-a*b/3)==-27*J,"twist-sensitive Jacobian isomorphism")
        Ip,Jp=invariants(V.mul(q,r))
        require(Ip==a*a+12*b and Jp==2*a*(36*b-a*a),"product invariants")
        require([Q(row["product_quartic_I"]),Q(row["product_quartic_J"])]==[Ip,Jp],"product invariant labels")
        require(list(map(Q,row["quotient_a_invariants"]))==[0,-2*a,0,res,0],"quotient coefficients")
        # Clear x^3 from the 2-isogeny identity. No numerical sampling.
        T=[b,a,Q(1)]
        left=V.mul(T,V.power([-b,Q(0),Q(1)],2))
        right=V.add(V.add(V.power(T,3),V.scale([Q(0)]+V.power(T,2),-2*a)),V.scale([Q(0),Q(0)]+T,res))
        require(left==right,"2-isogeny polynomial identity")
        # The fixed map X=9*x-6*a,Y=27*y gives the quartic-product Jacobian.
        require(81*(res-4*a*a/3)==-27*Ip and
                729*(-16*a**3/27+2*a*res/3)==-27*Jp,"quotient-to-product Jacobian isomorphism")
        require(row["rational_points_on_paired_base"]=="UNKNOWN" and row["complete_selmer_bound"]=="UNKNOWN","global arithmetic overclaim")
    return {"schema":"hill1152-paired-jacobians-independent-replay-v1","status":"PASS",
        "paired_genus_one_bases":15,"Jacobian_Q_models":15,"explicit_degree2_isogenies":15,
        "complete_selmer_spaces":0,"new_torsor_points":0,
        "boundary":"All polynomial identities and Q-isomorphism scales checked. Classical binary-quartic Jacobian theorem assumed. No global Selmer completeness, torsor solubility, rank or new rational point follows."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--certificate",type=Path,default=DEFAULT)
    p.add_argument("--output",type=Path);a=p.parse_args()
    result=verify(json.loads(a.certificate.read_text()))
    if a.output:
        with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps(result,sort_keys=True))

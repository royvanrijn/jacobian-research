#!/usr/bin/env python3
"""Independent character sums and Frobenius replay for compatible triples."""
import argparse
import base64
from pathlib import Path
import struct
import subprocess
import time
import numpy as np
import retrospective as r
import equal_class_picard as ex
from verify_mixed_character import Field

OUTPUT=r.OUT/"rank_jump_equal_class_verification_v1.json"


def field_check(case,p,saved,rec,old):
    row=r.read(ex.INPUT)["cases"][case]
    A,B=[r.mod(v,p) for v in row["model"][3:]]
    degree=saved["degree"];q=saved["q"]
    F=Field(p,degree,saved["modulus_ascending"])
    x=np.arange(q,dtype=np.int64);xx=F.mul(x,x);xxx=F.mul(xx,x)
    chi=F.power(x,(q-1)//2)
    assert set(map(int,np.unique(chi)))=={0,1,p-1}
    chi=np.where(chi==p-1,-1,chi)
    packed=base64.b64decode(saved["trace_values_i16_le_base64"])
    traces=np.array(struct.unpack("<"+"h"*q,packed),dtype=np.int64)
    direct=0
    if not old:
        y0=F.add(F.add(xxx,F.scale(x,A)),B)
        y1=F.add(F.scale(xx,2*A),F.scale(x,3*B))
        y2=F.add(F.scale(x,A*A),A*B%p)
        frobenius=F.power(x,p);seen=np.zeros(q,dtype=bool)
        assert len(np.unique(frobenius))==q
        for u in range(q):
            if seen[u]:continue
            orbit=[];v=u
            while not seen[v]:
                orbit.append(v);seen[v]=True;v=int(frobenius[v])
            assert v==u and degree%len(orbit)==0
            uu=F.mul(u,u);uuu=F.mul(uu,u)
            values=F.add(F.add(F.add(y0,F.mul(u,y1)),F.mul(uu,y2)),F.scale(uuu,-B*B))
            value=int(chi[values].sum())
            assert all(value==int(traces[v]) for v in orbit)
            direct+=1
        assert seen.all()
    xs=[r.mod(P[0],p) for P in row["generic_points"]]
    d=np.ones(q,dtype=np.int64)
    for a in xs:d=F.mul(d,F.add(1,F.scale(x,-a)))
    finite=int((chi[d]*traces).sum());k=-xs[0]*xs[1]*xs[2]%p
    values=F.add(F.add(F.add(xxx,F.scale(xx,2*A*k)),F.scale(x,A*A*k*k)),(-B*B*k**3)%p)
    infinity=int(chi[values].sum())
    expected=next(v for v in rec["fields"] if v["degree"]==degree)
    assert (finite,infinity,finite+infinity)==(expected["finite_sum"],expected["infinity_sum"],expected["residual_H2_trace"])
    return {"degree":degree,"base_parameters":q,"direct_orbit_character_sums":direct,
            "reused":bool(old),"trace_sha256":r.digest(packed),"finite_sum":finite,"infinity_sum":infinity}


def prime_check(case,p):
    counts=r.read(ex.RAW);report=r.read(ex.OUTPUT)
    raw=next(v for v in counts["records"] if (v["case"],v["p"])==(case,p))
    rec=next(v for v in report["reductions"] if (v["case"],v["p"])==(case,p))
    assert raw["status"]=="COUNTED"
    prior_path=r.OUT/"rank_jump_triple_character_verification_v1.json"
    prior=r.read(prior_path)
    assert prior["status"]=="PASS" and prior["counts_sha256"]==r.digest(ex.previous.RAW.read_bytes())
    old=next((v for v in r.read(ex.previous.RAW)["records"] if (v["case"],v["p"])==(case,p)),None)
    if old:
        assert raw["fields"]==old["fields"]
        certified=next(v for v in prior["records"] if (v["case"],v["p"])==(case,p))
    fields=[]
    for saved in raw["fields"]:
        checked=field_check(case,p,saved,rec,old)
        if old:
            previous=next(v for v in certified["fields"] if v["degree"]==saved["degree"])
            assert checked["trace_sha256"]==previous["trace_sha256"]
        fields.append(checked)
        r.write_new(ex.WORK/f"replay-case{case}-p{p}-degree{saved['degree']}.json",checked)
        print("PASS field",case,p,saved["degree"],flush=True)
    return {"case":case,"p":p,"fields":fields,"status":"PASS"}


def geometry_check():
    from sage.all import QQ,ZZ,PolynomialRing,companion_matrix,identity_matrix,matrix
    inputs=r.read(ex.INPUT);gate=r.read(ex.gate.OUTPUT);original=r.read(ex.previous.INPUT)
    assert inputs["component_certificate_sha256"]==r.digest(ex.gate.OUTPUT.read_bytes())
    assert inputs["original_triple_input_sha256"]==r.digest(ex.previous.INPUT.read_bytes())
    for row,old,proof in zip(inputs["cases"],original["cases"],gate["cases"]):
        assert row["model"]==old["model"] and row["id"]==old["id"]
        assert row["generic_points"]==[old["generic_points"][0]]+[v["point"] for v in proof["same_class_control"]["derived_points"]]
        assert abs(matrix(ZZ,row["generic_combinations"]).det())==4
        assert proof["same_class_control"]["finite_additive_component_gate"]=="PASS"
    R=PolynomialRing(QQ,"X");X=R.gen();checks=[]
    for rec in r.read(ex.OUTPUT)["reductions"]:
        p=rec["p"];pol=R(list(map(QQ,rec["polynomial_ascending"])))
        M=companion_matrix(pol)
        assert [int((M**i).trace()) for i in (1,2,3)]==rec["traces"]
        normpol=R(pol(p*X)/p**5);cycles=0
        for factor,mult in normpol.factor():
            if all(c.denominator()==1 for c in factor.list()) and PolynomialRing(ZZ,"X")(factor).is_cyclotomic():
                cycles+=int(factor.degree())*int(mult)
        assert 17+cycles==rec["reduction_geometric_Picard_rank"]
        if rec["status"]=="RHO_18_REDUCTION":
            root=next(e*p for e in (-1,1) if pol(e*p)==0)
            V=companion_matrix(pol//(X-root))/p
            determinant=(identity_matrix(QQ,4)-V**6).det()
            assert determinant==QQ(rec["normalized_transcendental_degree6"])
            signed=-determinant/QQ(p**6)
            assert int(ZZ(signed.numerator()*signed.denominator()).squarefree_part())==rec["NS_discriminant_squareclass"]
        checks.append({"case":rec["case"],"p":p,"PASS":True})
    return checks


def verify():
    records=[];started=time.monotonic()
    for raw in r.read(ex.RAW)["records"]:
        case,p=raw["case"],raw["p"];dest=ex.WORK/f"replay-case{case}-p{p}.json"
        if not dest.exists():
            log=ex.WORK/f"replay-case{case}-p{p}.log"
            with log.open("x") as handle:
                try:
                    result=subprocess.run(["python3",str(Path(__file__).resolve()),"prime","--case",str(case),
                                           "--prime",str(p),"--destination",str(dest)],cwd=r.ROOT,stdout=handle,stderr=handle,timeout=60)
                    assert result.returncode==0,log.read_text()
                except subprocess.TimeoutExpired:
                    r.write_new(dest,{"status":"UNKNOWN","case":case,"p":p,"reason":"60-second independent replay cap"})
        record=r.read(dest);records.append(record);print("checkpoint replay",case,p,record["status"],flush=True)
    assert all(v["status"]=="PASS" for v in records),"Independent count replay incomplete"
    checks=geometry_check()
    return {"schema":"rank-jump.equal-class-verification.v1","status":"PASS","records":records,
            "frobenius_checks":checks,"analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),
            "counts_sha256":r.digest(ex.RAW.read_bytes()),"verifier_sha256":r.digest(Path(__file__).read_bytes()),
            "field_arithmetic_sha256":r.digest((ex.HERE/"verify_mixed_character.py").read_bytes()),
            "prior_verification_sha256":r.digest((r.OUT/"rank_jump_triple_character_verification_v1.json").read_bytes()),
            "numpy":np.__version__,"elapsed_seconds":round(time.monotonic()-started,3)}


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["prime","verify"])
    parser.add_argument("--case",type=int);parser.add_argument("--prime",type=int)
    parser.add_argument("--destination",type=Path);args=parser.parse_args()
    result=prime_check(args.case,args.prime) if args.mode=="prime" else verify()
    if args.destination:r.write_new(args.destination,result)

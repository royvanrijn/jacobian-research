#!/usr/bin/env python3
"""Independent new trace arrays and all character weights, including infinity."""
import argparse
import base64
from pathlib import Path
import struct
import numpy as np
import retrospective as r
import triple_character as ex
from verify_mixed_character import Field


def verify():
    counts=r.read(ex.RAW);report=r.read(ex.OUTPUT);cases=r.read(ex.INPUT)["cases"]
    previous=r.read(ex.previous.COUNTS)
    prior=r.read(r.OUT/"rank_jump_mixed_character_verification_v1.json")["records"]
    records=[];new_checked=0;reused=0
    for raw in counts["records"]:
        assert raw["status"]=="COUNTED"
        case,p=raw["case"],raw["p"];row=cases[case]
        A,B=[r.mod(x,p) for x in row["model"][3:]]
        old=next((x for x in previous["reductions"] if x["case"]==case and x["p"]==p),None)
        if old:
            assert raw["fields"]==old["fields"]
            certified=next(x for x in prior if x["case"]==case and x["p"]==p)
            assert certified["status"]=="PASS" and certified["counts_sha256"]==r.digest(ex.previous.COUNTS.read_bytes())
        fields=[]
        for saved in raw["fields"]:
            degree=saved["degree"];q=saved["q"];F=Field(p,degree,saved["modulus_ascending"])
            x=np.arange(q,dtype=np.int64);xx=F.mul(x,x);xxx=F.mul(xx,x)
            chi=F.power(x,(q-1)//2);assert set(map(int,np.unique(chi)))=={0,1,p-1}
            chi=np.where(chi==p-1,-1,chi)
            packed=base64.b64decode(saved["trace_values_i16_le_base64"])
            traces=np.array(struct.unpack("<"+"h"*q,packed),dtype=np.int64)
            if old:
                field_certificate=next(x for x in certified["fields"] if x["degree"]==degree)
                assert r.digest(packed)==field_certificate["trace_sha256"]
                reused+=q
            else:
                y0=F.add(F.add(xxx,F.scale(x,A)),B)
                y1=F.add(F.scale(xx,2*A),F.scale(x,3*B))
                y2=F.add(F.scale(x,A*A),A*B%p)
                for u in range(q):
                    uu=F.mul(u,u);uuu=F.mul(uu,u)
                    values=F.add(F.add(F.add(y0,F.mul(u,y1)),F.mul(uu,y2)),F.scale(uuu,-B*B))
                    assert int(chi[values].sum())==int(traces[u])
                new_checked+=q
            weights=[]
            for rec in [v for v in report["reductions"] if v["case"]==case and v["p"]==p and "fields" in v]:
                mask=rec["mask"];xs=[r.mod(P[0],p) for i,P in enumerate(row["generic_points"]) if mask>>i&1]
                d=np.ones(q,dtype=np.int64)
                for a in xs:d=F.mul(d,F.add(1,F.scale(x,-a)))
                finite=int((chi[d]*traces).sum());infinity=0
                if len(xs)==3:
                    k=-xs[0]*xs[1]*xs[2]%p
                    values=F.add(F.add(F.add(xxx,F.scale(xx,2*A*k)),
                                      F.scale(x,A*A*k*k)),(-B*B*k**3)%p)
                    infinity=int(chi[values].sum())
                expected=next(v for v in rec["fields"] if v["degree"]==degree)
                assert (finite,infinity)==(expected["finite_sum"],expected["infinity_sum"])
                weights.append({"mask":mask,"finite_sum":finite,"infinity_sum":infinity})
            fields.append({"degree":degree,"base_parameters":q,"trace_sha256":r.digest(packed),"weights":weights})
        records.append({"case":case,"p":p,"fields":fields})
        print("PASS triple weights and traces",case,p,flush=True)
    out={"schema":"rank-jump.triple-character-verification.v1","status":"PASS",
         "new_directly_counted_base_parameters":new_checked,"previous_independently_verified_base_parameters":reused,
         "records":records,"counts_sha256":r.digest(ex.RAW.read_bytes()),"analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),
         "previous_verification_sha256":r.digest((r.OUT/"rank_jump_mixed_character_verification_v1.json").read_bytes()),
         "verifier_sha256":r.digest(Path(__file__).read_bytes()),
         "field_arithmetic_source_sha256":r.digest((ex.HERE/"verify_mixed_character.py").read_bytes()),
         "numpy":np.__version__}
    print("PASS",new_checked,"new and",reused,"reused base parameters")
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--destination",type=Path);args=p.parse_args()
    out=verify()
    if args.destination:r.write_new(args.destination,out)

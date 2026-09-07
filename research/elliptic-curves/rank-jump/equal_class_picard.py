#!/usr/bin/env python3
"""Frozen Picard test of the two equal-class triple controls."""
import argparse
import base64
from pathlib import Path
import struct
import subprocess
import retrospective as r
import triple_character as previous
import component_kummer_gate as gate

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"EQUAL_CLASS_PICARD_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_equal_class_inputs_v1.json"
RAW=r.OUT/"rank_jump_equal_class_counts_v1.json"
OUTPUT=r.OUT/"rank_jump_equal_class_picard_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-equal-class-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,INPUT,previous.RAW,Path(previous.__file__),gate.OUTPUT)}


def worker(case,p):
    current=r.read(INPUT)["cases"][case]
    old=r.read(previous.INPUT)["cases"][case]
    assert current["model"]==old["model"] and p==r.read(PROTOCOL)["limits"]["new_prime"]
    # Reuse the sealed untwisted counter on the identical anchor equation.
    # Only its local checkpoint destination changes in this process.
    previous.WORK=WORK
    out=previous.new_count(case,p)
    producer=out.pop("bindings")
    r.write_new(WORK/f"case{case}-p{p}.json",{**out,"bindings":bindings(),"producer_bindings":producer})


def capture():
    WORK.mkdir(parents=True,exist_ok=True);out=[]
    old=r.read(previous.RAW)
    for case,row in enumerate(r.read(INPUT)["cases"]):
        for p in r.read(PROTOCOL)["primes"][row["role"]]:
            assert previous.good(row,p,7)
            retained=next((x for x in old["records"] if x["case"]==case and x["p"]==p),None)
            if retained is not None:
                assert retained["status"]=="COUNTED"
                out.append({"case":case,"p":p,"status":"COUNTED","fields":retained["fields"],
                            "origin":"retained untwisted arrays","source_sha256":r.digest(previous.RAW.read_bytes())})
                continue
            path=WORK/f"case{case}-p{p}.json"
            if not path.exists():
                with (WORK/f"case{case}-p{p}.log").open("x") as log:
                    try:
                        result=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                               "--case",str(case),"--prime",str(p)],cwd=r.ROOT,
                                              stdout=log,stderr=log,timeout=40)
                        error=None if result.returncode==0 else "worker failure"
                    except subprocess.TimeoutExpired:error="40-second cap"
                if error:r.write_new(path,{"case":case,"p":p,"status":"UNKNOWN","bindings":bindings(),"reason":error,
                                          "transcript":(WORK/f"case{case}-p{p}.log").read_text()})
            record=r.read(path);assert record["bindings"]==bindings();out.append(record)
            print("checkpoint",case,p,record["status"],flush=True)
    r.write_new(RAW,{"schema":"rank-jump.equal-class-counts.v1","bindings":bindings(),"records":out})


def build(check=False):
    from sage.all import GF,PolynomialRing,EllipticCurve,prod
    data=r.read(RAW);assert data["bindings"]==bindings();rows=[];cases=r.read(INPUT)["cases"]
    for raw in data["records"]:
        case,p=raw["case"],raw["p"];row=cases[case]
        if raw["status"]!="COUNTED":
            rows.append({"case":case,"p":p,"status":"UNKNOWN"});continue
        assert previous.good(row,p,7)
        A,B=[r.mod(x,p) for x in row["model"][3:]];xs=[r.mod(P[0],p) for P in row["generic_points"]]
        fields=[]
        for saved in raw["fields"]:
            degree=saved["degree"];q=saved["q"];Fp=GF(p);R=PolynomialRing(Fp,"w")
            F=Fp if degree==1 else GF(q,name="w",modulus=R(saved["modulus_ascending"]))
            w=F(0) if degree==1 else F.gen()
            traces=struct.unpack("<"+"h"*q,base64.b64decode(saved["trace_values_i16_le_base64"]))
            total=0
            for i,value in enumerate(traces):
                u=sum(F((i//p**j)%p)*w**j for j in range(degree))
                d=prod(1-F(a)*u for a in xs)
                total+=(0 if d==0 else (1 if d.is_square() else -1))*value
            k=-prod(F(a) for a in xs)
            infinity=EllipticCurve(F,[0,2*F(A)*k,0,F(A*A)*k*k,-F(B*B)*k**3])
            tail=int(infinity.cardinality(algorithm="pari")-q-1)
            fields.append({"degree":degree,"finite_sum":total,"infinity_sum":tail,
                           "residual_H2_trace":total+tail})
        rec=previous.recover(p,[x["residual_H2_trace"] for x in fields])
        rows.append({"case":case,"p":p,"fields":fields,**rec})
    bounds=[]
    for case,row in enumerate(cases):
        relevant=[x for x in rows if x["case"]==case]
        rank18=[x for x in relevant if x["status"]=="RHO_18_REDUCTION"]
        witness=next(([a["p"],b["p"]] for i,a in enumerate(rank18) for b in rank18[i+1:]
                      if a["NS_discriminant_squareclass"]!=b["NS_discriminant_squareclass"]),None)
        upper=min([3]+[x["reduction_geometric_Picard_rank"]-17 for x in relevant if "reduction_geometric_Picard_rank" in x])
        if witness:upper=0
        bounds.append({"id":row["id"],"role":row["role"],"geometric_mixed_rank_interval":[0,int(upper)],
                       "arithmetic_mixed_rank_interval":[0,min(1,int(upper))],
                       "zero_rank_witness_primes":witness,
                       "full_base_arithmetic_rank_interval":[3,3+min(1,int(upper))],
                       "production_curve_rank":"UNKNOWN"})
    out={"schema":"rank-jump.equal-class-picard.v1","bindings":bindings(),"counts_sha256":r.digest(RAW.read_bytes()),
         "reductions":rows,"bounds":bounds,
         "boundary":"Ranks of the new function-field construction only; no rank upper bound on a production curve."}
    if check:assert out==r.read(OUTPUT);print("PASS equal-class Picard replay")
    else:r.write_new(OUTPUT,out)
    for row in rows:print(row["case"],row["p"],row.get("traces"),row["status"],row.get("NS_discriminant_squareclass"))
    for row in bounds:print(row)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["worker","capture","build","check"])
    p.add_argument("--case",type=int);p.add_argument("--prime",type=int);args=p.parse_args()
    if args.mode=="worker":worker(args.case,args.prime)
    elif args.mode=="capture":capture()
    else:build(args.mode=="check")

#!/usr/bin/env python3
"""Bounded reductions of the three residual branch-character K3 surfaces."""
import argparse
import base64
from pathlib import Path
import struct
import subprocess
import retrospective as r
import local_collision as lc

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"BRANCH_SPLIT_PICARD_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_branch_split_picard_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_branch_split_picard_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-branch-split-picard-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,lc.INPUT)}


def decode(row):
    raw=base64.b64decode(row["trace_values_i16_le_base64"],validate=True)
    assert len(raw)==2*row["q"]
    return list(struct.unpack("<"+"h"*row["q"],raw))


def character_sums(p,degree,roots):
    import numpy as np
    q=p**degree
    nonsquare=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
    legendre=np.array([0]+[1 if pow(a,(p-1)//2,p)==1 else -1 for a in range(1,p)],dtype=np.int64)
    xa=np.arange(q,dtype=np.int64)%p;xb=np.arange(q,dtype=np.int64)//p
    def mul(a,b,c,d):return ((a*c+nonsquare*b*d)%p,(a*d+b*c)%p)
    def chi(a,b):
        return legendre[a%p] if degree==1 else legendre[(a*a-nonsquare*b*b)%p]
    traces=[];weighted=[0,0,0]
    for u in range(q):
        ua,ub=u%p,u//p
        pa,pb=np.ones(q,dtype=np.int64),np.zeros(q,dtype=np.int64)
        for t in roots:
            ra=(t+ua*t*t)%p;rb=ub*t*t%p
            pa,pb=mul(pa,pb,(xa-ra)%p,(xb-rb)%p)
        value=int(chi(pa,pb).sum());traces.append(value)
        gamma=[((1-ua*t)%p,(-ub*t)%p) for t in roots]
        for k in range(3):
            i,j=[v for v in range(3) if v!=k]
            a,b=mul(*gamma[i],*gamma[j])
            weighted[k]+=int(chi(a,b))*value
        if u and u%512==0:print("checkpoint",p,degree,u,flush=True)
    packed=struct.pack("<"+"h"*q,*traces)
    return {"p":p,"degree":degree,"q":q,"quadratic_nonresidue":nonsquare,
            "field_encoding":"a+p*b; b=0 in degree one, otherwise omega^2=quadratic_nonresidue",
            "roots":roots,"trace_values_i16_le_base64":base64.b64encode(packed).decode(),
            "residual_H2_traces":weighted,
            "resolved_surface_point_counts":[q*q+19*q+1+t for t in weighted]}


def worker(p):
    anchor=r.read(lc.INPUT)["anchor"];A,B=map(int,anchor["short_model_ainvariants"][3:])
    assert B*(-4*A**3-27*B*B)%p
    roots=[x for x in range(p) if (x**3+A*x+B)%p==0];assert len(roots)==3
    rows=[]
    for degree in (1,2):
        path=WORK/f"p{p}-degree{degree}.json"
        row=character_sums(p,degree,roots)
        r.write_new(path,{"bindings":bindings(),**row});rows.append(row)
        print("completed",p,degree,row["residual_H2_traces"],flush=True)
    r.write_new(WORK/f"p{p}.json",{"bindings":bindings(),"p":p,"status":"COUNTED","fields":rows})


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for p in r.read(PROTOCOL)["primes"]:
        path=WORK/f"p{p}.json"
        if not path.exists():
            with (WORK/f"p{p}.log").open("x") as log:
                try:
                    proc=subprocess.run(["python3",str(Path(__file__).resolve()),"worker","--prime",str(p)],
                                        cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if proc.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:reason="30-second timeout"
                if reason:r.write_new(path,{"bindings":bindings(),"p":p,"status":"UNKNOWN","reason":reason,
                                           "transcript":(WORK/f"p{p}.log").read_text()})
        row=r.read(path);assert row["bindings"]==bindings();rows.append(row)
        print("checkpoint",p,row["status"],flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.branch-split-picard-inputs.v1","bindings":bindings(),"primes":rows})


def squarefree(n):
    assert n
    sign=-1 if n<0 else 1;n=abs(n);p=2
    while p*p<=n:
        while n%(p*p)==0:n//=p*p
        p+=1
    return sign*n


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings();rows=[]
    for raw in data["primes"]:
        if raw["status"]!="COUNTED":continue
        one,two=raw["fields"];p=raw["p"]
        for k,t1 in enumerate(one["residual_H2_traces"]):
            t2=two["residual_H2_traces"][k]
            options=[(e,t1-e*p) for e in (-1,1) if (t1-e*p)**2-p*p==t2]
            record={"p":p,"root":one["roots"][k],"t1":t1,"t2":t2,"status":"UNKNOWN"}
            if len(options)==1:
                e,t=options[0]
                assert abs(t)<=2*p
                record.update({"algebraic_eigenvalue":e*p,"quadratic_factor_ascending":[p*p,-t,1]})
                if t not in (-2*p,-p,0,p,2*p):
                    record.update({"status":"RHO_20_REDUCTION","geometric_picard_rank":20,
                                   "NS_discriminant_squareclass":squarefree(t*t-4*p*p)})
            rows.append(record)
    usable=[x for x in rows if x["status"]=="RHO_20_REDUCTION"]
    witness=next(([a,b] for i,a in enumerate(usable) for b in usable[i+1:]
                  if a["NS_discriminant_squareclass"]!=b["NS_discriminant_squareclass"]),None)
    out={"schema":"rank-jump.branch-split-picard.v1","bindings":bindings(),
         "input_sha256":r.digest(INPUT.read_bytes()),"reductions":rows,"two_place_witness":witness,
         "generic_geometric_Picard_rank":19 if witness else "UNKNOWN",
         "product_branch_twist_geometric_MW_rank":0 if witness else "UNKNOWN",
         "full_branch_cover_geometric_MW_rank":1 if witness else "UNKNOWN",
         "boundary":"Function-field statements only. No specialized curve rank or rational solubility follows from full branch splitting alone."}
    if check:assert r.read(OUTPUT)==out;print("PASS Picard-certificate replay")
    else:r.write_new(OUTPUT,out)
    for row in rows:print(row)
    print("generic Picard",out["generic_geometric_Picard_rank"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("worker","capture","build","check"))
    p.add_argument("--prime",type=int);args=p.parse_args()
    if args.mode=="worker":worker(args.prime)
    elif args.mode=="capture":capture()
    else:build(args.mode=="check")

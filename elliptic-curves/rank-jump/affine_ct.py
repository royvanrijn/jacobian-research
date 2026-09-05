#!/usr/bin/env python3
"""Four bounded Fisher pairings testing the new affine Selmer cosets."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import affine_selmer as af

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"AFFINE_CT_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_affine_ct_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-affine-ct-v1"
CAS=r.ROOT/"elliptic-curves/cas"


def bindings():
    paths=[PROTOCOL,af.INPUT,af.OUTPUT,lc.INPUT,Path(__file__),
           CAS/"research_runtime/sage_subspace.py",CAS/"research_runtime/sage_arithmetic.py",
           CAS/"research_runtime/fisher.py",CAS/"research_runtime/subspace.py"]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def setup(index):
    from sage.all import QQ,pari
    sys.path.insert(0,str(CAS))
    from research_runtime.store import FactStore
    from research_runtime.arithmetic import TwoTorsionContext
    from research_runtime.sage_arithmetic import SageArithmetic
    from research_runtime.subspace import GlobalSquareclasses
    from research_runtime.sage_subspace import SageSubspaceBackend
    protocol=r.read(PROTOCOL)
    pari.allocatemem(64000000,protocol["limits"]["pari_stack_bytes"],silent=True)
    pari.setrand(20260906)
    u,rad=protocol["pairs"][index]
    old=r.read(lc.INPUT)
    row=next(x for x in old["rows"] if int(x["parameter_u"])==u)
    matrix=next(x["matrix"] for x in old["ct"] if x["u"]==u)
    assert rad in [lc.lift(v,row["W_u_basis"]) for v in lc.orthogonal(map(r.pack,matrix),len(matrix))]
    loc=next(x for x in r.read(af.INPUT)["cases"] if x["u"]==u)
    affine=next(x for x in r.read(af.OUTPUT)["cases"] if x["u"]==u)
    assert affine["affine_solution"]["consistent"]
    mask=affine["affine_solution"]["particular_anchor_mask"]
    A,B=map(QQ,loc["anchor_model"][3:])
    coeff=[0,2*A*u,0,A+3*B*u+A*A*u*u,B+A*B*u*u-B*B*u**3]
    store=FactStore(WORK/f"facts-{index}")
    arithmetic=SageArithmetic(store)
    algebra=TwoTorsionContext((""+str(B),str(A),"0","1"))
    primes=[x["place"] for x in loc["local"] if x["place"]!="infinity"]
    context=arithmetic.prepare_congruent(list(map(str,coeff)),algebra,(0,1,u),
                                          factor_primes=primes,discover=True)
    arithmetic.field(algebra,factor_primes=[p["prime"] for p in old["anchor"]["base_discriminant_factorization"]],discover=True)
    backend=SageSubspaceBackend(arithmetic,context,None,
               local_candidate_cap=protocol["limits"]["local_witness_node_cap"],
               cover_policy="minimize-reduce")
    anchor=[backend.field_element(beta) for beta in loc["anchor_beta_coordinates"]]
    th=backend.theta
    # eta=kappa*gamma^2, so kappa is a smaller exact representative.
    first=1+u*th+u*u*(A+th*th)
    second=backend.K(1)
    for j,beta in enumerate(anchor):
        if mask>>j&1:first*=beta
        if rad>>j&1:second*=beta
    representatives=[backend.coordinates(first),backend.coordinates(second)]
    classes=GlobalSquareclasses(algebra.key,representatives,
              r.digest(af.OUTPUT.read_bytes())+":new-prime valuation plus inherited independence")
    return u,rad,mask,context,classes,backend,store


def compute(index,destination):
    from sage.all import QQ,pari
    u,rad,mask,context,classes,backend,store=setup(index)
    covers=[]
    for m in (1,2,3):
        print("COVER",u,rad,m,flush=True)
        cover=backend.cover(context,classes,m)
        backend.verify_cover(context,classes,m,cover)
        covers.append(cover)
    print("PAIR",u,rad,flush=True)
    result=backend._pair(classes,(1,2,3),covers)
    for term in result["local_terms"]:
        p=term["place"]
        if p!="infinity":
            assert int(pari.hilbert(QQ(covers[1]["quartic"][4]),QQ(term["gamma_value"]),p))==term["hilbert_symbol"]
    record={"schema":"rank-jump.affine-ct-pair.v1","bindings":bindings(),"pair_index":index,
            "u":u,"old_radical_mask":rad,"affine_anchor_correction":mask,
            "classes":list(classes.representatives),"context":context.record(),
            "covers":covers,"pair":result}
    # Portable proof objects contain the exact maps and Hilbert witnesses.
    r.write_new(destination,record)
    print("DONE",u,rad,result["value"],flush=True)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    results=[]
    for i,(u,rad) in enumerate(r.read(PROTOCOL)["pairs"]):
        path=WORK/f"pair-{i}.json"
        if not path.exists():
            try:
                process=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                     "--index",str(i),"--destination",str(path)],cwd=r.ROOT,
                     capture_output=True,text=True,timeout=r.read(PROTOCOL)["limits"]["seconds_per_pair"])
                log=process.stdout+process.stderr
            except subprocess.TimeoutExpired as exc:
                log=str(exc.stdout or "")+str(exc.stderr or "")+"\nTIMEOUT: UNKNOWN"
            with (WORK/f"pair-{i}.log").open("x") as f:f.write(log)
        if path.exists():
            result=r.read(path);assert result["bindings"]==bindings()
            results.append(result);print("checkpoint",u,rad,result["pair"]["value"],flush=True)
        else:
            results.append({"u":u,"old_radical_mask":rad,"status":"UNKNOWN",
                            "checkpoint_log":str((WORK/f"pair-{i}.log").relative_to(r.ROOT))})
            print("UNKNOWN",u,rad,flush=True)
    r.write_new(OUTPUT,{"schema":"rank-jump.affine-ct.v1","bindings":bindings(),"pairs":results})


def verify():
    from sage.all import QQ,pari
    data=r.read(OUTPUT);assert data["bindings"]==bindings()
    for record in data["pairs"]:
        if record.get("status")=="UNKNOWN":
            print("UNKNOWN",record["u"],record["old_radical_mask"]);continue
        u,rad,mask,context,classes,backend,_=setup(record["pair_index"])
        assert json.loads(json.dumps(context.record()))==record["context"]
        assert list(map(list,classes.representatives))==record["classes"]
        for m,c in zip((1,2,3),record["covers"]):
            backend.verify_cover(context,classes,m,c)
        result=backend._pair(classes,(1,2,3),record["covers"],retained=record["pair"])
        for t in result["local_terms"]:
            if t["place"]!="infinity":
                assert int(pari.hilbert(QQ(record["covers"][1]["quartic"][4]),QQ(t["gamma_value"]),t["place"]))==t["hilbert_symbol"]
        print("PASS exact cover and Fisher replay",u,rad,result["value"])


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","verify"))
    parser.add_argument("--index",type=int);parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="worker":compute(args.index,args.destination)
    elif args.mode=="capture":capture()
    else:verify()

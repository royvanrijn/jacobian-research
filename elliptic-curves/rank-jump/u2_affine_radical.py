#!/usr/bin/env python3
"""Complete one CT column and expose exactly two affine radical covers."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import affine_selmer as af
import affine_ct as old
from ct_variation import solve

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"U2_AFFINE_RADICAL_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_u2_affine_ct_inputs_v1.json"
REPORT=r.OUT/"rank_jump_u2_affine_radical_v1.json"
COVERS=r.OUT/"rank_jump_u2_affine_covers_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-u2-affine-radical-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (PROTOCOL,Path(__file__),HERE/"affine_ct.py",lc.INPUT,af.INPUT,af.OUTPUT,old.OUTPUT)}


def basis():
    inp=r.read(lc.INPUT)
    w=next(x["W_u_basis"] for x in inp["rows"] if int(x["parameter_u"])==2)
    B=next(x["matrix"] for x in inp["ct"] if x["u"]==2)
    return w,B


def setup(kind,index):
    from sage.all import QQ
    sys.path.insert(0,str(old.CAS))
    from research_runtime.subspace import GlobalSquareclasses
    # Process-local destination override; source file and prior store unchanged.
    old.WORK=WORK/"arithmetic"
    _,_,_,context,original,backend,store=old.setup(2)
    loc=next(x for x in r.read(af.INPUT)["cases"] if x["u"]==2)
    betas=[backend.field_element(b) for b in loc["anchor_beta_coordinates"]]
    if kind=="pair":
        j=r.read(PROTOCOL)["old_basis_indices"][index]
        mask=basis()[0][j]
        first=backend.field_element(original.representatives[0])
        second=backend.K(1)
        for i,b in enumerate(betas):
            if mask>>i&1:second*=b
    else:
        report=r.read(REPORT)
        mask=report["affine_anchor_masks"][0]
        A=QQ(loc["anchor_model"][3]);th=backend.theta
        first=1+2*th+4*(A+th*th)
        second=backend.K(1)
        for i,b in enumerate(betas):
            if mask>>i&1:first*=b
            if report["old_radical_mask"]>>i&1:second*=b
    classes=GlobalSquareclasses(original.algebra_key,
                    [backend.coordinates(first),backend.coordinates(second)],
                    r.digest(af.OUTPUT.read_bytes())+":independent affine class and nonzero inherited class")
    return context,classes,backend


def worker_pair(index,destination):
    from sage.all import QQ,pari
    context,classes,backend=setup("pair",index)
    covers=[]
    for mask in (1,2,3):
        print("COVER",index,mask,flush=True)
        cover=backend.cover(context,classes,mask)
        backend.verify_cover(context,classes,mask,cover)
        covers.append(cover)
    result=backend._pair(classes,(1,2,3),covers)
    for t in result["local_terms"]:
        if t["place"]!="infinity":
            assert int(pari.hilbert(QQ(covers[1]["quartic"][4]),QQ(t["gamma_value"]),t["place"]))==t["hilbert_symbol"]
    r.write_new(destination,{"index":index,"bindings":bindings(),
                 "old_basis_index":r.read(PROTOCOL)["old_basis_indices"][index],
                 "classes":list(classes.representatives),"context":context.record(),
                 "covers":covers,"pair":result})
    print("DONE",index,result["value"],flush=True)


def run_worker(mode,index,dest,seconds):
    if dest.exists():
        record=r.read(dest);assert record["bindings"]==bindings()
        return record
    result=subprocess.run(["sage","-python",str(Path(__file__).resolve()),mode,
             "--index",str(index),"--destination",str(dest)],cwd=r.ROOT,
             capture_output=True,text=True,timeout=seconds)
    with dest.with_suffix(".log").open("x") as f:f.write(result.stdout+result.stderr)
    if result.returncode or not dest.exists():
        raise RuntimeError("worker incomplete; see "+str(dest.with_suffix(".log")))
    return r.read(dest)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    rows=[]
    for i in range(12):
        row=run_worker("pair",i,WORK/f"pair-{i}.json",60)
        rows.append(row);print("checkpoint pair",i,row["pair"]["value"],flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.u2-affine-ct-inputs.v1","bindings":bindings(),"pairs":rows})


def build(check=False):
    inp=r.read(INPUT);assert inp["bindings"]==bindings()
    protocol=r.read(PROTOCOL);w,B=basis();n=len(w)
    rad=lc.orthogonal(map(r.pack,B),n)
    assert len(rad)==1 and lc.lift(rad[0],w)==protocol["existing_old_radical_mask"]
    omit=protocol["omitted_basis_index"]
    assert omit==rad[0].bit_length()-1
    oldpair=next(p for p in r.read(old.OUTPUT)["pairs"] if p["u"]==2)
    assert oldpair["old_radical_mask"]==protocol["existing_old_radical_mask"]
    assert oldpair["pair"]["value"]==0
    values=[None]*n
    for rec in inp["pairs"]:
        values[rec["old_basis_index"]]=rec["pair"]["value"]
    values[omit]=sum(values[i] for i in range(n) if i!=omit and rad[0]>>i&1)%2
    column=r.pack(values)
    correction=solve(list(map(r.pack,B)),column,n)
    anchor_mask=next(x for x in r.read(af.OUTPUT)["cases"] if x["u"]==2)["affine_solution"]["particular_anchor_mask"]
    first=anchor_mask^lc.lift(correction,w)
    second=first^protocol["existing_old_radical_mask"]
    extended=[r.pack(row)|(values[i]<<n) for i,row in enumerate(B)]+[column]
    radical=lc.orthogonal(extended,n+1)
    global_basis=w+[anchor_mask|(1<<20)]
    embedded=lc.canonical([lc.lift(v,global_basis) for v in radical])
    expected=lc.canonical([protocol["existing_old_radical_mask"],first|(1<<20)])
    assert embedded==expected and r.rank(extended)==12
    out={"schema":"rank-jump.u2-affine-radical.v1","bindings":bindings(),
         "pair_input_sha256":r.digest(INPUT.read_bytes()),"u":2,
         "old_basis":w,"new_CT_column":values,"omitted_column_index":omit,
         "old_radical_coordinate_mask":rad[0],"old_radical_mask":protocol["existing_old_radical_mask"],
         "CT_correction_old_basis_mask":correction,"affine_anchor_masks":[first,second],
         "extended_CT_matrix_packed_rows":extended,"embedded_radical_basis":embedded,
         "extended_CT_rank":12,"restricted_radical_dimension":2,
         "scope":"The two affine masks specify eta times their anchor products; solubility remains UNKNOWN."}
    if check:
        assert r.read(REPORT)==out
        print("PASS u2 affine radical linear replay")
    else:r.write_new(REPORT,out)
    print("column",values,"affine masks",first,second,flush=True)


def worker_cover(index,destination):
    context,classes,backend=setup("cover",index)
    mask=(1,3)[index]
    cover=backend.cover(context,classes,mask)
    backend.verify_cover(context,classes,mask,cover)
    r.write_new(destination,{"index":index,"bindings":bindings(),
          "global_affine_anchor_mask":r.read(REPORT)["affine_anchor_masks"][index],
          "class_mask":mask,"classes":list(classes.representatives),
          "context":context.record(),"cover":cover})


def worker_points(index,destination):
    from sage.all import QQ,pari
    context,classes,backend=setup("cover",index)
    record=r.read(WORK/f"cover-{index}.json")
    cover=record["cover"];mask=record["class_mask"]
    backend.verify_cover(context,classes,mask,cover)
    q=backend.R(list(map(QQ,cover["quartic"])));den=q.denominator()
    height=r.read(PROTOCOL)["limits"]["point_search_height"]
    raw=pari.hyperellratpoints(q*den**2,height)
    points=[(QQ(P[0]),QQ(1),QQ(P[1])/den) for P in raw]
    if q[4].is_square():
        points.append((QQ(1),QQ(0),q[4].sqrt()))
    recovered=[]
    for P in points:
        image=backend.point_from_cover(context,classes,mask,cover,P)
        recovered.append({"quartic_point":list(map(str,P)),"point_on_input_curve":image})
    r.write_new(destination,{"index":index,"bindings":bindings(),"height":height,
              "cover_sha256":r.digest((WORK/f"cover-{index}.json").read_bytes()),
              "point_count":len(points),"points":recovered,
              "status":"RATIONAL_POINTS_CERTIFIED" if points else "UNKNOWN_BOUNDED_MISS"})


def capture_covers():
    rows=[]
    for i in range(2):
        cover=run_worker("cover",i,WORK/f"cover-{i}.json",60)
        print("checkpoint cover",i,flush=True)
        points=run_worker("points",i,WORK/f"points-{i}.json",10)
        rows.append({"cover":cover,"search":points})
        print("checkpoint solubility",i,points["status"],flush=True)
    r.write_new(COVERS,{"schema":"rank-jump.u2-affine-covers.v1","bindings":bindings(),
                        "radical_report_sha256":r.digest(REPORT.read_bytes()),"rows":rows})


def verify():
    from sage.all import QQ,pari
    inp=r.read(INPUT);assert inp["bindings"]==bindings()
    for rec in inp["pairs"]:
        context,classes,backend=setup("pair",rec["index"])
        assert list(map(list,classes.representatives))==rec["classes"]
        assert json.loads(json.dumps(context.record()))==rec["context"]
        for mask,cover in zip((1,2,3),rec["covers"]):
            backend.verify_cover(context,classes,mask,cover)
        result=backend._pair(classes,(1,2,3),rec["covers"],retained=rec["pair"])
        for t in result["local_terms"]:
            if t["place"]!="infinity":
                assert int(pari.hilbert(QQ(rec["covers"][1]["quartic"][4]),QQ(t["gamma_value"]),t["place"]))==t["hilbert_symbol"]
        print("PASS pair",rec["index"],result["value"],flush=True)
    if COVERS.exists():
        data=r.read(COVERS);assert data["bindings"]==bindings()
        for row in data["rows"]:
            rec=row["cover"];context,classes,backend=setup("cover",rec["index"])
            assert list(map(list,classes.representatives))==rec["classes"]
            backend.verify_cover(context,classes,rec["class_mask"],rec["cover"])
            for P in row["search"]["points"]:
                assert backend.point_from_cover(context,classes,rec["class_mask"],rec["cover"],
                        P["quartic_point"])==P["point_on_input_curve"]
            print("PASS cover",rec["index"],row["search"]["status"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","pair","build","check","cover","points","capture-covers","verify"))
    parser.add_argument("--index",type=int);parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="pair":worker_pair(args.index,args.destination)
    elif args.mode=="cover":worker_cover(args.index,args.destination)
    elif args.mode=="points":worker_points(args.index,args.destination)
    elif args.mode=="capture":capture()
    elif args.mode=="capture-covers":capture_covers()
    elif args.mode=="verify":verify()
    else:build(args.mode=="check")

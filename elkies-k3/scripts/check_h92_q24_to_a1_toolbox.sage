#!/usr/bin/env sage -python
"""
Check/export the selected q24/orbit85 D12 -> A1 corridor while applying the
reusable q32-era diagnostics to the q24 marking.

No q32 frame is imported.
"""
import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, lcm, matrix, pari, vector

def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [h/"Documents"/"jacobian-research",h/"jacobian-research"]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
WORK=LOCAL/"q24-to-a1-toolbox"
FRAMES=WORK/"frames"
WORK.mkdir(parents=True,exist_ok=True)
FRAMES.mkdir(parents=True,exist_ok=True)

SIG=LOCAL/f"q24-orbit85-d12-signature-mod-{args.prime}.json"
if not SIG.exists():
    raise SystemExit(
        f"missing {SIG}; run extract_h92_q24_d12_modp_signature.sage first"
    )
sig=json.loads(SIG.read_text())
assert sig["status"]=="PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE"
assert sig["rr"]=={
    "ambient":56,
    "collision_rank":48,
    "post_collision":8,
    "resolved_rank":6,
    "kernel":2,
    "geometric_fibre_twist":-8,
}
assert (
    int(sig["quartic_degree"]),
    int(sig["child_root_rank"]),
    int(sig["child_root_det"]),
    int(sig["child_euler"]),
)==(4,12,4,24)

print(
    "Q24A1CHECK_ANCHOR|"
    "q24=orbit85|child=D12/MW5|RR=56,48,8,6,2|"
    "quartic=4|root=12,4|euler=24|status=PASS",
    flush=True,
)

VERIFY=SCRIPTS/"verify_h3_d13_to_mw17_path.sage"
REPLAY=WORK/"q24-historical-lattice-replay.json"
if not VERIFY.exists():
    raise SystemExit(f"missing {VERIFY}")

saved=list(sys.argv)
scope={"__name__":"__main__","__file__":str(VERIFY)}
buf=io.StringIO()
try:
    sys.argv=[str(VERIFY),"--output",str(REPLAY)]
    with contextlib.redirect_stdout(buf):
        exec(compile(VERIFY.read_text(),str(VERIFY),"exec"),scope)
finally:
    sys.argv=saved

for line in buf.getvalue().splitlines():
    if line.startswith("H3D13MW17|"):
        print(line,flush=True)

path=json.loads(REPLAY.read_text())
assert path["status"]=="PASS_H3_D13_TO_MW17_LATTICE_PATH"
assert len(path["steps"])==11

expected=[
    (1,24,85,"D12",5,(12,264,4)),
    (2,6,42,"A11",6,(11,132,12)),
    (3,8,922,"A5+A5",7,(10,60,36)),
    (4,4,472,"A3+A3+A3",8,(9,36,64)),
    (5,4,323,"A2+A2+A3",10,(7,24,36)),
    (6,4,207,"A1+A1+A1+A1+A1",12,(5,10,32)),
    (7,4,52,"A1+A1+A1+A1",13,(4,8,16)),
    (8,4,114,"A1+A1+A1",14,(3,6,8)),
    (9,4,498,"A1+A1",15,(2,4,4)),
    (10,4,981,"A1",16,(1,2,2)),
]

def save_frame(path,frame,meta):
    lines=[f"# {k} = {v}" for k,v in meta.items()]
    lines += [" ".join(map(str,row)) for row in frame.rows()]
    path.write_text("\n".join(lines)+"\n")

records=[]
current_parent=None
d12_frame=None
d12_target=None

for idx,q,orbit,ade,mw,root_data in expected:
    step=path["steps"][idx-1]
    assert (
        int(step["step"]),int(step["q"]),int(step["orbit_index"]),
        step["ade"],int(step["mw_rank"]),tuple(step["root_data"])
    )==(idx,q,orbit,ade,mw,root_data)

    artifact=ROOT/step["artifact"]
    if not artifact.exists():
        raise SystemExit(f"missing selected neighbor artifact {artifact}")
    data=json.loads(artifact.read_text())
    assert data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    hits=[r for r in data["neighbors"] if int(r["orbit_index"])==orbit]
    assert len(hits)==1
    rec=hits[0]
    assert int(rec["q"])==q
    assert int(rec["old_fiber_degree"])==2
    assert rec["child_ade"]==ade
    assert int(rec["child_mw_rank"])==mw
    assert tuple(rec["child_root_data"])==root_data

    child=matrix(ZZ,rec["child_root_adapted_frame"])
    assert child.dimensions()==(17,17)
    assert child.det()==948

    frame_path=FRAMES/f"step{idx:02d}-{ade.replace('+','_')}-mw{mw}.txt"
    save_frame(frame_path,child,{
        "q":q,
        "orbit":orbit,
        "ADE":ade,
        "MW":mw,
        "source_artifact":step["artifact"],
    })

    if current_parent is not None:
        rr=current_parent["root_rank"]
        parent=current_parent["frame"]
        cartan=parent[:rr,:rr]
        coupling=parent[:rr,rr:]
        tail=parent[rr:,rr:]
        height=tail-coupling.transpose()*cartan.inverse()*coupling

        witness=vector(ZZ,rec["witness"])
        mwv=vector(ZZ,rec["mw_projection"])
        labels=vector(ZZ,rec["dominant_labels"])
        rootcoords=vector(ZZ,witness[:rr])
        assert vector(ZZ,witness[rr:])==mwv
        assert cartan*rootcoords+coupling*mwv==labels

        root_norm=QQ(labels*cartan.inverse()*labels)
        mw_norm=QQ(mwv*height*mwv)
        assert root_norm+mw_norm==QQ(2*q)

        entry={
            "step":idx,
            "source":current_parent["ade"],
            "child":ade,
            "q":q,
            "orbit":orbit,
            "root_rank":rr,
            "mw_rank_parent":17-rr,
            "mw_projection":list(map(int,mwv)),
            "dominant_labels":list(map(int,labels)),
            "root_norm":str(root_norm),
            "mw_norm":str(mw_norm),
            "total_norm":str(root_norm+mw_norm),
            "witness":list(map(int,witness)),
        }
        records.append(entry)

        print(
            "Q24A1CHECK_TARGET|"
            f"step={idx}|source={current_parent['ade']}|child={ade}|"
            f"q={q}|orbit={orbit}|root_norm={root_norm}|mw_norm={mw_norm}|"
            f"labels_L1={sum(abs(int(v)) for v in labels)}|"
            "status=PASS_SELECTED_WITNESS",
            flush=True,
        )

        if idx==2:
            d12_frame=parent
            d12_target=(mwv,rec,entry)

    current_parent={
        "frame":child,
        "root_rank":root_data[0],
        "ade":ade,
        "mw":mw,
    }

assert current_parent["ade"]=="A1"
assert current_parent["root_rank"]==1
print(
    "Q24A1CHECK_LATTICE|steps=10|final=A1/MW16|"
    "all_old_fibre_degree=2|status=PASS_Q24_TO_A1_SELECTED_SUFFIX",
    flush=True,
)

assert d12_frame is not None and d12_target is not None
G=d12_frame
root=G[:12,:12]
coupling=G[:12,12:]
tail=G[12:,12:]
H=tail-coupling.transpose()*root.inverse()*coupling
assert H.dimensions()==(5,5)

def class_order(dual):
    order=ZZ(1)
    for value in dual:
        order=lcm(order,ZZ(QQ(value).denominator()))
    return order

def mod2(x):
    x=QQ(x)
    return x-2*(x/2).floor()

def correction_for(z):
    z=vector(ZZ,z)
    base=vector(ZZ,[0]*12+list(z))
    pair=vector(QQ,base*G[:,:12])
    dual=pair*root.inverse()
    order=class_order(dual)
    raw=QQ(dual*root*dual)
    if order==1:
        return QQ(0),order
    candidates=[QQ(0),QQ(1),QQ(3)]
    candidates=[c for c in candidates if mod2(c)==mod2(raw)]
    h=QQ(z*H*z)
    valid=[]
    for c in candidates:
        po=(h+c-4)/2
        if po in ZZ and po>=0:
            valid.append((ZZ(po),c))
    if not valid:
        return None,order
    valid.sort()
    return valid[0][1],order

scale=ZZ(1)
for v in H.list():
    scale=lcm(scale,ZZ(QQ(v).denominator()))
IH=(scale*H).change_ring(ZZ)
qf=pari(IH).qfminim(ZZ(scale*12))
cols=matrix(ZZ,qf[2]).columns()

rows=[]
seen=set()
for col in cols:
    for sign in (1,-1):
        z=sign*vector(ZZ,col)
        key=tuple(map(int,z))
        if key in seen:
            continue
        seen.add(key)
        h=QQ(z*H*z)
        if h>12:
            continue
        corr,order=correction_for(z)
        if corr is None:
            continue
        po=(h+corr-4)/2
        if po in ZZ and po>=0:
            rows.append((ZZ(po),h,corr,order,z))

rows.sort(key=lambda t:(t[0],t[1],sum(abs(int(x)) for x in t[4]),tuple(t[4])))
zero=[z for po,h,c,o,z in rows if po==0]
zero_span_rank=matrix(ZZ,[list(z) for z in zero]).rank() if zero else 0
L=matrix(ZZ,[list(z) for z in zero]).row_module() if zero else None

target_mw=d12_target[0]
target_corr,target_order=correction_for(target_mw)
target_h=QQ(target_mw*H*target_mw)
target_po=None if target_corr is None else (target_h+target_corr-4)/2
target_in_zero_span=bool(L is not None and target_mw in L)

counts={}
for po,h,c,o,z in rows:
    counts[int(po)]=counts.get(int(po),0)+1

print(
    "Q24A1CHECK_D12_SECTION_PROFILE|"
    f"zero_pole_count={len(zero)}|zero_pole_span_rank={zero_span_rank}|"
    f"target_mw={','.join(map(str,target_mw))}|target_height={target_h}|"
    f"target_correction={target_corr}|target_PdotO={target_po}|"
    f"target_in_zero_pole_Zspan={int(target_in_zero_span)}|"
    "status=PASS_Q32_STYLE_D12_PROFILE",
    flush=True,
)

if target_in_zero_span:
    next_method="DIRECT_ZERO_POLE_SECTION_SPAN"
elif target_po is not None and target_po<=2:
    next_method="LOW_POLE_POLYNOMIAL_SECTION_ENUMERATION_THEN_ACTUAL_TWIST"
else:
    next_method="GENERIC_DIVISORIAL_RR_WITH_TARGETED_Q6_WITNESS"

print(
    "Q24A1CHECK_FRONTIER|"
    "source=D12/MW5|target=A11/MW6|q=6|orbit=42|"
    f"method={next_method}|status=ACTIONABLE",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-to-a1-toolbox.v1",
    "status":"PASS_Q24_TO_A1_TOOLBOX_CHECK",
    "prime":int(args.prime),
    "anchor_signature":str(SIG.relative_to(ROOT)),
    "lattice_replay":str(REPLAY.relative_to(ROOT)),
    "stop":"A1/MW16",
    "selected_suffix":[
        {
            "step":idx,
            "q":q,
            "orbit":orbit,
            "ade":ade,
            "mw_rank":mw,
            "root_data":list(root_data),
        }
        for idx,q,orbit,ade,mw,root_data in expected
    ],
    "target_decompositions":records,
    "d12_to_a11_frontier":{
        "q":6,
        "orbit":42,
        "mw_height_gram":[[str(v) for v in row] for row in H.rows()],
        "short_section_counts_by_PdotO":counts,
        "zero_pole_count":len(zero),
        "zero_pole_span_rank":int(zero_span_rank),
        "target_mw":list(map(int,target_mw)),
        "target_height":str(target_h),
        "target_correction":None if target_corr is None else str(target_corr),
        "target_P_dot_O":None if target_po is None else str(target_po),
        "target_in_zero_pole_integral_span":target_in_zero_span,
        "recommended_next_method":next_method,
    },
    "exported_frames":str(FRAMES.relative_to(ROOT)),
    "proof_boundary":(
        "The q24 D12 anchor is explicit over GF(p); D12->A1 is an exact "
        "selected integral/nef lattice corridor. This checker does not claim "
        "equation-level execution of the nine downstream pencils."
    ),
}

OUT=(args.output.resolve() if args.output else WORK/f"q24-to-a1-toolbox-p{args.prime}.json")
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24A1CHECK_RESULT|anchor_modp=PASS|lattice_to_A1=PASS|"
    f"frontier_method={next_method}|status=PASS_Q24_TO_A1_TOOLBOX_CHECK",
    flush=True,
)

#!/usr/bin/env python3
"""
Walk the q24 equation-side D12 frame natively to A1 using targeted exact
degree-two Weyl-neighbor searches.

At each stage the next ADE/root signature is known in advance.  We do not
reuse q32 frames or historical orbit numbers: the search is rerun in the
actual q24-derived child frame and accepts only the requested root data.
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT=Path("/Users/royvanrijn/Documents/jacobian-research")
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3/q24-native-suffix"
SEARCH=SCRIPTS/"search_root_adapted_weyl_neighbors_targeted.sage"
EXPORT=SCRIPTS/"export_h92_q24_native_d12_frame.sage"

STEPS=[
    # name, q, parent root rank, target root data, ADE, MW
    ("a11",      6,12,(11,132,12),"A11",6),
    ("2a5",      8,11,(10,60,36),"A5+A5",7),
    ("3a3",      4,10,(9,36,64),"A3+A3+A3",8),
    ("a3-2a2",   4, 9,(7,24,36),"A2+A2+A3",10),
    ("5a1",      4, 7,(5,10,32),"A1+A1+A1+A1+A1",12),
    ("4a1",      4, 5,(4,8,16),"A1+A1+A1+A1",13),
    ("3a1",      4, 4,(3,6,8),"A1+A1+A1",14),
    ("2a1",      4, 3,(2,4,4),"A1+A1",15),
    ("a1",       4, 2,(1,2,2),"A1",16),
]

def run(cmd):
    print("+"," ".join(map(str,cmd)),flush=True)
    subprocess.run(cmd,cwd=str(ROOT),check=True)

def write_frame(path,record):
    path.write_text(
        f"# q24-native child {record['child_ade']}/MW{record['child_mw_rank']}\n"
        f"# q = {record['q']}\n"
        f"# discovered orbit index = {record['orbit_index']}\n"
        + "\n".join(
            " ".join(map(str,row))
            for row in record["child_root_adapted_frame"]
        )
        + "\n"
    )

if not SEARCH.exists():
    raise SystemExit(f"missing targeted search engine {SEARCH}")

LOCAL.mkdir(parents=True,exist_ok=True)

# Always regenerate the anchor from D24eq so stale q32/historical frames cannot
# silently enter the run.
run(["sage","-python",str(EXPORT)])
current=LOCAL/"step00-d12-mw5-frame.txt"
if not current.exists():
    raise SystemExit("q24 native D12 frame was not produced")

ledger=[]
for idx,(slug,q,parent_rr,target,ade,mw) in enumerate(STEPS,1):
    out=LOCAL/f"step{idx:02d}-{slug}-search.json"
    frames=LOCAL/f"step{idx:02d}-{slug}-frames"

    cmd=[
        "sage","-python",str(SEARCH),
        "--frame",str(current),
        "--root-rank",str(parent_rr),
        "--q",str(q),
        "--degree","2",
        "--adapt-mw-at-least",str(mw),
        "--rank-growth-only",
        "--stop-after-first-growth",
        "--stream-first-growth",
        "--target-root-rank",str(target[0]),
        "--target-root-count",str(target[1]),
        "--target-root-determinant",str(target[2]),
        "--output",str(out),
        "--frames-dir",str(frames),
        "--stream-progress-every","250",
    ]
    run(cmd)

    data=json.loads(out.read_text())
    if data.get("status")!="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS":
        raise SystemExit(f"stage {slug}: unexpected search status {data.get('status')}")

    hits=[
        r for r in data.get("neighbors",[])
        if tuple(r.get("child_root_data",()))==target
        and r.get("child_ade")==ade
        and int(r.get("child_mw_rank",-1))==mw
        and int(r.get("q",-1))==q
    ]
    if not hits:
        raise SystemExit(
            f"stage {slug}: targeted search produced no "
            f"{ade}/MW{mw} child with root_data={target}"
        )

    # Streaming target mode should normally produce one.  If there are several,
    # use the deterministic first record and retain all hits in the ledger.
    rec=hits[0]
    child=LOCAL/f"step{idx:02d}-{slug}-mw{mw}-frame.txt"
    write_frame(child,rec)

    ledger.append({
        "step":idx,
        "source_frame":str(current.relative_to(ROOT)),
        "q":q,
        "old_fibre_degree":2,
        "discovered_orbit_index":int(rec["orbit_index"]),
        "target_ade":ade,
        "target_root_data":list(target),
        "target_mw_rank":mw,
        "mw_projection":rec["mw_projection"],
        "dominant_labels":rec["dominant_labels"],
        "witness":rec["witness"],
        "child_frame":str(child.relative_to(ROOT)),
        "search_artifact":str(out.relative_to(ROOT)),
        "matching_hits":len(hits),
    })

    print(
        "Q24NATIVE_SUFFIX|"
        f"step={idx}|q={q}|orbit={rec['orbit_index']}|"
        f"ADE={ade}|root_data={target[0]},{target[1]},{target[2]}|"
        f"MW={mw}|status=PASS",
        flush=True,
    )
    current=child

result={
    "schema":"elkies-k3.h3-q24-native-d12-to-a1.v1",
    "status":"PASS_Q24_NATIVE_D12_TO_A1",
    "anchor":"artifacts/local/elkies-k3/q24-native-suffix/step00-d12-mw5-frame.json",
    "route":"D12-q6-A11-q8-2A5-q4-3A3-q4-A3+2A2-q4-5A1-q4-4A1-q4-3A1-q4-2A1-q4-A1",
    "all_old_fibre_degree":2,
    "steps":ledger,
    "final_frame":str(current.relative_to(ROOT)),
    "final_ade":"A1",
    "final_mw_rank":16,
    "proof_boundary":(
        "This is a fresh exact lattice search starting from the equation-side "
        "q24 D24eq child. It proves the native degree-two corridor to A1. "
        "It does not yet execute the downstream equation-level pencils."
    ),
}
outfile=LOCAL/"q24-native-d12-to-a1.json"
outfile.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")

print(f"OUTPUT|{outfile}",flush=True)
print(
    "Q24NATIVE_SUFFIX_RESULT|steps=9|final=A1/MW16|"
    "status=PASS_Q24_NATIVE_D12_TO_A1",
    flush=True,
)

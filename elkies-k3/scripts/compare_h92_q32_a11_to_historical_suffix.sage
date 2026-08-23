
#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import QQ, ZZ, lcm, matrix, pari

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

NEW=LOCAL/"route-scout/q32-d12-next-q4q6q8.json"
HIST=GEN/"elkies-k3-h3-d12-o85-q6-degree2.json"
OUT=LOCAL/"route-scout/q32-a11-historical-splice.json"

for p in (NEW,HIST):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")

new=json.loads(NEW.read_text())
hist=json.loads(HIST.read_text())
assert new["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert hist["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

nrec=next(
    r for r in new["neighbors"]
    if int(r["q"])==6 and int(r["orbit_index"])==12
    and r["child_ade"]=="A11" and int(r["child_mw_rank"])==6
)
hrec=next(
    r for r in hist["neighbors"]
    if int(r["orbit_index"])==42
)

N=matrix(ZZ,nrec["child_root_adapted_frame"])
H=matrix(ZZ,hrec["child_root_adapted_frame"])
assert N.dimensions()==H.dimensions()==(17,17)
assert N.det()==H.det()==948
assert tuple(nrec["child_root_data"])==(11,132,12)
assert tuple(hrec["child_root_data"])==(11,132,12)

NH=matrix(QQ,nrec["child_mw_height"])
HH=matrix(QQ,hrec["child_mw_height"])
assert NH.dimensions()==HH.dimensions()==(6,6)

exact_frame=(N==H)
exact_height=(NH==HH)

# Compare rational MW lattices by clearing denominators.
scale=ZZ(1)
for v in list(NH.list())+list(HH.list()):
    scale=lcm(scale,ZZ(QQ(v).denominator()))
Ni=(scale*NH).change_ring(ZZ)
Hi=(scale*HH).change_ring(ZZ)

mw_iso=None
mw_iso_matrix=None
try:
    qiso=pari(Ni).qfisom(pari(Hi))
    if str(qiso)!="0":
        mw_iso=True
        mw_iso_matrix=matrix(ZZ,qiso)
        assert mw_iso_matrix*Ni*mw_iso_matrix.transpose()==Hi
    else:
        mw_iso=False
except Exception as exc:
    mw_iso=f"ERROR:{type(exc).__name__}:{exc}"

# Full-frame isometry.  Often deterministic adaptation makes this exact;
# otherwise PARI can still prove the integral lattices are isometric.
full_iso=None
full_iso_matrix=None
if exact_frame:
    full_iso=True
    full_iso_matrix=matrix.identity(ZZ,17)
else:
    try:
        qiso=pari(N).qfisom(pari(H))
        if str(qiso)!="0":
            full_iso=True
            full_iso_matrix=matrix(ZZ,qiso)
            assert full_iso_matrix*N*full_iso_matrix.transpose()==H
        else:
            full_iso=False
    except Exception as exc:
        full_iso=f"ERROR:{type(exc).__name__}:{exc}"

# Compare the q8 shells directly as an additional splice test.  If the frames
# are exactly equal, the historical orbit922 witness should have norm 16 in N
# and can be replayed verbatim.
A11Q8=GEN/"elkies-k3-h3-a11-middle-q8-degree2.json"
historical_q8_replay=False
historical_q8_orbit=None
if A11Q8.exists():
    q8=json.loads(A11Q8.read_text())
    h8=next(r for r in q8["neighbors"] if int(r["orbit_index"])==922)
    historical_q8_orbit=h8
    w=matrix(ZZ,1,17,h8["witness"]).row(0)
    historical_q8_replay=(w*N*w==16)

print(
    "Q32A11SPLICE|"
    f"new_orbit=12|historical_orbit=42|"
    f"frame_equal={int(exact_frame)}|height_equal={int(exact_height)}|"
    f"mw_height_isometric={mw_iso}|full_frame_isometric={full_iso}|"
    f"historical_q8_witness_replays={int(historical_q8_replay)}|"
    "status=PASS_COMPARISON",
    flush=True,
)

if full_iso is True:
    status=(
        "PASS_LITERAL_SUFFIX_SPLICE"
        if exact_frame and historical_q8_replay
        else "PASS_SUFFIX_SPLICE_BY_ISOMETRY"
    )
else:
    status="DISTINCT_A11_FRAME_REQUIRES_FRESH_SUFFIX"

print(
    "Q32A11SPLICE_RESULT|"
    f"status={status}",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q32-a11-historical-splice.v1",
    "status":status,
    "new":{
        "q":6,"orbit":12,
        "mw_projection":nrec["mw_projection"],
        "dominant_labels":nrec["dominant_labels"],
        "witness":nrec["witness"],
        "frame":nrec["child_root_adapted_frame"],
        "mw_height":nrec["child_mw_height"],
    },
    "historical":{
        "q":6,"orbit":42,
        "mw_projection":hrec["mw_projection"],
        "dominant_labels":hrec["dominant_labels"],
        "witness":hrec["witness"],
        "frame":hrec["child_root_adapted_frame"],
        "mw_height":hrec["child_mw_height"],
    },
    "comparison":{
        "frame_equal":exact_frame,
        "height_equal":exact_height,
        "mw_height_isometric":mw_iso,
        "full_frame_isometric":full_iso,
        "mw_isometry":None if mw_iso_matrix is None else [
            list(map(int,row)) for row in mw_iso_matrix.rows()
        ],
        "full_frame_isometry":None if full_iso_matrix is None else [
            list(map(int,row)) for row in full_iso_matrix.rows()
        ],
        "historical_q8_orbit922_witness_replays":historical_q8_replay,
    },
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)

#!/usr/bin/env sage -python
"""
Export the actual q24/orbit85 D12 child frame directly from the equation-side
q24 divisor D24eq.

This uses the exact NS class produced by close_h92_q8_q24_by_q6_translation
and the generic exact neighbor engine.  No q32 D12 marking and no historical
orbit85 child frame are imported into the construction.
"""
import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, lcm, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-native-suffix"
OUTFRAME=OUTDIR/"step00-d12-mw5-frame.txt"
OUTMETA=OUTDIR/"step00-d12-mw5-frame.json"

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"

for path in (ENGINE,CLOSE):
    if not path.exists():
        raise SystemExit(f"missing {path}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

saved=list(sys.argv)
scope={"__name__":"__embedded_q24_close__","__file__":str(CLOSE)}
buf=io.StringIO()
try:
    sys.argv=[str(CLOSE)]
    with contextlib.redirect_stdout(buf):
        exec(compile(CLOSE.read_text(),str(CLOSE),"exec"),scope)
finally:
    sys.argv=saved

need=("ns","D24eq","F8eq")
missing=[k for k in need if k not in scope]
if missing:
    raise SystemExit("q24 close scope missing: "+",".join(missing))

ns=matrix(ZZ,scope["ns"])
D=vector(ZZ,scope["D24eq"])
Fold=vector(ZZ,scope["F8eq"])

assert D*ns*D==0
assert D*ns*Fold==2

split=primitive_hyperbolic_split(ns,D)
raw=matrix(ZZ,split["child_frame"])
minimized=minimize_child_frame(raw)
frame=matrix(ZZ,minimized["frame"])
root_data=tuple(map(int,minimized["root_data"]))

assert frame.dimensions()==(17,17)
assert frame.det()==948
assert root_data==(12,264,4),root_data

rr=12
R=frame[:rr,:rr]
C=frame[:rr,rr:]
T=frame[rr:,rr:]
H=T-C.transpose()*R.inverse()*C
assert H.dimensions()==(5,5)
assert H.det()==QQ(237)

OUTDIR.mkdir(parents=True,exist_ok=True)
OUTFRAME.write_text(
    "# q24 equation-side D12/MW5 child\n"
    "# source = close_h92_q8_q24_by_q6_translation.sage::D24eq\n"
    "# old-fibre degree = 2\n"
    + "\n".join(" ".join(map(str,row)) for row in frame.rows())
    + "\n"
)

payload={
    "schema":"elkies-k3.h3-q24-native-d12-frame.v1",
    "status":"PASS_Q24_NATIVE_D12_FRAME",
    "source_divisor":"D24eq",
    "old_fibre_degree":2,
    "root_data":list(root_data),
    "mw_rank":5,
    "mw_height":[[str(v) for v in row] for row in H.rows()],
    "mw_height_det":str(H.det()),
    "frame":str(OUTFRAME.relative_to(ROOT)),
    "split_transport":[[int(v) for v in row] for row in matrix(ZZ,split["transport"]).rows()],
}
OUTMETA.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q24NATIVE_D12|"
    f"root_data={root_data[0]},{root_data[1]},{root_data[2]}|MW=5|"
    f"height_det={H.det()}|frame={OUTFRAME.relative_to(ROOT)}|"
    "status=PASS_Q24_NATIVE_D12_FRAME",
    flush=True,
)
print(f"OUTPUT|{OUTMETA}",flush=True)

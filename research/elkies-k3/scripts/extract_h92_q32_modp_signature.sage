#!/usr/bin/env sage -python
import argparse, contextlib, io, json, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
COMP=SCRIPTS/"compile_h92_q32_d12_from_divval_modp.sage"

parser=argparse.ArgumentParser()
parser.add_argument("--prime",type=int,required=True)
args=parser.parse_args()

saved=list(sys.argv)
scope={"__name__":"__embedded__"}
buf=io.StringIO()
try:
    sys.argv=[str(COMP),"--prime",str(args.prime)]
    with contextlib.redirect_stdout(buf):
        exec(compile(COMP.read_text(),str(COMP),"exec"),scope)
finally:
    sys.argv=saved

need=("ambient_kernel","jacA","jacB","quartic","VR","VF","root_rank","root_det")
missing=[k for k in need if k not in scope]
if missing:
    raise SystemExit("compiler scope missing "+",".join(missing))

plane=scope["ambient_kernel"].echelon_form()
assert plane.dimensions()==(2,56)
VR=scope["VR"]
VF=scope["VF"]

def norm_rf(v):
    v=VF(v)
    n=VR(v.numerator())
    d=VR(v.denominator())
    lc=d.leading_coefficient()
    n/=lc
    d/=lc
    return {
        "num_degree":int(n.degree()),
        "den_degree":int(d.degree()),
        "num":[int(x) for x in n.list()],
        "den":[int(x) for x in d.list()],
    }

quartic=scope["quartic"]
payload={
    "schema":"elkies-k3.h3-q32-modp-signature.v1",
    "status":"PASS_Q32_MODP_SIGNATURE",
    "prime":int(args.prime),
    "plane_rref_2x56":[[int(v) for v in row] for row in plane.rows()],
    "plane_pivots":[int(x) for x in plane.pivots()],
    "quartic_degree":int(quartic.degree()),
    "quartic_coefficients":[norm_rf(quartic[i]) for i in range(quartic.degree()+1)],
    "jacobian_A":norm_rf(scope["jacA"]),
    "jacobian_B":norm_rf(scope["jacB"]),
    "child_root_rank":int(scope["root_rank"]),
    "child_root_det":int(scope["root_det"]),
}
out=LOCAL/f"q32-signature-mod-{args.prime}.json"
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(
    "Q32SIGNATURE|"
    f"prime={args.prime}|pivots={','.join(map(str,payload['plane_pivots']))}|"
    f"quartic={payload['quartic_degree']}|"
    f"A={payload['jacobian_A']['num_degree']}/{payload['jacobian_A']['den_degree']}|"
    f"B={payload['jacobian_B']['num_degree']}/{payload['jacobian_B']['den_degree']}|"
    f"root={payload['child_root_rank']},{payload['child_root_det']}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{out}",flush=True)

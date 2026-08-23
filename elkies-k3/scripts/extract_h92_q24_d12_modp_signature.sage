#!/usr/bin/env sage -python
"""
Extract the explicit modular D12 child produced by the q24 resolved RR compiler.
This is the q24/orbit85 analogue of the old q32 signature extractor.
"""
import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
    ]
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
COMP=SCRIPTS/"probe_h92_q24_d12_component_valuation_rr_modp.sage"

if not COMP.exists():
    raise SystemExit(f"missing q24 compiler {COMP}")

saved=list(sys.argv)
scope={"__name__":"__embedded_q24_signature__"}
buf=io.StringIO()
try:
    sys.argv=[str(COMP),"--prime",str(args.prime)]
    with contextlib.redirect_stdout(buf):
        exec(compile(COMP.read_text(),str(COMP),"exec"),scope)
finally:
    sys.argv=saved

for line in buf.getvalue().splitlines():
    if (
        line.startswith("Q24GEOM_")
        or line.startswith("Q24DIVVAL_RR")
        or line.startswith("Q24DIVVAL_QUARTIC")
        or line.startswith("Q24DIVVAL_CHILD")
        or line.startswith("Q24DIVVAL_RESULT")
    ):
        print(line,flush=True)

need=(
    "final58","jacA","jacB","quartic","VR","VF",
    "root_rank","root_det","euler","final_dimension",
    "resolved_rank","ambient","collision_rank","post_dim",
    "terminal_status",
)
missing=[name for name in need if name not in scope]
if missing:
    raise SystemExit("q24 compiler scope missing: "+",".join(missing))

plane=scope["final58"].echelon_form()
assert plane.dimensions()==(2,56)
assert int(scope["final_dimension"])==2
assert int(scope["resolved_rank"])==6
assert int(scope["collision_rank"])==48
assert int(scope["post_dim"])==8
assert int(scope["root_rank"])==12
assert int(scope["root_det"])==4
assert int(scope["euler"])==24
quartic=scope["quartic"]
assert int(quartic.degree())==4
compiler_terminal_status=str(scope["terminal_status"])
if compiler_terminal_status == "PASS_H3_Q24_EFFECTIVE_D13_D12_MODP":
    signature_status="PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE"
else:
    signature_status="CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE"

VR=scope["VR"]
VF=scope["VF"]

def norm_rf(value):
    value=VF(value)
    n=VR(value.numerator())
    d=VR(value.denominator())
    lc=d.leading_coefficient()
    n/=lc
    d/=lc
    return {
        "num_degree":int(n.degree()),
        "den_degree":int(d.degree()),
        "num":[int(v) for v in n.list()],
        "den":[int(v) for v in d.list()],
    }

payload={
    "schema":"elkies-k3.h3-q24-orbit85-d12-modp-signature.v1",
    "status":signature_status,
    "component_compiler_status":compiler_terminal_status,
    "prime":int(args.prime),
    "source_neighbor":{
        "q":24,
        "orbit":85,
        "source":"D13/MW4",
        "child":"D12/MW5",
    },
    "rr":{
        "ambient":56,
        "collision_rank":48,
        "post_collision":8,
        "resolved_rank":6,
        "kernel":2,
        "geometric_fibre_twist":-8,
    },
    "plane_rref_2x56":[[int(v) for v in row] for row in plane.rows()],
    "plane_pivots":[int(v) for v in plane.pivots()],
    "quartic_degree":4,
    "quartic_coefficients":[norm_rf(quartic[i]) for i in range(5)],
    "jacobian_A":norm_rf(scope["jacA"]),
    "jacobian_B":norm_rf(scope["jacB"]),
    "child_root_rank":12,
    "child_root_det":4,
    "child_euler":24,
    "proof_boundary":(
        "Explicit modular q24/orbit85 resolved-RR signature. "
        "The embedded component compiler status is recorded separately; "
        "characteristic-zero replay remains separate."
    ),
}

OUT=(args.output.resolve() if args.output else
     LOCAL/f"q24-orbit85-d12-signature-mod-{args.prime}.json")
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q24D12SIGNATURE|"
    f"prime={args.prime}|plane=2x56|pivots={','.join(map(str,payload['plane_pivots']))}|"
    "RR=56,48,8,6,2|quartic=4|root=12,4|euler=24|"
    f"compiler_status={compiler_terminal_status}|status={signature_status}",
    flush=True,
)
print(f"OUTPUT|{OUT}",flush=True)

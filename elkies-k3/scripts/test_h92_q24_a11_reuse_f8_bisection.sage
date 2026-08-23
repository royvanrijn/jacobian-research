#!/usr/bin/env sage -python
"""
Test whether a q24-native D12 -> A11 divisor has the same generic degree-two
class as the old q8 fibre F8.

If D_A11 - F8 is only q24-fibre + D12-root vertical support, then no new D12
section is needed: the generic degree-two linear system is already represented
by the old q8 pencil, and the A11 hop reduces to imposing the target vertical
D12 component profile on that known bisection pencil.
"""
import argparse
import contextlib
import io
import json
import subprocess
import sys
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift/bisection-reuse"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"
SEARCH=SCRIPTS/"search_root_adapted_weyl_neighbors.sage"
for path in (ENGINE,CLOSE,SEARCH):
    if not path.exists():
        raise SystemExit(f"missing {path}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

# Exact q24 divisor and old q8 fibre in H3 equation NS.
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
D24=vector(ZZ,scope["D24eq"])
F8=vector(ZZ,scope["F8eq"])
assert D24*ns*D24==0
assert D24*ns*F8==2

# Build q24-native D12 frame and retain full parent transport.
split=primitive_hyperbolic_split(ns,D24)
raw=matrix(ZZ,split["child_frame"])
mini=minimize_child_frame(raw)
G=matrix(ZZ,mini["frame"])
A=matrix(ZZ,mini["basis"])
assert tuple(map(int,mini["root_data"]))==(12,264,4)

B=block_diagonal_matrix(identity_matrix(ZZ,2),A)*matrix(ZZ,split["transport"])
U=matrix(ZZ,((0,1),(1,0)))
assert B*ns*B.transpose()==block_diagonal_matrix(U,-G)

def child_coords(parent):
    q=vector(QQ,parent)*B.inverse()
    assert all(v in ZZ for v in q)
    return vector(ZZ,[ZZ(v) for v in q])

f8c=child_coords(F8)
assert f8c[1]==2

R=G[:12,:12]
C=G[:12,12:]
T=G[12:,12:]
H=T-C.transpose()*R.inverse()*C
assert H.dimensions()==(5,5)

# In a root-adapted frame, the MW quotient coordinate of a full NS vector is
# simply the last five positive-frame coordinates.
f8_mw=vector(ZZ,f8c[-5:])

OUTDIR.mkdir(parents=True,exist_ok=True)
FRAME=OUTDIR/"q24-d12-frame.txt"
FRAME.write_text(
    "# q24-native D12/MW5 frame for F8-bisection reuse test\n"
    + "\n".join(" ".join(map(str,row)) for row in G.rows())
    + "\n"
)

SEARCHOUT=OUTDIR/"q24-d12-q6-all.json"
FRAMES=OUTDIR/"frames"
cmd=[
    "sage","-python",str(SEARCH),
    "--frame",str(FRAME),
    "--root-rank","12",
    "--q","6",
    "--degree","2",
    "--adapt-mw-at-least","6",
    "--rank-growth-only",
    "--output",str(SEARCHOUT),
    "--frames-dir",str(FRAMES),
]
print("+"," ".join(cmd),flush=True)
subprocess.run(cmd,cwd=str(ROOT),check=True)

data=json.loads(SEARCHOUT.read_text())
assert data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
a11=[
    rec for rec in data["neighbors"]
    if tuple(rec["child_root_data"])==(11,132,12)
    and rec["child_ade"]=="A11"
    and int(rec["child_mw_rank"])==6
]
assert a11

records=[]
reuse=[]
for rec in a11:
    D11=vector(ZZ,rec["fiber"])
    assert D11[1]==2

    delta=D11-f8c
    # Same generic degree-two class iff MW quotient vanishes.
    delta_mw=vector(ZZ,delta[-5:])
    same_generic=(delta_mw==0)

    vertical_root=None
    fibre_shift=None
    reconstructed=False
    if same_generic:
        # delta has b=0. It must be k*F + root-lattice vector.
        assert delta[1]==0
        fibre_shift=ZZ(delta[0])
        vertical_root=vector(ZZ,delta[2:14])
        assert all(v==0 for v in delta[14:])
        recon=(
            fibre_shift*vector(ZZ,[1,0]+[0]*17)
            + vector(ZZ,[0,0]+list(vertical_root)+[0]*5)
        )
        reconstructed=(recon==delta)
        assert reconstructed

    row={
        "orbit_index":int(rec["orbit_index"]),
        "A11_mw_projection":rec["mw_projection"],
        "A11_dominant_labels":rec["dominant_labels"],
        "A11_fibre":rec["fiber"],
        "F8_child_coordinates":list(map(int,f8c)),
        "F8_mw_projection":list(map(int,f8_mw)),
        "delta_mw_projection":list(map(int,delta_mw)),
        "same_generic_degree_two_class_as_F8":same_generic,
        "q24_fibre_shift":None if fibre_shift is None else int(fibre_shift),
        "D12_root_vertical_coefficients":(
            None if vertical_root is None else list(map(int,vertical_root))
        ),
    }
    records.append(row)
    if same_generic:
        reuse.append(row)

    print(
        "Q24A11_F8REUSE|"
        f"orbit={rec['orbit_index']}|"
        f"F8_mw={','.join(map(str,f8_mw))}|"
        f"A11_mw={','.join(map(str,rec['mw_projection']))}|"
        f"delta_mw={','.join(map(str,delta_mw))}|"
        f"same_generic={int(same_generic)}|"
        f"fibre_shift={fibre_shift}|status={'PASS_REUSE' if same_generic else 'NO'}",
        flush=True,
    )

# Also test affine shifts by the already-explicit A0-R3 section direction if
# a prior explicit-curve audit is available.  A degree-two divisor class can
# be translated by a Pic^0 section; this tells us how far each A11 target is
# from the known F8 bisection in MW quotient.
AUDIT=LOCAL/"q24-downstream-lift/explicit-curves-a11-span-p100003.json"
known_direction=None
if AUDIT.exists():
    audit=json.loads(AUDIT.read_text())
    known_direction=vector(ZZ,audit["A0_mw_projection"])
    for row in records:
        d=vector(ZZ,row["delta_mw_projection"])
        scalar=None
        for n in range(-12,13):
            if d==n*known_direction:
                scalar=n
                break
        row["delta_is_scalar_A0_minus_R3"]=scalar is not None
        row["delta_A0_minus_R3_scalar"]=scalar

payload={
    "schema":"elkies-k3.h3-q24-a11-f8-bisection-reuse.v1",
    "status":(
        "PASS_Q24_A11_REUSES_F8_GENERIC_BISECTION"
        if reuse else
        "Q24_A11_NEEDS_NONTRIVIAL_PIC0_TRANSLATION_OF_F8"
    ),
    "prime":int(args.prime),
    "F8_child_coordinates":list(map(int,f8c)),
    "F8_mw_projection":list(map(int,f8_mw)),
    "known_A0_minus_R3_direction":(
        None if known_direction is None else list(map(int,known_direction))
    ),
    "A11_candidates":records,
    "reusable_candidates":reuse,
    "interpretation":(
        "same_generic_degree_two_class_as_F8 means D_A11-F8 is only q24 fibre "
        "plus D12 root support. Then the old q8 pencil already supplies the "
        "generic two-dimensional degree-two system; only vertical D12 "
        "trivialization/conditions differ."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/"q24-a11-f8-bisection-reuse.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24A11_F8REUSE_RESULT|"
    f"A11={len(records)}|reusable={len(reuse)}|status={payload['status']}",
    flush=True,
)

#!/usr/bin/env sage -python
"""
Lift every q24-native D12 --q6--> A11 neighbour back into the exact q8/D13
equation frame.

This avoids the missing second D12 section.  For each A11 fibre we compute:
  * its exact class in the original H3/q8 NS frame;
  * old-q8 degree;
  * D13 root/MW coordinates in the anchored q8 equation frame;
  * the equivalent marked q8 section MW vector and predicted P.O;
  * exact vertical correction relative to (degree-1) O8 + P.

The cheapest candidate is exported for direct equation construction on the
already-certified q8/D13 Weierstrass model.
"""
import argparse
import contextlib
import io
import json
import subprocess
import sys
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix,
    lcm, matrix, vector
)

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
OUTDIR=LOCAL/"q24-downstream-lift/q8-backlift"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"
SEARCH=SCRIPTS/"search_root_adapted_weyl_neighbors.sage"

for path in (ENGINE,CLOSE,SEARCH):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

exec(compile(ENGINE.read_text(),str(ENGINE),"exec"))

# -------------------------------------------------------------------------
# 1. Load the exact q8 equation frame and q24 divisor.
# -------------------------------------------------------------------------
saved=list(sys.argv)
cl={"__name__":"__embedded_q24_close__","__file__":str(CLOSE)}
buf=io.StringIO()
try:
    sys.argv=[str(CLOSE)]
    with contextlib.redirect_stdout(buf):
        exec(compile(CLOSE.read_text(),str(CLOSE),"exec"),cl)
finally:
    sys.argv=saved

need=("ns","D24eq","F8eq","O8","Badapt","adapted","H13")
missing=[k for k in need if k not in cl]
if missing:
    raise SystemExit("q24 close scope missing: "+",".join(missing))

ns=matrix(ZZ,cl["ns"])
D24=vector(ZZ,cl["D24eq"])
F8=vector(ZZ,cl["F8eq"])
O8=vector(ZZ,cl["O8"])
Badapt=matrix(ZZ,cl["Badapt"])
G13=matrix(ZZ,cl["adapted"])
H13=matrix(QQ,cl["H13"])

assert D24*ns*D24==0
assert D24*ns*F8==2
assert O8*ns*F8==1
assert G13.dimensions()==(17,17)
assert H13.dimensions()==(4,4)
assert abs(Badapt.det())==1

# -------------------------------------------------------------------------
# 2. Rebuild the exact q24-native D12 frame + full parent transport.
# -------------------------------------------------------------------------
split=primitive_hyperbolic_split(ns,D24)
raw=matrix(ZZ,split["child_frame"])
mini=minimize_child_frame(raw)
G12=matrix(ZZ,mini["frame"])
A12=matrix(ZZ,mini["basis"])
assert tuple(map(int,mini["root_data"]))==(12,264,4)
assert G12.det()==948

B24=block_diagonal_matrix(identity_matrix(ZZ,2),A12)*matrix(ZZ,split["transport"])
U2=matrix(ZZ,((0,1),(1,0)))
assert B24*ns*B24.transpose()==block_diagonal_matrix(U2,-G12)

OUTDIR.mkdir(parents=True,exist_ok=True)
FRAME=OUTDIR/"q24-native-d12-frame.txt"
FRAME.write_text(
    "# q24-native D12/MW5 frame for exhaustive A11 back-lift\n"
    + "\n".join(" ".join(map(str,row)) for row in G12.rows())
    + "\n"
)

# -------------------------------------------------------------------------
# 3. Exhaust all q6 neighbours and keep A11/MW6.
# -------------------------------------------------------------------------
SEARCHOUT=OUTDIR/"q24-native-d12-q6-all.json"
FRAMES=OUTDIR/"q24-native-d12-q6-frames"
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
if not a11:
    raise SystemExit("no q24-native q6 -> A11 neighbours found")

# -------------------------------------------------------------------------
# 4. Helpers for D13 section profile.
# -------------------------------------------------------------------------
R13=G13[:13,:13]
C13=G13[:13,13:]
T13=G13[13:,13:]
assert T13-C13.transpose()*R13.inverse()*C13==H13

def class_order(dual):
    o=ZZ(1)
    for x in dual:
        o=lcm(o,ZZ(QQ(x).denominator()))
    return o

def d13_correction(z):
    z=vector(ZZ,z)
    base=vector(ZZ,[0]*13+list(z))
    dual=vector(QQ,base*G13[:,:13])*R13.inverse()
    order=class_order(dual)
    expected={ZZ(1):QQ(0),ZZ(2):QQ(1),ZZ(4):QQ(13)/4}
    if order not in expected:
        return order,None
    corr=expected[order]
    raw=QQ(dual*R13*dual)
    mod2=lambda x: QQ(x)-2*(QQ(x)/2).floor()
    if mod2(raw)!=mod2(corr):
        raise ArithmeticError("D13 correction parity mismatch")
    return order,corr

def minimal_section_frame(z,corr):
    z=vector(ZZ,z)
    h=QQ(z*H13*z)
    target=h+corr
    base=vector(ZZ,[0]*13+list(z))
    dual=vector(QQ,base*G13[:,:13])*R13.inverse()

    L=IntegralLattice(R13)
    it=L.enumerate_close_vectors(-dual)
    for _ in range(100000):
        shift=vector(ZZ,next(it))
        pf=base+vector(ZZ,list(shift)+[0]*4)
        norm=QQ(pf*G13*pf)
        if norm==target:
            return pf
        if norm>target:
            break
    raise ArithmeticError("could not recover minimal D13 section frame")

# -------------------------------------------------------------------------
# 5. Lift every A11 fibre back through q24 into q8/D13.
# -------------------------------------------------------------------------
records=[]
for rec in a11:
    fibre12=vector(ZZ,rec["fiber"])
    assert fibre12[1]==2

    # B24 rows are D12 child NS basis vectors in the q8/H3 parent.
    D11=vector(ZZ,fibre12*B24)
    assert D11*ns*D11==0

    degree=ZZ(D11*ns*F8)
    dO=ZZ(D11*ns*O8)

    q8coords=vector(QQ,D11)*Badapt.inverse()
    if not all(v in ZZ for v in q8coords):
        raise ArithmeticError("A11 back-lift not integral in q8 adapted frame")
    q8coords=vector(ZZ,[ZZ(v) for v in q8coords])
    assert q8coords[1]==degree

    z=vector(ZZ,q8coords[-4:])
    height=QQ(z*H13*z)
    order,corr=d13_correction(z)

    po=None
    pframe=None
    vertical_root=None
    fibre_twist=None
    if corr is not None:
        poQ=(height+corr-4)/2
        if poQ in ZZ and poQ>=0:
            po=ZZ(poQ)
            pframe=minimal_section_frame(z,corr)

            # Generic restriction in q8 adapted NS:
            # (degree-1)O + P + V + kF.
            # O = (-1,1,0), P=(po+1,1,pframe)
            # target D coords are q8coords.
            Ocoords=vector(ZZ,[-1,1]+[0]*17)
            Pcoords=vector(ZZ,[po+1,1]+list(pframe))
            residual=q8coords-(degree-1)*Ocoords-Pcoords

            # Residual must be vertical: MW tail zero.
            assert all(v==0 for v in residual[2+13:])
            fibre_twist=ZZ(residual[0])
            vertical_root=vector(ZZ,residual[2:2+13])

            # Regression.
            reconstructed=(
                (degree-1)*Ocoords
                +Pcoords
                +vector(ZZ,[0,0]+list(vertical_root)+[0]*4)
                +fibre_twist*vector(ZZ,[1,0]+[0]*17)
            )
            assert reconstructed==q8coords

    row={
        "orbit_index":int(rec["orbit_index"]),
        "native_d12_mw_projection":rec["mw_projection"],
        "native_d12_dominant_labels":rec["dominant_labels"],
        "native_d12_witness":rec["witness"],
        "q8_degree":int(degree),
        "q8_D_dot_O":int(dO),
        "q8_full_coordinates":list(map(int,q8coords)),
        "q8_horizontal_mw_projection":list(map(int,z)),
        "q8_horizontal_height":str(height),
        "q8_component_class_order":int(order),
        "q8_local_correction":None if corr is None else str(corr),
        "q8_marked_section_P_dot_O":None if po is None else int(po),
        "q8_marked_section_frame":None if pframe is None else list(map(int,pframe)),
        "q8_vertical_root_coefficients":(
            None if vertical_root is None else list(map(int,vertical_root))
        ),
        "q8_fibre_twist":None if fibre_twist is None else int(fibre_twist),
        "formula":"D=(q8_degree-1)O8+P+V+kF8",
    }
    records.append(row)

    print(
        "Q24A11BACKLIFT|"
        f"orbit={rec['orbit_index']}|q8_degree={degree}|DdotO={dO}|"
        f"mw={','.join(map(str,z))}|height={height}|"
        f"class_order={order}|corr={corr}|PdotO={po}|"
        f"fibre_twist={fibre_twist}|status=PASS",
        flush=True,
    )

# Prefer low q8 degree, then low marked-section collision degree, then small MW.
def score(row):
    po=row["q8_marked_section_P_dot_O"]
    return (
        row["q8_degree"],
        10**9 if po is None else po,
        sum(abs(int(x)) for x in row["q8_horizontal_mw_projection"]),
        row["orbit_index"],
    )

records.sort(key=score)
best=records[0]

print(
    "Q24A11BACKLIFT_BEST|"
    f"orbit={best['orbit_index']}|q8_degree={best['q8_degree']}|"
    f"mw={','.join(map(str,best['q8_horizontal_mw_projection']))}|"
    f"PdotO={best['q8_marked_section_P_dot_O']}|"
    f"fibre_twist={best['q8_fibre_twist']}|status=SELECTED",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-a11-backlift-to-q8.v1",
    "status":"PASS_Q24_A11_BACKLIFT_TO_Q8",
    "prime":int(args.prime),
    "source":"q24-native D12 derived directly from D24eq",
    "A11_candidate_count":len(records),
    "candidates":records,
    "selected":best,
    "q8_frame":{
        "root_data":[13,312,4],
        "mw_rank":4,
        "old_fibre":"F8eq",
        "zero":"O8",
    },
    "proof_boundary":(
        "Every A11 fibre is lifted exactly through the q24 neighbour to the "
        "anchored q8/D13 NS frame. The horizontal MW projection, local "
        "correction, P.O and D=(d-1)O+P+V+kF decomposition are exact lattice "
        "data. This does not yet provide equation coordinates for the marked "
        "q8 section P when that section is not already explicit."
    ),
}
OUT=args.output.resolve() if args.output else OUTDIR/"q24-a11-backlift-to-q8.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24A11BACKLIFT_RESULT|"
    f"candidates={len(records)}|selected_orbit={best['orbit_index']}|"
    f"q8_degree={best['q8_degree']}|PdotO={best['q8_marked_section_P_dot_O']}|"
    "status=PASS_Q24_A11_BACKLIFT_TO_Q8",
    flush=True,
)

#!/usr/bin/env sage -python
"""
Close q24-native A1 -> pinned R17 by transporting the already-certified
historical final q6 witness through an exact integral isometry of the A1
positive frames.

This is NOT a historical-frame splice in the route: the output is an explicit
q6 witness in native-A1 coordinates.  The historical A1 frame is used only as
a coordinate chart from which an already-certified isotropic class is pulled
back by a checked determinant-one isometry.

Certificates:
  B * G_native_A1 * B^t = G_historical_A1
  native_witness = historical_witness * B
  native fibre = (3,2,native_witness), primitive, square 0, old-fibre degree 2
  the transported full NS map carries that fibre to the historical final q6
  fibre and then exactly to pinned R17.
"""

import json
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector
)

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3" / "q24-native-suffix"
GEN = ROOT / "artifacts" / "generated-results"

SUFFIX = LOCAL / "q24-native-d12-to-a1.json"
HIST_A1_SEARCH = GEN / "elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json"
HIST_FINAL_SEARCH = GEN / "elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json"
PINNED_PATH = ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt"
ENDPOINT_ISO_PATH = (
    ROOT / "elkies-k3" / "data" / "fibrations"
    / "h3_rootless_mw17_to_pinned_rank17_isometry.txt"
)
OUT = LOCAL / "q24-native-a1-to-pinned-r17-transported.json"

U2 = matrix(ZZ, ((0,1),(1,0)))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(M):
    return [[int(v) for v in row] for row in M.rows()]


def find_basis_map(source, target):
    """
    Return B with B*source*B^t = target.

    PARI qfisom(source,target) convention:
        source = S^t * target * S.
    Therefore B=S^{-t}.  We still test several equivalent orientations and
    accept only an exact integral determinant-one identity.
    """
    trials = []
    for G,H in ((source,target),(target,source)):
        raw = pari(G).qfisom(pari(H))
        if str(raw) == "0":
            continue
        S = matrix(ZZ,raw)
        trials += [S,S.transpose()]
        try:
            Si=S.inverse()
            if Si.change_ring(ZZ)==Si:
                Si=Si.change_ring(ZZ)
                trials += [Si,Si.transpose()]
        except Exception:
            pass

    seen=set()
    for B in trials:
        key=tuple(B.list())
        if key in seen:
            continue
        seen.add(key)
        if B.dimensions()!=(17,17):
            continue
        if B*source*B.transpose()==target:
            assert abs(ZZ(B.det()))==1
            return B
    return None


for path in (
    SUFFIX, HIST_A1_SEARCH, HIST_FINAL_SEARCH,
    PINNED_PATH, ENDPOINT_ISO_PATH,
):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

suffix=json.loads(SUFFIX.read_text())
assert suffix["status"]=="PASS_Q24_NATIVE_D12_TO_A1"
native_path=ROOT/suffix["final_frame"]
if not native_path.exists():
    raise SystemExit(f"missing native A1 frame: {native_path}")
N=load_gram(native_path)
assert N.det()==948

# Historical A1 is the selected orbit981 child of the 2A1 parent.
a1data=json.loads(HIST_A1_SEARCH.read_text())
assert a1data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
a1hits=[
    r for r in a1data["neighbors"]
    if int(r["orbit_index"])==981
    and r["child_ade"]=="A1"
    and int(r["child_mw_rank"])==16
]
assert len(a1hits)==1
a1rec=a1hits[0]
H=matrix(ZZ,a1rec["child_root_adapted_frame"])
assert H.det()==948
assert tuple(a1rec["child_root_data"])==(1,2,2)

# Confirm final-q6 artifact literally uses this historical A1 frame.
final=json.loads(HIST_FINAL_SEARCH.read_text())
assert final["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
declared=ROOT/final["frame"]
assert load_gram(declared)==H
fhits=[
    r for r in final["neighbors"]
    if int(r["orbit_index"])==2247
    and int(r["q"])==6
    and int(r["old_fiber_degree"])==2
    and tuple(r["child_root_data"])==(0,0,1)
    and r["child_ade"]=="rootless"
    and int(r["child_mw_rank"])==17
]
assert len(fhits)==1
frec=fhits[0]

print(
    "Q24NATIVE_R17_TRANSPORT_INPUT|"
    f"native={native_path.relative_to(ROOT)}|"
    "historical=A1/orbit981|historical_final_q6=orbit2247|status=PASS",
    flush=True,
)

# Exact A1-frame isometry.
B=find_basis_map(N,H)
if B is None:
    print(
        "Q24NATIVE_R17_A1_ISOMETRY|isometric=0|"
        "status=NATIVE_A1_DISTINCT_FROM_HISTORICAL_A1",
        flush=True,
    )
    raise SystemExit(
        "Native and historical A1 positive frames are not integrally isometric. "
        "Then the only remaining route is a larger fresh q6 shell search."
    )

assert B*N*B.transpose()==H
assert abs(ZZ(B.det()))==1
print(
    "Q24NATIVE_R17_A1_ISOMETRY|"
    f"isometric=1|det={B.det()}|status=PASS",
    flush=True,
)

# B gives historical basis vectors in native coordinates:
# H = B*N*B^t.  Hence historical row coordinates w_H pull back as w_N=w_H*B.
wH=vector(ZZ,frec["witness"])
wN=vector(ZZ,wH*B)
assert wH*H*wH==12
assert wN*N*wN==12

fH=vector(ZZ,frec["fiber"])
fN=vector(ZZ,[3,2]+list(wN))
assert fH==vector(ZZ,[3,2]+list(wH))

NSN=block_diagonal_matrix(U2,-N)
NSH=block_diagonal_matrix(U2,-H)
fullB=block_diagonal_matrix(identity_matrix(ZZ,2),B)
assert fullB*NSN*fullB.transpose()==NSH
assert fN*fullB.transpose()==fH
assert fN*NSN*fN==0
oldF=vector(ZZ,[1,0]+[0]*17)
assert fN*NSN*oldF==2

from sage.all import gcd
assert gcd(tuple(fN))==1

print(
    "Q24NATIVE_R17_NATIVE_Q6|"
    f"historical_orbit=2247|native_witness_norm={wN*N*wN}|"
    "primitive=1|old_fibre_degree=2|"
    f"mw_projection_native={','.join(map(str,wN[1:])) if False else 'transported'}|"
    "status=PASS_EXPLICIT_NATIVE_Q6_FIBRE",
    flush=True,
)

# Historical q6 transition: H A1 -> historical rootless child.
raw=matrix(ZZ,frec["child_frame"])
adapt=matrix(ZZ,frec["child_root_adapted_basis"])
hist_child=matrix(ZZ,frec["child_root_adapted_frame"])
neighbor=matrix(ZZ,frec["neighbor_basis"])
assert adapt*raw*adapt.transpose()==hist_child
assert ZZ(pari(hist_child).qfminim(2)[0])==0
assert hist_child.det()==948

Thist=block_diagonal_matrix(identity_matrix(ZZ,2),adapt)*neighbor
NShistchild=block_diagonal_matrix(U2,-hist_child)
assert Thist*NSH*Thist.transpose()==NShistchild
assert abs(ZZ(Thist.det()))==1

# Stored endpoint isometry C satisfies C^t * pinned * C = historical endpoint.
pinned=load_gram(PINNED_PATH)
C=load_gram(ENDPOINT_ISO_PATH)
assert C.dimensions()==pinned.dimensions()==(17,17)
assert C.transpose()*pinned*C==hist_child
assert C.det()==1

# E maps historical endpoint NS coordinates to pinned NS coordinates.
# diag(I,C^t) maps pinned -> historical endpoint, so invert it.
pinned_to_hist=block_diagonal_matrix(identity_matrix(ZZ,2),C.transpose())
assert pinned_to_hist*block_diagonal_matrix(U2,-pinned)*pinned_to_hist.transpose()==NShistchild
E=pinned_to_hist.inverse().change_ring(ZZ)
assert E*NShistchild*E.transpose()==block_diagonal_matrix(U2,-pinned)
assert abs(ZZ(E.det()))==1

# Native A1 -> historical A1 -> historical rootless -> pinned R17.
Ttotal=E*Thist*fullB
assert (
    Ttotal*NSN*Ttotal.transpose()
    == block_diagonal_matrix(U2,-pinned)
)
assert abs(ZZ(Ttotal.det()))==1

# Check fibre semantics across first transport and historical q6 step.
assert fN*fullB.transpose()==fH

print(
    "Q24NATIVE_R17_ENDPOINT|"
    "historical_child_rootless=1|historical_child_det=948|"
    "pinned_endpoint_match=1|"
    f"total_transport_det={Ttotal.det()}|"
    "status=PASS_PINNED_R17",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-native-a1-to-pinned-r17-transport.v1",
    "status":"PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    "native_parent":{
        "frame":str(native_path.relative_to(ROOT)),
        "ade":"A1","root_data":[1,2,2],"mw_rank":16,
    },
    "historical_coordinate_chart":{
        "a1_source_artifact":str(HIST_A1_SEARCH.relative_to(ROOT)),
        "historical_a1_orbit":981,
        "native_to_historical_positive_frame_basis_map":rows(B),
        "basis_map_identity":"B * G_native * B^t = G_historical",
    },
    "native_q6":{
        "q":6,
        "factor_order":[3,2],
        "old_fibre_degree":2,
        "historical_source_orbit":2247,
        "historical_witness":list(map(int,wH)),
        "native_witness":list(map(int,wN)),
        "native_fibre":list(map(int,fN)),
        "native_witness_norm":12,
        "primitive":True,
    },
    "endpoint":{
        "historical_rootless_frame":rows(hist_child),
        "pinned_rank17_frame":str(PINNED_PATH.relative_to(ROOT)),
        "historical_endpoint_isometry_file":str(ENDPOINT_ISO_PATH.relative_to(ROOT)),
        "native_a1_to_pinned_r17_ns_transport":rows(Ttotal),
        "transport_determinant":int(Ttotal.det()),
    },
    "proof_boundary":(
        "Exact lattice/NS certificate. The historical A1 frame is used only "
        "as an integrally isometric coordinate chart to pull back the already "
        "certified final q6 fibre. The resulting q6 fibre is explicitly given "
        "and checked in the q24-native A1 coordinates, and the composite "
        "determinant-one NS transport lands on the pinned recovered R17 frame. "
        "This does not execute the characteristic-zero A1->R17 Weierstrass pencil."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24NATIVE_R17_RESULT|"
    "native_q6=1|rootless=1|pinned_R17=1|"
    "status=PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    flush=True,
)

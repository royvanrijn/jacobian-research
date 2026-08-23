#!/usr/bin/env sage -python
"""
Pull the pinned R17 fibration itself back into the q24-native A1 coordinates.

This avoids shell search entirely.  We compose the exact 19-dimensional
basis transports

    pinned R17 -> H3 source -> q24-native D12 -> ... -> native A1

and invert the composite.  The first row of the resulting A1->pinned basis
map is exactly the pinned R17 fibre class expressed in native-A1 coordinates.

The output therefore determines, rather than guesses/searches, the direct
neighbor parameters (a,b), q=a*b, and old-fibre degree b from native A1 to
the pinned R17 fibration.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, pari, vector
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL_ROOT = ROOT / "artifacts" / "local" / "elkies-k3"
LOCAL = LOCAL_ROOT / "q24-native-suffix"
GEN = ROOT / "artifacts" / "generated-results"

ENGINE = SCRIPTS / "exact_neighbor_engine.sage"
CLOSE = SCRIPTS / "close_h92_q8_q24_by_q6_translation.sage"
SUFFIX = LOCAL / "q24-native-d12-to-a1.json"
REVERSE = GEN / "elkies-k3-rank17-to-h3-reverse-transport.json"
PINNED_PATH = ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt"
OUT = LOCAL / "q24-native-a1-pulled-pinned-r17.json"

U2 = matrix(ZZ, ((0,1),(1,0)))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(M):
    return [[int(x) for x in row] for row in M.rows()]


def ns_of(frame):
    return block_diagonal_matrix(U2, -frame)


for p in (ENGINE, CLOSE, SUFFIX, REVERSE, PINNED_PATH):
    if not p.exists():
        raise SystemExit(f"missing prerequisite: {p}")

# -------------------------------------------------------------------------
# 1. Reconstruct the exact H3-source -> native D12 transport.
# -------------------------------------------------------------------------

engine_scope = {"__name__":"__q24_pullback_engine__","__file__":str(ENGINE)}
exec(compile(ENGINE.read_text(), str(ENGINE), "exec"), engine_scope)
primitive_hyperbolic_split = engine_scope["primitive_hyperbolic_split"]
minimize_child_frame = engine_scope["minimize_child_frame"]

saved = list(sys.argv)
scope = {"__name__":"__q24_pullback_close__","__file__":str(CLOSE)}
buf = io.StringIO()
try:
    sys.argv = [str(CLOSE)]
    with contextlib.redirect_stdout(buf):
        exec(compile(CLOSE.read_text(), str(CLOSE), "exec"), scope)
finally:
    sys.argv = saved

for name in ("ns","D24eq"):
    if name not in scope:
        raise SystemExit(f"close script did not expose {name}")

G_H3 = matrix(ZZ, scope["ns"])
D24 = vector(ZZ, scope["D24eq"])
assert G_H3.dimensions() == (19,19)
assert D24 * G_H3 * D24 == 0

split = primitive_hyperbolic_split(G_H3, D24)
raw_d12 = matrix(ZZ, split["child_frame"])
mini = minimize_child_frame(raw_d12)
A0 = matrix(ZZ, mini["basis"])
D12 = matrix(ZZ, mini["frame"])
assert tuple(map(int, mini["root_data"])) == (12,264,4)

T = block_diagonal_matrix(identity_matrix(ZZ,2), A0) * matrix(ZZ, split["transport"])
assert abs(ZZ(T.det())) == 1
assert T * G_H3 * T.transpose() == ns_of(D12)

print(
    "Q24R17_PULLBACK_D12|"
    "source=H3_NS|q24=D24eq|child=D12/MW5|"
    f"transport_det={T.det()}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Compose the actual local q24-native D12 -> A1 transitions.
# -------------------------------------------------------------------------

suffix = json.loads(SUFFIX.read_text())
assert suffix["status"] == "PASS_Q24_NATIVE_D12_TO_A1"
assert len(suffix["steps"]) == 9

current = D12
for item in suffix["steps"]:
    artifact = ROOT / item["search_artifact"]
    if not artifact.exists():
        raise SystemExit(f"missing native suffix search artifact: {artifact}")
    data = json.loads(artifact.read_text())
    assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

    # Strong source-frame identity: no historical frame may enter silently.
    source_path = ROOT / data["frame"]
    assert load_gram(source_path) == current

    orbit = int(item["discovered_orbit_index"])
    hits = [r for r in data["neighbors"] if int(r["orbit_index"]) == orbit]
    assert len(hits) == 1
    rec = hits[0]

    raw = matrix(ZZ, rec["child_frame"])
    adapt = matrix(ZZ, rec["child_root_adapted_basis"])
    child = matrix(ZZ, rec["child_root_adapted_frame"])
    neighbor = matrix(ZZ, rec["neighbor_basis"])
    assert adapt * raw * adapt.transpose() == child

    stepT = block_diagonal_matrix(identity_matrix(ZZ,2), adapt) * neighbor
    assert abs(ZZ(stepT.det())) == 1
    assert stepT * ns_of(current) * stepT.transpose() == ns_of(child)

    T = stepT * T
    current = child

    print(
        "Q24R17_PULLBACK_SUFFIX|"
        f"step={item['step']}|q={item['q']}|"
        f"orbit={orbit}|child={item['target_ade']}|"
        f"MW={item['target_mw_rank']}|status=PASS",
        flush=True,
    )

A1 = current
a1_path = ROOT / suffix["final_frame"]
assert load_gram(a1_path) == A1
assert tuple(map(int, suffix["steps"][-1]["target_root_data"])) == (1,2,2)
assert T * G_H3 * T.transpose() == ns_of(A1)

print(
    "Q24R17_PULLBACK_A1|"
    f"transport_det={T.det()}|root_data=1,2,2|MW=16|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 3. Exact pinned-R17 -> H3 source map from the canonical reverse ledger.
# -------------------------------------------------------------------------

rev = json.loads(REVERSE.read_text())
assert rev["status"] == "PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT"
R = matrix(ZZ, rev["complete_reverse_transport"]["pinned_rank17_to_h3_basis_R"])
assert R.dimensions() == (19,19)
assert abs(ZZ(R.det())) == 1

pinned = load_gram(PINNED_PATH)
G_PIN = ns_of(pinned)

assert R * G_PIN * R.transpose() == G_H3

# C gives the native-A1 basis in pinned-R17 coordinates:
#   G_A1 = C * G_PIN * C^t.
C = T * R
assert abs(ZZ(C.det())) == 1
assert C * G_PIN * C.transpose() == ns_of(A1)

# Therefore C^{-1} gives the pinned basis in native-A1 coordinates.
P = C.inverse()
assert P.change_ring(ZZ) == P
P = P.change_ring(ZZ)
assert abs(ZZ(P.det())) == 1
assert P * ns_of(A1) * P.transpose() == G_PIN

# Normalize sign of the pulled pinned U simultaneously if needed.
D = vector(ZZ, P.row(0))
M = vector(ZZ, P.row(1))
Fold = vector(ZZ, [1,0] + [0]*17)
degree = ZZ(D * ns_of(A1) * Fold)
if degree < 0:
    signfix = block_diagonal_matrix(matrix(ZZ,[[-1,0],[0,-1]]), identity_matrix(ZZ,17))
    P = signfix * P
    D = vector(ZZ, P.row(0))
    M = vector(ZZ, P.row(1))
    degree = ZZ(D * ns_of(A1) * Fold)

assert D * ns_of(A1) * D == 0
assert M * ns_of(A1) * M == 0
assert D * ns_of(A1) * M == 1
assert gcd(tuple(D)) == 1
assert degree > 0

a = ZZ(D[0])
b = ZZ(D[1])
w = vector(ZZ, D[2:])
q = a*b

assert b == degree
assert q > 0
assert w * A1 * w == 2*q
assert D == vector(ZZ, [a,b] + list(w))

# The complete split is already explicit: P itself sends native A1 NS
# exactly to U + (-pinned R17), with D as its first basis vector.
assert P * ns_of(A1) * P.transpose() == G_PIN
assert P.row(0) == D
assert ZZ(pari(pinned).qfminim(2)[0]) == 0
assert pinned.det() == 948

print(
    "Q24R17_PULLBACK_TARGET|"
    f"a={a}|b={b}|q={q}|old_fibre_degree={degree}|"
    f"witness_norm={w*A1*w}|primitive=1|"
    "child=rootless/MW17|pinned_R17=1|status=PASS",
    flush=True,
)

classification = (
    "EXACT_SELECTED_Q6_DEGREE2"
    if q == 6 and degree == 2
    else "EXACT_DEGREE2_DIFFERENT_Q"
    if degree == 2
    else "EXACT_DIRECT_R17_FIBRATION_NOT_DEGREE2"
)

payload = {
    "schema":"elkies-k3.h3-q24-native-a1-pinned-r17-pullback.v1",
    "status":"PASS_Q24_NATIVE_A1_PINNED_R17_PULLBACK",
    "classification":classification,
    "native_a1":{
        "frame":str(a1_path.relative_to(ROOT)),
        "root_data":[1,2,2],
        "mw_rank":16,
    },
    "pinned_r17":{
        "frame":str(PINNED_PATH.relative_to(ROOT)),
        "root_data":[0,0,1],
        "mw_rank":17,
    },
    "direct_pulled_fibre":{
        "a":int(a),
        "b":int(b),
        "q":int(q),
        "old_fibre_degree":int(degree),
        "witness":list(map(int,w)),
        "fiber":list(map(int,D)),
        "mate":list(map(int,M)),
        "witness_norm":int(w*A1*w),
        "primitive":True,
    },
    "transports":{
        "h3_to_native_a1_basis":rows(T),
        "pinned_r17_to_h3_basis":rows(R),
        "native_a1_basis_in_pinned_r17":rows(C),
        "pinned_r17_basis_in_native_a1":rows(P),
        "determinant":int(P.det()),
    },
    "proof_boundary":(
        "Exact integral NS transport. The pinned recovered R17 fibre is pulled "
        "through the certified pinned-R17->H3 map and the freshly replayed "
        "q24-native H3->A1 transports. Thus the direct fibre class, its "
        "(a,b), q, and old-fibre degree are determined exactly rather than "
        "found by bounded shell search. The transform P itself splits the "
        "native A1 NS as U plus the pinned R17 positive frame. This is a "
        "lattice/marking certificate, not characteristic-zero equation execution."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24R17_PULLBACK_RESULT|"
    f"classification={classification}|q={q}|degree={degree}|"
    "pinned_R17=1|status=PASS_Q24_NATIVE_A1_PINNED_R17_PULLBACK",
    flush=True,
)

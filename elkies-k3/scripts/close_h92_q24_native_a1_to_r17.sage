#!/usr/bin/env sage -python
"""
Close the final q24-native H3 lattice step

    native A1/MW16 --q6, old-fibre degree 2--> rootless/MW17 = pinned R17.

Crucially, the parent is the A1 frame actually produced by
run_h92_q24_native_suffix_to_a1.py.  No historical A1 frame is substituted.

The script:
  1. runs a fresh targeted q6 rootless search from that native A1 frame;
  2. verifies the exact primitive isotropic witness and U-neighbour;
  3. verifies the child is rootless of determinant 948;
  4. finds and checks an integral determinant-one isometry to
     data/lattice/rank17_gram.txt;
  5. exports the full 19-dimensional native-A1 -> pinned-R17 transport.
"""

import json
import math
import subprocess
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, pari, vector
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3" / "q24-native-suffix"
SEARCH = SCRIPTS / "search_root_adapted_weyl_neighbors_targeted.sage"
SUFFIX = LOCAL / "q24-native-d12-to-a1.json"
PINNED_PATH = ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt"

OUT = LOCAL / "q24-native-a1-to-r17.json"
SEARCH_OUT = LOCAL / "step10-rootless-q6-search.json"
FRAMES = LOCAL / "step10-rootless-q6-frames"
MW_CACHE = LOCAL / "step10-rootless-q6-mw-vectors-cap10000.json"

U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(v) for v in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(M):
    return [[int(v) for v in row] for row in M.rows()]


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], cwd=ROOT, check=True)


def normalize_qfisom(child, pinned):
    """
    PARI/Sage wrappers have historically exposed qfisom with slightly
    different orientation conventions.  Accept only an explicitly checked
    integral map Q satisfying Q*child*Q^t = pinned.
    """
    raw = pari(child).qfisom(pari(pinned))
    if str(raw) == "0":
        return None

    C = matrix(ZZ, raw)
    candidates = [C, C.transpose()]

    try:
        Ci = C.inverse()
        if Ci.change_ring(ZZ) == Ci:
            Ci = Ci.change_ring(ZZ)
            candidates.extend([Ci, Ci.transpose()])
    except Exception:
        pass

    seen = set()
    for Q in candidates:
        key = tuple(Q.list())
        if key in seen:
            continue
        seen.add(key)
        if Q.dimensions() != (17, 17):
            continue
        if Q * child * Q.transpose() == pinned:
            assert abs(ZZ(Q.det())) == 1
            return Q
    return None


for path in (SEARCH, SUFFIX, PINNED_PATH):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

suffix = json.loads(SUFFIX.read_text())
assert suffix["status"] == "PASS_Q24_NATIVE_D12_TO_A1"
assert suffix["final_ade"] == "A1"
assert int(suffix["final_mw_rank"]) == 16

a1_path = ROOT / suffix["final_frame"]
if not a1_path.exists():
    raise SystemExit(f"native A1 frame missing: {a1_path}")

A1 = load_gram(a1_path)
assert A1.dimensions() == (17, 17)
assert A1.det() == 948

# Parent root data must literally be A1.
qf = pari(A1).qfminim(2)
root_count = ZZ(qf[0])
root_rank = matrix(ZZ, qf[2]).rank() if root_count else 0
assert (root_rank, root_count) == (1, 2)

print(
    "Q24NATIVE_R17_PARENT|"
    f"frame={a1_path.relative_to(ROOT)}|root_data=1,2,2|MW=16|"
    "status=PASS_NATIVE_A1",
    flush=True,
)

# The historical final search needed only a bounded MW-vector storage sample
# to find a q6 rootless existence witness.  This is sufficient for an exact
# positive certificate: once found, every lattice identity below is checked
# without relying on search exhaustiveness.
cmd = [
    "sage", "-python", SEARCH,
    "--frame", a1_path,
    "--root-rank", "1",
    "--q", "6",
    "--degree", "2",
    "--adapt-mw-at-least", "17",
    "--rank-growth-only",
    "--stop-after-first-growth",
    "--stream-first-growth",
    "--target-root-rank", "0",
    "--target-root-count", "0",
    "--target-root-determinant", "1",
    "--mw-vector-cap", "10000",
    "--mw-vectors-cache", MW_CACHE,
    "--output", SEARCH_OUT,
    "--frames-dir", FRAMES,
    "--stream-progress-every", "250",
]
run(cmd)

data = json.loads(SEARCH_OUT.read_text())
assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert data["input_root_data"] == [1, 2, 2]
assert int(data["input_mw_rank"]) == 16

hits = [
    r for r in data.get("neighbors", [])
    if int(r["q"]) == 6
    and int(r["old_fiber_degree"]) == 2
    and tuple(r["child_root_data"]) == (0, 0, 1)
    and r["child_ade"] == "rootless"
    and int(r["child_mw_rank"]) == 17
]

if not hits:
    summary = data.get("summaries", [{}])[-1]
    print(
        "Q24NATIVE_R17_SEARCH|"
        f"screened={summary.get('screened_orbits')}|"
        f"mw_vectors={summary.get('mw_projection_representatives')}|"
        "status=NO_ROOTLESS_IN_CAP10000",
        flush=True,
    )
    raise SystemExit(
        "No rootless q6 hit in the 10k MW-vector sample. "
        "This is not a mathematical obstruction; rerun with a larger cap."
    )

rec = hits[0]
w = vector(ZZ, rec["witness"])
fiber = vector(ZZ, rec["fiber"])

assert w * A1 * w == 12
assert fiber == vector(ZZ, [3, 2] + list(w))
assert gcd(tuple(fiber)) == 1

NS_A1 = block_diagonal_matrix(U2, -A1)
assert fiber * NS_A1 * fiber == 0
old_fiber = vector(ZZ, [1, 0] + [0] * 17)
assert fiber * NS_A1 * old_fiber == 2

raw_child = matrix(ZZ, rec["child_frame"])
adapt = matrix(ZZ, rec["child_root_adapted_basis"])
child = matrix(ZZ, rec["child_root_adapted_frame"])
neighbor_basis = matrix(ZZ, rec["neighbor_basis"])

assert raw_child.dimensions() == child.dimensions() == (17, 17)
assert raw_child.det() == child.det() == 948
assert adapt * raw_child * adapt.transpose() == child
assert abs(ZZ(adapt.det())) == 1
assert abs(ZZ(neighbor_basis.det())) == 1

# Rootless must be checked independently, not trusted from the JSON label.
cqf = pari(child).qfminim(2)
assert ZZ(cqf[0]) == 0

transition = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), adapt)
    * neighbor_basis
)
assert abs(ZZ(transition.det())) == 1

NS_child = block_diagonal_matrix(U2, -child)
assert transition * NS_A1 * transition.transpose() == NS_child

print(
    "Q24NATIVE_R17_Q6|"
    f"orbit={rec['orbit_index']}|q=6|ab=3,2|"
    f"witness_norm={w*A1*w}|primitive=1|old_fibre_degree=2|"
    "root_data=0,0,1|MW=17|status=PASS_NATIVE_ROOTLESS_Q6",
    flush=True,
)

pinned = load_gram(PINNED_PATH)
assert pinned.dimensions() == (17, 17)
assert pinned.det() == 948

Q = normalize_qfisom(child, pinned)
if Q is None:
    raise SystemExit(
        "Fresh native q6 child is rootless/det948, but qfisom did not "
        "produce an explicitly verified integral isometry to pinned R17."
    )

assert Q * child * Q.transpose() == pinned
assert abs(ZZ(Q.det())) == 1

endpoint_map = block_diagonal_matrix(identity_matrix(ZZ, 2), Q)
NS_pinned = block_diagonal_matrix(U2, -pinned)
assert endpoint_map * NS_child * endpoint_map.transpose() == NS_pinned

native_a1_to_pinned = endpoint_map * transition
assert abs(ZZ(native_a1_to_pinned.det())) == 1
assert (
    native_a1_to_pinned
    * NS_A1
    * native_a1_to_pinned.transpose()
    == NS_pinned
)

print(
    "Q24NATIVE_R17_ISOMETRY|"
    f"det={Q.det()}|"
    "Q_child_Qt_equals_pinned=1|"
    f"ns_transport_det={native_a1_to_pinned.det()}|"
    "status=PASS_PINNED_R17",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-q24-native-a1-to-pinned-r17.v1",
    "status": "PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    "parent": {
        "frame": str(a1_path.relative_to(ROOT)),
        "root_data": [1, 2, 2],
        "ade": "A1",
        "mw_rank": 16,
    },
    "neighbor": {
        "q": 6,
        "factor_order": [3, 2],
        "old_fibre_degree": 2,
        "search_is_exhaustive": False,
        "mw_vector_cap": 10000,
        "discovered_orbit_index": int(rec["orbit_index"]),
        "mw_projection": rec["mw_projection"],
        "dominant_labels": rec["dominant_labels"],
        "witness": rec["witness"],
        "fiber": rec["fiber"],
        "child_root_data": [0, 0, 1],
        "child_ade": "rootless",
        "child_mw_rank": 17,
        "search_artifact": str(SEARCH_OUT.relative_to(ROOT)),
    },
    "rootless_child_frame": rows(child),
    "pinned_rank17_frame": str(PINNED_PATH.relative_to(ROOT)),
    "rootless_child_to_pinned_r17_isometry": rows(Q),
    "native_a1_to_rootless_ns_transport": rows(transition),
    "native_a1_to_pinned_r17_ns_transport": rows(native_a1_to_pinned),
    "proof_boundary": (
        "Exact positive lattice/NS certificate from the q24-native A1 frame "
        "to the pinned recovered R17 frame. The bounded search is used only "
        "to discover a witness; the selected q6 neighbour, rootlessness, "
        "primitivity, determinant-one transports, and endpoint isometry are "
        "all checked exactly. This does not execute the characteristic-zero "
        "A1->R17 Weierstrass pencil."
    ),
}

OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24NATIVE_R17_RESULT|"
    f"orbit={rec['orbit_index']}|final=rootless/MW17|"
    "pinned_R17=1|status=PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    flush=True,
)

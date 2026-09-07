#!/usr/bin/env sage -python
"""Rank all q4 exits from the actual equation-marked orbit12 A5+A5 child."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-o12-q4-degree2-all.json"
O12_CERT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-marked-target-neighbor-ranking.json"
INPUTS = (NEIGHBORS, O12_CERT, CROSSOVERS, MANIFEST, FINGERPRINT)
TARGETS = ("pinned_R17", "q25_mw7", "q25_mw4", "mw3_a5_d4_2a2_a1", "mw2_e6_d4_2a2_a1")
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


neighbors = json.loads(NEIGHBORS.read_text())
o12 = json.loads(O12_CERT.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert o12["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"

parent_frame_path = ROOT / neighbors["frame"]
parent = load_matrix(parent_frame_path)
ns = block_diagonal_matrix(U2, -parent)
equation_a11_to_parent = matrix(ZZ, o12["transport"]["parent_to_child_basis"])
equation_in_parent = equation_a11_to_parent.inverse().change_ring(ZZ)
target_in_equation = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"]
    if item["state"] == "equation_A11"
}
target_in_parent = {
    name: fibre * equation_in_parent for name, fibre in target_in_equation.items()
}

records = []
for raw in neighbors["neighbors"]:
    fibre = vector(ZZ, raw["fiber"])
    marked = {
        target: int(target_in_parent[target] * ns * fibre)
        for target in TARGETS
    }
    assert all(value >= 0 for value in marked.values())
    records.append({
        "candidate_id": {
            "q": int(raw["q"]),
            "old_fibre_degree": int(raw["old_fiber_degree"]),
            "orbit_index": int(raw["orbit_index"]),
        },
        "child": {
            "ade": raw["child_ade"],
            "mw_rank": int(raw["child_mw_rank"]),
            "root_data": raw["child_root_data"],
        },
        "marked_target_degrees": marked,
        "coordinate_growth_max": int(max(abs(value) for value in fibre)),
        "dominant_labels": raw["dominant_labels"],
        "mw_projection": raw["mw_projection"],
    })

rankings = {
    target: sorted(
        records,
        key=lambda item: (
            item["marked_target_degrees"][target],
            item["coordinate_growth_max"],
            item["candidate_id"]["orbit_index"],
        ),
    )[:100]
    for target in TARGETS
}

# Locate the currently certified historical q4 edge in this equation-marked
# parent by transporting its class and reducing to the dominant Weyl chamber.
historical_a11_in_equation_frame = matrix(
    ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]
)
historical_a11_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2), historical_a11_in_equation_frame
)
a11_to_historical_child = matrix(ZZ, manifest["forward_steps"][2]["transition"])
historical_child_in_equation = a11_to_historical_child * historical_a11_in_equation
historical_child_in_parent = historical_child_in_equation * equation_in_parent
historical_zero_in_equation = vector(
    ZZ, historical_child_in_equation.row(1) - historical_child_in_equation.row(0)
)
explicit_degree_one_curves = {
    "old_A11_component_9": vector(
        ZZ, [0, 0] + [-ZZ(index == 9) for index in range(17)]
    ),
    "old_A11_affine": vector(ZZ, [1, 0] + [1] * 11 + [0] * 6),
}
historical_zero_matches = [
    name for name, curve in explicit_degree_one_curves.items()
    if curve == historical_zero_in_equation
]
current_route_fibre_historical = vector(
    ZZ, manifest["forward_steps"][3]["new_fibre_in_parent"]
)
current_route_fibre_parent = current_route_fibre_historical * historical_child_in_parent
assert current_route_fibre_parent * ns * current_route_fibre_parent == 0

# Sage may return the same mutable vector when the base ring is unchanged.
dominant = vector(ZZ, list(current_route_fibre_parent))
for unused in range(10000):
    labels = list(dominant[2:] * parent[:, :10])
    negative = [index for index, value in enumerate(labels) if value < 0]
    if not negative:
        break
    index = negative[0]
    dominant[2 + index] -= labels[index]
else:
    raise RuntimeError("Weyl chamber reduction did not terminate")
current_route_marked_degrees = {
    target: int(target_in_parent[target] * ns * current_route_fibre_parent)
    for target in TARGETS
}
assert all(value >= 0 for value in current_route_marked_degrees.values())
current_route_equation_zero_q = int(
    current_route_fibre_parent[0] * current_route_fibre_parent[1]
)

payload = {
    "schema": "elkies-k3.h3-a5a5-marked-target-neighbor-ranking.v1",
    "status": "PASS_EXACT_A5A5_MARKED_TARGET_NEIGHBOR_RANKING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(parent_frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (parent_frame_path,)
        },
    },
    "candidate_count": len(records),
    "target_fibres_in_parent": {
        name: [int(value) for value in fibre] for name, fibre in target_in_parent.items()
    },
    "rankings_top_100": rankings,
    "current_route": {
        "historical_orbit": 472,
        "historical_q": 4,
        "equation_zero_q": current_route_equation_zero_q,
        "historical_zero_in_equation_A11": [int(value) for value in historical_zero_in_equation],
        "historical_zero_explicit_curve_matches": historical_zero_matches,
        "fibre_in_equation_marked_parent": [int(value) for value in current_route_fibre_parent],
        "dominant_fibre": [int(value) for value in dominant],
        "dominant_labels": [int(value) for value in dominant[2:] * parent[:, :10]],
        "marked_target_degrees": current_route_marked_degrees,
        "child": {"ade": "A3+A3+A3", "mw_rank": 8, "root_data": [9, 36, 64]},
    },
    "proof_boundary": (
        "Exact exhaustive q4 marked-degree ranking and exact identification of "
        "the current historical edge in the equation child chamber. Full section-"
        "wall nefness and equation-cost scoring remain separate certification gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in TARGETS:
    item = rankings[target][0]
    print(
        "A5A5TARGET|target={}|orbit={}|degree={}|child={}/MW{}|growth={}".format(
            target, item["candidate_id"]["orbit_index"],
            item["marked_target_degrees"][target], item["child"]["ade"],
            item["child"]["mw_rank"], item["coordinate_growth_max"],
        ),
        flush=True,
    )
print(
    "A5A5TARGET|current_historical_orbit=472|historical_q=4|equation_zero_q={}|pinned_degree={}".format(
        current_route_equation_zero_q,
        current_route_marked_degrees["pinned_R17"],
    ),
    flush=True,
)
print(f"A5A5TARGET|status={payload['status']}|output={OUTPUT.resolve()}", flush=True)

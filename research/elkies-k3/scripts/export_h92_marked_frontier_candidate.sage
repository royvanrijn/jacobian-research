#!/usr/bin/env sage -python
"""Promote one exact ranked neighbour to a reusable fully marked beam state."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--frontier", type=Path, required=True)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--q", type=int, required=True)
parser.add_argument("--degree", type=int, required=True)
parser.add_argument("--orbit", type=int, required=True)
parser.add_argument("--edge-operational-score", type=int, required=True)
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

frontier_path = args.frontier.resolve()
source_path = args.source.resolve()
frame_output = args.frame_output.resolve()
output = args.output.resolve()
frontier = json.loads(frontier_path.read_text())
source = json.loads(source_path.read_text())
assert frontier["status"] == "PASS_EXACT_MARKED_ROOT_ADAPTED_FRONTIER_RANKING"
assert source["status"] in {
    "PASS_EXACT_A11_EQUATION_MARKING",
    "PASS_EXACT_D13_ZERO_CHILD_MARKING",
    "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
}
selected = next(
    item for item in frontier["ranked_candidates"]
    if item["candidate_id"] == {
        "q": args.q,
        "old_fibre_degree": args.degree,
        "orbit_index": args.orbit,
    }
)
raw = selected["source_neighbor_record"]
raw_fibre = vector(ZZ, raw["fiber"])
child = matrix(ZZ, raw["child_root_adapted_frame"])
neighbor = matrix(ZZ, raw["neighbor_basis"])
adaptation = matrix(ZZ, raw["child_root_adapted_basis"])
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * neighbor
inverse = transition.inverse().change_ring(ZZ)

source_frame_path = ROOT / source["frame_output"]
source_frame = matrix(
    ZZ,
    [[ZZ(x) for x in line.split()] for line in source_frame_path.read_text().splitlines()
     if line.strip() and not line.lstrip().startswith("#")],
)
g_source = block_diagonal_matrix(U2, -source_frame)
g_child = block_diagonal_matrix(U2, -child)
assert abs(transition.det()) == 1
assert transition * g_source * transition.transpose() == g_child
target_key = (
    "target_fibres_in_root_adapted_hub"
    if "target_fibres_in_root_adapted_hub" in source
    else "target_fibres_in_child"
)
targets = {
    name: vector(ZZ, value) * inverse
    for name, value in source[target_key].items()
}
assert all(value * g_child * value == 0 for value in targets.values())
explicit_curves = {
    name: vector(ZZ, value) * inverse
    for name, value in source.get("equation_explicit_curves_in_child", {}).items()
}
assert all(value * g_child * value == -2 for value in explicit_curves.values())

frame_output.write_text("\n".join(" ".join(map(str, row)) for row in child.rows()) + "\n")
inputs = (frontier_path, source_path, source_frame_path)
payload = {
    "schema": "elkies-k3.h3-marked-beam-state.v1",
    "status": "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "hub": "{}_q{}_orbit{}".format(
        source.get("hub", source.get("source_hub", "state")), args.q, args.orbit
    ),
    "candidate_id": selected["candidate_id"],
    "child": selected["child"],
    "root_data": selected["child"]["root_data"],
    "frame_output": str(frame_output.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_output.read_bytes()).hexdigest(),
    "source_to_root_adapted_hub_basis": [[int(x) for x in row] for row in transition.rows()],
    "root_adapted_hub_to_source_basis": [[int(x) for x in row] for row in inverse.rows()],
    "prefix_operational_score": int(
        source.get("prefix_operational_score", 0) + args.edge_operational_score
    ),
    "edge_operational_score": args.edge_operational_score,
    "edge_nef_audit": {
        "component_pairings": selected["component_pairings"],
        "affine_pairings": selected["affine_pairings"],
        "minimum_section_intersection": selected["minimum_section_intersection"],
        # Older exact frontier artifacts predate the redundant P.O field.
        # In the marked U convention O=e-f, so P.O=P[0]-P[1].
        "P_dot_O": int(selected.get("P_dot_O", raw_fibre[0] - raw_fibre[1])),
    },
    "target_fibres_in_root_adapted_hub": {
        name: [int(x) for x in value] for name, value in targets.items()
    },
    "equation_explicit_curves_in_child": {
        name: [int(x) for x in value] for name, value in explicit_curves.items()
    },
    "proof_boundary": (
        "Exact primitive nef neighbour, marked U, child root data, and determinant-one "
        "bidirectional source transport. The supplied operational score is a planning "
        "weight and must be recomputed by a final route certifier."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "BEAMSTATE|hub={}|root={}|prefix_operational={}|det={}|status={}|output={}".format(
        payload["hub"], ",".join(map(str, payload["root_data"])),
        payload["prefix_operational_score"], int(transition.det()), payload["status"], output,
    ), flush=True,
)

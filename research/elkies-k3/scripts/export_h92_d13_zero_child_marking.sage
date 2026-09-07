#!/usr/bin/env sage -python
"""Export a selected explicit-zero D13 neighbour as a fully marked state."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--search", type=Path, required=True)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--q", type=int, required=True)
parser.add_argument("--degree", type=int, required=True)
parser.add_argument("--orbit", type=int, required=True)
parser.add_argument("--zero", required=True)
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

search_path = args.search.resolve()
source_path = args.source.resolve()
frame_output = args.frame_output.resolve()
output = args.output.resolve()
search = json.loads(search_path.read_text())
source = json.loads(source_path.read_text())
assert search["status"] == "PASS_EXACT_D13_ZERO_CHANGING_D12_PRESENTATION_SEARCH"
assert source["status"] in {
    "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "PASS_EXACT_D13_ZERO_LOOP_RETURNED_MARKING",
}
selected = next(
    item for item in search["ranked_presentations"]
    if item["first_edge_candidate_id"] == {
        "q": args.q,
        "old_fibre_degree": args.degree,
        "orbit_index": args.orbit,
    } and item["explicit_zero_curve"] == args.zero
)

child = matrix(ZZ, selected["explicit_child_frame"])
transition = matrix(ZZ, selected["source_to_explicit_child_basis"])
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
targets = {
    name: vector(ZZ, value) * inverse
    for name, value in source["target_fibres_in_root_adapted_hub"].items()
}
assert all(value * g_child * value == 0 for value in targets.values())

frame_output.write_text("\n".join(" ".join(map(str, row)) for row in child.rows()) + "\n")
prefix_raw = int(source.get("prefix_raw_score", 0) + selected["first_edge_score"])
prefix_operational = int(
    source.get("prefix_operational_score", 0) + max(500, selected["first_edge_score"])
)
inputs = (search_path, source_path, source_frame_path)
payload = {
    "schema": "elkies-k3.h3-d13-explicit-zero-child-marking.v1",
    "status": "PASS_EXACT_D13_ZERO_CHILD_MARKING",
    "hub": "{}_q{}_orbit{}_{}".format(
        source.get("hub", "d13"), args.q, args.orbit, args.zero
    ),
    "candidate_id": selected["first_edge_candidate_id"],
    "explicit_zero_curve": selected["explicit_zero_curve"],
    "root_data": list(map(int, selected["explicit_child_root_data"])),
    "frame_output": str(frame_output.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_output.read_bytes()).hexdigest(),
    "source_to_root_adapted_hub_basis": [[int(x) for x in row] for row in transition.rows()],
    "root_adapted_hub_to_source_basis": [[int(x) for x in row] for row in inverse.rows()],
    "prefix_raw_score": prefix_raw,
    "prefix_operational_score": prefix_operational,
    "inherited_explicit_curve_degrees": selected["inherited_explicit_curve_degrees"]["first_edge"],
    "inherited_explicit_curve_names": selected["inherited_explicit_curve_degrees"]["names"],
    "target_fibres_in_root_adapted_hub": {
        name: [int(x) for x in value] for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact marked U, explicit zero, nef first edge, and determinant-one "
        "bidirectional source transport. New outgoing edges need separate nef certificates."
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
    "D13ZEROCHILD|hub={}|root={}|prefix_operational={}|det={}|status={}|output={}".format(
        payload["hub"], ",".join(map(str, payload["root_data"])), prefix_operational,
        int(transition.det()), payload["status"], output,
    ), flush=True,
)

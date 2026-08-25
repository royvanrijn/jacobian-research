#!/usr/bin/env sage -python
"""Export an exact D13 state returned by a selected zero-changing loop."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--search", type=Path,
                    default=GENERATED / "elkies-k3-h3-d13-zero-changing-d12-presentations.json")
parser.add_argument("--source", type=Path,
                    default=GENERATED / "elkies-k3-h3-equation-d13-marking.json")
parser.add_argument("--q", type=int, default=4)
parser.add_argument("--degree", type=int, default=2)
parser.add_argument("--orbit", type=int, default=11)
parser.add_argument("--zero", default="old_D13_component_5")
parser.add_argument("--frame-output", type=Path,
                    default=GENERATED / "elkies-k3-h3-d13-q4o11-returned-frame.txt")
parser.add_argument("--output", type=Path,
                    default=GENERATED / "elkies-k3-h3-d13-q4o11-returned-marking.json")
args = parser.parse_args()
SEARCH = args.search.resolve()
SOURCE = args.source.resolve()
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


search = json.loads(SEARCH.read_text())
source = json.loads(SOURCE.read_text())
assert search["status"] == "PASS_EXACT_D13_ZERO_CHANGING_D12_PRESENTATION_SEARCH"
assert source["status"] in {
    "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "PASS_EXACT_D13_ZERO_LOOP_RETURNED_MARKING",
}
selected = next(row for row in search["ranked_presentations"]
                if row["first_edge_candidate_id"] ==
                {"q": args.q, "old_fibre_degree": args.degree, "orbit_index": args.orbit}
                and row["explicit_zero_curve"] == args.zero)
returned = matrix(ZZ, selected["returned_frame"])
transition = matrix(ZZ, selected["source_to_returned_D13_basis"])
inverse = transition.inverse().change_ring(ZZ)
source_frame = load_matrix(ROOT / source["frame_output"])
g_source = block_diagonal_matrix(U2, -source_frame)
g_returned = block_diagonal_matrix(U2, -returned)
assert abs(transition.det()) == 1
assert transition * g_source * transition.transpose() == g_returned

targets = {
    name: vector(ZZ, value) * inverse
    for name, value in source["target_fibres_in_root_adapted_hub"].items()
}
assert targets["current_0_D12"] == vector(ZZ, selected["exit_D12_fibre_in_returned_D13"])
assert all(value * g_returned * value == 0 for value in targets.values())

FRAME_OUTPUT.write_text("\n".join(" ".join(map(str, row)) for row in returned.rows()) + "\n")
inputs = (SEARCH, SOURCE, ROOT / source["frame_output"])
payload = {
    "schema": "elkies-k3.h3-d13-zero-loop-returned-marking.v1",
    "status": "PASS_EXACT_D13_ZERO_LOOP_RETURNED_MARKING",
    "hub": "{}_q{}_orbit{}_{}".format(source.get("hub", "d13"), args.q, args.orbit, args.zero),
    "candidate_id": selected["first_edge_candidate_id"],
    "explicit_zero_curve": selected["explicit_zero_curve"],
    "prefix_raw_score": int(source.get("prefix_raw_score", 0)
                            + selected["first_edge_score"] + selected["return_score"]),
    "prefix_operational_score": int(source.get("prefix_operational_score", 0)
                                    + sum(max(500, int(selected[key]))
                                          for key in ("first_edge_score", "return_score"))),
    "root_data": [13, 312, 4],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_D13_to_root_adapted_hub_basis": rows(transition),
    "root_adapted_hub_to_equation_D13_basis": rows(inverse),
    "target_fibres_in_root_adapted_hub": {
        name: [int(x) for x in value] for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact returned D13 frame, marked U, full bidirectional equation-D13 transport, "
        "and transported route targets. Scores are compiler-planning estimates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("D13RETURNMARK|prefix_raw={}|prefix_operational={}|det={}|status={}|output={}".format(
    payload["prefix_raw_score"], payload["prefix_operational_score"], transition.det(),
    payload["status"], OUTPUT,
), flush=True)

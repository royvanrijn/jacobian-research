#!/usr/bin/env sage -python
"""Export the returned E8+E6 marking of the promoted first-q8 zero loop."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SEARCH = GENERATED / "elkies-k3-h3-first-q8-zero-changing-d13-presentations.json"
SOURCE = GENERATED / "elkies-k3-h3-first-q8-source-marking.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-first-q8-q4o11-c1-returned-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-first-q8-q4o11-c1-returned-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


search = json.loads(SEARCH.read_text())
source = json.loads(SOURCE.read_text())
assert search["status"] == "PASS_EXACT_E8E6_ZERO_CHANGING_D13_PRESENTATION_SEARCH"
assert source["status"] == "PASS_EXACT_FIRST_Q8_SOURCE_MARKING"
selected = next(row for row in search["ranked_presentations"]
                if row["first_edge_candidate_id"] ==
                {"q": 4, "old_fibre_degree": 2, "orbit_index": 11}
                and row["explicit_zero_curve"] == "old_E8E6_component_1")
returned = matrix(ZZ, selected["returned_frame"])
transition = matrix(ZZ, selected["source_to_returned_E8E6_basis"])
inverse = transition.inverse().change_ring(ZZ)
source_frame = load_matrix(ROOT / source["frame_output"])
g_source = block_diagonal_matrix(U2, -source_frame)
g_returned = block_diagonal_matrix(U2, -returned)
assert abs(transition.det()) == 1
assert transition * g_source * transition.transpose() == g_returned

targets = {name: vector(ZZ, value) * inverse
           for name, value in source["target_fibres_in_root_adapted_hub"].items()}
assert targets["equation_D13"] == vector(ZZ, selected["exit_D13_fibre_in_returned_E8E6"])
assert all(value * g_returned * value == 0 for value in targets.values())

FRAME_OUTPUT.write_text("\n".join(" ".join(map(str, row)) for row in returned.rows()) + "\n")
inputs = (SEARCH, SOURCE, ROOT / source["frame_output"])
payload = {
    "schema": "elkies-k3.h3-first-q8-zero-loop-returned-marking.v1",
    "status": "PASS_EXACT_FIRST_Q8_ZERO_LOOP_RETURNED_MARKING",
    "hub": "first_q8_E8_plus_E6_q4_orbit11_old_component_1",
    "candidate_id": selected["first_edge_candidate_id"],
    "explicit_zero_curve": selected["explicit_zero_curve"],
    "prefix_raw_score": int(selected["first_edge_score"] + selected["return_score"]),
    "prefix_operational_score": int(sum(max(500, int(selected[key]))
                                          for key in ("first_edge_score", "return_score"))),
    "root_data": [14, 312, 3],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_E8E6_to_root_adapted_hub_basis": rows(transition),
    "root_adapted_hub_to_equation_E8E6_basis": rows(inverse),
    "target_fibres_in_root_adapted_hub": {
        name: [int(x) for x in value] for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact returned E8+E6 frame, marked U, full bidirectional source transport, "
        "and transported equation-D13 target. Scores are compiler-planning estimates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("FIRSTQ8RETURNMARK|prefix_raw={}|prefix_operational={}|det={}|status={}|output={}".format(
    payload["prefix_raw_score"], payload["prefix_operational_score"], transition.det(),
    payload["status"], OUTPUT,
), flush=True)

#!/usr/bin/env sage -python
"""Export a returned 2A5 zero state from the exact zero-loop search."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
LOOPS = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
A11_MARKING = GENERATED / "elkies-k3-h3-current_A11-marked-frame.json"
SOURCE_ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q", type=int, required=True)
parser.add_argument("--orbit", type=int, required=True)
parser.add_argument("--zero", required=True)
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


def entries(value):
    return [int(x) for x in vector(ZZ, value)]


loops = json.loads(LOOPS.read_text())
a11 = json.loads(A11_MARKING.read_text())
matches = [
    row for row in loops["ranked_loops"]
    if row["first_edge_candidate_id"] == {
        "q": args.q, "old_fibre_degree": 2, "orbit_index": args.orbit,
    } and row["explicit_zero_curve"] == args.zero
]
assert len(matches) == 1
selected = matches[0]
returned_frame = matrix(ZZ, selected["returned_A5A5_frame"])
equation_to_explicit = matrix(ZZ, selected["equation_A11_to_explicit_child_basis"])
explicit_to_returned = matrix(ZZ, selected["return_transition"])
equation_to_returned = explicit_to_returned * equation_to_explicit
returned_to_equation = equation_to_returned.inverse().change_ring(ZZ)
assert abs(equation_to_returned.det()) == 1

source_zero = json.loads(SOURCE_ZERO.read_text())["selected"]
source_zero_frame = matrix(ZZ, source_zero["frame"])
g_source_zero = block_diagonal_matrix(U2, -source_zero_frame)
equation_to_source_zero = matrix(ZZ, source_zero["equation_A11_to_explicit_zero_basis"])
source_zero_to_equation = equation_to_source_zero.inverse().change_ring(ZZ)
g_equation = source_zero_to_equation * g_source_zero * source_zero_to_equation.transpose()
g_returned = block_diagonal_matrix(U2, -returned_frame)
assert equation_to_returned * g_equation * equation_to_returned.transpose() == g_returned

manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
cumulative = identity_matrix(ZZ, 19)
targets_equation = {"current_A5A5": vector(ZZ, historical_in_equation.row(0))}
for index, step in enumerate(manifest["forward_steps"]):
    if index < 2:
        continue
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    if step["child"] == "2A5/MW7":
        targets_equation["current_A5A5"] = vector(ZZ, (cumulative * historical_in_equation).row(0))
    elif step["child"] == "3A3/MW8":
        targets_equation["current_3A3"] = vector(ZZ, (cumulative * historical_in_equation).row(0))
        break
crossovers = json.loads(CROSSOVERS.read_text())
targets_equation["pinned_R17"] = vector(ZZ, next(
    item["target_fibre_in_state"] for item in crossovers["records"]
    if item["state"] == "equation_A11" and item["target"] == "pinned_R17"
))
targets = {name: value * returned_to_equation for name, value in targets_equation.items()}
assert targets["current_A5A5"] == vector(ZZ, [1, 0] + [0] * 17)
assert targets["current_3A3"] == vector(ZZ, selected["exit_3A3_fibre_in_returned_A5A5"])

FRAME_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
FRAME_OUTPUT.write_text("\n".join(" ".join(map(str, row)) for row in returned_frame.rows()) + "\n")
inputs = (LOOPS, A11_MARKING, SOURCE_ZERO, MANIFEST, FINGERPRINT, CROSSOVERS)
payload = {
    "schema": "elkies-k3.h3-a5a5-zero-loop-returned-marking.v1",
    "status": "PASS_EXACT_A5A5_ZERO_LOOP_RETURNED_MARKING",
    "hub": f"a5a5_q{args.q}_orbit{args.orbit}_{args.zero}",
    "candidate_id": selected["first_edge_candidate_id"],
    "explicit_zero_curve": args.zero,
    "prefix_equation_cost_score": int(
        selected["first_edge_equation_cost_score"] + selected["return_equation_cost_score"]
    ),
    "root_data": [10, 60, 36],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": rows(equation_to_returned),
    "root_adapted_hub_to_equation_A11_basis": rows(returned_to_equation),
    "target_fibres_in_root_adapted_hub": {name: entries(value) for name, value in targets.items()},
    "proof_boundary": (
        "The returned 2A5 frame, full equation-A11 transport, and all marked suffix fibres are exact. "
        "The equation-cost prefix remains a planning estimate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in inputs},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("A5RETURNMARK|q={}|orbit={}|zero={}|prefix={}|det={}|status={}|output={}".format(
    args.q, args.orbit, args.zero, payload["prefix_equation_cost_score"],
    equation_to_returned.det(), payload["status"], OUTPUT,
), flush=True)

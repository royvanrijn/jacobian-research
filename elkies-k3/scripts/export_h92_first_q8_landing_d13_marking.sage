#!/usr/bin/env sage -python
"""Export the changed-zero D13 landing of the promoted first-q8 detour."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CERTIFICATE = GENERATED / "elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json"
SEARCH = GENERATED / "elkies-k3-h3-first-q8-zero-changing-d13-presentations.json"
EQUATION_MARKING = GENERATED / "elkies-k3-h3-equation-d13-marking.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-first-q8-q4o11-landing-d13-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-first-q8-q4o11-landing-d13-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


certificate = json.loads(CERTIFICATE.read_text())
search = json.loads(SEARCH.read_text())
equation = json.loads(EQUATION_MARKING.read_text())
assert certificate["status"] == "PASS_EXACT_PROMOTED_FIRST_Q8_TRIPLE_Q4_ROUTE_TO_PINNED_R17"
assert equation["status"] == "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
selected = search["ranked_presentations"][0]
landing_frame = matrix(ZZ, selected["exit_child_frame"])
g_landing = block_diagonal_matrix(U2, -landing_frame)
equation_in_landing = matrix(
    ZZ, certificate["equation_D13_identification"]["canonical_equation_D13_basis_in_landing"]
)
landing_in_equation = equation_in_landing.inverse().change_ring(ZZ)
assert abs(equation_in_landing.det()) == 1

targets = {
    name: vector(ZZ, value) * equation_in_landing
    for name, value in equation["target_fibres_in_root_adapted_hub"].items()
}
assert targets["equation_D13"] == vector(ZZ, [1, 0] + [0] * 17)
assert all(value * g_landing * value == 0 for value in targets.values())

FRAME_OUTPUT.write_text("\n".join(" ".join(map(str, row)) for row in landing_frame.rows()) + "\n")
inputs = (CERTIFICATE, SEARCH, EQUATION_MARKING)
payload = {
    "schema": "elkies-k3.h3-first-q8-landing-d13-marking.v1",
    "status": "PASS_EXACT_FIRST_Q8_LANDING_D13_MARKING",
    "hub": "first_q8_q4o11_changed_zero_landing_D13",
    "root_data": [13, 312, 4],
    "prefix_raw_score": int(selected["total_equation_cost_score"]),
    "prefix_operational_score": int(
        selected["inherited_explicit_equation_cost"]["operational_total_score"]
    ),
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_D13_to_root_adapted_hub_basis": rows(landing_in_equation),
    "root_adapted_hub_to_equation_D13_basis": rows(equation_in_landing),
    "target_fibres_in_root_adapted_hub": {
        name: [int(x) for x in value] for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact changed-zero D13 landing, marked U, full bidirectional current-equation-D13 "
        "transport, and all transported current-suffix and pinned targets."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("FIRSTQ8LANDINGD13|prefix={}|d12_degree={}|det={}|status={}|output={}".format(
    payload["prefix_operational_score"], targets["current_0_D12"][1],
    landing_in_equation.det(), payload["status"], OUTPUT,
), flush=True)

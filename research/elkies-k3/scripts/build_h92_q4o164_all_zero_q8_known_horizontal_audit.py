#!/usr/bin/env python3
"""Aggregate the bounded q8 known-horizontal scans for all effective q4/o164 zeros.

status: ACTIVE_SEARCH
claim: no known-horizontal hit in the recorded cap-10000 q8 shells
inputs: eight effective-zero markings, neighbor shells, and equation-cost filters
outputs: generated all-zero bounded negative audit
"""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
OUTPUT = GEN / "elkies-k3-h3-q4o164-all-effective-zeros-q8-known-horizontal-audit.json"
ZEROS = {
    "first_I6_affine_component": "first-affine",
    "old_A11_component_1": "old_a11_component_1",
    "old_A11_component_2": "old_a11_component_2",
    "old_A11_component_3": "old_a11_component_3",
    "old_A11_component_4": "old_a11_component_4",
    "old_A11_component_7": "old_a11_component_7",
    "old_A11_component_8": "c8",
    "old_A11_component_10": "old_a11_component_10",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


records = []
inputs = []
for zero, tag in ZEROS.items():
    if tag == "c8":
        neighbor = GEN / "elkies-k3-h3-q4o164-c8-q8d2-cap10000-growth-neighbors.json"
        score = GEN / "elkies-k3-h3-q4o164-c8-q8d2-cap10000-known-horizontal-equation-cost.json"
    else:
        neighbor = GEN / f"elkies-k3-h3-q4o164-{tag}-q8d2-cap10000-neighbors.json"
        score = GEN / f"elkies-k3-h3-q4o164-{tag}-q8d2-cap10000-known-horizontal-equation-cost.json"
    ndata = json.loads(neighbor.read_text())
    sdata = json.loads(score.read_text())
    assert ndata["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert sdata["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
    assert sdata["candidate_count"] == sdata["retained_count"] == 0
    summary = next(item for item in ndata["summaries"] if item["q"] == 8)
    assert summary["mw_vector_cap"] == summary["mw_projection_representatives"] == 10000
    records.append({
        "zero": zero,
        "neighbor_artifact": str(neighbor.relative_to(ROOT)),
        "primitive_q8_candidates": int(summary["primitive_neighbors"]),
        "known_horizontal_candidates": 0,
        "score_artifact": str(score.relative_to(ROOT)),
    })
    inputs.extend((neighbor, score))

payload = {
    "schema": "elkies-k3.h3-q4o164-all-zero-q8-known-horizontal-audit.v1",
    "status": "PASS_BOUNDED_Q4O164_ALL_ZERO_Q8_NO_KNOWN_HORIZONTAL",
    "effective_zero_count": len(records),
    "q": 8,
    "old_fibre_degree": 2,
    "mw_vector_cap_per_zero": 10000,
    "records": records,
    "total_primitive_candidates_scored": sum(
        item["primitive_q8_candidates"] for item in records
    ),
    "total_known_horizontal_candidates": 0,
    "construction_consequence": (
        "Changing among the eight certified effective origins does not expose a q8 "
        "horizontal in the already explicit section subgroup within any recorded shell. "
        "Keep q8/o376 and construct a new rational direction; the current preferred "
        "no-Groebner carrier is the exact inherited degree-seven P1 divisor followed by "
        "fibrewise Abel reduction."
    ),
    "proof_boundary": (
        "Each individual subgroup membership test is exact for its stored neighbor shell, "
        "but the Mordell-Weil enumeration was capped at 10,000 representatives per zero. "
        "This is a bounded changed-origin obstruction, not a global non-existence theorem "
        "for q8 neighbors or rational construction words."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ALLZEROQ8|zeros={}|primitive={}|known=0|status={}|output={}".format(
        len(records), payload["total_primitive_candidates_scored"],
        payload["status"], OUTPUT,
    )
)

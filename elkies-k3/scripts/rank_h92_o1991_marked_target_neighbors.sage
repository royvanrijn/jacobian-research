#!/usr/bin/env sage -python
"""Rank the complete orbit1991 shell by exact marked R17-hub degrees."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-o1991-explicit-zero-degree2-neighbors.json"
SCORES = GENERATED / "elkies-k3-h3-o1991-explicit-zero-equation-cost-neighbors.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
OUTPUT = GENERATED / "elkies-k3-h3-o1991-marked-target-neighbor-ranking.json"
INPUTS = (NEIGHBORS, SCORES, CROSSOVERS)
TARGETS = ("pinned_R17", "q25_mw7", "q25_mw4", "mw3_a5_d4_2a2_a1", "mw2_e6_d4_2a2_a1")


def candidate_id(item):
    return (int(item["q"]), int(item["old_fiber_degree"]), int(item["orbit_index"]))


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


neighbors = json.loads(NEIGHBORS.read_text())
scores = json.loads(SCORES.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert scores["status"] == "PASS_EXACT_O1991_EXPLICIT_ZERO_EQUATION_COST_SCORING"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"

score_by_id = {
    (
        int(item["candidate_id"]["q"]),
        int(item["candidate_id"]["old_fibre_degree"]),
        int(item["candidate_id"]["orbit_index"]),
    ): item
    for item in scores["retained_candidates"]
}
target_in_parent = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"]
    if item["state"] == "q8_orbit1991_explicit_zero"
}
assert set(target_in_parent) == set(TARGETS)

records = []
not_root_adapted = 0
for raw in neighbors["neighbors"]:
    if "child_root_adapted_basis" not in raw:
        not_root_adapted += 1
        continue
    key = candidate_id(raw)
    transition = block_diagonal_matrix(
        identity_matrix(ZZ, 2), matrix(ZZ, raw["child_root_adapted_basis"])
    ) * matrix(ZZ, raw["neighbor_basis"])
    inverse = transition.inverse().change_ring(ZZ)
    marked = {}
    for target in TARGETS:
        target_fibre = target_in_parent[target] * inverse
        marked[target] = {
            "degree": int(target_fibre[1]),
            "fibre_max_abs_coordinate": int(max(abs(value) for value in target_fibre)),
            "fibre_in_child": entries(target_fibre),
        }
    score = score_by_id.get(key)
    records.append({
        "candidate_id": {"q": key[0], "old_fibre_degree": key[1], "orbit_index": key[2]},
        "child": {
            "ade": raw["child_ade"],
            "mw_rank": int(raw["child_mw_rank"]),
            "root_data": [int(value) for value in raw["child_root_data"]],
        },
        "marked_targets": marked,
        "equation_cost": None if score is None else {
            "score": int(score["equation_cost_score"]),
            "P_dot_O": int(score["horizontal"]["P_dot_O"]),
            "expected_RR_ambient": int(score["expected_RR_ambient"]),
            "explicit_degree_zero_curves": score["explicit_degree_zero_curves"],
            "explicit_degree_one_curves": score["explicit_degree_one_curves"],
            "declared_curve_and_affine_nef_gate": score["declared_curve_and_affine_nef_gate"],
        },
    })

rankings = {
    target: sorted(
        records,
        key=lambda item: (
            item["marked_targets"][target]["degree"],
            item["marked_targets"][target]["fibre_max_abs_coordinate"],
        ),
    )[:100]
    for target in TARGETS
}

with_equation_cost = [item for item in records if item["equation_cost"] is not None]
pareto = []
for item in sorted(
    with_equation_cost,
    key=lambda value: (
        value["marked_targets"]["pinned_R17"]["degree"],
        value["equation_cost"]["score"],
    ),
):
    degree = item["marked_targets"]["pinned_R17"]["degree"]
    cost = item["equation_cost"]["score"]
    if not any(
        other["marked_targets"]["pinned_R17"]["degree"] <= degree
        and other["equation_cost"]["score"] <= cost
        and (
            other["marked_targets"]["pinned_R17"]["degree"] < degree
            or other["equation_cost"]["score"] < cost
        )
        for other in with_equation_cost
    ):
        pareto.append(item)

payload = {
    "schema": "elkies-k3.h3-o1991-marked-target-neighbor-ranking.v1",
    "status": "PASS_EXACT_MARKED_TARGET_NEIGHBOR_RANKING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "candidate_count": len(records),
    "excluded_not_root_adapted_count": not_root_adapted,
    "rankings_top_100": rankings,
    "equation_cost_retained_count": len(with_equation_cost),
    "pinned_degree_equation_cost_pareto_front": pareto,
    "proof_boundary": (
        "All marked degrees and transports used in ranking are exact. Candidate "
        "shell records have component-dominance/root-adaptation checks but require "
        "a separate full section-wall nef certificate before route use."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in TARGETS:
    item = rankings[target][0]
    print(
        "O1991TARGET|target={}|best={},{},{}|degree={}|child={}/MW{}|equation_cost={}".format(
            target,
            item["candidate_id"]["q"], item["candidate_id"]["old_fibre_degree"],
            item["candidate_id"]["orbit_index"], item["marked_targets"][target]["degree"],
            item["child"]["ade"], item["child"]["mw_rank"],
            None if item["equation_cost"] is None else item["equation_cost"]["score"],
        ),
        flush=True,
    )
print(f"O1991TARGET|pareto={len(pareto)}|status={payload['status']}|output={OUTPUT.resolve()}", flush=True)

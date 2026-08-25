#!/usr/bin/env sage -python
"""Rank the exhaustive equation-side A11 q8 shell by marked R17 distances."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
SCORES = GENERATED / "elkies-k3-h3-a11-equation-cost-neighbors-all.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
OUTPUT = GENERATED / "elkies-k3-h3-a11-marked-target-neighbor-ranking.json"
INPUTS = (NEIGHBORS, SCORES, CROSSOVERS)
TARGETS = (
    "pinned_R17",
    "q25_mw7",
    "q25_mw4",
    "mw3_a5_d4_2a2_a1",
    "mw2_e6_d4_2a2_a1",
    "mw2_a5_a4_2a3_semistable",
)
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
scores = json.loads(SCORES.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert scores["status"] == "PASS_EXACT_A11_EQUATION_COST_SCORING"
assert scores["retained_count"] == scores["candidate_count"] == 2333
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"

frame_path = ROOT / neighbors["frame"]
frame = load_matrix(frame_path)
ns = block_diagonal_matrix(U2, -frame)
target_fibres = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"]
    if item["state"] == "equation_A11"
}
assert set(target_fibres) == set(TARGETS)
score_by_orbit = {
    int(item["candidate_id"]["orbit_index"]): item
    for item in scores["retained_candidates"]
}

records = []
for raw in neighbors["neighbors"]:
    orbit = int(raw["orbit_index"])
    fibre = vector(ZZ, raw["fiber"])
    score = score_by_orbit[orbit]
    marked = {
        target: int(target_fibres[target] * ns * fibre)
        for target in TARGETS
    }
    assert all(value >= 0 for value in marked.values())
    records.append({
        "candidate_id": {"q": 8, "old_fibre_degree": 2, "orbit_index": orbit},
        "child": score["child"],
        "marked_target_degrees": marked,
        "equation_cost_score": int(score["equation_cost_score"]),
        "P_dot_O": int(score["horizontal"]["P_dot_O"]),
        "expected_RR_ambient": int(score["expected_RR_ambient"]),
        "construction_tier": int(score["target_coset_mod_exact_sections"]["construction_tier"]),
        "total_explicit_degree_zero_count": int(
            score["explicit_curve_degrees"]["total_explicit_degree_zero_count"]
        ),
        "total_explicit_degree_one_count": int(
            score["explicit_curve_degrees"]["total_explicit_degree_one_count"]
        ),
        "declared_curve_nef_gate": score["declared_curve_nef_gate"],
    })

rankings = {
    target: sorted(
        records,
        key=lambda item: (
            item["marked_target_degrees"][target],
            item["equation_cost_score"],
            item["candidate_id"]["orbit_index"],
        ),
    )[:100]
    for target in TARGETS
}

pareto = []
for item in records:
    degree = item["marked_target_degrees"]["pinned_R17"]
    cost = item["equation_cost_score"]
    if not any(
        other["marked_target_degrees"]["pinned_R17"] <= degree
        and other["equation_cost_score"] <= cost
        and (
            other["marked_target_degrees"]["pinned_R17"] < degree
            or other["equation_cost_score"] < cost
        )
        for other in records
    ):
        pareto.append(item)
pareto.sort(key=lambda item: item["marked_target_degrees"]["pinned_R17"])

payload = {
    "schema": "elkies-k3.h3-a11-marked-target-neighbor-ranking.v1",
    "status": "PASS_EXACT_A11_MARKED_TARGET_NEIGHBOR_RANKING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (frame_path,)
        },
    },
    "candidate_count": len(records),
    "rankings_top_100": rankings,
    "pinned_degree_equation_cost_pareto_front": pareto,
    "named_orbits": {
        str(orbit): next(item for item in records if item["candidate_id"]["orbit_index"] == orbit)
        for orbit in (12, 829, 849, 1991, 2162)
    },
    "proof_boundary": (
        "Marked target degrees and equation-side curve costs are exact/planning "
        "data respectively. Component-dominance is inherited from the exhaustive "
        "shell; every promoted candidate still requires a full section-wall nef, "
        "marked-U, root, and bidirectional unimodular transport certificate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in TARGETS:
    item = rankings[target][0]
    print(
        "A11TARGET|target={}|orbit={}|degree={}|cost={}|PO={}|RR={}|child={}/MW{}".format(
            target, item["candidate_id"]["orbit_index"],
            item["marked_target_degrees"][target], item["equation_cost_score"],
            item["P_dot_O"], item["expected_RR_ambient"],
            item["child"]["ade"], item["child"]["mw_rank"],
        ),
        flush=True,
    )
for orbit in (12, 2162):
    item = payload["named_orbits"][str(orbit)]
    print(
        "A11TARGET|named_orbit={}|pinned_degree={}|cost={}|PO={}|RR={}".format(
            orbit, item["marked_target_degrees"]["pinned_R17"],
            item["equation_cost_score"], item["P_dot_O"], item["expected_RR_ambient"],
        ),
        flush=True,
    )
print(f"A11TARGET|pareto={len(pareto)}|status={payload['status']}|output={OUTPUT.resolve()}", flush=True)

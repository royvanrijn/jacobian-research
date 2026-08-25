#!/usr/bin/env sage -python
"""Rank compiler-gated orbit12 exits against the pinned semistable MW2 hub."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
GATE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-explicit-curve-gate.json"
SCORES = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-semistable-hub-survivor-ranking.json"
INPUTS = (GATE, SCORES, CROSSOVERS, ZERO, FRAME)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


gate = json.loads(GATE.read_text())
scores = json.loads(SCORES.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
zero = json.loads(ZERO.read_text())
parent = load_matrix(FRAME)
g_parent = block_diagonal_matrix(U2, -parent)
equation_to_parent = matrix(ZZ, zero["selected"]["equation_A11_to_explicit_zero_basis"])
equation_in_parent = equation_to_parent.inverse().change_ring(ZZ)
target_equation = vector(ZZ, next(
    item["target_fibre_in_state"] for item in crossovers["records"]
    if item["state"] == "equation_A11" and item["target"] == "mw2_a5_a4_2a3_semistable"
))
target_parent = target_equation * equation_in_parent
cost_by_id = {
    (item["candidate_id"]["q"], item["candidate_id"]["orbit_index"]): item
    for item in scores["ranked_candidates"]
}
records = []
for item in gate["survivors"]:
    raw = item["source_neighbor_record"]
    key = (int(raw["q"]), int(raw["orbit_index"]))
    score = cost_by_id[key]
    fibre = vector(ZZ, raw["fiber"])
    records.append({
        "candidate_id": item["candidate_id"],
        "child": item["child"],
        "semistable_hub_degree": int(fibre * g_parent * target_parent),
        "pinned_R17_degree": int(item["marked_target_degrees"]["pinned_R17"]),
        "equation_cost_score": int(score["equation_cost_score"]),
        "P_dot_O": int(score["horizontal"]["P_dot_O"]),
        "expected_RR_ambient": int(score["expected_RR_ambient"]),
        "explicit_degree_zero_count": len(item["explicit_degree_zero_curves"]),
        "explicit_degree_one_count": len(item["explicit_degree_one_curves"]),
    })
assert all(item["semistable_hub_degree"] >= 0 for item in records)
records.sort(key=lambda item: (item["semistable_hub_degree"], item["equation_cost_score"], item["candidate_id"]["q"], item["candidate_id"]["orbit_index"]))
pareto = []
for item in records:
    if not any(
        other["semistable_hub_degree"] <= item["semistable_hub_degree"]
        and other["equation_cost_score"] <= item["equation_cost_score"]
        and (other["semistable_hub_degree"] < item["semistable_hub_degree"] or other["equation_cost_score"] < item["equation_cost_score"])
        for other in records
    ):
        pareto.append(item)
payload = {
    "schema": "elkies-k3.h3-a5a5-semistable-hub-survivor-ranking.v1",
    "status": "PASS_EXACT_A5A5_SEMISTABLE_HUB_SURVIVOR_RANKING",
    "candidate_count": len(records),
    "top_100": records[:100],
    "degree_cost_pareto_front": pareto,
    "proof_boundary": "Exact marked degrees for candidates already passing the explicit-curve/parent-affine gate; full all-section/root/transport certification remains separate.",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0]
print("A5A5SEMHUB|candidates={}|best=q{}o{}|degree={}|cost={}|PO={}|RR={}|child={}/MW{}|pareto={}|status={}".format(
    len(records), best["candidate_id"]["q"], best["candidate_id"]["orbit_index"], best["semistable_hub_degree"],
    best["equation_cost_score"], best["P_dot_O"], best["expected_RR_ambient"], best["child"]["ade"],
    best["child"]["mw_rank"], len(pareto), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT}", flush=True)

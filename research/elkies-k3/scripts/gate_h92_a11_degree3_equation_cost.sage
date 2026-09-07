#!/usr/bin/env sage -python
"""Apply the exact all-section nef gate to cost-ranked A11 degree-3 exits."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-all.json"
COST = GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-equation-cost.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
O12 = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
OUTPUT = GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-full-nef-equation-cost.json"
INPUTS = (NEIGHBORS, COST, CROSSOVERS, O12)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


neighbors = json.loads(NEIGHBORS.read_text())
cost = json.loads(COST.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
o12 = json.loads(O12.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert cost["status"] == "PASS_EXACT_A11_EQUATION_COST_SCORING"
assert cost["candidate_count"] == cost["retained_count"] == len(neighbors["neighbors"])
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert o12["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"

frame_path = ROOT / neighbors["frame"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)
full_lattice = IntegralLattice(frame)
raw_by_id = {
    (int(raw["q"]), int(raw["old_fiber_degree"]), int(raw["orbit_index"])): raw
    for raw in neighbors["neighbors"]
}
targets = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}
targets["orbit12"] = vector(
    ZZ, matrix(ZZ, o12["selected"]["equation_A11_to_explicit_zero_basis"]).row(0)
)

survivors = []
affine_rejected = 0
section_rejected = 0
o385_audit = None
for score in cost["retained_candidates"]:
    identifier = score["candidate_id"]
    key = (identifier["q"], identifier["old_fibre_degree"], identifier["orbit_index"])
    raw = raw_by_id[key]
    degree = ZZ(raw["old_fiber_degree"])
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairing = degree - sum(labels)  # A11 highest root is (1,...,1).
    if affine_pairing < 0:
        affine_rejected += 1
        continue
    witness = vector(ZZ, raw["witness"])
    center = vector(QQ, witness) / degree
    closest = vector(ZZ, next(full_lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    nef = minimum_section >= 0
    audit = {
        "component_pairings": entries(labels),
        "affine_pairings": [int(affine_pairing)],
        "closest_section_distance": str(distance),
        "closest_section_vector": entries(closest),
        "minimum_section_intersection": str(minimum_section),
        "nef": bool(nef),
    }
    if key == (6, 3, 385):
        o385_audit = audit
    if not nef:
        section_rejected += 1
        continue
    fibre = vector(ZZ, raw["fiber"])
    marked = {name: int(fibre * g * target) for name, target in targets.items()}
    assert min(marked.values()) >= 0
    survivors.append({
        "candidate_id": identifier,
        "child": score["child"],
        "equation_cost_score": score["equation_cost_score"],
        "horizontal": score["horizontal"],
        "expected_RR_ambient": score["expected_RR_ambient"],
        "explicit_curve_degrees": score["explicit_curve_degrees"],
        "target_coset_mod_exact_sections": score["target_coset_mod_exact_sections"],
        "equation_cost_terms": score["equation_cost_terms"],
        "marked_target_degrees": marked,
        "nef_audit": audit,
        "source_neighbor_record": raw,
    })

assert o385_audit is not None and not o385_audit["nef"]
# The source cost artifact is already sorted by compiler score and stable
# tie-breaks, so filtering preserves the exact compiler ordering.
orbit12_cost_pareto = []
best_cost_seen = None
for item in sorted(
    survivors,
    key=lambda row: (
        row["marked_target_degrees"]["orbit12"], row["equation_cost_score"]
    ),
):
    if best_cost_seen is None or item["equation_cost_score"] < best_cost_seen:
        orbit12_cost_pareto.append(item)
        best_cost_seen = item["equation_cost_score"]
payload = {
    "schema": "elkies-k3.h3-a11-degree3-full-nef-equation-cost.v1",
    "status": "PASS_EXACT_A11_DEGREE3_FULL_NEF_EQUATION_COST_GATE",
    "input_candidate_count": len(neighbors["neighbors"]),
    "declared_curve_nef_count": sum(
        item["declared_curve_nef_gate"] == "PASS" for item in cost["retained_candidates"]
    ),
    "affine_rejected_count": affine_rejected,
    "section_rejected_count": section_rejected,
    "full_nef_candidate_count": len(survivors),
    "rejected_cost_leader_q6d3o385": o385_audit,
    "best_candidate": survivors[0] if survivors else None,
    "ranked_candidates_top_200": survivors[:200],
    "closest_to_pinned_top_100": sorted(
        survivors,
        key=lambda item: (
            item["marked_target_degrees"]["pinned_R17"],
            item["equation_cost_score"],
        ),
    )[:100],
    "closest_to_orbit12_top_100": sorted(
        survivors,
        key=lambda item: (
            item["marked_target_degrees"]["orbit12"],
            item["equation_cost_score"],
        ),
    )[:100],
    "orbit12_degree_equation_cost_pareto_front": orbit12_cost_pareto,
    "proof_boundary": (
        "Exact old-A11 component and all-section nef gate over the exhaustive "
        "degree-3 q6/q9/q12 shells, joined to exact equation-curve cost data and "
        "marked target degrees. A selected child still needs a standalone full "
        "marked-U/root/bidirectional-transport certificate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (frame_path,)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if survivors:
    best = survivors[0]
    closest = min(survivors, key=lambda item: item["marked_target_degrees"]["pinned_R17"])
    print(
        "A11D3GATE|inputs={}|nef={}|best=q{}d{}o{}|cost={}|PO={}|RR={}|pinned={}|closest=q{}d{}o{}:{}|o385_min={}|status={}".format(
            len(neighbors["neighbors"]), len(survivors), best["candidate_id"]["q"],
            best["candidate_id"]["old_fibre_degree"], best["candidate_id"]["orbit_index"],
            best["equation_cost_score"], best["horizontal"]["P_dot_O"],
            best["expected_RR_ambient"], best["marked_target_degrees"]["pinned_R17"],
            closest["candidate_id"]["q"], closest["candidate_id"]["old_fibre_degree"],
            closest["candidate_id"]["orbit_index"], closest["marked_target_degrees"]["pinned_R17"],
            o385_audit["minimum_section_intersection"], payload["status"],
        ),
        flush=True,
    )
else:
    print(
        "A11D3GATE|inputs={}|nef=0|o385_min={}|status={}".format(
            len(neighbors["neighbors"]), o385_audit["minimum_section_intersection"], payload["status"]
        ),
        flush=True,
    )
print(f"OUTPUT|{OUTPUT}", flush=True)

#!/usr/bin/env sage -python
"""Score exact-curve-nef q4/q6 exits from the explicit-zero orbit12 child."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--gate", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-explicit-curve-gate.json")
parser.add_argument("--neighbors", type=Path, default=GENERATED / "elkies-k3-h3-a11-o12-explicit-zero-q4q6-degree2-all.json")
parser.add_argument("--root-rank", type=int, default=10)
parser.add_argument("--output", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json")
args = parser.parse_args()
GATE = args.gate.resolve()
NEIGHBORS = args.neighbors.resolve()
OUTPUT = args.output.resolve()
INPUTS = (GATE, NEIGHBORS)


def load_matrix(path):
    return matrix(ZZ, [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [unseen.pop()]
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if cartan[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def highest_roots(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    result = []
    for component in components(cartan):
        candidates = [item for item in roots if all(value >= 0 for value in item)
                      and all(index in component or item[index] == 0 for index in range(cartan.nrows()))]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


gate = json.loads(GATE.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
assert gate["status"] == "PASS_EXACT_A5A5_Q6Q8_EXPLICIT_CURVE_GATE"
parent = load_matrix(ROOT / neighbors["frame"])
root_rank = args.root_rank
root = parent[:root_rank, :root_rank]
coupling = parent[:root_rank, root_rank:]
height = parent[root_rank:, root_rank:] - coupling.transpose() * root.inverse() * coupling
lattice = IntegralLattice(root)
root_edges = [(i, j) for i in range(root_rank) for j in range(i + 1, root_rank) if root[i, j] == -1]


def vertical_layers(coefficients):
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    previous = 0
    total = 0
    for level in sorted(set(value for value in magnitudes if value)):
        active = set(index for index, value in enumerate(magnitudes) if value >= level)
        count = 0
        while active:
            count += 1
            todo = [active.pop()]
            while todo:
                node = todo.pop()
                for left, right in root_edges:
                    other = right if left == node else left if right == node else None
                    if other in active:
                        active.remove(other)
                        todo.append(other)
        total += (level - previous) * count
        previous = level
    return int(total)


profile_cache = {}


def horizontal_profile(z):
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * parent[:, :root_rank]) * root.inverse()
    iterator = lattice.enumerate_close_vectors(-dual)
    minimum = None
    choices = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
        norm = QQ(lifted * parent * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            choices.append((ZZ(pole), vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))))
    assert choices
    choices.sort(key=lambda item: (item[0], tuple(item[1])))
    profile_cache[key] = choices[0]
    return choices[0]


highest = highest_roots(root)
records = []
for item in gate["survivors"]:
    raw = item["source_neighbor_record"]
    fibre = vector(ZZ, raw["fiber"])
    labels = vector(ZZ, raw["dominant_labels"])
    affine = [int(ZZ(raw["old_fiber_degree"]) - top * labels) for top in highest]
    negative_affine = [index for index, value in enumerate(affine) if value < 0]
    pole, section = horizontal_profile(raw["mw_projection"])
    residual = fibre - vector(ZZ, [-1, 1] + [0] * 17) - section
    vertical = vector(ZZ, residual[2:2 + root_rank])
    layers = vertical_layers(vertical)
    rr = 2 + 2 * int(pole) + layers
    zeros = item["explicit_degree_zero_curves"]
    ones = item["explicit_degree_one_curves"]
    root_adapted = "child_root_adapted_basis" in raw
    terms = {
        "parent_affine_non_nef_penalty": 1000000 if negative_affine else 0,
        "root_not_adapted_penalty": 1000000 if not root_adapted else 0,
        "P_dot_O": 900 * int(pole),
        "horizontal_degree": 250 * int(raw["old_fiber_degree"]),
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * layers,
        "vertical_support": 25 * sum(value != 0 for value in vertical),
        "child_root_count": int(raw["child_root_data"][1]),
        "coordinate_growth": item["coordinate_growth_max"],
        "no_explicit_degree_one_curve": 4000 if not ones else 0,
        "explicit_degree_one_credit": -500 * min(len(ones), 6),
        "explicit_degree_zero_credit": -100 * min(len(zeros), 12),
    }
    records.append({
        "candidate_id": item["candidate_id"],
        "child": item["child"],
        "root_adapted": root_adapted,
        "full_declared_nef_gate": "PASS" if not negative_affine and root_adapted else "REJECT",
        "parent_affine_component_pairings": affine,
        "negative_parent_affine_components": negative_affine,
        "horizontal": {"P_dot_O": int(pole), "section": entries(section), "vertical": entries(vertical), "vertical_layers": layers, "vertical_support": sum(value != 0 for value in vertical)},
        "expected_RR_ambient": rr,
        "explicit_degree_zero_curves": zeros,
        "explicit_degree_one_curves": ones,
        "marked_target_degrees": item["marked_target_degrees"],
        "coordinate_growth_max": item["coordinate_growth_max"],
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
    })

records.sort(key=lambda item: (item["equation_cost_score"], item["marked_target_degrees"]["pinned_R17"], item["candidate_id"]["q"], item["candidate_id"]["orbit_index"]))
passes = [item for item in records if item["full_declared_nef_gate"] == "PASS"]
pareto = []
for item in passes:
    if not any(other["equation_cost_score"] <= item["equation_cost_score"]
               and other["marked_target_degrees"]["pinned_R17"] <= item["marked_target_degrees"]["pinned_R17"]
               and (other["equation_cost_score"] < item["equation_cost_score"]
                    or other["marked_target_degrees"]["pinned_R17"] < item["marked_target_degrees"]["pinned_R17"])
               for other in passes):
        pareto.append(item)
pareto.sort(key=lambda item: item["marked_target_degrees"]["pinned_R17"])
payload = {
    "schema": "elkies-k3.h3-a5a5-explicit-zero-q4q6-equation-cost.v1",
    "status": "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "explicit_curve_gate_survivors": len(records),
    "full_declared_nef_gate_passes": len(passes),
    "ranked_candidates": records,
    "pinned_degree_equation_cost_pareto_front": pareto,
    "proof_boundary": "Exact explicit-curve/affine intersections and lattice profiles; full section-wall nefness and endpoint certification remain separate.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = passes[0]
print("A5A5EXCOST|survivors={}|passes={}|best=q{}o{}|score={}|PO={}|RR={}|pinned={}|deg0={}|deg1={}|child={}/MW{}|pareto={}|status={}".format(
    len(records), len(passes), best["candidate_id"]["q"], best["candidate_id"]["orbit_index"], best["equation_cost_score"],
    best["horizontal"]["P_dot_O"], best["expected_RR_ambient"], best["marked_target_degrees"]["pinned_R17"],
    len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]), best["child"]["ade"], best["child"]["mw_rank"], len(pareto), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)

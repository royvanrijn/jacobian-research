#!/usr/bin/env sage -python
"""Score exact low-q neighbours of the pinned two-polynomial MW2 hub."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, default=GENERATED / "elkies-k3-h3-semistable-mw2-q4q6-degree2-all.json")
parser.add_argument("--output", type=Path, default=GENERATED / "elkies-k3-h3-semistable-mw2-q4q6-equation-cost.json")
args = parser.parse_args()
NEIGHBORS = args.neighbors.resolve()
OUTPUT = args.output.resolve()
MARKING = GENERATED / "elkies-k3-h3-semistable-mw2-equation-marking.json"
INPUTS = (NEIGHBORS, MARKING)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


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


neighbors = json.loads(NEIGHBORS.read_text())
marking = json.loads(MARKING.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert marking["status"] == "PASS_EXACT_SEMISTABLE_MW2_EQUATION_MARKING"
frame_path = ROOT / marking["frame_output"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)
root_rank = 15
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
height = frame[root_rank:, root_rank:] - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
full_lattice = IntegralLattice(frame)
highest = highest_roots(root)
root_edges = [(i, j) for i in range(root_rank) for j in range(i + 1, root_rank) if root[i, j] == -1]


def horizontal_profile(z):
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    choices = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            choices.append((ZZ(pole), vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))))
    assert choices
    return min(choices, key=lambda item: (item[0], tuple(item[1])))


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


# Recover an actual optimal two-polynomial section basis in this root-adapted
# marking, rather than assuming that the LLL tail vectors already have it.
height_scale = lcm(entry.denominator() for entry in height.list())
short = pari((height_scale * height).change_ring(ZZ)).qfminim(4 * height_scale)
mw_vectors = []
for column in matrix(ZZ, short[2]).columns():
    for sign in (1, -1):
        z = sign * vector(ZZ, column)
        if z != 0 and z not in mw_vectors:
            mw_vectors.append(z)
polynomial = []
for z in mw_vectors:
    pole, section = horizontal_profile(z)
    if pole == 0:
        polynomial.append((z, section))
pairs = []
for left_z, left in polynomial:
    for right_z, right in polynomial:
        if abs(matrix(ZZ, [left_z, right_z]).det()) == 1 and left * g * right == 1:
            pairs.append((tuple(left), tuple(right), left_z, right_z, left, right))
assert pairs
unused_left_tuple, unused_right_tuple, p1_z, p2_z, p1, p2 = min(pairs)

old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
simple = [vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)]) for node in range(root_rank)]
affines = [vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank)) for top in highest]
curves = [old_zero] + simple + affines + [p1, p2]
curve_names = ["old_zero"] + [f"component_{i}" for i in range(root_rank)] + [f"affine_{i}" for i in range(len(affines))] + ["polynomial_P1", "polynomial_P2"]
assert all(curve * g * curve == -2 for curve in curves)

targets = {name: vector(ZZ, value) for name, value in marking["target_fibres_in_root_adapted_semistable"].items()}
records = []
for raw in neighbors["neighbors"]:
    fibre = vector(ZZ, raw["fiber"])
    old_fibre_degree = ZZ(raw["old_fiber_degree"])
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairings = [ZZ(old_fibre_degree - top * labels) for top in highest]
    if min(affine_pairings) < 0:
        continue
    witness = vector(ZZ, raw["witness"])
    center = vector(QQ, witness) / old_fibre_degree
    closest = vector(ZZ, next(full_lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = old_fibre_degree * (distance - 2) / 2
    if minimum_section < 0:
        continue
    degrees = [int(curve * g * fibre) for curve in curves]
    assert min(degrees) >= 0
    zeros = [curve_names[index] for index, value in enumerate(degrees) if value == 0]
    ones = [curve_names[index] for index, value in enumerate(degrees) if value == 1]
    if old_fibre_degree == 2:
        pole, section = horizontal_profile(raw["mw_projection"])
        residual = fibre - old_zero - section
        vertical = vector(ZZ, residual[2:2 + root_rank])
        layers = vertical_layers(vertical)
        horizontal = {
            "mode": "exact_degree2_decomposition",
            "P_dot_O": int(pole),
            "section": entries(section),
            "vertical": entries(vertical),
            "vertical_layers": layers,
        }
    else:
        section_pairs = []
        one_indices = [index for index, value in enumerate(degrees) if value == 1]
        for left_position, left in enumerate(one_indices):
            for right in one_indices[left_position + 1:]:
                section_pairs.append((int(curves[left] * g * curves[right]), left, right))
        if section_pairs:
            pole, left, right = min(section_pairs)
            assert pole >= 0
            section_pair = [curve_names[left], curve_names[right]]
        else:
            # No explicit pair is available.  The numerical placeholder is
            # accompanied by the large missing-curve penalty below and is not
            # presented as an exact section intersection.
            pole = 10
            section_pair = []
        vertical = labels
        layers = vertical_layers(vertical)
        horizontal = {
            "mode": "explicit_degree_one_pair_proxy",
            "P_dot_O": int(pole) if section_pairs else None,
            "section_pair": section_pair,
            "vertical": entries(vertical),
            "vertical_layers": layers,
        }
    rr = 2 + 2 * int(pole) + layers
    marked = {name: int(fibre * g * target) for name, target in targets.items()}
    terms = {
        "P_dot_O": 900 * int(pole),
        "horizontal_degree": 250 * int(old_fibre_degree),
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * layers,
        "vertical_support": 25 * sum(value != 0 for value in vertical),
        "coordinate_growth": int(max(abs(value) for value in fibre)),
        "no_explicit_degree_one_curve": 4000 if not ones else 0,
        "explicit_degree_one_credit": -500 * min(len(ones), 6),
        "explicit_degree_zero_credit": -100 * min(len(zeros), 12),
    }
    records.append({
        "candidate_id": {"q": int(raw["q"]), "old_fibre_degree": int(old_fibre_degree), "orbit_index": int(raw["orbit_index"])},
        "child": {"ade": raw["child_ade"], "mw_rank": int(raw["child_mw_rank"]), "root_data": raw["child_root_data"]},
        "component_pairings": entries(labels),
        "affine_pairings": entries(affine_pairings),
        "closest_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "horizontal": horizontal,
        "expected_RR_ambient": rr,
        "explicit_curve_degrees": dict(zip(curve_names, degrees)),
        "explicit_degree_zero_curves": zeros,
        "explicit_degree_one_curves": ones,
        "marked_target_degrees": marked,
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
        "source_neighbor_record": raw,
    })

records.sort(key=lambda item: (item["equation_cost_score"], item["marked_target_degrees"]["orbit12_fibre"], item["candidate_id"]["q"], item["candidate_id"]["orbit_index"]))
pareto = []
for item in records:
    if not any(other["equation_cost_score"] <= item["equation_cost_score"]
               and other["marked_target_degrees"]["orbit12_fibre"] <= item["marked_target_degrees"]["orbit12_fibre"]
               and (other["equation_cost_score"] < item["equation_cost_score"]
                    or other["marked_target_degrees"]["orbit12_fibre"] < item["marked_target_degrees"]["orbit12_fibre"])
               for other in records):
        pareto.append(item)
pareto.sort(key=lambda item: item["marked_target_degrees"]["orbit12_fibre"])
payload = {
    "schema": "elkies-k3.h3-semistable-mw2-frontier-cost.v1",
    "status": "PASS_EXACT_SEMISTABLE_MW2_FRONTIER_COST_SCORING",
    "polynomial_section_basis": {"P1_MW": entries(p1_z), "P1": entries(p1), "P2_MW": entries(p2_z), "P2": entries(p2), "pair_intersection": int(p1 * g * p2)},
    "input_candidate_count": len(neighbors["neighbors"]),
    "full_nef_candidate_count": len(records),
    "ranked_candidates": records,
    "orbit12_degree_equation_cost_pareto_front": pareto,
    "proof_boundary": "Exact component/all-section nef and marked-degree scoring with an explicit optimal polynomial section basis. Selected candidates still require a separate full marked-U/root/transport certificate.",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS + (frame_path,)}},
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if not records:
    print(
        "SEMFRONT|inputs={}|nef=0|status={}".format(
            len(neighbors["neighbors"]), payload["status"]
        ),
        flush=True,
    )
    print(f"OUTPUT|{OUTPUT}", flush=True)
    raise SystemExit(0)
best = records[0]
closest = min(records, key=lambda item: item["marked_target_degrees"]["orbit12_fibre"])
print("SEMFRONT|inputs={}|nef={}|best=q{}o{}|cost={}|PO={}|RR={}|orbit12={}|deg0={}|deg1={}|closest=q{}o{}:{}|pareto={}|status={}".format(
    len(neighbors["neighbors"]), len(records), best["candidate_id"]["q"], best["candidate_id"]["orbit_index"], best["equation_cost_score"],
    best["horizontal"]["P_dot_O"], best["expected_RR_ambient"], best["marked_target_degrees"]["orbit12_fibre"],
    len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]), closest["candidate_id"]["q"],
    closest["candidate_id"]["orbit_index"], closest["marked_target_degrees"]["orbit12_fibre"], len(pareto), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT}", flush=True)

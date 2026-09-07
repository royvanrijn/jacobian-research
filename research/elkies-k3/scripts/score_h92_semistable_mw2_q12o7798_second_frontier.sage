#!/usr/bin/env sage -python
"""Score a second exact beam layer beyond semistable q12/orbit7798.

The cost inventory uses the same already-explicit curves as the semistable
MW2 equation frame, transported through the exact determinant-one first hop.
Thus degree-zero/one credits refer to curves that are genuinely available to
an equation compiler, not merely short abstract classes in the child frame.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--neighbors",
    type=Path,
    default=GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-q4q6-degree2-all.json",
)
parser.add_argument(
    "--checkpoint",
    type=Path,
    default=GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-equation-marking.json",
)
parser.add_argument("--first-q", type=int, default=12)
parser.add_argument("--first-orbit", type=int, default=7798)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-q4q6-equation-cost.json",
)
args = parser.parse_args()
NEIGHBORS = args.neighbors.resolve()
OUTPUT = args.output.resolve()
SEM_MARKING = GENERATED / "elkies-k3-h3-semistable-mw2-equation-marking.json"
SEM_SCORE = GENERATED / "elkies-k3-h3-semistable-mw2-q8q10q12-equation-cost.json"
CHECKPOINT = args.checkpoint.resolve()
INPUTS = (NEIGHBORS, SEM_MARKING, SEM_SCORE, CHECKPOINT)
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
        candidates = [
            item for item in roots
            if all(value >= 0 for value in item)
            and all(index in component or item[index] == 0 for index in range(cartan.nrows()))
        ]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


neighbors = json.loads(NEIGHBORS.read_text())
sem_marking = json.loads(SEM_MARKING.read_text())
sem_score = json.loads(SEM_SCORE.read_text())
checkpoint = json.loads(CHECKPOINT.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert sem_marking["status"] == "PASS_EXACT_SEMISTABLE_MW2_EQUATION_MARKING"
assert sem_score["status"] == "PASS_EXACT_SEMISTABLE_MW2_FRONTIER_COST_SCORING"
assert checkpoint["status"] == "PASS_EXACT_SEMISTABLE_FRONTIER_EQUATION_MARKING_SEARCH_CHECKPOINT"

frame_path = ROOT / checkpoint["frame_output"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)
root_rank = checkpoint["child"]["root_data"][0]
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
height = frame[root_rank:, root_rank:] - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
full_lattice = IntegralLattice(frame)
highest = highest_roots(root)
root_edges = [
    (i, j) for i in range(root_rank) for j in range(i + 1, root_rank)
    if root[i, j] == -1
]


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


# Rebuild exactly the explicit semistable curve inventory used by the first
# scorer, then carry it to this child.  This includes its two polynomial MW
# generators and all old fibre components/affine components.
sem_frame_path = ROOT / sem_marking["frame_output"]
sem_frame = load_matrix(sem_frame_path)
g_sem = block_diagonal_matrix(U2, -sem_frame)
sem_root_rank = 15
sem_root = sem_frame[:sem_root_rank, :sem_root_rank]
sem_highest = highest_roots(sem_root)
old_zero_sem = vector(ZZ, [-1, 1] + [0] * 17)
simple_sem = [
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(sem_root_rank)
]
affines_sem = [
    vector(ZZ, [1, 0] + list(top) + [0] * (17 - sem_root_rank))
    for top in sem_highest
]
polynomial_data = sem_score["polynomial_section_basis"]
p1_sem = vector(ZZ, polynomial_data["P1"])
p2_sem = vector(ZZ, polynomial_data["P2"])
curves_sem = [old_zero_sem] + simple_sem + affines_sem + [p1_sem, p2_sem]
curve_names = (
    ["old_zero"]
    + [f"component_{index}" for index in range(sem_root_rank)]
    + [f"affine_{index}" for index in range(len(affines_sem))]
    + ["polynomial_P1", "polynomial_P2"]
)
transition = matrix(ZZ, checkpoint["root_adapted_semistable_to_child_basis"])
transition_inverse = transition.inverse().change_ring(ZZ)
curves = [curve * transition_inverse for curve in curves_sem]
assert all(curve in ZZ**19 and curve * g * curve == -2 for curve in curves)
first_fibre_sem = vector(ZZ, transition.row(0))
first_degrees = [int(curve * g_sem * first_fibre_sem) for curve in curves_sem]
first_record = next(
    item for item in sem_score["ranked_candidates"]
    if item["candidate_id"]["q"] == args.first_q
    and item["candidate_id"]["orbit_index"] == args.first_orbit
)
assert dict(zip(curve_names, first_degrees)) == first_record["explicit_curve_degrees"]

targets = {
    name: vector(ZZ, value) for name, value in checkpoint["target_fibres_in_child"].items()
}
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
records = []
explicit_curve_negative_rejected = 0
for raw in neighbors["neighbors"]:
    fibre = vector(ZZ, raw["fiber"])
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairings = [ZZ(2 - top * labels) for top in highest]
    if min(affine_pairings) < 0:
        continue
    witness = vector(ZZ, raw["witness"])
    center = vector(QQ, witness) / 2
    closest = vector(ZZ, next(full_lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = distance - 2
    if minimum_section < 0:
        continue
    degrees = [int(curve * g * fibre) for curve in curves]
    # Current components and sections do not exhaust effective multisections.
    # Every transported curve here is already explicit on the equation model,
    # so a negative pairing is an exact additional non-nef obstruction.
    if min(degrees) < 0:
        explicit_curve_negative_rejected += 1
        continue
    pole, section = horizontal_profile(raw["mw_projection"])
    residual = fibre - old_zero - section
    vertical = vector(ZZ, residual[2:2 + root_rank])
    layers = vertical_layers(vertical)
    rr = 2 + 2 * int(pole) + layers
    zeros = [curve_names[index] for index, value in enumerate(degrees) if value == 0]
    ones = [curve_names[index] for index, value in enumerate(degrees) if value == 1]
    marked = {name: int(fibre * g * target) for name, target in targets.items()}
    terms = {
        "P_dot_O": 900 * int(pole),
        "horizontal_degree": 500,
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * layers,
        "vertical_support": 25 * sum(value != 0 for value in vertical),
        "coordinate_growth": int(max(abs(value) for value in fibre)),
        "no_explicit_degree_one_curve": 4000 if not ones else 0,
        "explicit_degree_one_credit": -500 * min(len(ones), 6),
        "explicit_degree_zero_credit": -100 * min(len(zeros), 12),
    }
    records.append({
        "candidate_id": {
            "q": int(raw["q"]),
            "old_fibre_degree": 2,
            "orbit_index": int(raw["orbit_index"]),
        },
        "child": {
            "ade": raw["child_ade"],
            "mw_rank": int(raw["child_mw_rank"]),
            "root_data": raw["child_root_data"],
        },
        "component_pairings": entries(labels),
        "affine_pairings": entries(affine_pairings),
        "closest_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "horizontal": {
            "P_dot_O": int(pole),
            "section": entries(section),
            "vertical": entries(vertical),
            "vertical_layers": layers,
        },
        "expected_RR_ambient": rr,
        "explicit_curve_degrees": dict(zip(curve_names, degrees)),
        "explicit_degree_zero_curves": zeros,
        "explicit_degree_one_curves": ones,
        "marked_target_degrees": marked,
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
        "source_neighbor_record": raw,
    })

records.sort(
    key=lambda item: (
        item["equation_cost_score"],
        item["marked_target_degrees"]["orbit12_fibre"],
        item["candidate_id"]["q"],
        item["candidate_id"]["orbit_index"],
    )
)
pareto = []
for item in records:
    if not any(
        other["equation_cost_score"] <= item["equation_cost_score"]
        and other["marked_target_degrees"]["orbit12_fibre"] <= item["marked_target_degrees"]["orbit12_fibre"]
        and (
            other["equation_cost_score"] < item["equation_cost_score"]
            or other["marked_target_degrees"]["orbit12_fibre"] < item["marked_target_degrees"]["orbit12_fibre"]
        )
        for other in records
    ):
        pareto.append(item)
pareto.sort(key=lambda item: item["marked_target_degrees"]["orbit12_fibre"])

payload = {
    "schema": "elkies-k3.h3-semistable-mw2-second-frontier-cost.v1",
    "status": "PASS_EXACT_SEMISTABLE_SECOND_FRONTIER_COST_SCORING",
    "first_hop": {
        "candidate_id": first_record["candidate_id"],
        "child": first_record["child"],
        "equation_cost_score": first_record["equation_cost_score"],
        "marked_target_degrees": first_record["marked_target_degrees"],
    },
    "input_candidate_count": len(neighbors["neighbors"]),
    "full_nef_candidate_count": len(records),
    "explicit_curve_negative_rejected_count": explicit_curve_negative_rejected,
    "ranked_candidates": records,
    "orbit12_degree_equation_cost_pareto_front": pareto,
    "proof_boundary": (
        "Exact component/all-section nef and marked-degree scoring for the second layer, "
        "using transported already-explicit semistable curves. Selected candidates still "
        "require a separate full marked-U/root/transport certificate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [
            str(frame_path.relative_to(ROOT)), str(sem_frame_path.relative_to(ROOT))
        ],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (frame_path, sem_frame_path)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0]
closest = min(records, key=lambda item: item["marked_target_degrees"]["orbit12_fibre"])
print(
    "SEM2FRONT|inputs={}|nef={}|best=q{}o{}|cost={}|PO={}|RR={}|orbit12={}|deg0={}|deg1={}|closest=q{}o{}:{}|pareto={}|status={}".format(
        len(neighbors["neighbors"]), len(records), best["candidate_id"]["q"],
        best["candidate_id"]["orbit_index"], best["equation_cost_score"],
        best["horizontal"]["P_dot_O"], best["expected_RR_ambient"],
        best["marked_target_degrees"]["orbit12_fibre"],
        len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]),
        closest["candidate_id"]["q"], closest["candidate_id"]["orbit_index"],
        closest["marked_target_degrees"]["orbit12_fibre"], len(pareto), payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)

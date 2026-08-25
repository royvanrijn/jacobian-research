#!/usr/bin/env sage-python
"""Score an exact marked neighbour frontier by inherited equation cost.

This is the reusable beam scorer for a root-adapted source whose marking
contains ``equation_explicit_curves_in_child``.  It prices every candidate by
the cheapest exact root lift of its horizontal MW direction, resolved-RR size,
vertical complexity, coordinate growth, and degree-zero/one inherited curves.
If ``--nef-frontier`` is supplied, only candidates that passed that artifact's
component, affine, and all-section nef gates are retained.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, required=True)
parser.add_argument("--marking", type=Path, required=True)
parser.add_argument("--nef-frontier", type=Path)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--retain", type=int, default=300)
args = parser.parse_args()
if args.retain <= 0:
    parser.error("--retain must be positive")

NEIGHBORS = args.neighbors.resolve()
MARKING = args.marking.resolve()
NEF = args.nef_frontier.resolve() if args.nef_frontier else None
OUTPUT = args.output.resolve()


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


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
        candidates = [
            item for item in roots
            if all(value >= 0 for value in item)
            and all(index in component or item[index] == 0 for index in range(cartan.nrows()))
        ]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


def graph_edges(cartan):
    return [
        (left, right)
        for left in range(cartan.nrows())
        for right in range(left + 1, cartan.nrows())
        if cartan[left, right] == -1
    ]


def connected_count(edges, active):
    active = set(active)
    count = 0
    while active:
        count += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    return count


def vertical_layers(coefficients, edges):
    magnitudes = [abs(ZZ(item)) for item in coefficients]
    total = previous = ZZ(0)
    for level in sorted(set(item for item in magnitudes if item)):
        active = [index for index, item in enumerate(magnitudes) if item >= level]
        total += (level - previous) * connected_count(edges, active)
        previous = level
    return int(total)


def candidate_key(record):
    candidate = record.get("candidate_id", record)
    return (
        int(candidate["q"]),
        int(candidate.get("old_fibre_degree", candidate.get("old_fiber_degree"))),
        int(candidate["orbit_index"]),
    )


neighbors = json.loads(NEIGHBORS.read_text())
marking = json.loads(MARKING.read_text())
assert neighbors["status"] in {
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS",
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS_TARGET_FILTERED",
}
assert marking["status"] in {
    "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING",
    "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING",
    "PASS_EXACT_Q4O208_CHILD_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING",
    "PASS_EXACT_Q4O1584_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_CORRECTED_A3_2A2_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_Q8_5A1_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_AN_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_A2_EFFECTIVE_ZERO_MARKING",
}
assert "equation_explicit_curves_in_child" in marking

allowed = None
nef = None
marked_targets_by_candidate = {}
if NEF:
    nef = json.loads(NEF.read_text())
    assert nef["status"] == "PASS_EXACT_MARKED_ROOT_ADAPTED_FRONTIER_RANKING"
    allowed = {candidate_key(item) for item in nef["ranked_candidates"]}
    marked_targets_by_candidate = {
        candidate_key(item): item["marked_target_degrees"]
        for item in nef["ranked_candidates"]
    }

frame_path = ROOT / marking["frame_output"]
frame = load_matrix(frame_path)
root_rank = int(marking["root_data"][0])
mw_rank = 17 - root_rank
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
edges = graph_edges(root)
highest = highest_roots(root)
g_parent = block_diagonal_matrix(U2, -frame)

explicit_curves = {
    name: vector(ZZ, value)
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
assert all(curve * g_parent * curve == -2 for curve in explicit_curves.values())
explicit_sections = [curve for curve in explicit_curves.values() if curve[1] == 1]
section_mw = [vector(ZZ, curve[-mw_rank:]) for curve in explicit_sections if any(curve[-mw_rank:])]
known_section_lattice = (
    matrix(ZZ, section_mw).row_module()
    if section_mw else matrix(ZZ, 0, mw_rank).row_module()
)

profile_cache = {}


def section_profiles(z):
    z = vector(ZZ, z)
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    horizontal_height = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    result = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * mw_rank)
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            result.append((lifted, horizontal_height, norm - horizontal_height, ZZ(pole)))
    assert result
    profile_cache[key] = result
    return result


records = []
for raw in neighbors["neighbors"]:
    key = candidate_key(raw)
    if allowed is not None and key not in allowed:
        continue
    fibre = vector(ZZ, raw["fiber"])
    z = vector(ZZ, raw["mw_projection"])
    degree = int(raw["old_fiber_degree"])
    profiles = []
    for lifted, horizontal_height, correction, pole in section_profiles(z):
        section = vector(ZZ, [pole + 1, 1] + list(lifted))
        residual = fibre - (degree - 1) * vector(ZZ, [-1, 1] + [0] * 17) - section
        assert residual[1] == 0 and not any(residual[2 + root_rank:])
        vertical = vector(ZZ, residual[2:2 + root_rank])
        profiles.append({
            "height": str(horizontal_height),
            "local_correction": str(correction),
            "P_dot_O": int(pole),
            "section": entries(section),
            "vertical": entries(vertical),
            "fibre_twist": int(residual[0]),
            "vertical_support": sum(item != 0 for item in vertical),
            "vertical_L1": int(sum(abs(item) for item in vertical)),
            "vertical_layers": vertical_layers(vertical, edges),
        })
    profiles.sort(key=lambda item: (
        item["P_dot_O"], item["vertical_layers"], item["vertical_support"],
        item["vertical_L1"], tuple(item["section"]),
    ))
    horizontal = profiles[0]

    curve_degrees = {name: int(curve * g_parent * fibre) for name, curve in explicit_curves.items()}
    negative_curves = sorted(name for name, value in curve_degrees.items() if value < 0)
    explicit_degree_zero = sorted(name for name, value in curve_degrees.items() if value == 0)
    explicit_degree_one = sorted(name for name, value in curve_degrees.items() if value == 1)
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairings = [int(degree - top * labels) for top in highest]
    negative_affine = [index for index, value in enumerate(affine_pairings) if value < 0]
    declared_nef = not negative_curves and not negative_affine

    in_known = bool(known_section_lattice.rank() and z in known_section_lattice)
    rank_gap = 0 if in_known else mw_rank - int(known_section_lattice.rank())
    rr = 2 + 2 * horizontal["P_dot_O"] + horizontal["vertical_layers"]
    child_root_count = int(raw["child_root_data"][1])
    coordinate_growth = max(abs(int(item)) for item in fibre)
    terms = {
        "declared_non_nef_penalty": 1000000 if not declared_nef else 0,
        "unspanned_horizontal_rank_gap": 5000 * rank_gap,
        "P_dot_O": 900 * horizontal["P_dot_O"],
        "horizontal_degree": 250 * degree,
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * horizontal["vertical_layers"],
        "vertical_support": 25 * horizontal["vertical_support"],
        "child_root_count": child_root_count,
        "coordinate_growth": coordinate_growth,
        "no_explicit_degree_one_curve": 3000 if not explicit_degree_one else 0,
        "explicit_degree_one_credit": -500 * min(len(explicit_degree_one), 4),
        "explicit_degree_zero_credit": -100 * min(len(explicit_degree_zero), 8),
    }
    targets = marked_targets_by_candidate.get(key, {})
    records.append({
        "candidate_id": {"q": key[0], "old_fibre_degree": key[1], "orbit_index": key[2]},
        "declared_curve_and_affine_nef_gate": "PASS" if declared_nef else "REJECT",
        "child": {
            "ade": raw["child_ade"], "mw_rank": int(raw["child_mw_rank"]),
            "root_data": raw["child_root_data"],
        },
        "fibre": entries(fibre),
        "mw_projection": entries(z),
        "horizontal": horizontal,
        "expected_RR_ambient": rr,
        "exact_curve_degrees": curve_degrees,
        "negative_exact_curves": negative_curves,
        "explicit_degree_zero_curves": explicit_degree_zero,
        "explicit_degree_one_curves": explicit_degree_one,
        "affine_component_pairings": affine_pairings,
        "negative_affine_components": negative_affine,
        "known_section_subgroup": {
            "rank": int(known_section_lattice.rank()), "target_in_subgroup": in_known,
            "unspanned_rank_gap": rank_gap,
        },
        "marked_target_degrees": targets,
        "equation_cost_terms": terms,
        "equation_cost_score": sum(terms.values()),
        "source_neighbor_record": raw,
    })

records.sort(key=lambda item: (
    item["equation_cost_score"], item["horizontal"]["P_dot_O"],
    -len(item["explicit_degree_one_curves"]), item["expected_RR_ambient"],
    item["candidate_id"]["q"], item["candidate_id"]["orbit_index"],
))
retained = records[:args.retain]
inputs = [NEIGHBORS, MARKING, frame_path] + ([NEF] if NEF else [])


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


payload = {
    "schema": "elkies-k3.h3-marked-frontier-equation-cost.v1",
    "status": "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING",
    "source_hub": marking.get("source_hub", marking.get("hub")),
    "prefix_operational_score": marking.get("prefix_operational_score"),
    "candidate_count": len(records),
    "retained_count": len(retained),
    "known_explicit_section_subgroup_rank": int(known_section_lattice.rank()),
    "best_candidate": retained[0] if retained else None,
    "retained_candidates": retained,
    "proof_boundary": (
        "Exact inherited-curve intersections, horizontal closest-root lifts, and, when supplied, "
        "exact component/affine/all-section nef filtering. RR dimensions and weighted totals are "
        "planning estimates; a selected path still requires composed marked-U certification."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if retained:
    best = retained[0]
    print(
        "MARKEDCOST|hub={}|candidates={}|best=q{}o{}|child={}/MW{}|PO={}|RR={}|deg0={}|deg1={}|known_rank={}|score={}|combined={}|status={}".format(
            payload["source_hub"], len(records), best["candidate_id"]["q"],
            best["candidate_id"]["orbit_index"], best["child"]["ade"], best["child"]["mw_rank"],
            best["horizontal"]["P_dot_O"], best["expected_RR_ambient"],
            len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]),
            payload["known_explicit_section_subgroup_rank"], best["equation_cost_score"],
            (payload["prefix_operational_score"] or 0) + best["equation_cost_score"], payload["status"],
        ), flush=True,
    )
else:
    print(f"MARKEDCOST|hub={payload['source_hub']}|candidates=0|status={payload['status']}", flush=True)
print(f"OUTPUT|{OUTPUT}", flush=True)

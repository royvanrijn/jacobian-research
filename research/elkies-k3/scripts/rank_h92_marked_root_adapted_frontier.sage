#!/usr/bin/env sage -python
"""Apply exact component/section nef gates and rank a marked reverse frontier."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, required=True)
parser.add_argument("--marking", type=Path, required=True)
parser.add_argument("--target", default="orbit12")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--candidate-filter", type=Path,
    help="optional compact explicit-curve gate; only its exact survivor IDs are tested",
)
parser.add_argument(
    "--retain", type=int,
    help="store only this many ranked records while retaining exact total counts",
)
args = parser.parse_args()
if args.retain is not None and args.retain <= 0:
    parser.error("--retain must be positive")
NEIGHBORS = args.neighbors.resolve()
MARKING = args.marking.resolve()
OUTPUT = args.output.resolve()
CANDIDATE_FILTER = args.candidate_filter.resolve() if args.candidate_filter else None
INPUTS = (NEIGHBORS, MARKING) + ((CANDIDATE_FILTER,) if CANDIDATE_FILTER else ())
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
marking = json.loads(MARKING.read_text())
assert neighbors["status"] in {
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS",
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS_TARGET_FILTERED",
    "PASS_ROOT_ADAPTED_WEYL_SELECTED_NEIGHBORS",
}
assert marking["status"] in {
    "PASS_EXACT_A11_EQUATION_MARKING",
    "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "PASS_EXACT_A11_DEGREE3_CANDIDATE_LATTICE_CERTIFICATE",
    "PASS_EXACT_PINNED_R17_Q4_TARGETED_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_PINNED_R17_TARGETED_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING",
    "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_A5A5_ZERO_LOOP_RETURNED_MARKING",
    "PASS_EXACT_D13_ZERO_LOOP_RETURNED_MARKING",
    "PASS_EXACT_D12_ZERO_LOOP_RETURNED_MARKING",
    "PASS_EXACT_D13_ZERO_CHILD_MARKING",
    "PASS_EXACT_FIRST_Q8_SOURCE_MARKING",
    "PASS_EXACT_FIRST_Q8_ZERO_LOOP_RETURNED_MARKING",
    "PASS_EXACT_FIRST_Q8_LANDING_D13_MARKING",
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
frame_path = ROOT / marking["frame_output"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)
root_data = marking["root_data"] if "root_data" in marking else marking["child"]["root_data"]
root_rank = int(root_data[0])
root = frame[:root_rank, :root_rank]
highest = highest_roots(root)
full_lattice = IntegralLattice(frame)
target_key = (
    "target_fibres_in_root_adapted_hub"
    if "target_fibres_in_root_adapted_hub" in marking
    else "target_fibres_in_child"
)
targets = {
    name: vector(ZZ, value)
    for name, value in marking[target_key].items()
}
assert args.target in targets

allowed = None
if CANDIDATE_FILTER:
    candidate_filter = json.loads(CANDIDATE_FILTER.read_text())
    assert candidate_filter["status"] == "PASS_EXACT_MARKED_EXPLICIT_CURVE_GATE"
    allowed = {
        (
            int(item["candidate_id"]["q"]),
            int(item["candidate_id"]["old_fibre_degree"]),
            int(item["candidate_id"]["orbit_index"]),
        )
        for item in candidate_filter["survivors"]
    }

records = []
affine_rejected = 0
section_rejected = 0
prefiltered_candidate_count = 0
for raw in neighbors["neighbors"]:
    raw_key = (
        int(raw["q"]), int(raw["old_fiber_degree"]), int(raw["orbit_index"])
    )
    if allowed is not None and raw_key not in allowed:
        continue
    prefiltered_candidate_count += 1
    degree = ZZ(raw["old_fiber_degree"])
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairings = [ZZ(degree - top * labels) for top in highest]
    if affine_pairings and min(affine_pairings) < 0:
        affine_rejected += 1
        continue
    witness = vector(ZZ, raw["witness"])
    center = vector(QQ, witness) / degree
    closest = vector(ZZ, next(full_lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    if minimum_section < 0:
        section_rejected += 1
        continue
    fibre = vector(ZZ, raw["fiber"])
    marked = {name: int(fibre * g * target) for name, target in targets.items()}
    assert min(marked.values()) >= 0
    component_degrees = list(map(int, labels)) + list(map(int, affine_pairings))
    records.append({
        "candidate_id": {
            "q": int(raw["q"]),
            "old_fibre_degree": int(degree),
            "orbit_index": int(raw["orbit_index"]),
        },
        "child": {
            "ade": raw["child_ade"],
            "mw_rank": int(raw["child_mw_rank"]),
            "root_data": raw["child_root_data"],
        },
        "component_pairings": list(map(int, labels)),
        "affine_pairings": list(map(int, affine_pairings)),
        "component_degree_zero_count": component_degrees.count(0),
        "component_degree_one_count": component_degrees.count(1),
        "closest_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "P_dot_O": int(fibre[0] - fibre[1]),
        "marked_target_degrees": marked,
        "coordinate_growth_max": int(max(abs(value) for value in fibre)),
        "source_neighbor_record": raw,
    })

records.sort(
    key=lambda item: (
        item["marked_target_degrees"][args.target],
        -item["component_degree_one_count"],
        -item["component_degree_zero_count"],
        item["candidate_id"]["q"],
        item["candidate_id"]["orbit_index"],
    )
)
payload = {
    "schema": "elkies-k3.h3-marked-root-adapted-frontier-ranking.v1",
    "status": "PASS_EXACT_MARKED_ROOT_ADAPTED_FRONTIER_RANKING",
    "hub": marking.get("hub", marking.get("source_hub", "equation_A11_degree3_candidate")),
    "ranking_target": args.target,
    "input_candidate_count": len(neighbors["neighbors"]),
    "prefiltered_candidate_count": prefiltered_candidate_count,
    "affine_rejected_count": affine_rejected,
    "section_rejected_count": section_rejected,
    "full_nef_candidate_count": len(records),
    "retained_candidate_count": min(len(records), args.retain) if args.retain else len(records),
    "ranked_candidates": records[:args.retain] if args.retain else records,
    "proof_boundary": (
        "Exact current component and all-section nef gates, marked target degrees, "
        "and source-to-child unimodular data inherited from the exhaustive search. "
        "A selected continuation still needs a composed full marked route certificate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in INPUTS + (frame_path,)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if records:
    best = records[0]
    print(
        "MARKEDFRONT|hub={}|inputs={}|nef={}|best=q{}d{}o{}|target={}:{}|child={}/MW{}|deg0={}|deg1={}|status={}".format(
            marking.get("hub", marking.get("source_hub", "equation_A11_degree3_candidate")), len(neighbors["neighbors"]), len(records),
            best["candidate_id"]["q"], best["candidate_id"]["old_fibre_degree"],
            best["candidate_id"]["orbit_index"], args.target,
            best["marked_target_degrees"][args.target], best["child"]["ade"],
            best["child"]["mw_rank"], best["component_degree_zero_count"],
            best["component_degree_one_count"], payload["status"],
        ),
        flush=True,
    )
else:
    print(
        "MARKEDFRONT|hub={}|inputs={}|nef=0|status={}".format(
            marking.get("hub", marking.get("source_hub", "equation_A11_degree3_candidate")), len(neighbors["neighbors"]), payload["status"]
        ),
        flush=True,
    )
print(f"OUTPUT|{OUTPUT}", flush=True)

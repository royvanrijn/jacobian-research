#!/usr/bin/env sage
"""Certify a minimum-pole MW basis for the NS0024 A6+A4+A3 source.

The calculation is entirely in the exact root-adapted positive frame.  It
enumerates every norm-four physical section, proves that their MW images have
rank three, and exhibits a unimodular MW basis after adjoining one norm-six
section.  Thus an equation chart needs at least one section with P.O=1, and
the emitted basis has the minimum pole multiset (0,0,0,1).
"""

import argparse
import json
from itertools import combinations
from pathlib import Path

from sage.all import QQ, ZZ, matrix, pari, vector


ROOT = Path.cwd().resolve()
FRAME = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-r13-root-adapted.txt"
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in FRAME.read_text().splitlines()
    if line and not line.startswith("#")
])
source = json.loads(SOURCE.read_text())
height = matrix(QQ, source["source"]["mw_height_gram"])
assert frame.nrows() == 17 and height.nrows() == 4

# The deterministic simple-root order in the root-adapted frame consists of
# A4, A3, A6.  These chain orders were checked directly from the Cartan block.
chains = ((7, 12, 10, 9, 8, 11), (1, 0, 2, 3), (5, 4, 6))
orders = (7, 5, 4)


def physical_vectors(norm):
    half = matrix(ZZ, pari(frame).qfminim(norm)[2].sage()).columns()
    result = []
    for column in half:
        if column * frame * column != norm:
            continue
        for sign in (1, -1):
            candidate = sign * vector(ZZ, column)
            pairings = candidate * frame
            if all(pairings[index] >= 0 for index in range(13)) and all(
                pairings[index] in (0, 1) for index in range(13)
            ):
                result.append(candidate)
    return [vector(ZZ, item) for item in sorted({tuple(item) for item in result})]


def profile(section):
    pairings = section * frame
    answer = []
    for chain in chains:
        met = [index for index, root in enumerate(chain, 1) if pairings[root] == 1]
        assert len(met) <= 1
        answer.append(met[0] if met else 0)
    return tuple(answer)


polynomial = physical_vectors(4)
simple_pole = physical_vectors(6)
polynomial_images = matrix(ZZ, [list(item[13:]) for item in polynomial])
assert polynomial_images.rank() == 3

minimum_bases = []
for polynomial_part in combinations(polynomial, 3):
    for pole_section in simple_pole:
        candidate = polynomial_part + (pole_section,)
        quotient = matrix(ZZ, [list(item[13:]) for item in candidate])
        if abs(quotient.det()) != 1:
            continue
        candidate_profiles = tuple(profile(item) for item in candidate)
        node_count = sum(label != 0 for item in candidate_profiles for label in item)
        pair_intersection_sum = 0
        for left in range(4):
            for right in range(left):
                kl = (candidate[left] * frame * candidate[left] - 2) // 2
                kr = (candidate[right] * frame * candidate[right] - 2) // 2
                pair_intersection_sum += kl + kr - candidate[left] * frame * candidate[right]
        minimum_bases.append((
            -node_count,
            pair_intersection_sum,
            candidate_profiles,
            tuple(tuple(item) for item in candidate),
        ))
minimum_bases.sort()
assert minimum_bases
best_score = minimum_bases[0][:2]
best_bases = [item for item in minimum_bases if item[:2] == best_score]

chosen = (
    vector(ZZ, [0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,-1,0]),
    vector(ZZ, [1,1,1,0,0,0,0,1,1,1,1,1,1,0,-1,0,0]),
    vector(ZZ, [1,1,1,1,0,1,0,1,2,2,2,1,2,-1,0,0,0]),
    vector(ZZ, [1,1,1,1,0,1,0,1,1,1,1,0,1,0,0,0,-1]),
)
assert tuple(tuple(item) for item in chosen) == best_bases[0][3]
assert all(item in polynomial for item in chosen[:3])
assert chosen[3] in simple_pole
quotient_basis = matrix(ZZ, [list(item[13:]) for item in chosen])
assert abs(quotient_basis.det()) == 1

pole_orders = tuple((item * frame * item - 4) // 2 for item in chosen)
assert pole_orders == (0, 0, 0, 1)
profiles = tuple(profile(item) for item in chosen)
assert profiles == ((1,0,0), (2,1,3), (2,1,1), (1,1,1))
intersection = matrix(ZZ, 4, 4)
for left in range(4):
    for right in range(4):
        kl = (chosen[left] * frame * chosen[left] - 2) // 2
        kr = (chosen[right] * frame * chosen[right] - 2) // 2
        intersection[left, right] = kl + kr - chosen[left] * frame * chosen[right]
assert intersection == matrix(ZZ, [
    [-2,1,2,1], [1,-2,0,1], [2,0,-2,1], [1,1,1,-2]
])


def fibre_correction(left, right):
    correction = QQ(0)
    for order, left_label, right_label in zip(orders, profiles[left], profiles[right]):
        correction += QQ(min(left_label, right_label) * (order - max(left_label, right_label))) / order
    return correction


shioda = matrix(QQ, 4, 4)
for left in range(4):
    for right in range(4):
        shioda[left, right] = (
            2 + pole_orders[left] + pole_orders[right]
            - intersection[left, right] - fibre_correction(left, right)
        )
assert shioda == quotient_basis * height * quotient_basis.transpose()
assert shioda.det() == QQ(95) / 14


def optimal_basis_record(entry):
    candidate = tuple(vector(ZZ, item) for item in entry[3])
    candidate_profiles = tuple(profile(item) for item in candidate)
    candidate_quotient = matrix(ZZ, [list(item[13:]) for item in candidate])
    candidate_poles = tuple((item * frame * item - 4) // 2 for item in candidate)
    candidate_intersection = matrix(ZZ, 4, 4)
    for left in range(4):
        for right in range(4):
            kl = (candidate[left] * frame * candidate[left] - 2) // 2
            kr = (candidate[right] * frame * candidate[right] - 2) // 2
            candidate_intersection[left, right] = (
                kl + kr - candidate[left] * frame * candidate[right]
            )
    target_q4_coordinate = (-1, 0, 0, 0)
    return {
        "profiles_I7_I5_I4": [list(map(int, row)) for row in candidate_profiles],
        "component_depth_sum": int(sum(
            min(label, order - label)
            for row in candidate_profiles
            for label, order in zip(row, orders)
            if label
        )),
        "pole_orders": list(map(int, candidate_poles)),
        "quotient_coordinates": rows(candidate_quotient),
        "quotient_determinant": int(candidate_quotient.det()),
        "section_intersection_gram": rows(candidate_intersection),
        "contains_q4_orbit1_projection": any(
            tuple(map(int, row)) == target_q4_coordinate
            for row in candidate_quotient.rows()
        ),
        "q4_orbit1_basis_index": next((
            index
            for index, row in enumerate(candidate_quotient.rows(), 1)
            if tuple(map(int, row)) == target_q4_coordinate
        ), None),
    }


optimal_basis_records = [optimal_basis_record(entry) for entry in best_bases]
q4_basis_frontier = []
for entry in minimum_bases:
    record = optimal_basis_record(entry)
    if not record["contains_q4_orbit1_projection"]:
        continue
    node_count = sum(
        label != 0 for row in record["profiles_I7_I5_I4"] for label in row
    )
    pair_intersection_sum = sum(
        record["section_intersection_gram"][left][right]
        for left in range(4) for right in range(left)
    )
    record["labelled_node_incidences"] = int(node_count)
    record["total_pair_intersection"] = int(pair_intersection_sum)
    q4_basis_frontier.append(record)
q4_basis_frontier.sort(key=lambda record: (
    -record["component_depth_sum"],
    record["total_pair_intersection"],
    -record["labelled_node_incidences"],
    record["profiles_I7_I5_I4"],
    record["quotient_coordinates"],
))
assert q4_basis_frontier
resolved_chart_recommendation = q4_basis_frontier[0]

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-mw4-minimum-basis.v1",
    "status": "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS",
    "inputs": {
        str(FRAME.relative_to(ROOT)): source["source"]["root_adapted_gram"],
        str(SOURCE.relative_to(ROOT)): source["source"]["mw_height_gram"],
    },
    "enumeration": {
        "physical_norm_four_section_count": len(polynomial),
        "norm_four_mw_image_rank": int(polynomial_images.rank()),
        "physical_norm_six_section_count": len(simple_pole),
        "minimum_total_P_dot_O_for_a_basis": int(1),
        "minimum_pole_bases_tested": len(minimum_bases),
        "equation_complexity_objective": "maximize labelled singular-node incidences, then minimize total pair intersection",
        "maximum_labelled_singular_node_incidences": int(-best_score[0]),
        "minimum_total_pair_intersection_at_maximum_incidences": int(best_score[1]),
        "number_of_bases_at_optimum": len(best_bases),
        "optimal_bases": optimal_basis_records,
        "q4_containing_basis_count": len(q4_basis_frontier),
        "resolved_component_depth_recommendation": resolved_chart_recommendation,
        "completeness": "PARI qfminim exhausts all frame vectors through norm six; the physical chamber inequalities are exact.",
    },
    "basis": [
        {
            "name": f"P{index + 1}",
            "frame_vector": list(map(int, item)),
            "mw_quotient_coordinates": list(map(int, item[13:])),
            "P_dot_O": int(pole_orders[index]),
            "components_I7_I5_I4": list(map(int, profiles[index])),
        }
        for index, item in enumerate(chosen)
    ],
    "mw_quotient_basis_determinant": int(quotient_basis.det()),
    "section_intersection_gram": rows(intersection),
    "shioda_height_gram": [[str(entry) for entry in row] for row in shioda.rows()],
    "shioda_height_determinant": str(shioda.det()),
    "proof_boundary": {
        "proved": "Exact minimum-pole basis and its component/intersection/height marking in the NS0024 source lattice.",
        "not_proved": "This lattice certificate alone does not construct a finite-field or characteristic-zero equation family.",
    },
    "reproduce": "sage elkies-k3/scripts/certify_lattice_foundry_ns0024_mw4_basis.sage",
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if args.output.read_text() != serialized:
        raise SystemExit("minimum-basis artifact is stale")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)

print(
    "NS0024MW4BASIS|status=PASS"
    f"|norm4={len(polynomial)}|norm4_mw_rank={polynomial_images.rank()}"
    f"|norm6={len(simple_pole)}|pole_profile=0,0,0,1|det={shioda.det()}"
)

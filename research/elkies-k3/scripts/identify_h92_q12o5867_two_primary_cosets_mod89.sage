#!/usr/bin/env sage -python
"""Identify the three missing p=89 MW cosets by literal doubling.

The historical two-profile shell supplies the rank-13, index-eight current
subgroup.  For each lattice class occurring among the Abel-trace alternatives
of the complete sixteen-profile shell, express twice its marked MW vector in
that subgroup.  Retain a modular section/class pairing only when the literal
function-field identity [2]Q=P holds.  This is a classifier/certificate pass
over the already bounded shell, not a new section search.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, MixedIntegerLinearProgram, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--current-shell", type=Path,
    default=LOCAL / "q12o5867-p0-shell-all-records-mod89.json",
)
parser.add_argument(
    "--current-classifier", type=Path,
    default=LOCAL / "q12o5867-p0-shell-lattice-classification-mod89.json",
)
parser.add_argument(
    "--complete-shell", type=Path,
    default=LOCAL / "q12o5867-p0-shell-all-profiles-mod89.json",
)
parser.add_argument(
    "--complete-classifier", type=Path,
    default=LOCAL / "q12o5867-p0-shell-lattice-classification-all-profiles-mod89.json",
)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q12o5867-two-primary-cosets-mod89.json",
)
args = parser.parse_args()
started = time.monotonic()


def resolved(path):
    return path if path.is_absolute() else ROOT / path


paths = {
    "q8": Q8,
    "current_shell": resolved(args.current_shell),
    "current_classifier": resolved(args.current_classifier),
    "complete_shell": resolved(args.complete_shell),
    "complete_classifier": resolved(args.complete_classifier),
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
current_shell = json.loads(paths["current_shell"].read_text())
current_classifier = json.loads(paths["current_classifier"].read_text())
complete_shell = json.loads(paths["complete_shell"].read_text())
complete_classifier = json.loads(paths["complete_classifier"].read_text())
prime = ZZ(complete_shell["prime"])
assert prime == 89
assert current_shell["prime"] == complete_shell["prime"]
assert current_classifier["prime"] == complete_classifier["prime"] == int(prime)

F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()


def reduce_qq(value):
    value = QQ(value)
    assert value.denominator() % prime
    return F(value.numerator())/F(value.denominator())


A = K(R([reduce_qq(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]]))
B = K(R([reduce_qq(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]]))


def section_point(record):
    point = (
        K(R(record["x_coefficients_low_to_high"])),
        K(R(record["y_coefficients_low_to_high"])),
    )
    assert point[1]**2 == point[0]**3+A*point[0]+B
    return point


def point_neg(point):
    return None if point is None else (point[0], -point[1])


def point_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+A)/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    answer = (x3, slope*(x1-x3)-y1)
    assert answer[1]**2 == answer[0]**3+A*answer[0]+B
    return answer


def point_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return point_mul(-coefficient, point_neg(point))
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = point_add(answer, addend)
        addend = point_add(addend, addend)
        coefficient >>= 1
    return answer


def point_word(coefficients, points):
    answer = None
    for coefficient, point in zip(coefficients, points):
        if coefficient:
            answer = point_add(answer, point_mul(coefficient, point))
    return answer


def function_record(value):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    return {
        "numerator_coefficients_low_to_high": list(map(int, numerator.list())),
        "denominator_coefficients_low_to_high": list(map(int, denominator.list())),
        "numerator_degree": int(numerator.degree()),
        "denominator_degree": int(denominator.degree()),
    }


classes = complete_classifier["lattice_shell"]["classes"]
class_vectors = {
    int(row["class_index"]): vector(ZZ, row["current_4A1_mw"])
    for row in classes
}

# Use only literal one-candidate classifications from the historical shell.
# Pairwise refinements are deliberately excluded because their anchor chain is
# exactly what this certificate is auditing.
current_by_class = {}
for row in current_classifier["polynomial_shell"]["records"]:
    alternatives = row["profile_compatible_lattice_class_indices"]
    if row["ordinary_coefficient_jacobian_rank"] == 12 and len(alternatives) == 1:
        current_by_class.setdefault(int(alternatives[0]), int(row["shell_index"]))

current_class_indices = sorted(current_by_class)
assert len(current_class_indices) == 45
current_vectors = [class_vectors[index] for index in current_class_indices]
current_matrix = matrix(ZZ, current_vectors)
assert current_matrix.rank() == 13
current_module = current_matrix.row_module(ZZ)
current_hnf = current_module.basis_matrix()
assert abs(current_hnf.det()) == 8
assert current_hnf.elementary_divisors() == [1]*10+[2, 2, 2]
current_points = [
    section_point(current_shell["all_records"][current_by_class[index]])
    for index in current_class_indices
]


def quotient_key(mw_vector):
    remainder = vector(ZZ, mw_vector)
    for row_index in range(10):
        assert current_hnf[row_index, row_index] == 1
        remainder -= remainder[row_index]*current_hnf[row_index]
    assert all(remainder[index] == 0 for index in range(10))
    return tuple(int(remainder[index] % 2) for index in range(10, 13))


assert {quotient_key(vector(ZZ, [1 if i == j else 0 for i in range(13)]))
        for j in range(10, 13)} == {(1, 0, 0), (0, 1, 0), (0, 0, 1)}


def minimum_l1_word(target):
    problem = MixedIntegerLinearProgram(maximization=False, solver="GLPK")
    positive = problem.new_variable(integer=True, nonnegative=True)
    negative = problem.new_variable(integer=True, nonnegative=True)
    for coordinate in range(13):
        problem.add_constraint(
            sum(
                (positive[i]-negative[i])*current_vectors[i][coordinate]
                for i in range(len(current_vectors))
            ) == target[coordinate]
        )
    problem.set_objective(sum(
        positive[i]+negative[i] for i in range(len(current_vectors))
    ))
    problem.solve()
    coefficients = vector(ZZ, [
        ZZ(round(problem.get_values(positive[i])))
        -ZZ(round(problem.get_values(negative[i])))
        for i in range(len(current_vectors))
    ])
    assert sum(coefficients[i]*current_vectors[i]
               for i in range(len(current_vectors))) == target
    return coefficients


hnf_coordinate_doubled_targets = []
for coordinate in range(10, 13):
    half_vector = vector(ZZ, [1 if index == coordinate else 0 for index in range(13)])
    doubled_vector = 2*half_vector
    coefficients = minimum_l1_word(doubled_vector)
    target_point = point_word(coefficients, current_points)
    hnf_coordinate_doubled_targets.append({
        "quotient_generator_current_4A1_mw": list(map(int, half_vector)),
        "quotient_key": list(quotient_key(half_vector)),
        "doubled_target_current_4A1_mw": list(map(int, doubled_vector)),
        "minimum_l1_current_section_word": [
            {"coefficient": int(coefficient), "physical_class_index": class_id,
             "shell_index": current_by_class[class_id]}
            for coefficient, class_id in zip(coefficients, current_class_indices)
            if coefficient
        ],
        "minimum_l1_norm": int(sum(abs(value) for value in coefficients)),
        "target_section_mod89": {
            "x": function_record(target_point[0]),
            "y": function_record(target_point[1]),
        },
    })


candidate_class_indices = [
    int(row["class_index"])
    for row in classes
    if quotient_key(class_vectors[int(row["class_index"])]) != (0, 0, 0)
]

target_cache = {}
for class_index in candidate_class_indices:
    doubled_vector = 2*class_vectors[class_index]
    assert doubled_vector in current_module
    coefficients = minimum_l1_word(doubled_vector)
    target_cache[class_index] = {
        "coefficients": coefficients,
        "point": point_word(coefficients, current_points),
    }

abstract_target_candidates = []
for class_index in candidate_class_indices:
    coefficients = target_cache[class_index]["coefficients"]
    target_point = target_cache[class_index]["point"]
    abstract_target_candidates.append({
        "physical_class_index": class_index,
        "quotient_key": list(quotient_key(class_vectors[class_index])),
        "half_target_current_4A1_mw": list(map(int, class_vectors[class_index])),
        "doubled_target_current_4A1_mw": list(map(int, 2*class_vectors[class_index])),
        "minimum_l1_current_section_word": [
            {"coefficient": int(coefficient), "physical_class_index": class_id,
             "shell_index": current_by_class[class_id]}
            for coefficient, class_id in zip(coefficients, current_class_indices)
            if coefficient
        ],
        "minimum_l1_norm": int(sum(abs(value) for value in coefficients)),
        "target_section_mod89": {
            "x": function_record(target_point[0]),
            "y": function_record(target_point[1]),
        },
    })

abstract_target_candidates.sort(key=lambda row: (
    row["minimum_l1_norm"], tuple(row["quotient_key"]),
    row["physical_class_index"],
))
canonical_doubled_targets = []
selected_span = matrix(GF(2), 0, 3, []).row_module()
for row in abstract_target_candidates:
    key_vector = vector(GF(2), row["quotient_key"])
    if key_vector in selected_span:
        continue
    canonical_doubled_targets.append(row)
    selected_span = matrix(
        GF(2), [item["quotient_key"] for item in canonical_doubled_targets]
    ).row_module()
    if selected_span.dimension() == 3:
        break
assert len(canonical_doubled_targets) == 3
short_target_candidates_by_coset = {}
for row in abstract_target_candidates:
    key = "".join(map(str, row["quotient_key"]))
    if len(short_target_candidates_by_coset.setdefault(key, [])) < 10:
        short_target_candidates_by_coset[key].append(row)

shell_rows = {
    int(row["shell_index"]): row
    for row in complete_classifier["polynomial_shell"]["records"]
}
doubled_shell_points = {}
for shell_index in shell_rows:
    doubled_point = point_mul(
        2, section_point(complete_shell["all_records"][shell_index])
    )
    doubled_shell_points.setdefault(doubled_point, []).append(shell_index)

verified = []
for class_index in candidate_class_indices:
    for shell_index in doubled_shell_points.get(target_cache[class_index]["point"], []):
        row = shell_rows[shell_index]
        coefficients = target_cache[class_index]["coefficients"]
        word = [
            {"coefficient": int(coefficient), "physical_class_index": class_id,
             "shell_index": current_by_class[class_id]}
            for coefficient, class_id in zip(coefficients, current_class_indices)
            if coefficient
        ]
        verified.append({
            "physical_class_index": class_index,
            "shell_index": shell_index,
            "quotient_key": list(quotient_key(class_vectors[class_index])),
            "equation_component_profile": row["equation_component_profile"],
            "current_component_pairings": classes[class_index]["current_component_pairings"],
            "doubled_target_current_4A1_mw": list(map(int, 2*class_vectors[class_index])),
            "minimum_l1_current_section_word": word,
            "minimum_l1_norm": int(sum(abs(value) for value in coefficients)),
            "literal_function_field_doubling_verified": True,
            "abel_trace_matching_lattice_class": (
                class_index in row["trace_matching_lattice_class_indices"]
            ),
            "profile_compatible_under_selected_transport": (
                class_index in row["profile_compatible_lattice_class_indices"]
            ),
        })

verified.sort(key=lambda row: (
    tuple(row["quotient_key"]), row["minimum_l1_norm"],
    row["physical_class_index"], row["shell_index"],
))
by_coset = {}
for row in verified:
    key = "".join(map(str, row["quotient_key"]))
    by_coset.setdefault(key, []).append(row)

nonzero_cosets = [f"{value:03b}" for value in range(1, 8)]
print("literal doubling coset counts:", {key: len(value) for key, value in by_coset.items()})
assert set(by_coset) == {"010", "101", "111"}
literal_shell_span = current_module + matrix(ZZ, [
    class_vectors[row["physical_class_index"]] for row in verified
]).row_module(ZZ)
assert current_module.index_in(literal_shell_span) == 4
assert abs(literal_shell_span.basis_matrix().det()) == 2
output = {
    "schema": "q12o5867-two-primary-cosets-mod89-v1",
    "status": "PASS_EXACT_MOD89_COMPLETE_PROFILE_SHELL_ADDS_TWO_OF_THREE_TWO_DIRECTIONS",
    "prime": int(prime),
    "inputs": {name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
               for name, path in paths.items()},
    "current_subgroup": {
        "strict_unique_section_count": len(current_class_indices),
        "rank": int(current_matrix.rank()),
        "smith_diagonal": list(map(int, current_hnf.elementary_divisors())),
        "index": int(abs(current_hnf.det())),
        "row_hnf": [list(map(int, row)) for row in current_hnf.rows()],
        "physical_class_indices": current_class_indices,
    },
    "quotient": {
        "abstract_group": "(Z/2Z)^3",
        "coordinate_convention": (
            "reduce current_4A1_mw by the first ten unit pivots of row_hnf; "
            "read the final three residues modulo 2"
        ),
        "all_literal_doubling_verified_representatives_by_coset": by_coset,
        "canonical_doubled_targets": canonical_doubled_targets,
        "shortest_doubled_target_candidates_by_coset": short_target_candidates_by_coset,
        "hnf_coordinate_doubled_targets": hnf_coordinate_doubled_targets,
        "literal_complete_profile_shell_coset_keys": sorted(by_coset),
        "literal_complete_profile_shell_span_index_in_saturated_lattice": 2,
        "absent_nonzero_coset_keys": sorted(set(nonzero_cosets)-set(by_coset)),
    },
    "method": {
        "boundary": "the pre-existing complete sixteen-profile p=89 shell only",
        "class_filter": "every record in the declared complete profile shell",
        "acceptance": "literal [2]Q equality in GF(89)(u)",
        "current_word_optimization": "exact minimum L1 integer program over 45 strict-unique current sections",
    },
    "elapsed_seconds": time.monotonic()-started,
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
print(f"wrote {args.output}")
print(f"verified pairings: {len(verified)}")
print("literal complete-profile cosets:", sorted(by_coset))
for row in canonical_doubled_targets:
    print(
        "coset", "".join(map(str, row["quotient_key"])),
        "class", row["physical_class_index"],
        "L1", row["minimum_l1_norm"],
        "denominator degree", row["target_section_mod89"]["x"]["denominator_degree"],
    )

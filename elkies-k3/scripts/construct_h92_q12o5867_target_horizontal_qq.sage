#!/usr/bin/env sage
"""Construct the q12/o5867 target section by exact q8 group law.

This is deliberately a group-law-only gate.  It consumes the two regular
replacement lifts, the already certified Q2 and Q3 sections, and three exact
sections obtained by restricting the old vertical components.  It does not
perform a resolved Riemann--Roch or Jacobian calculation.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, vector


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corrections",
        type=Path,
        default=ROOT / "artifacts/local/elkies-k3/q8o376-correction-sections-qq-helper.json",
    )
    parser.add_argument(
        "--component9",
        type=Path,
        default=ROOT / "artifacts/local/elkies-k3/q8o376-old-a11-component9-section-qq.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "artifacts/local/elkies-k3/q12o5867-target-horizontal-qq.json",
    )
    return parser.parse_args()


args = parse_args()
paths = {
    "surface": ROOT / "artifacts/local/elkies-k3/q4o164-q8o376-smooth-rr-qq.json",
    "replacement": ROOT / "artifacts/local/elkies-k3/q12o5867-replacement-word-seeds-qq.json",
    "named": ROOT / "artifacts/local/elkies-k3/q12o5867-abel-trace-named-seeds-qq.json",
    "Q3": ROOT / "artifacts/local/elkies-k3/q12o5867-degree1-compiler-branch-qq.json",
    "shell89": ROOT / "artifacts/local/elkies-k3/q12o5867-p0-shell-all-records-mod89.json",
    "classifier": ROOT / "artifacts/local/elkies-k3/q12o5867-p0-shell-lattice-classification-mod89.json",
    "corrections": args.corrections.resolve(),
    "component9": args.component9.resolve(),
    "marking": ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json",
    "frame": ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt",
}
for path in paths.values():
    if not path.is_file():
        raise FileNotFoundError(path)

surface = read_json(paths["surface"])
replacement = read_json(paths["replacement"])
named = read_json(paths["named"])
q3_artifact = read_json(paths["Q3"])
shell89 = read_json(paths["shell89"])
classifier = read_json(paths["classifier"])
correction_artifact = read_json(paths["corrections"])
component9_artifact = read_json(paths["component9"])

R = PolynomialRing(QQ, "v")
v = R.gen()
K = R.fraction_field()
A = R(surface["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface["child"]["minimal_B_coefficients_low_to_high"])


def coordinate(record, name):
    if name in record and isinstance(record[name], dict):
        nested = record[name]
        if "numerator_coefficients_low_to_high" in nested:
            return K(R(nested["numerator_coefficients_low_to_high"])) / K(
                R(nested["denominator_coefficients_low_to_high"])
            )
    coefficient_key = f"{name}_coefficients_low_to_high"
    if coefficient_key in record:
        return K(R(record[coefficient_key]))
    numerator_key = f"{name}_numerator_coefficients_low_to_high"
    denominator_key = f"{name}_denominator_coefficients_low_to_high"
    if numerator_key in record and denominator_key in record:
        return K(R(record[numerator_key])) / K(R(record[denominator_key]))
    raise KeyError(f"no {name}-coordinate in section record")


def point(record):
    answer = coordinate(record, "x"), coordinate(record, "y")
    assert answer[1] ** 2 == answer[0] ** 3 + K(A) * answer[0] + K(B)
    return answer


def find_named_record(payload, name):
    if isinstance(payload, dict):
        if name in payload:
            candidate = payload[name]
            if isinstance(candidate, dict):
                return candidate.get("section", candidate)
        if payload.get("name") == name:
            return payload.get("section", payload)
        for value in payload.values():
            try:
                return find_named_record(value, name)
            except KeyError:
                pass
    elif isinstance(payload, list):
        for value in payload:
            try:
                return find_named_record(value, name)
            except KeyError:
                pass
    raise KeyError(name)


def neg(P):
    return None if P is None else (P[0], -P[1])


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1 ** 2 + K(A)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope ** 2 - x1 - x2
    y3 = slope * (x1 - x3) - y1
    assert y3 ** 2 == x3 ** 3 + K(A) * x3 + K(B)
    return x3, y3


def polynomial_coefficients(poly):
    return [str(value) for value in R(poly).list()]


def rational_record(value, prefix):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    return {
        f"{prefix}_numerator_coefficients_low_to_high": polynomial_coefficients(numerator),
        f"{prefix}_denominator_coefficients_low_to_high": polynomial_coefficients(denominator),
        f"degrees_{prefix}_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def maximum_rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(answer, abs(value.numerator()).nbits(), value.denominator().nbits())
    return int(answer)


def finite_pole_intersection(P):
    x, y = P
    dx = R(x.denominator())
    dy = R(y.denominator())
    total = ZZ(0)
    for factor, exponent in dx.factor():
        y_exponent = dy.valuation(factor)
        assert exponent % 2 == 0 and y_exponent == 3 * (exponent // 2)
        total += factor.degree() * (exponent // 2)
    return total


def intersection_with_zero(P):
    x, y = P
    finite = finite_pole_intersection(P)
    infinity_x = max(ZZ(0), ZZ(R(x.numerator()).degree() - R(x.denominator()).degree() - 4))
    infinity_y = max(ZZ(0), ZZ(R(y.numerator()).degree() - R(y.denominator()).degree() - 6))
    assert infinity_x % 2 == 0 and infinity_y % 3 == 0
    assert infinity_x // 2 == infinity_y // 3
    return int(finite + infinity_x // 2), int(finite), int(infinity_x // 2)


Zx = PolynomialRing(QQ, "z")
z = Zx.gen()
reducible_supports = []
reducible_nodes = []
for fibre in surface["child"]["finite_reducible_fibres"]:
    factor = R(fibre["factor"].replace("u", "v"))
    support = -factor[0] / factor[1]
    cubic = z ** 3 + A(support) * z + B(support)
    repeated = cubic.gcd(cubic.derivative())
    reducible_supports.append(support)
    reducible_nodes.append(-repeated[0] / repeated[1])
infinity_cubic = z ** 3 + A[A.degree()] * z + B[B.degree()]
infinity_repeated = infinity_cubic.gcd(infinity_cubic.derivative())
infinity_node = -infinity_repeated[0] / infinity_repeated[1]


def value_at_infinity(function, weight):
    function = K(function)
    numerator = R(function.numerator())
    denominator = R(function.denominator())
    difference = numerator.degree() - denominator.degree()
    if difference < weight:
        return QQ.zero()
    if difference == weight:
        return numerator.leading_coefficient() / denominator.leading_coefficient()
    return None


def component_profile(P):
    x, y = P
    answer = []
    for support, node in zip(reducible_supports, reducible_nodes):
        if x.denominator()(support) == 0 or y.denominator()(support) == 0:
            answer.append(0)
        else:
            answer.append(int(x(support) == node and y(support) == 0))
    infinity_x = value_at_infinity(x, 4)
    infinity_y = value_at_infinity(y, 6)
    answer.append(int(infinity_x == infinity_node and infinity_y == 0))
    return answer


def section_height(P):
    p_dot_o = intersection_with_zero(P)[0]
    return QQ(4 + 2 * p_dot_o) - QQ(sum(component_profile(P))) / 2


def height_pairing(P, Q):
    return (section_height(P) + section_height(Q) - section_height(add(P, neg(Q)))) / 2


section_records = {
    "class499_shell206": replacement["sections"]["replacement_class499_shell206"]["section"],
    "Q2_class500_shell116": named["sections"]["Q2_unique_rank12_seed"]["section"],
    "Q3_class69_shell90": q3_artifact["section"],
    "class511_shell172": replacement["sections"]["replacement_class511_shell172"]["section"],
}
correction_names = [
    "first_old_I6_I4_missing_component",
    "old_A11_component_5",
    "old_A11_component_9",
]
for name in correction_names[:2]:
    section_records[name] = find_named_record(correction_artifact, name)
section_records["old_A11_component_9"] = component9_artifact["section"]

points = {name: point(record) for name, record in section_records.items()}
orientation_anchor_names = [
    "class499_shell206", "Q2_class500_shell116",
    "Q3_class69_shell90", "class511_shell172",
]
orientation_diagnostics = {
    name: {
        "component_profile": component_profile(points[name]),
        "height": str(section_height(points[name])),
        "pairings_with_499_500_69_511": [
            str(height_pairing(points[name], points[anchor]))
            for anchor in orientation_anchor_names
        ],
    }
    for name in correction_names
}
orientation_diagnostics["Q2_Q3_anchor_gate"] = {
    "Q2_profile": component_profile(points["Q2_class500_shell116"]),
    "Q2_height": str(section_height(points["Q2_class500_shell116"])),
    "Q3_profile": component_profile(points["Q3_class69_shell90"]),
    "Q3_height": str(section_height(points["Q3_class69_shell90"])),
    "pairing": str(height_pairing(points["Q2_class500_shell116"], points["Q3_class69_shell90"])),
}
print(f"Q12O5867TARGETQQ|orientation_diagnostics={orientation_diagnostics}")
start = time.monotonic()
word = [
    (1, "class499_shell206"),
    (1, "Q2_class500_shell116"),
    (1, "Q3_class69_shell90"),
    (1, "class511_shell172"),
    (-1, "first_old_I6_I4_missing_component"),
    (1, "old_A11_component_5"),
    (-1, "old_A11_component_9"),
]
target = None
for coefficient, name in word:
    assert abs(coefficient) == 1
    target = add(target, points[name] if coefficient == 1 else neg(points[name]))
assert target is not None
assert target[1] ** 2 == target[0] ** 3 + K(A) * target[0] + K(B)

# The exact lattice identity is checked independently of the equation addition.
class_indices = [499, 500, 69, 511, 489, 487, 460]
classes = {row["class_index"]: row for row in classifier["lattice_shell"]["classes"]}
marking = read_json(paths["marking"])
frame_rows = [
    [ZZ(value) for value in line.split()]
    for line in paths["frame"].read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -matrix(ZZ, frame_rows))
source_in_basis = matrix(ZZ, marking["source_in_basis"])


def lattice_section(class_index):
    return vector(ZZ, classes[class_index]["q4o164_parent_curve"]) * source_in_basis


def lattice_height_pairing(left_index, right_index):
    left = classes[left_index]
    right = classes[right_index]
    local = sum(
        a * b for a, b in zip(left["current_component_pairings"], right["current_component_pairings"])
    )
    return QQ(2) - lattice_section(left_index) * gram * lattice_section(right_index) - QQ(local) / 2


anchor_indices = [499, 500, 69, 511]
equation_helper_classes = {}
for name in correction_names:
    # The equation finite-I2 order interchanges the first two marked A1 roots.
    profile = orientation_diagnostics[name]["component_profile"]
    marked_profile = [profile[1], profile[0], profile[2], profile[3]]
    equation_pairings = [QQ(value) for value in orientation_diagnostics[name]["pairings_with_499_500_69_511"]]
    matches = [
        index for index, row in classes.items()
        if row["q4o164_parent_degree"] == 0
        and row["current_component_pairings"] == marked_profile
        and lattice_height_pairing(index, index) == section_height(points[name])
        and [lattice_height_pairing(index, anchor) for anchor in anchor_indices] == equation_pairings
    ]
    assert len(matches) == 1
    equation_helper_classes[name] = matches[0]
assert equation_helper_classes == {
    "first_old_I6_I4_missing_component": 489,
    "old_A11_component_5": 933,
    "old_A11_component_9": 913,
}

lattice_sum = []
for coordinate_index in range(13):
    lattice_sum.append(
        sum(classes[index]["current_4A1_mw"][coordinate_index] for index in class_indices[:4])
        - classes[489]["current_4A1_mw"][coordinate_index]
        - classes[487]["current_4A1_mw"][coordinate_index]
        + classes[460]["current_4A1_mw"][coordinate_index]
    )
expected_lattice_target = [-11, 2, 13, 11, -6, -41, 6, 5, -14, 0, -3, -1, -2]
assert lattice_sum == expected_lattice_target
helper_lattice_sum = []
for coordinate_index in range(13):
    helper_lattice_sum.append(
        sum(classes[index]["current_4A1_mw"][coordinate_index] for index in class_indices[:4])
        - classes[489]["current_4A1_mw"][coordinate_index]
        + classes[933]["current_4A1_mw"][coordinate_index]
        - classes[913]["current_4A1_mw"][coordinate_index]
    )
assert helper_lattice_sum == lattice_sum

# Match the four replacement reductions to their independently enumerated p=89
# shell records.  The vertical correction sections need not belong to that
# bounded dominant-parent shell, so reduce their exact coordinates directly.
Fp = GF(89)
Rp = PolynomialRing(Fp, "v")
vp = Rp.gen()
Kp = Rp.fraction_field()
Ap = Rp([Fp(QQ(value).numerator()) / Fp(QQ(value).denominator()) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
Bp = Rp([Fp(QQ(value).numerator()) / Fp(QQ(value).denominator()) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])


def reduce_poly(poly):
    return Rp([Fp(value.numerator()) / Fp(value.denominator()) for value in R(poly).list()])


def reduce_coordinate(value):
    value = K(value)
    return Kp(reduce_poly(value.numerator())) / Kp(reduce_poly(value.denominator()))


def reduce_point(P):
    answer = reduce_coordinate(P[0]), reduce_coordinate(P[1])
    assert answer[1] ** 2 == answer[0] ** 3 + Kp(Ap) * answer[0] + Kp(Bp)
    return answer


def add_modp(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1 ** 2 + Kp(Ap)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope ** 2 - x1 - x2
    y3 = slope * (x1 - x3) - y1
    assert y3 ** 2 == x3 ** 3 + Kp(Ap) * x3 + Kp(Bp)
    return x3, y3


shell_points = {
    index: (Kp(Rp(row["x_coefficients_low_to_high"])), Kp(Rp(row["y_coefficients_low_to_high"])))
    for index, row in enumerate(shell89["all_records"])
}
main_shells = {
    "class499_shell206": 206,
    "Q2_class500_shell116": 116,
    "Q3_class69_shell90": 90,
    "class511_shell172": 172,
}
for name, shell_index in main_shells.items():
    assert reduce_point(points[name]) == shell_points[shell_index]

correction_shells = {}
for name in correction_names:
    reduced = reduce_point(points[name])
    matches = [index for index in [16, 29, 124, 160] if shell_points[index] == reduced]
    correction_shells[name] = matches[0] if len(matches) == 1 else None

modular_target = None
modular_word = [
    (1, shell_points[206]), (1, shell_points[116]),
    (1, shell_points[90]), (1, shell_points[172]),
    (-1, reduce_point(points["first_old_I6_I4_missing_component"])),
    (1, reduce_point(points["old_A11_component_5"])),
    (-1, reduce_point(points["old_A11_component_9"])),
]
for coefficient, P in modular_word:
    modular_target = add_modp(modular_target, P if coefficient == 1 else (P[0], -P[1]))
assert modular_target == reduce_point(target)

p_dot_o, finite_p_dot_o, infinity_p_dot_o = intersection_with_zero(target)
assert p_dot_o == 10, (p_dot_o, finite_p_dot_o, infinity_p_dot_o)
target_component_profile = component_profile(target)
target_height = section_height(target)
assert target_component_profile == [1, 1, 1, 1]
assert target_height == 22
coordinate_payload = {}
coordinate_payload.update(rational_record(target[0], "x"))
coordinate_payload.update(rational_record(target[1], "y"))
coordinate_payload.update({
    "exact_weierstrass_identity": True,
    "P_dot_O": p_dot_o,
    "finite_P_dot_O": finite_p_dot_o,
    "infinity_P_dot_O": infinity_p_dot_o,
    "equation_component_profile": target_component_profile,
    "height": str(target_height),
    "maximum_rational_bits": maximum_rational_bits(
        list(R(target[0].numerator())) + list(R(target[0].denominator()))
        + list(R(target[1].numerator())) + list(R(target[1].denominator()))
    ),
})

payload = {
    "schema": "h92-q12o5867-target-horizontal-qq-v1",
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/construct_h92_q12o5867_target_horizontal_qq.sage "
        "--output artifacts/local/elkies-k3/q12o5867-target-horizontal-qq.json"
    ),
    "inputs": {name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for name, path in paths.items()},
    "method": {
        "construction": "exact affine group law over QQ(v)",
        "word": [{"coefficient": coefficient, "section": name} for coefficient, name in word],
        "no_groebner_or_elimination": True,
        "stopped_before_resolved_RR_or_Jacobian": True,
    },
    "target": coordinate_payload,
    "mod89": {
        "exact_reduction_matches_mixed_shell_and_exact_correction_group_sum": True,
        "dominant_replacement_shell_indices": main_shells,
        "correction_shell_indices": correction_shells,
        "correction_reductions_checked_on_child_equation": True,
    },
    "lattice": {
        "exact_MW_word_identity": True,
        "named_class_word": "+499+500+69+511-489-487+460",
        "equation_helper_class_word": "+499+500+69+511-489+933-913",
        "equation_helper_classification": equation_helper_classes,
        "helper_and_named_correction_words_have_identical_MW_tail": True,
        "helper_orientation_diagnostics": orientation_diagnostics,
        "current_4A1_mw_target": lattice_sum,
    },
    "runtime_seconds": time.monotonic() - start,
    "proof_boundary": "Exact section/equation, mod-89, P.O, and marked MW-word gate only; no resolved RR or Jacobian compiler is attempted here.",
    "status": "PASS_EXACT_QQ_Q12O5867_TARGET_HORIZONTAL_SECTION",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
with open(args.output, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
print(
    "Q12O5867TARGETQQ|"
    f"degrees_x={coordinate_payload['degrees_x_numerator_denominator']}|"
    f"degrees_y={coordinate_payload['degrees_y_numerator_denominator']}|"
    f"P.O={p_dot_o}|correction_shells={correction_shells}|"
    f"runtime={payload['runtime_seconds']:.3f}|status={payload['status']}|output={args.output.resolve()}"
)

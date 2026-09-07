#!/usr/bin/env sage -python
"""Replay the three short duplication targets exactly over QQ(u)."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
TARGETS = LOCAL / "q12o5867-two-primary-cosets-mod89.json"
SUPPORTS = (
    LOCAL / "q12o5867-two-primary-target-support-sections-v2-qq.json",
    LOCAL / "q12o5867-support-class170-shell32-qq.json",
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q12o5867-three-target-halvings-qq.json",
)
args = parser.parse_args()
output = args.output if args.output.is_absolute() else ROOT / args.output
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
targets = json.loads(TARGETS.read_text())
R = PolynomialRing(QQ, "u")
K = R.fraction_field()
S = PolynomialRing(K, "X")
X = S.gen()
A = K(R(q8["child"]["minimal_A_coefficients_low_to_high"]))
B = K(R(q8["child"]["minimal_B_coefficients_low_to_high"]))


def section_point(section):
    point = (
        K(R(section["x_coefficients_low_to_high"])),
        K(R(section["y_coefficients_low_to_high"])),
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


def function_record(value):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    values = list(numerator.list())+list(denominator.list())
    return {
        "numerator_coefficients_low_to_high": list(map(str, numerator.list())),
        "denominator_coefficients_low_to_high": list(map(str, denominator.list())),
        "numerator_degree": int(numerator.degree()),
        "denominator_degree": int(denominator.degree()),
        "maximum_coefficient_bits": max(
            max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
            for value in values
        ),
    }


support_points = {}
support_metadata = {}
for path in SUPPORTS:
    data = json.loads(path.read_text())
    for key, record in data["sections"].items():
        class_index = int(record["lattice_class_index"])
        assert class_index not in support_points
        support_points[class_index] = section_point(record["section"])
        support_metadata[class_index] = {
            "artifact": str(path.relative_to(ROOT)),
            "key": key,
            "selection": record["selection"],
            "full_NS_name_certified_by_lift_artifact": record["full_NS_name_certified"],
        }

attempts = []
for target in targets["quotient"]["canonical_doubled_targets"]:
    target_started = time.monotonic()
    point = None
    support_records = []
    for term in target["minimum_l1_current_section_word"]:
        class_index = int(term["physical_class_index"])
        coefficient = ZZ(term["coefficient"])
        assert class_index in support_points
        point = point_add(point, point_mul(coefficient, support_points[class_index]))
        support_records.append({
            "coefficient": int(coefficient),
            "supplied_physical_class_index": class_index,
            **support_metadata[class_index],
        })
    assert point is not None
    x_target, y_target = point
    assert y_target**2 == x_target**3+A*x_target+B
    duplication = (
        X**4-2*A*X**2-8*B*X+A**2
        -4*x_target*(X**3+A*X+B)
    )
    factorization = duplication.factor()
    factor_degrees = []
    verified_halves = []
    for factor, multiplicity in factorization:
        factor_degrees.extend([int(factor.degree())]*int(multiplicity))
        if factor.degree() != 1:
            continue
        x_half = -factor[0]/factor[1]
        rhs = x_half**3+A*x_half+B
        if not rhs.is_square():
            continue
        y_half = rhs.sqrt()
        for signed_y in (y_half, -y_half):
            candidate = (x_half, signed_y)
            if point_add(candidate, candidate) == point:
                verified_halves.append({
                    "x": function_record(x_half),
                    "y": function_record(signed_y),
                    "literal_curve_substitution": True,
                    "literal_doubling_verified": True,
                })
    factor_degrees.sort()
    attempts.append({
        "declared_mod89_half_target_class_index": target["physical_class_index"],
        "declared_quotient_key": target["quotient_key"],
        "declared_doubled_target_current_4A1_mw": target[
            "doubled_target_current_4A1_mw"
        ],
        "exact_QQ_target_support_word": support_records,
        "exact_QQ_target_section": {
            "x": function_record(x_target),
            "y": function_record(y_target),
            "literal_curve_substitution": True,
        },
        "duplication_polynomial_degree": int(duplication.degree()),
        "irreducible_factor_degrees_with_multiplicity": factor_degrees,
        "linear_factor_count": factor_degrees.count(1),
        "verified_rational_halves": verified_halves,
        "rational_half_exists_over_QQ_u": bool(verified_halves),
        "elapsed_seconds": time.monotonic()-target_started,
    })
    print(
        "class", target["physical_class_index"],
        "factors", factor_degrees,
        "halves", len(verified_halves),
        "seconds", attempts[-1]["elapsed_seconds"],
    )

payload = {
    "schema": "q12o5867-three-target-halvings-qq-v1",
    "status": "PASS_EXACT_THREE_DECLARED_DUPLICATION_FACTORIZATIONS_OVER_QQ_U",
    "inputs": {
        "q8": {"path": str(Q8.relative_to(ROOT)), "sha256": sha256(Q8)},
        "targets": {"path": str(TARGETS.relative_to(ROOT)), "sha256": sha256(TARGETS)},
        "support_artifacts": [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for path in SUPPORTS
        ],
    },
    "attempt_count": len(attempts),
    "attempts": attempts,
    "proof_boundary": (
        "Exactly three short current-section words are tested over QQ(u). The modular "
        "class labels are not promoted beyond their separately recorded cross-prime audit."
    ),
    "elapsed_seconds": time.monotonic()-started,
}
assert len(attempts) == 3
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(f"wrote {output}")

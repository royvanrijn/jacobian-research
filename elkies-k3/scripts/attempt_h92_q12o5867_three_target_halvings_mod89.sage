#!/usr/bin/env sage -python
"""Factor the three declared p=89 duplication equations over GF(89)(u)."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
TARGETS = LOCAL / "q12o5867-two-primary-cosets-mod89.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--targets", type=Path, default=TARGETS)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q12o5867-three-target-halvings-mod89.json",
)
args = parser.parse_args()
targets_path = args.targets if args.targets.is_absolute() else ROOT / args.targets
output = args.output if args.output.is_absolute() else ROOT / args.output
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
targets = json.loads(targets_path.read_text())
prime = ZZ(targets["prime"])
assert prime == 89
F = GF(prime)
R = PolynomialRing(F, "u")
K = R.fraction_field()
S = PolynomialRing(K, "X")
X = S.gen()


def reduce_qq(value):
    value = QQ(value)
    assert value.denominator() % prime
    return F(value.numerator())/F(value.denominator())


def read_function(record):
    return K(R(record["numerator_coefficients_low_to_high"])) / K(
        R(record["denominator_coefficients_low_to_high"])
    )


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


A = K(R([reduce_qq(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]]))
B = K(R([reduce_qq(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]]))


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
    return x3, slope*(x1-x3)-y1


def point_double(point):
    return point_add(point, point)


attempts = []
for target in targets["quotient"]["canonical_doubled_targets"]:
    target_started = time.monotonic()
    target_point = (
        read_function(target["target_section_mod89"]["x"]),
        read_function(target["target_section_mod89"]["y"]),
    )
    assert target_point[1]**2 == target_point[0]**3+A*target_point[0]+B
    x_target = target_point[0]
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
            if point_double(candidate) == target_point:
                verified_halves.append({
                    "x": function_record(x_half),
                    "y": function_record(signed_y),
                    "literal_doubling_verified": True,
                })
    factor_degrees.sort()
    attempts.append({
        "physical_half_target_class_index": target["physical_class_index"],
        "quotient_key": target["quotient_key"],
        "doubled_target_current_4A1_mw": target["doubled_target_current_4A1_mw"],
        "target_current_section_word": target["minimum_l1_current_section_word"],
        "target_section_mod89": target["target_section_mod89"],
        "duplication_polynomial_degree": int(duplication.degree()),
        "irreducible_factor_degrees_with_multiplicity": factor_degrees,
        "linear_factor_count": factor_degrees.count(1),
        "verified_rational_halves": verified_halves,
        "rational_half_exists_over_GF89_u": bool(verified_halves),
        "elapsed_seconds": time.monotonic()-target_started,
    })

payload = {
    "schema": "q12o5867-three-target-halvings-mod89-v1",
    "status": "PASS_EXACT_THREE_DECLARED_DUPLICATION_FACTORIZATIONS_OVER_GF89_U",
    "prime": int(prime),
    "inputs": {
        "q8": {"path": str(Q8.relative_to(ROOT)), "sha256": sha256(Q8)},
        "targets": {"path": str(targets_path.relative_to(ROOT)),
                    "sha256": sha256(targets_path)},
    },
    "attempt_count": len(attempts),
    "attempts": attempts,
    "proof_boundary": (
        "Exactly the three canonical independent doubled targets are factored. "
        "No unrelated section or parent-fibre enumeration is performed."
    ),
    "elapsed_seconds": time.monotonic()-started,
}
assert len(attempts) == 3
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(f"wrote {output}")
for row in attempts:
    print(
        "coset", "".join(map(str, row["quotient_key"])),
        "class", row["physical_half_target_class_index"],
        "factors", row["irreducible_factor_degrees_with_multiplicity"],
        "halves", len(row["verified_rational_halves"]),
    )

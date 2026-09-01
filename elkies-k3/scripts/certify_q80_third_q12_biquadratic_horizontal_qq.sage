#!/usr/bin/env sage -python
"""Compose the exact u=-2 third-q12 horizontal in its biquadratic algebra."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, QQ, ZZ, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--operands", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--surface", type=Path,
    default=RESULTS / "q80-fixed-u-minus2-p19-height-shell-with-po1.json",
)
parser.add_argument(
    "--held-out-surface", type=Path,
    default=RESULTS / "q80-fixed-u-minus2-p71-heldout-good-reduction.json",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-horizontal-qq.json",
)
args = parser.parse_args()
for name in ("operands", "surface", "held_out_surface", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


started = time.monotonic()
operands = json.loads(args.operands.read_text())
surface = json.loads(args.surface.read_text())
held_out = json.loads(args.held_out_surface.read_text())
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("biquadratic operands are not exact-certified")
if surface.get("schema") != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected exact surface schema")
if held_out.get("schema") != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected held-out surface schema")

q1 = read_rational(operands["biquadratic_field"]["q1"])
q2 = read_rational(operands["biquadratic_field"]["q2"])
operand1, operand2 = operands["operands"]
x1_coefficients = list(map(read_rational, operand1["x_coefficients_low_to_high"]))
x2_coefficients = list(map(read_rational, operand2["x_coefficients_low_to_high"]))
if x1_coefficients[-1] != q1 or x2_coefficients[-1] != q2:
    raise ArithmeticError("operand leading coefficients do not equal l^2")

W_ring = PolynomialRing(QQ, "W")
W = W_ring.gen()
equation = surface["parameters"][0]["exact_equations"]["second_q4"]
A = W_ring(equation["A_coefficients_low_to_high"])
B = W_ring(equation["B_coefficients_low_to_high"])
x1 = W_ring(x1_coefficients)
x2 = W_ring(x2_coefficients)


def y_multiplier(x_value, q_value):
    square = x_value**3 + A * x_value + B
    quotient = square / q_value
    coefficients = [QQ.zero() for unused in range(7)]
    coefficients[6] = q_value
    for degree in range(11, 5, -1):
        index = degree - 6
        partial = W_ring(coefficients)
        coefficients[index] = (
            quotient[degree] - (partial**2)[degree]
        ) / (2 * coefficients[6])
    answer = W_ring(coefficients)
    if answer**2 != quotient:
        raise ArithmeticError("closure operand y multiplier is not an exact square root")
    return answer


u1 = y_multiplier(x1, q1)
u2 = y_multiplier(x2, q2)
function_field = W_ring.fraction_field()

# Basis order is 1,a,b,ab with a^2=q1 and b^2=q2.  The addition formula
# needs only rational-function division because x2-x1 is rational.
def bq_add(left, right):
    return tuple(left[index] + right[index] for index in range(4))


def bq_neg(value):
    return tuple(-entry for entry in value)


def bq_sub(left, right):
    return bq_add(left, bq_neg(right))


def bq_mul(left, right):
    x0, x1c, x2c, x3c = left
    y0, y1c, y2c, y3c = right
    return (
        x0*y0 + q1*x1c*y1c + q2*x2c*y2c + q1*q2*x3c*y3c,
        x0*y1c + x1c*y0 + q2*(x2c*y3c + x3c*y2c),
        x0*y2c + x2c*y0 + q1*(x1c*y3c + x3c*y1c),
        x0*y3c + x3c*y0 + x1c*y2c + x2c*y1c,
    )


zero = function_field.zero()
X1 = (function_field(x1), zero, zero, zero)
X2 = (function_field(x2), zero, zero, zero)
Y1 = (zero, function_field(u1), zero, zero)
Y2 = (zero, zero, function_field(u2), zero)
denominator = function_field(x2 - x1)
slope = (
    zero,
    -function_field(u1) / denominator,
    function_field(u2) / denominator,
    zero,
)
X_horizontal = bq_sub(bq_sub(bq_mul(slope, slope), X1), X2)
Y_horizontal = bq_sub(bq_mul(slope, bq_sub(X1, X_horizontal)), Y1)
if X_horizontal[1] or X_horizontal[2] or Y_horizontal[0] or Y_horizontal[3]:
    raise ArithmeticError("composed horizontal has the wrong biquadratic support")

curve_rhs = bq_add(
    bq_add(bq_mul(bq_mul(X_horizontal, X_horizontal), X_horizontal), tuple(function_field(A)*entry for entry in X_horizontal)),
    (function_field(B), zero, zero, zero),
)
if bq_mul(Y_horizontal, Y_horizontal) != curve_rhs:
    raise ArithmeticError("exact biquadratic horizontal fails the Weierstrass equation")


def rational_function_record(value):
    value = function_field(value)
    return {
        "numerator_coefficients_low_to_high": [str(coefficient) for coefficient in value.numerator().list()],
        "denominator_coefficients_low_to_high": [str(coefficient) for coefficient in value.denominator().list()],
        "degrees_numerator_denominator": [
            int(value.numerator().degree()), int(value.denominator().degree())
        ],
    }


def reduce_rational(value, finite):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ZeroDivisionError("exact horizontal denominator vanishes at held-out prime")
    return finite(value.numerator()) / denominator


def reduce_function(value, finite, target_ring):
    value = function_field(value)
    numerator = target_ring([reduce_rational(coefficient, finite) for coefficient in value.numerator().list()])
    denominator = target_ring([reduce_rational(coefficient, finite) for coefficient in value.denominator().list()])
    return target_ring.fraction_field()(numerator / denominator)


held_parameter = held_out["parameters"][0]
held_modular = held_parameter["modular"][0]
if held_parameter["u"] != "-2" or held_modular.get("status") != "PASS_GOOD_REDUCTION_AUDIT":
    raise ValueError("held-out input is not the audited u=-2 good reduction")
prime = int(held_modular["prime"])
if prime != 71:
    raise ValueError("this certificate reserves p=71 as the held-out prime")
prime_field = GF(prime)
finite = GF(prime**2, "r")
base = PolynomialRing(finite, "W")
held_equation = held_parameter["exact_equations"]["second_q4"]
held_A = base([reduce_rational(value, finite) for value in held_equation["A_coefficients_low_to_high"]])
held_B = base([reduce_rational(value, finite) for value in held_equation["B_coefficients_low_to_high"]])
held_function = base.fraction_field()
held_curve = EllipticCurve(held_function, [0, 0, 0, held_function(held_A), held_function(held_B)])
a0 = finite(reduce_rational(q1, finite)).sqrt()
b0 = finite(reduce_rational(q2, finite)).sqrt()


def instantiate(value, a_value, b_value):
    components = [reduce_function(entry, finite, base) for entry in value]
    return held_function(
        components[0]
        + a_value*components[1]
        + b_value*components[2]
        + a_value*b_value*components[3]
    )


def po_from_x(point):
    x_value = point[0]
    twice_intersection = max(
        x_value.denominator().degree(),
        x_value.numerator().degree() - 4,
    )
    if twice_intersection < 0 or twice_intersection % 2:
        raise ArithmeticError("invalid held-out section pole divisor")
    return int(twice_intersection // 2)


delta = 4 * held_A**3 + 27 * held_B**2
star_factor = next(factor.monic() for factor, exponent in delta.factor() if int(exponent) == 7)
star_root = -star_factor[0] / star_factor[1]
node_ring = PolynomialRing(finite, "Xnode")
Xnode = node_ring.gen()
node_cubic = Xnode**3 + held_A(star_root)*Xnode + held_B(star_root)
singular_roots = node_cubic.gcd(node_cubic.derivative()).roots(multiplicities=False)
if len(singular_roots) != 1:
    raise ArithmeticError("held-out I1* node is not unique")
singular_x = singular_roots[0]

held_out_candidates = []
unsigned_keys = set()
for sign_a in (1, -1):
    for sign_b in (1, -1):
        a_value = sign_a * a0
        b_value = sign_b * b0
        point = held_curve(
            instantiate(X_horizontal, a_value, b_value),
            instantiate(Y_horizontal, a_value, b_value),
        )
        fourth = 4 * point
        eighth = 2 * fourth
        height4 = QQ(4 + 2*po_from_x(fourth)) / 16
        height8 = QQ(4 + 2*po_from_x(eighth)) / 64
        x_value, y_value = point[0], point[1]
        target = (
            po_from_x(point) == 2
            and height4 == height8 == 8
            and x_value.denominator().degree() == 4
            and y_value.denominator().degree() == 6
            and x_value.numerator().degree() - x_value.denominator().degree() == 4
            and y_value.numerator().degree() - y_value.denominator().degree() == 6
            and x_value.denominator()(star_root)
            and x_value(star_root) != singular_x
        )
        key = min(str(point), str(-point))
        unsigned_keys.add(key)
        held_out_candidates.append({
            "sign_a": sign_a,
            "sign_b": sign_b,
            "P_dot_O": po_from_x(point),
            "height_from_fourth": str(height4),
            "height_from_eighth": str(height8),
            "finite_I1star_identity": bool(
                x_value.denominator()(star_root) and x_value(star_root) != singular_x
            ),
            "target_profile": bool(target),
        })
if not all(record["target_profile"] for record in held_out_candidates):
    raise ArithmeticError("an exact biquadratic conjugate fails the held-out target profile")
if len(unsigned_keys) != 2:
    raise ArithmeticError("held-out conjugates do not form two unsigned classes")

output = {
    "schema": "elkies-k3.q80-third-q12-biquadratic-horizontal-qq.v1",
    "status": "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_HORIZONTAL_HELDOUT_P71",
    "specialization": {"u": "-2"},
    "field": {
        "basis": ["1", "a", "b", "a*b"],
        "relations": [f"a^2=({q1})", f"b^2=({q2})"],
        "degree": 4,
        "horizontal_x_subfield": "QQ(a*b)",
        "horizontal_x_subfield_square_class": "q1*q2",
    },
    "horizontal": {
        "construction": "P1+P2",
        "x_basis_coefficients": [rational_function_record(value) for value in X_horizontal],
        "y_basis_coefficients": [rational_function_record(value) for value in Y_horizontal],
        "support": {
            "x": ["1", "a*b"],
            "y": ["a", "b"],
        },
        "exact_weierstrass_substitution": True,
    },
    "held_out_p71": {
        "surface_path": str(args.held_out_surface.relative_to(ROOT)),
        "surface_sha256": sha256(args.held_out_surface),
        "q1_character": -1 if not prime_field(reduce_rational(q1, prime_field)).is_square() else 1,
        "q2_character": -1 if not prime_field(reduce_rational(q2, prime_field)).is_square() else 1,
        "q1q2_character": -1 if not prime_field(reduce_rational(q1*q2, prime_field)).is_square() else 1,
        "signed_conjugates_tested": len(held_out_candidates),
        "unsigned_classes": len(unsigned_keys),
        "candidates": held_out_candidates,
        "all_target_profile": True,
    },
    "inputs": {
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
        "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "exact characteristic-zero sum of the two closure operands in a degree-four biquadratic algebra",
            "literal characteristic-zero Weierstrass substitution",
            "x-coordinate descent to the third quadratic subfield QQ(sqrt(q1*q2))",
            "held-out p=71 replay of all four signed conjugates as two unsigned target-profile horizontals",
        ],
        "not_proved": [
            "the characteristic-zero resolved pencil or child Jacobian",
            "characteristic-zero birational maps, fibre marking, or Mordell--Weil rank",
        ],
    },
    "runtime_seconds": time.monotonic() - started,
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_q80_third_q12_biquadratic_horizontal_qq.sage "
        f"--output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12BIQUADRATICHORIZONTAL|field_degree=4|x_subfield_degree=2|"
    f"heldout=71|signed={len(held_out_candidates)}|unsigned={len(unsigned_keys)}|"
    "status=PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_HORIZONTAL_HELDOUT_P71",
    flush=True,
)

#!/usr/bin/env sage -python
"""Compress the exact q4/orbit164 equation by a weighted base rescaling.

The expanded q4/orbit164 model is arithmetically large because its finite
reducible-fibre supports are written in a badly scaled base coordinate.  The
two finite I2 supports have a small rational ratio.  Put the second one at
``t=1`` and normalize the nodal x-coordinate at infinity to ``3``.  This is
an exact Weierstrass isomorphism; no elimination or Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
POINTING = LOCAL / "q4o164-c8-equation-marking-qq.json"
OUTPUT = LOCAL / "q4o164-compact-weierstrass-qq.json"

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def polynomial_bits(polynomial):
    return max(coefficient_bits(value) for value in polynomial)


def polynomial_record(polynomial):
    return [str(value) for value in polynomial.list()]


model = json.loads(MODEL.read_text())
pointing = json.loads(POINTING.read_text())
assert model["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"

RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
RT = PolynomialRing(QQ, "t")
t = RT.gen()
KT = RT.fraction_field()

A_old = RU([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B_old = RU([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
records = model["child"]["finite_reducible_fibres"]
assert [record["kodaira"] for record in records] == ["I2", "I2", "I4"]
supports_old = [QQ(record["support"]) for record in records]
assert supports_old[2] == 0

# U = c*t puts the second finite I2 at t=1.  At infinity the old nodal
# cubic has repeated root r=3*m^2.  With s=c^2*m, the weighted change
# x_old=s^2*x_new, y_old=s^3*y_new preserves a polynomial K3 model and
# makes the new infinity cubic x^3-27*x+54.
c = supports_old[1]
r_infinity = -3 * B_old[12] / (2 * A_old[8])
assert A_old[8] == -3 * r_infinity**2
assert B_old[12] == 2 * r_infinity**3
assert (r_infinity / 3).is_square()
m = (r_infinity / 3).sqrt()
s = c**2 * m

A = RT(A_old(c * t)) / s**4
B = RT(B_old(c * t)) / s**6
assert A.degree() == 8 and B.degree() == 12
assert A[8] == -27 and B[12] == 54

Delta_old = -16 * (4 * A_old**3 + 27 * B_old**2)
Delta = -16 * (4 * A**3 + 27 * B**2)
assert Delta == RT(Delta_old(c * t)) / s**12
assert Delta.degree() == 20

supports = [value / c for value in supports_old]
assert supports[1:] == [1, 0]
for support, order in zip(supports, (2, 2, 4)):
    shifted = Delta(t + support)
    assert shifted.valuation() == order
nodal_factor = Delta
for support, order in zip(supports, (2, 2, 4)):
    nodal_factor //= (t - support)**order
assert nodal_factor.degree() == 12
assert nodal_factor.gcd(nodal_factor.derivative()).degree() == 0
assert all(nodal_factor(support) for support in supports)


def load_rational(record):
    numerator = RU([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = RU([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return KU(numerator) / KU(denominator)


def transform_rational(function, weight):
    numerator = RT(function.numerator()(c * t))
    denominator = RT(function.denominator()(c * t))
    return KT(numerator) / (KT(denominator) * s**weight)


def rational_record(function):
    numerator = RT(function.numerator())
    denominator = RT(function.denominator())
    return {
        "numerator_coefficients_low_to_high": polynomial_record(numerator),
        "denominator_coefficients_low_to_high": polynomial_record(denominator),
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
        "maximum_rational_bits": max(polynomial_bits(numerator), polynomial_bits(denominator)),
    }


old_point = pointing["opposite_constant_support_section"]
x_point = transform_rational(load_rational(old_point["x"]), 2)
y_point = transform_rational(load_rational(old_point["y"]), 3)
assert y_point**2 == x_point**3 + A * x_point + B

payload = {
    "schema": "elkies-k3.q4o164-compact-weierstrass-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION",
    "compact_model": {
        "equation": "y^2 = x^3 + A(t)*x + B(t)",
        "A_coefficients_low_to_high": polynomial_record(A),
        "B_coefficients_low_to_high": polynomial_record(B),
        "degrees_A_B_Delta": [8, 12, 20],
        "maximum_A_rational_bits": polynomial_bits(A),
        "maximum_B_rational_bits": polynomial_bits(B),
        "reducible_fibres": [
            {"kodaira": record["kodaira"], "support": str(support)}
            for record, support in zip(records, supports)
        ] + [{"kodaira": "I4", "support": "infinity"}],
        "finite_nodal_factor_degree": 12,
        "finite_nodal_factor_squarefree": True,
        "ADE": "2A3+2A1",
        "MW_rank_if_rho19": 9,
    },
    "exact_coordinate_change": {
        "base": "U = c*t",
        "x": "x_old = s^2*x",
        "y": "y_old = s^3*y",
        "c": str(c),
        "m": str(m),
        "s": str(s),
        "relations": ["s=c^2*m", "r_infinity=3*m^2"],
        "old_infinity_repeated_root": str(r_infinity),
        "new_infinity_repeated_root": "3",
    },
    "small_moduli": {
        "second_I2_support": "1",
        "first_I2_support": str(supports[0]),
        "first_I2_support_expected": "25281/168246841",
        "support_identity": supports[0] == QQ(25281) / QQ(168246841),
    },
    "transported_exact_section": {
        "source_label": "opposite_constant_support_section",
        "x": rational_record(x_point),
        "y": rational_record(y_point),
        "exact_compact_model_identity": True,
    },
    "coefficient_growth": {
        "old_maximum_A_B_rational_bits": int(model["child"]["maximum_A_B_rational_bits"]),
        "new_maximum_A_rational_bits": polynomial_bits(A),
        "new_maximum_B_rational_bits": polynomial_bits(B),
        "transported_section_x_rational_bits": rational_record(x_point)["maximum_rational_bits"],
        "transported_section_y_rational_bits": rational_record(y_point)["maximum_rational_bits"],
    },
    "method": {
        "large_Groebner_required": False,
        "polynomial_factorization_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This is an exact QQ Weierstrass isomorphism of the already-certified q4/orbit164 "
        "model, with the discriminant transformation, reducible-fibre orders, residual "
        "squarefreeness, and one transported exact section checked literally. It does not "
        "yet recover the promoted q8/orbit376 horizontal section."
    ),
    "next_required": (
        "Extend the exact rank-8 pair-node section subgroup by one zero- or one-node "
        "integral direction, then synthesize the q8/orbit376 horizontal by group law."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, POINTING)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, POINTING)
        },
    },
}
assert payload["small_moduli"]["support_identity"]
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164COMPACT|A_bits={}|B_bits={}|section_bits={}/{}|status={}|output={}".format(
        payload["coefficient_growth"]["new_maximum_A_rational_bits"],
        payload["coefficient_growth"]["new_maximum_B_rational_bits"],
        payload["coefficient_growth"]["transported_section_x_rational_bits"],
        payload["coefficient_growth"]["transported_section_y_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)

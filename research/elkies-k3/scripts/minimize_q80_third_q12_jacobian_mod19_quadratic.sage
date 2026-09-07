#!/usr/bin/env sage -python
"""Minimize and classify the interpolated p=19 third-q12 Jacobian."""

import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"
EXPECTED = "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(INPUT.read_text())
if payload.get("status") != EXPECTED:
    raise ValueError("interpolated Jacobian is not certified")

base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
polynomial_ring = PolynomialRing(finite, "V")
V = polynomial_ring.gen()
function_field = polynomial_ring.fraction_field()


def element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


def polynomial(coordinates):
    return polynomial_ring([element(value) for value in coordinates])


def rational(record):
    return function_field(
        polynomial(record["numerator_coefficients_low_to_high_1_r"])
        / polynomial(record["denominator_coefficients_low_to_high_1_r"])
    )


a1, a2, a3, a4, a6 = [rational(payload["weierstrass"][name]) for name in ("a1", "a2", "a3", "a4", "a6")]
b2 = a1**2 + 4 * a2
b4 = a1 * a3 + 2 * a4
b6 = a3**2 + 4 * a6
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
delta = rational(payload["discriminant"])
if c4**3 - c6**2 != finite(1728) * delta:
    raise ArithmeticError("long Weierstrass invariant identity failed")

# The Laurent gauge has one pole divisor d with ord_d(Delta)=-18.  The
# admissible integral scaling is q=d^2, producing the missing I6 fibre.
denominator_factorization = tuple(delta.denominator().factor())
if len(denominator_factorization) != 1 or int(denominator_factorization[0][1]) != 18:
    raise ArithmeticError("unexpected Laurent-gauge discriminant denominator")
d = denominator_factorization[0][0].monic()
q = d**2
A = -c4 / 48
B = -c6 / 864
A_minimal = function_field(q**4 * A)
B_minimal = function_field(q**6 * B)
delta_minimal = function_field(q**12 * delta)
if A_minimal.denominator() != 1 or B_minimal.denominator() != 1:
    raise ArithmeticError("gauge-pole scaling did not produce an integral model")
A_minimal = A_minimal.numerator()
B_minimal = B_minimal.numerator()
if (A_minimal.degree(), B_minimal.degree()) != (8, 12):
    raise ArithmeticError("minimal K3 coefficient degrees are not 8 and 12")
expected_delta = -finite(16) * (4 * A_minimal**3 + 27 * B_minimal**2)
if delta_minimal.denominator() != 1 or delta_minimal.numerator() != expected_delta:
    raise ArithmeticError("minimal discriminant replay failed")
delta_minimal = delta_minimal.numerator()
if delta_minimal.degree() != 24:
    raise ArithmeticError("minimal K3 discriminant does not have degree 24")

fibres = []
root_counts = {2: 0, 4: 0, 6: 0}
for factor, exponent in delta_minimal.factor():
    exponent = int(exponent)
    if c4.numerator() % factor == 0:
        # Test the minimal c4, since the original rational c4 has a pole at d.
        c4_minimal = -finite(48) * A_minimal
        if c4_minimal % factor == 0:
            raise ArithmeticError("nonmultiplicative finite fibre encountered")
    degree = int(factor.degree())
    kind = f"I{exponent}"
    fibres.append(
        {
            "factor": str(factor.monic()),
            "degree": degree,
            "ord_discriminant": exponent,
            "kodaira": kind,
            "root": None if exponent == 1 else f"A{exponent - 1}",
        }
    )
    if exponent in root_counts:
        root_counts[exponent] += degree

if root_counts != {2: 3, 4: 1, 6: 1}:
    raise ArithmeticError(f"wrong reducible-fibre multiplicities {root_counts}")
if sum(record["degree"] * record["ord_discriminant"] for record in fibres) != 24:
    raise ArithmeticError("finite fibre Euler sum is not 24")


def coordinates(value):
    result = list(finite(value).list())
    result += [base_finite.zero()] * (2 - len(result))
    return [int(result[0]), int(result[1])]


def polynomial_record(value):
    return [coordinates(coefficient) for coefficient in value.list()]


output = {
    "schema": "elkies-k3.q80-third-q12-jacobian-minimal-modp2.v1",
    "status": "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "minimal_short_weierstrass": {
        "equation": "Y^2=X^3+A(V)X+B(V)",
        "A_coefficients_low_to_high_1_r": polynomial_record(A_minimal),
        "B_coefficients_low_to_high_1_r": polynomial_record(B_minimal),
        "degrees_A_B": [8, 12],
        "discriminant_coefficients_low_to_high_1_r": polynomial_record(delta_minimal),
        "discriminant_factorization": [
            [str(factor.monic()), int(exponent)]
            for factor, exponent in delta_minimal.factor()
        ],
    },
    "long_to_minimal_map": {
        "gauge_pole": str(d),
        "q": str(q),
        "X_min": "q^2*(x_long+b2/12)",
        "Y_min": "q^3*(y_long+(a1*x_long+a3)/2)",
        "inverse_x_long": "X_min/q^2-b2/12",
        "inverse_y_long": "Y_min/q^3-(a1*x_long+a3)/2",
        "literal_invariant_replay": True,
    },
    "fibres": {
        "finite": fibres,
        "infinity": "smooth",
        "configuration": "I6+I4+3I2+8I1",
        "root_lattice": "A5+A3+3A1",
        "root_rank": 11,
        "euler_sum": 24,
    },
    "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "explicit integral short Weierstrass K3 model over GF(19^2)(V)",
            "exact long/short scaling and inverse",
            "minimal multiplicative fibre configuration I6+I4+3I2+8I1",
            "root lattice A5+A3+3A1",
        ],
        "not_proved": [
            "transported old-component and zero-section marking",
            "generic forward/inverse maps from the resolved plane model",
            "characteristic-zero MW rank or a second-prime alignment",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/minimize_q80_third_q12_jacobian_mod19_quadratic.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12MINIMAL|degrees=8,12|Delta_degree=24|"
    "fibres=I6+I4+3I2+8I1|roots=A5+A3+3A1|"
    "status=PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC",
    flush=True,
)

#!/usr/bin/env sage -python
"""Hensel-lift isolated polynomial-section residues of the H92 q=6 child.

This consumes an explicitly bounded modular polynomial-section search.  For
each signed solution, it tests the 13-by-12 coefficient Jacobian and, when it
has full column rank, performs deterministic coefficientwise Hensel lifting.
The output is p-adic data only.  It is not a rational reconstruction, an
identification of the MW basis, or a q=8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-43-iv-singular.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-hensel-43.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pad(value, length):
    return [value.parent()(value[index]) if index <= value.degree() else value.parent()(0)
            for index in range(length)] if value else [value.parent()(0)] * length


def coefficients_mod(ring, values):
    answer = []
    for value in values:
        value = QQ(value)
        denominator = ring(ZZ(value.denominator()))
        if not denominator:
            raise ValueError("prime divides a child-model coefficient denominator")
        answer.append(ring(ZZ(value.numerator())) / denominator)
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--precision", type=int, default=64)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.child = args.child.resolve()
args.input = args.input.resolve()
args.output = args.output.resolve()

if args.precision < 1:
    raise ValueError("precision must be positive")
child = json.loads(args.child.read_text())
search = json.loads(args.input.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert search["status"] == "EXPERIMENTAL_EXHAUSTIVE_MODULAR_ANSATZ"
prime = ZZ(search["prime"])
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("input prime must be an odd prime different from 3")
assert search["inputs"]["child_jacobian_sha256"] == digest(args.child)

field = GF(prime)
ring = PolynomialRing(field, "T")
T = ring.gen()
coefficient_a = ring(coefficients_mod(
    field, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
))
coefficient_b = ring(coefficients_mod(
    field, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
))

def section_polynomial(entry, name):
    return ring([field(value) for value in entry[name]])


def coefficient_jacobian(x_value, y_value):
    derivative_x = -3 * x_value**2 - coefficient_a
    derivative_y = 2 * y_value
    columns = []
    for power in range(5):
        columns.append(pad(T**power * derivative_x, 13))
    for power in range(7):
        columns.append(pad(T**power * derivative_y, 13))
    return matrix(field, columns).transpose()


def lift(entry):
    x_mod_p = section_polynomial(entry, "x_coefficients_low_to_high")
    y_mod_p = section_polynomial(entry, "y_coefficients_low_to_high")
    assert y_mod_p**2 == x_mod_p**3 + coefficient_a * x_mod_p + coefficient_b
    jacobian = coefficient_jacobian(x_mod_p, y_mod_p)
    jacobian_rank = jacobian.rank()
    record = {
        "x_mod_p": [int(value) for value in x_mod_p.list()],
        "y_mod_p": [int(value) for value in y_mod_p.list()],
        "jacobian_rank": int(jacobian_rank),
    }
    if jacobian_rank != 12:
        record["lift_status"] = "NONISOLATED_OR_SINGULAR_MOD_P"
        return record

    x_coefficients = [ZZ(value) for value in pad(x_mod_p, 5)]
    y_coefficients = [ZZ(value) for value in pad(y_mod_p, 7)]
    for level in range(1, args.precision):
        modulus = prime ** (level + 1)
        coefficient_ring = ZZ.quotient(modulus)
        polynomial_ring = PolynomialRing(coefficient_ring, "T")
        t = polynomial_ring.gen()
        a_lift = polynomial_ring(coefficients_mod(
            coefficient_ring,
            child["minimal_short_weierstrass"]["A_coefficients_low_to_high"],
        ))
        b_lift = polynomial_ring(coefficients_mod(
            coefficient_ring,
            child["minimal_short_weierstrass"]["B_coefficients_low_to_high"],
        ))
        x_lift = polynomial_ring(x_coefficients)
        y_lift = polynomial_ring(y_coefficients)
        residual = y_lift**2 - x_lift**3 - a_lift * x_lift - b_lift
        rhs = []
        for index in range(13):
            coefficient = ZZ(residual[index]) if index <= residual.degree() else ZZ(0)
            assert coefficient % prime**level == 0
            rhs.append(field(-(coefficient // prime**level)))
        correction = jacobian.solve_right(vector(field, rhs))
        x_coefficients = [
            coefficient + prime**level * ZZ(correction[index])
            for index, coefficient in enumerate(x_coefficients)
        ]
        y_coefficients = [
            coefficient + prime**level * ZZ(correction[5 + index])
            for index, coefficient in enumerate(y_coefficients)
        ]
    modulus = prime**args.precision
    record.update({
        "lift_status": "PASS_UNIQUE_HENSEL_LIFT",
        "precision": int(args.precision),
        "modulus": str(modulus),
        "x_coefficients_mod_p_to_precision": [str(value % modulus) for value in x_coefficients],
        "y_coefficients_mod_p_to_precision": [str(value % modulus) for value in y_coefficients],
    })
    return record


records = [lift(entry) for entry in search["sections"]]
payload = {
    "schema": "elkies-k3.h92-q6-child-polynomial-sections-hensel.v1",
    "status": "EXPERIMENTAL_P_ADIC_POLYNOMIAL_SECTION_LIFTS",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "modular_search": {"path": str(args.input.relative_to(ROOT)), "sha256": digest(args.input)},
    },
    "prime": int(prime),
    "precision": int(args.precision),
    "records": records,
    "boundary": (
        "These are isolated p-adic lifts of the supplied finite-field polynomial "
        "solutions only. No rational reconstruction, characteristic-zero section, "
        "MW identification, q8 pencil, bisection, collision, or rank claim is made."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDPOLYHENSEL|prime={}|precision={}|records={}|isolated={}|"
    "status=EXPERIMENTAL_P_ADIC_POLYNOMIAL_SECTION_LIFTS".format(
        prime, args.precision, len(records),
        sum(record["lift_status"] == "PASS_UNIQUE_HENSEL_LIFT" for record in records),
    ),
    flush=True,
)

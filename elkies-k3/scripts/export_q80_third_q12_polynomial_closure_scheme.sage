#!/usr/bin/env sage
"""Export the algebraic-closure polynomial-section gate for the Q80 q12 seed.

At ``u=-2`` modulo 19 the eight rational polynomial section pairs generate a
rank-three height lattice.  The exact MW5 embedding certificate has one
orientation in which the desired height-eight horizontal ``H`` satisfies

    (H - P4).O = 0,

where ``P4`` is the fourth ordered polynomial point.  Consequently the hard
denominator-two target can be recovered from the algebraic-closure
``P.O=0`` section scheme, followed by one exact group-law addition.

This script exports the exact-degree-four polynomial chart.  It uses

    X = l^2 W^4 + x3 W^3 + ... + x0,
    Y = l^3 W^6 + y5 W^5 + ... + y0,

solves ``y5,...,y0`` recursively from the top six coefficients, and writes
the six residual equations plus ``sat*l-1`` for msolve.  The output is a
finite-field search gate, not yet a section or characteristic-zero lift.
"""

import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT / "artifacts/generated-results/"
    "q80-third-q12-um2-p19-height-shell-complete.json"
)
LATTICE = ROOT / "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"
SYSTEM = (
    ROOT / "artifacts/local/elkies-k3/"
    "q80-third-q12-um2-p19-polynomial-closure.ms"
)
OUTPUT = (
    ROOT / "artifacts/generated-results/"
    "q80-third-q12-um2-p19-polynomial-closure-scheme.json"
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(INPUT.read_text())
lattice = json.loads(LATTICE.read_text())
record = payload["parameters"][0]
if record["u"] != "-2":
    raise ValueError("expected the pinned u=-2 specialization")
modular = record["modular"][0]
prime = int(modular["prime"])
if prime != 19 or modular["status"] != "PASS_GOOD_REDUCTION_MODULAR_SHELL":
    raise ValueError("expected the pinned good p=19 reduction")
embedding = lattice["polynomial_rank3_subgroup"]["embeddings"][0]
if embedding["ordered_P_dot_O_of_H_minus_eight_polynomial_points"][3] != 0:
    raise ValueError("the selected embedding does not have (H-P4).O=0")

finite = GF(prime)
base = PolynomialRing(finite, "W")
W = base.gen()


def reduce_rational(value):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ZeroDivisionError("coefficient denominator vanishes modulo p")
    return finite(value.numerator()) / denominator


equation = record["exact_equations"]["second_q4"]
A = base([reduce_rational(value) for value in equation["A_coefficients_low_to_high"]])
B = base([reduce_rational(value) for value in equation["B_coefficients_low_to_high"]])
delta = 4*A**3 + 27*B**2
signature = sorted((factor.degree(), int(exponent)) for factor, exponent in delta.factor())
if signature != [(1, 7), (2, 1), (6, 1)]:
    raise ArithmeticError(f"unexpected p=19 reduction signature: {signature}")

names = ("l", "x0", "x1", "x2", "x3", "sat")
scheme = PolynomialRing(finite, names=names, order="degrevlex")
l, x0, x1, x2, x3, sat = scheme.gens()
fraction = scheme.fraction_field()
polynomial = PolynomialRing(fraction, "W_section")
W_section = polynomial.gen()


def lift(poly):
    return polynomial([fraction(value) for value in poly.list()])


X = l**2*W_section**4 + x3*W_section**3 + x2*W_section**2 + x1*W_section + x0
square = X**3 + lift(A)*X + lift(B)
y = [fraction.zero() for unused in range(7)]
y[6] = fraction(l**3)
for degree in range(11, 5, -1):
    index = degree-6
    partial = sum(y[j]*W_section**j for j in range(7))
    known = (partial**2)[degree]
    y[index] = (square[degree]-known)/(2*y[6])
Y = sum(y[j]*W_section**j for j in range(7))
identity = Y**2-square
if any(identity[index] for index in range(6, 13)):
    raise ArithmeticError("top-down polynomial-section recursion did not close")
residual = [scheme(identity[index].numerator()) for index in range(6)]
residual.append(sat*l-1)
if any(not value for value in residual):
    raise ArithmeticError("unexpected zero equation in polynomial-section chart")

SYSTEM.parent.mkdir(parents=True, exist_ok=True)
with SYSTEM.open("w") as stream:
    stream.write(",".join(names)+"\n")
    stream.write(str(prime)+"\n")
    for index, value in enumerate(residual):
        stream.write(str(value).replace("**", "^"))
        stream.write(",\n" if index+1 < len(residual) else "\n")

output = {
    "schema": "elkies-k3.q80-third-q12-polynomial-closure-scheme.v1",
    "status": "PASS_EXACT_POLYNOMIAL_CLOSURE_SCHEME_EXPORTED",
    "prime": prime,
    "specialization": "u=-2",
    "chart": {
        "P_dot_O": 0,
        "x_degree": 4,
        "y_degree": 6,
        "variables": list(names),
        "equations": len(residual),
        "term_counts": [len(value.dict()) for value in residual],
        "leading_parameterization": "x4=l^2,y6=l^3,l!=0",
    },
    "target_reduction": {
        "embedding_index_one_based": 1,
        "known_polynomial_point_index_one_based": 4,
        "identity": "Q=H-P4 has Q.O=0; recover H=P4+Q",
        "ordered_target_intersection_fingerprint": embedding[
            "ordered_P_dot_O_of_H_minus_eight_polynomial_points"
        ],
    },
    "system": {
        "path": str(SYSTEM.relative_to(ROOT)),
        "sha256": sha256(SYSTEM),
        "size_bytes": SYSTEM.stat().st_size,
    },
    "inputs": {
        "modular_shell": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
        "height_lattice": {"path": str(LATTICE.relative_to(ROOT)), "sha256": sha256(LATTICE)},
    },
    "claim_boundary": {
        "proved": [
            "exact reduced p=19 polynomial-section scheme",
            "target-to-polynomial reduction from the complete MW5 embedding certificate",
        ],
        "not_proved": [
            "a non-rational algebraic polynomial section",
            "identification of H=P4+Q",
            "a third-q12 equation",
            "characteristic-zero lifting",
        ],
    },
    "reproduce": "sage elkies-k3/scripts/export_q80_third_q12_polynomial_closure_scheme.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True, default=int)+"\n")
print(
    "Q80POLYCLOSURE|prime={}|variables={}|equations={}|terms={}|bytes={}|"
    "status=PASS_EXACT_POLYNOMIAL_CLOSURE_SCHEME_EXPORTED".format(
        prime, len(names), len(residual), tuple(len(value.dict()) for value in residual),
        SYSTEM.stat().st_size,
    ),
    flush=True,
)

#!/usr/bin/env sage -python
"""Export the Q80 polynomial-section closure scheme at any certified (u,p)."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, prod, valuation


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--surface", type=Path, required=True)
parser.add_argument("--system", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.surface = args.surface.resolve()
args.system = args.system.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


surface = json.loads(args.surface.read_text())
if surface.get("schema") != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected surface schema")
if len(surface["parameters"]) != 1 or len(surface["parameters"][0].get("modular", [])) != 1:
    raise ValueError("producer requires exactly one (u,p) specialization")
parameter = surface["parameters"][0]
modular = parameter["modular"][0]
prime = int(modular["prime"])
if modular.get("status") not in (
    "PASS_GOOD_REDUCTION_AUDIT",
    "PASS_GOOD_REDUCTION_MODULAR_SHELL",
):
    raise ValueError("input does not certify good reduction")

finite = GF(prime)
base = PolynomialRing(finite, "W")
W = base.gen()


def reduce_rational(value):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ZeroDivisionError("coefficient denominator vanishes modulo p")
    return finite(value.numerator()) / denominator


equation = parameter["exact_equations"]["second_q4"]
A = base([reduce_rational(value) for value in equation["A_coefficients_low_to_high"]])
B = base([reduce_rational(value) for value in equation["B_coefficients_low_to_high"]])
delta = 4 * A**3 + 27 * B**2
factorization = tuple(delta.factor())
star_factors = [
    factor.monic()
    for factor, exponent in factorization
    if factor.degree() == 1 and int(exponent) == 7
]
residual_factors = [factor.monic() for factor, exponent in factorization if int(exponent) == 1]
if (
    len(star_factors) != 1
    or any(int(exponent) not in (1, 7) for factor, exponent in factorization)
    or sum(factor.degree() for factor in residual_factors) != 8
    or not prod(residual_factors, base.one()).is_squarefree()
):
    raise ArithmeticError("input reduction no longer has I1* plus eight simple residual roots")
star = star_factors[0]
if (valuation(A, star), valuation(B, star), valuation(delta, star)) != (2, 3, 7):
    raise ArithmeticError("finite exponent-seven place is not minimal I1*")

names = ("l", "x0", "x1", "x2", "x3", "sat")
scheme = PolynomialRing(finite, names=names, order="degrevlex")
l, x0, x1, x2, x3, sat = scheme.gens()
fraction = scheme.fraction_field()
polynomial = PolynomialRing(fraction, "W_section")
W_section = polynomial.gen()


def lift(poly):
    return polynomial([fraction(value) for value in poly.list()])


X = l**2 * W_section**4 + x3 * W_section**3 + x2 * W_section**2 + x1 * W_section + x0
square = X**3 + lift(A) * X + lift(B)
y = [fraction.zero() for unused in range(7)]
y[6] = fraction(l**3)
for degree in range(11, 5, -1):
    index = degree - 6
    partial = sum(y[j] * W_section**j for j in range(7))
    y[index] = (square[degree] - (partial**2)[degree]) / (2 * y[6])
Y = sum(y[j] * W_section**j for j in range(7))
identity = Y**2 - square
if any(identity[index] for index in range(6, 13)):
    raise ArithmeticError("top-down section recursion did not close")
residual = [scheme(identity[index].numerator()) for index in range(6)]
residual.append(sat * l - 1)
if any(not value for value in residual):
    raise ArithmeticError("unexpected zero closure equation")

args.system.parent.mkdir(parents=True, exist_ok=True)
with args.system.open("w") as stream:
    stream.write(",".join(names) + "\n")
    stream.write(str(prime) + "\n")
    for index, value in enumerate(residual):
        stream.write(str(value).replace("**", "^"))
        stream.write(",\n" if index + 1 < len(residual) else "\n")

output = {
    "schema": "elkies-k3.q80-third-q12-polynomial-closure-producer-modp.v1",
    "status": "PASS_EXACT_POLYNOMIAL_CLOSURE_PRODUCER_EXPORTED",
    "specialization": {"u": parameter["u"], "prime": prime},
    "good_reduction": {
        "factor_degrees_exponents": [
            [int(factor.degree()), int(exponent)] for factor, exponent in factorization
        ],
        "finite_I1star_factor": str(star),
        "residual_degree": 8,
        "residual_squarefree": True,
    },
    "chart": {
        "P_dot_O": 0,
        "x_degree": 4,
        "y_degree": 6,
        "variables": list(names),
        "equations": len(residual),
        "term_counts": [len(value.dict()) for value in residual],
        "leading_parameterization": "x4=l^2,y6=l^3,l!=0",
    },
    "system": {
        "path": str(args.system.relative_to(ROOT)),
        "sha256": sha256(args.system),
        "size_bytes": args.system.stat().st_size,
    },
    "input": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "arbitrary certified good-prime polynomial-section closure scheme",
            "no rational polynomial-shell enumeration is required",
        ],
        "not_proved": [
            "existence or selection of a target quadratic horizontal",
            "a child equation or characteristic-zero lift",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/produce_q80_third_q12_polynomial_closure_modp.sage "
        f"--surface {args.surface} --system {args.system} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80POLYCLOSUREPRODUCER|u={parameter['u']}|prime={prime}|variables=6|equations=7|"
    "status=PASS_EXACT_POLYNOMIAL_CLOSURE_PRODUCER_EXPORTED",
    flush=True,
)

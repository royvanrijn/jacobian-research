#!/usr/bin/env sage-python
"""Certify the obstruction in the simplest two-point A4 Tate slice.

status: ACTIVE_PROOF
claim: exact obstruction only for the declared r=s=1 two-point Tate slice

Start from

    y^2+a1*x*y+a3*y=x^3+a2*x^2,
    a1=1+A*u, h=H*u, k=B*u+C*u^2,
    a2=k, a3=h*(1-h-a1)+k.

The points P=(0,0) and Q=(h,h^2) are automatic.  On B != 0, imposing
ord_u(Delta)>=5 gives H=-B or H=B-A.  In both branches the residual
discriminant contains

    (B+(A*B-B^2+C)*u)^2.

Avoiding the resulting extra reducible fibre (generically I2) forces
C=B^2-A*B.  The two points then
satisfy P+2Q=0 on H=-B and 3P+2Q=0 on H=B-A.  Thus this normalized slice
cannot give two independent marked directions in the pure A4 stratum.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT/"artifacts/generated-results"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

parameter_ring = PolynomialRing(QQ, names=("A", "B", "C", "H"))
coefficient_field = parameter_ring.fraction_field()
A, B, C, H = coefficient_field.gens()
u_ring = PolynomialRing(coefficient_field, "u")
u = u_ring.gen()

a1 = 1+A*u
h = H*u
kappa = B*u+C*u**2
a2 = kappa
a3 = h*(1-h-a1)+kappa
b2 = a1**2+4*a2
b4 = a1*a3
b6 = a3**2
b8 = a2*a3**2
discriminant = -b2**2*b8-8*b4**3-27*b6**2+9*b2*b4*b6

if any(discriminant[degree] != 0 for degree in range(4)):
    raise ArithmeticError("the normalized Tate slice lost its automatic u^4")
expected_u4 = -B**2*(H+B)*(A+H-B)
if discriminant[4] != expected_u4:
    raise ArithmeticError("the two I5 branches changed")


def substitute_coefficients(polynomial, substitutions):
    return u_ring(
        sum(
            coefficient_field(polynomial[degree]).subs(substitutions)*u**degree
            for degree in range(polynomial.degree()+1)
        )
    )


branch_data = []
for branch_name, H_value, relation in (
    ("H=-B", -B, (1, 2)),
    ("H=B-A", B-A, (3, 2)),
):
    substitutions = {H: H_value}
    branch_discriminant = substitute_coefficients(discriminant, substitutions)
    residual = branch_discriminant//u**5
    linear = B+(A*B-B**2+C)*u
    quotient, remainder = residual.quo_rem(linear**2)
    if remainder != 0:
        raise ArithmeticError(f"forced square factor changed on {branch_name}")

    pure_substitutions = {H: H_value, C: B**2-A*B}
    pure_a1 = substitute_coefficients(a1, pure_substitutions)
    pure_a2 = substitute_coefficients(a2, pure_substitutions)
    pure_a3 = substitute_coefficients(a3, pure_substitutions)
    pure_h = substitute_coefficients(h, pure_substitutions)
    function_field = u_ring.fraction_field()
    curve = EllipticCurve(
        function_field,
        [pure_a1, pure_a2, pure_a3, 0, 0],
    )
    P = curve(0, 0)
    Q = curve(pure_h, pure_h**2)
    if relation[0]*P+relation[1]*Q != curve(0):
        raise ArithmeticError(f"marked-point relation changed on {branch_name}")
    branch_data.append(
        {
            "branch": branch_name,
            "forced_residual_square_factor": "(B+(A*B-B^2+C)*u)^2",
            "pure_A4_condition": "C=B^2-A*B",
            "marked_point_relation": f"{relation[0]}*P+{relation[1]}*Q=0",
            "residual_quotient_degree": int(quotient.degree()),
        }
    )

summary = {
    "status": "PASS_EXACT_SCOPED_OBSTRUCTION",
    "claim_boundary": (
        "Exact only for the normalized r=s=1 two-point Tate slice with "
        "a1=1+A*u, h=H*u, and kappa=B*u+C*u^2 on B nonzero."
    ),
    "automatic_marked_points": ["P=(0,0)", "Q=(h,h^2)"],
    "discriminant_u4_coefficient": "-B^2*(H+B)*(A+H-B)",
    "branches": branch_data,
    "conclusion": (
        "A pure A4 surface in this slice makes the two marked points dependent; "
        "keeping them potentially independent leaves a forced repeated "
        "discriminant root, generically an extra I2 fibre."
    ),
}

GENERATED.mkdir(parents=True, exist_ok=True)
output_path = GENERATED/"elkies-k3-a4-two-point-tate-slice-obstruction-v1.json"
serialized = json.dumps(summary, indent=2, sort_keys=True)+"\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale or missing generated summary: {output_path}")
else:
    output_path.write_text(serialized)

print("A4_TWO_POINT_TATE_SLICE|branches=2|status=PASS_EXACT_SCOPED_OBSTRUCTION")
print(output_path)

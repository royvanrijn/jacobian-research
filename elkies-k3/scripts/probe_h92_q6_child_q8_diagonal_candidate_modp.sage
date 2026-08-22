#!/usr/bin/env sage -python
"""Fast modular branch-degree screen for the diagonal finite-q-module pencil.

This is the reduction of the characteristic-zero diagonal candidate
``V=f_IV*(m-p)/(h*f_II^2*f_IV^3)``.  It rejects the candidate if the
quadratic x-discriminant has a branch degree other than four over a good
finite field.  Passing is only a necessary condition for the q8 pencil.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime divides an input denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--v", type=int, default=1)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
A = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = field(polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
)
sy = field(polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
)
assert sy**2 == sx**3 + field(A)*sx + field(B)
h = monic_power_root(ring(sx.denominator()), 2)
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).list())
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).list())
M = ii**2 * iv**3
p = -sy/sx
m = p + field(h*M/iv) * finite(args.v)
x_ring = PolynomialRing(field, "x")
x = x_ring.gen()
y = x_ring(m) * (x-x_ring(sx)) - x_ring(sy)
relation = y**2-x**3-x_ring(A)*x-x_ring(B)
quadratic, remainder = relation.quo_rem(x-x_ring(sx))
assert remainder == 0 and quadratic.degree() == 2
a, b, c = quadratic[2], quadratic[1], quadratic[0]
discriminant = field(b**2-4*a*c)
numerator = ring(discriminant.numerator())
denominator = ring(discriminant.denominator())
odd_degree = sum(
    factor.degree()
    for polynomial_value in (numerator, denominator)
    for factor, multiplicity in polynomial_value.squarefree_decomposition()
    if multiplicity % 2
)
infinity = (denominator.degree()-numerator.degree()) % 2
branch_degree = odd_degree + infinity
print(
    "H92Q6CHILDQ8DIAGONALMODP|prime={}|V={}|branch_degree={}|"
    "finite_odd_degree={}|infinity_branch={}|status={}".format(
        args.prime, args.v, branch_degree, odd_degree, infinity,
        "PASS_NECESSARY_GENUS_ONE" if branch_degree == 4 else "REJECTED_BRANCH_DEGREE",
    ),
    flush=True,
)

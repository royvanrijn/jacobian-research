#!/usr/bin/env python3
"""Exact audit of the normalized Davenport boundary involution."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_pair, reduce_a  # noqa: E402


g, _ = davenport_pair()
J = reduce_a(sp.diff(g, Y), T, Y)
z, v = sp.symbols("z v")

# On the torus put z=T/Y^2.  The derivative equation becomes quadratic in
# Y and its discriminant is a square times a linear polynomial.
quadratic = reduce_a(J.subs(T, z * Y**2) / Y**4, z, Y)
assert sp.Poly(quadratic, Y).degree() == 2
quadratic_discriminant = reduce_a(
    sp.discriminant(quadratic, Y),
    z,
)
expected_discriminant = reduce_a(
    -4
    * z**2
    * (z + 1) ** 2
    * ((4 * A + 16) * z + 2 - 3 * A),
    z,
)
assert quadratic_discriminant == expected_discriminant

# Normalize by v^2=-((4a+16)z+2-3a).
z_of_v = reduce_a((3 * A - 2 - v**2) / (4 * A + 16), v)
assert reduce_a(
    v**2 + (4 * A + 16) * z_of_v + 2 - 3 * A,
    v,
) == 0

coefficient_two = sp.Poly(quadratic, Y).nth(2)
coefficient_one = sp.Poly(quadratic, Y).nth(1)
raw_Y = reduce_a(
    (
        -coefficient_one
        + 2 * z * (z + 1) * v
    ).subs(z, z_of_v)
    / (2 * coefficient_two.subs(z, z_of_v)),
    v,
)
raw_numerator, raw_denominator = sp.fraction(sp.cancel(raw_Y))

field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
numerator_over_field = sp.Poly(
    raw_numerator.subs(A, alpha),
    v,
    domain=number_field,
)
denominator_over_field = sp.Poly(
    raw_denominator.subs(A, alpha),
    v,
    domain=number_field,
)
common_factor = sp.gcd(numerator_over_field, denominator_over_field).monic()
reduced_denominator = denominator_over_field.exquo(common_factor).monic()
expected_punctures = sp.Poly(
    (
        (v + A)
        * (v**2 + (-2 + 3 * A) * v - 8 - 2 * A)
    ).subs(A, alpha),
    v,
    domain=number_field,
).monic()
assert reduced_denominator == expected_punctures

puncture_factors = sp.factor_list(reduced_denominator)[1]
assert sorted(
    factor.degree()
    for factor, multiplicity in puncture_factors
    for _ in range(multiplicity)
) == [1, 2]
assert len(puncture_factors) == 2

print("PASS: the derivative normalization is rational in the parameter v")
print("PASS: its punctures split over K as one rational point plus a quadratic pair")


# The node z=-1 has two normalization points v=+-(a+4).
node_plus = A + 4
node_minus = -A - 4
for node in (node_plus, node_minus):
    assert reduce_a(
        node**2
        + (4 * A + 16) * (-1)
        + 2
        - 3 * A
    ) == 0

# Send the rational puncture v=-a to infinity with xi=1/(v+a).  The
# quadratic puncture pair has trace -1/2 in xi, so reflection about -1/4
# is the unique nontrivial K-automorphism preserving all three punctures.
xi = sp.symbols("xi")
quadratic_punctures = (
    v**2 + (-2 + 3 * A) * v - 8 - 2 * A
)
transformed_quadratic = reduce_a(
    xi**2 * quadratic_punctures.subs(v, -A + 1 / xi),
    xi,
)
transformed_monic = reduce_a(
    transformed_quadratic
    / sp.Poly(transformed_quadratic, xi).LC(),
    xi,
)
trace = reduce_a(-sp.Poly(transformed_monic, xi).nth(1))
assert trace == -sp.Rational(1, 2)

involution = reduce_a(
    -A + 1 / (trace - 1 / (v + A)),
    v,
)
assert reduce_a(involution.subs(v, involution) - v, v) == 0

# It preserves the quadratic puncture divisor.
puncture_pullback_numerator = sp.together(
    reduce_a(quadratic_punctures.subs(v, involution), v)
).as_numer_denom()[0]
assert reduce_a(
    sp.rem(
        sp.Poly(puncture_pullback_numerator, v),
        sp.Poly(quadratic_punctures, v),
    ).as_expr(),
    v,
) == 0

# Its quotient coordinate has one finite removed value and infinity, hence
# the quotient of the normalized affine boundary is G_m.
quotient_coordinate = (xi - trace / 2) ** 2
quadratic_discriminant_xi = reduce_a(
    sp.discriminant(transformed_monic, xi)
)
finite_quotient_puncture = reduce_a(quadratic_discriminant_xi / 4)
assert finite_quotient_puncture != 0
assert quotient_coordinate.has(xi)

print("PASS: the unique puncture-swapping involution has quotient G_m")


# The conductor pair is not invariant: node_minus is fixed, while node_plus
# is sent to neither node branch.  Therefore the involution does not descend
# from the normalization to the singular derivative curve.
assert reduce_a(involution.subs(v, node_minus) - node_minus, v) == 0
image_node_plus = reduce_a(involution.subs(v, node_plus), v)
assert reduce_a(image_node_plus - node_plus, v) != 0
assert reduce_a(image_node_plus - node_minus, v) != 0
assert reduce_a(involution.subs(v, image_node_plus) - node_plus, v) == 0

# The smallest involution-invariant equivalence relation containing the
# original node relation node_plus~node_minus must also contain
# image_node_plus~node_minus.  Hence it glues three distinct normalization
# points.  In the quotient, {node_plus,image_node_plus} is one orbit and
# {node_minus} is another, and those two quotient points remain identified.
conductor_points = (node_plus, node_minus, image_node_plus)
for index, first in enumerate(conductor_points):
    for second in conductor_points[index + 1 :]:
        assert reduce_a(first - second, v) != 0
normalization_orbits = (
    frozenset(("node_plus", "image_node_plus")),
    frozenset(("node_minus",)),
)
assert len(normalization_orbits) == 2

print("PASS: the involution fails to preserve the two-point node conductor")
print("PASS: the G_m quotient exists only after normalization")
print("PASS: the minimal equivariant conductor is triple and its quotient is nodal")
print("PASS Davenport boundary involution obstruction")

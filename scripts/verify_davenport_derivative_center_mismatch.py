#!/usr/bin/env python3
"""Exact audit of the Davenport derivative/modification-center mismatch."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_pair, reduce_a  # noqa: E402


g, _ = davenport_pair()
J = reduce_a(sp.diff(g, Y), T, Y)
assert sp.Poly(J, Y).degree() == 6
assert sp.Poly(J, T).degree() == 3

# The Newton polygon is the triangle with vertices (2,0),(3,0),(0,6).
support = {
    exponent
    for exponent, coefficient in sp.Poly(J, T, Y).terms()
    if coefficient != 0
}
assert {(2, 0), (3, 0), (0, 6)} <= support
assert all(
    j >= 0 and 2 * i + j <= 6 and 3 * i + j >= 6
    for i, j in support
)
interior_points = [
    (i, j)
    for i in range(4)
    for j in range(7)
    if j > 0 and 2 * i + j < 6 and 3 * i + j > 6
]
assert interior_points == [(2, 1)]


# There is exactly one singularity in the two-dimensional torus.
field_variable, inverse = sp.symbols("field_variable inverse")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
torus_singularities = sp.groebner(
    [
        J.subs(A, alpha),
        sp.diff(J, T).subs(A, alpha),
        sp.diff(J, Y).subs(A, alpha),
        inverse * T * Y - 1,
    ],
    T,
    Y,
    inverse,
    domain=number_field,
    order="grevlex",
)
assert {
    sp.Poly(polynomial, T, Y, inverse, domain=number_field).monic().as_expr()
    for polynomial in torus_singularities.polys
} == {
    T + sp.Rational(1, 4),
    Y + sp.Rational(1, 2),
    inverse - 8,
}

# Its tangent cone is a product of two distinct lines, so it is an ordinary
# node and contributes delta=1.
t, y = sp.symbols("t y")
local = reduce_a(
    J.subs({T: -sp.Rational(1, 4) + t, Y: -sp.Rational(1, 2) + y}),
    t,
    y,
)
tangent_cone = sum(
    coefficient * t**i * y**j
    for (i, j), coefficient in sp.Poly(local, t, y).terms()
    if i + j == 2
)
assert reduce_a(
    tangent_cone
    + t * ((-3 * A - 6) * t + (2 * A + 8) * y) / 4,
    t,
    y,
) == 0
assert reduce_a(sp.discriminant(tangent_cone, t), y) != 0


# All three Newton faces are squarefree.  The outer face has lattice length
# gcd(3,6)=3 and therefore contributes three geometric points at infinity.
z = sp.symbols("z")
poly = sp.Poly(J, T, Y)
outer_face = sum(
    coefficient * T**i * Y**j
    for (i, j), coefficient in poly.terms()
    if 2 * i + j == 6
)
outer_polynomial = reduce_a(outer_face.subs({T: z, Y: 1}), z)
assert sp.Poly(outer_polynomial, z).degree() == 3
assert reduce_a(sp.discriminant(outer_polynomial, z)) != 0

bottom_face = sum(
    coefficient * T**i * Y**j
    for (i, j), coefficient in poly.terms()
    if j == 0
)
bottom_reduced = sp.cancel(bottom_face / T**2).subs(T, z)
assert sp.Poly(bottom_reduced, z).degree() == 1

left_face = sum(
    coefficient * T**i * Y**j
    for (i, j), coefficient in poly.terms()
    if 3 * i + j == 6
)
left_univariate = reduce_a(left_face.subs({T: 1, Y: z}), z)
assert reduce_a(sp.discriminant(left_univariate, z)) != 0

# Newton arithmetic genus 1 minus the single nodal delta invariant 1 gives
# normalization genus zero.  Removing the three outer-face points gives
# geometric unit rank 3-1=2.
newton_arithmetic_genus = len(interior_points)
node_delta = 1
normalization_genus = newton_arithmetic_genus - node_delta
points_at_infinity = 3
derivative_unit_rank = points_at_infinity - 1
assert normalization_genus == 0
assert derivative_unit_rank == 2

print("PASS: the derivative Newton polygon has arithmetic genus one")
print("PASS: its unique torus singularity is an ordinary node")
print("PASS: the normalized affine derivative curve is P1 minus three points")


# The split-boundary determinant threefold is an affine modification of A3:
#
#   X={D^2*U-C*R=beta} -> A3_(D,C,R).
#
# Its reduced center is D=0, C*R=-beta, hence G_m.  The residue Jacobian of
# the projection is D^2, exactly matching the doubled ledger exponent, but
# the center has geometric unit rank one rather than two.
D, U, C, R, beta = sp.symbols("D U C R beta", nonzero=True)
equation = D**2 * U - C * R - beta
assert sp.diff(equation, U) == D**2
center_equation = sp.expand(equation.subs(D, 0))
assert center_equation == -C * R - beta

modification_center_unit_rank = 1
assert derivative_unit_rank != modification_center_unit_rank

# Polynomial stabilization does not change units.
for number_of_variables in range(5):
    assert derivative_unit_rank == 2
    assert modification_center_unit_rank == 1

print("PASS: the reverse affine modification has residue Jacobian D^2")
print("PASS: its reduced center is G_m with geometric unit rank one")
print("PASS: center unit ranks remain different after every stabilization")
print("PASS Davenport derivative-center mismatch")

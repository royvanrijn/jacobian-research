#!/usr/bin/env python3
"""Exact audit of the Davenport tangent-mark normalization curve."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_pair, reduce_a  # noqa: E402


g, _ = davenport_pair()
r, v, slope = sp.symbols("r v slope")

tangent_mark = reduce_a(
    (g.subs(Y, r) - g.subs(Y, 0) - r * sp.diff(g, Y).subs(Y, 0)) / r**2,
    T,
    r,
)
mark_poly = sp.Poly(tangent_mark, T, domain=sp.QQ.frac_field(A, r))
assert mark_poly.degree() == 2
assert sp.Poly(tangent_mark, r, domain=sp.QQ.frac_field(A, T)).degree() == 5

A2 = reduce_a(mark_poly.nth(2), r)
A1 = reduce_a(mark_poly.nth(1), r)
A0 = reduce_a(mark_poly.nth(0), r)
Q = reduce_a(sp.discriminant(tangent_mark, T) / r**4, r)
expected_Q = reduce_a(
    ((5 - A) * r**2 + (-2 * A - 6) * r + 7 * A - 7) / 7,
    r,
)
assert sp.factor(Q - expected_Q) == 0
assert reduce_a(Q.subs(r, 0) - (1 + A) ** 2) == 0
assert reduce_a(sp.discriminant(Q, r), r) != 0

# Completing the quadratic gives a birational map from the marking curve
# (away from r=0) to the smooth conic v^2=Q(r).
T_from_conic = r**2 * (-(1 + A) * (r + 1) + v) / (2 * A2)
quadratic_remainder = sp.cancel(
    A2 * T_from_conic**2 + A1 * T_from_conic + A0
)
assert sp.factor(
    reduce_a(
        quadratic_remainder
        - r**4 * (v**2 - Q) / (4 * A2),
        r,
        v,
    )
) == 0

# Lines v=(1+a)+slope*r through one rational point parametrize the conic.
r_parameter = -2 * (
    7 * A * slope + A + 7 * slope + 3
) / (A + 7 * slope**2 - 5)
v_parameter = 1 + A + slope * r_parameter
assert sp.factor(
    reduce_a(v_parameter**2 - Q.subs(r, r_parameter), slope)
) == 0
assert reduce_a(sp.diff(r_parameter, slope), slope) != 0

# The two distinct K-rational points above r=0 are removed by the nonzero
# tangent-mark condition.  On the normalization P^1 this leaves at least
# P^1\\{P_+,P_-}; its coordinate ring has a nonconstant unit and admits no
# nonconstant morphism from A^1.
point_plus = 1 + A
point_minus = -1 - A
assert reduce_a(point_plus**2 - Q.subs(r, 0)) == 0
assert reduce_a(point_minus**2 - Q.subs(r, 0)) == 0
assert reduce_a(point_plus - point_minus) != 0

# Allow both tangent points to move, but require them to be affine-linear in
# T.  The six coefficient equations have unit Groebner basis over K, so this
# entire polynomial escape hatch is empty.
q = sp.symbols("q")
q0, q1, r0, r1 = sp.symbols("q0 q1 r0 r1")
moving_tangent = reduce_a(
    (
        g.subs(Y, r)
        - g.subs(Y, q)
        - sp.diff(g, Y).subs(Y, q) * (r - q)
    )
    / (r - q) ** 2,
    T,
    q,
    r,
)
affine_substitution = {
    q: q0 + q1 * T,
    r: r0 + r1 * T,
}
coefficient_polynomial = sp.Poly(
    sp.expand(moving_tangent.subs(affine_substitution)),
    T,
)
coefficient_equations = [
    reduce_a(
        coefficient_polynomial.nth(degree),
        q0,
        q1,
        r0,
        r1,
    )
    for degree in range(coefficient_polynomial.degree() + 1)
]
field_variable = sp.symbols("field_variable")
coefficient_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = coefficient_field.ext
groebner = sp.groebner(
    [equation.subs(A, alpha) for equation in coefficient_equations],
    q0,
    q1,
    r0,
    r1,
    domain=coefficient_field,
    order="grevlex",
)
assert len(groebner.polys) == 1
assert groebner.polys[0].as_expr() == 1

print("PASS: the degree-five tangent-mark equation is quadratic in T")
print("PASS: its normalization is the explicit smooth rational conic v^2=Q(r)")
print("PASS: the conic has a rational slope parametrization")
print("PASS: the nonzero-mark open removes two distinct rational points")
print("PASS: the present relative weighted family has no nonconstant A1 reparametrization")
print("PASS: no tangent pair has both marks affine-linear in T")

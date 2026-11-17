#!/usr/bin/env python3
"""Exact audit of the two-center Davenport weighted-gluing obstruction."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    A,
    T,
    Y,
    Z,
    conjugate_a,
    davenport_pair,
    reduce_a,
)


g, h = davenport_pair()
s, W = sp.symbols("s W")

k_four = -A - 4
k_two = (9 * A - 48) / 23
mu = (5 + 2 * A) / 3


def reduced_primitive(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    first: sp.Expr,
    proportionality: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    difference = proportionality * first
    primitive = reduce_a(
        polynomial.subs(
            {
                variable: first + difference * W,
                T: -s**2,
            }
        )
        - polynomial.subs({variable: first, T: -s**2})
        - difference
        * sp.diff(polynomial, variable).subs(
            {variable: first, T: -s**2}
        )
        * W,
        s,
        W,
    )
    reduced = reduce_a(primitive / s**6, s, W)
    assert sp.Poly(reduced, W).degree() == 7
    return difference, reduced


data = []
for label, polynomial, variable, conjugate in (
    ("g", g, Y, False),
    ("h", h, Z, True),
):
    k4 = conjugate_a(k_four, s) if conjugate else k_four
    k2 = conjugate_a(k_two, s) if conjugate else k_two
    scale = conjugate_a(mu, s) if conjugate else mu

    d4, H4 = reduced_primitive(polynomial, variable, s, k4)
    d2, H2 = reduced_primitive(polynomial, variable, scale * s, k2)
    assert reduce_a(d4 - d2, s) == 0

    translation = reduce_a((s - scale * s) / d2, s)
    transition_difference = reduce_a(
        H4 - H2.subs(W, W + translation),
        s,
        W,
    )
    assert reduce_a(
        sp.diff(transition_difference, W, 2),
        s,
        W,
    ) == 0

    slope_shift = reduce_a(sp.diff(transition_difference, W), s)
    intercept_shift = reduce_a(
        transition_difference.subs(W, 0),
        s,
    )
    assert sp.Poly(slope_shift, s).degree() == 1
    assert sp.Poly(intercept_shift, s).degree() == 1
    assert reduce_a(
        sp.resultant(slope_shift, intercept_shift, s)
    ) != 0
    data.append(
        (
            label,
            translation,
            slope_shift,
            intercept_shift,
        )
    )

print("PASS: both Davenport overlaps have affine slope/intercept cocycles")
print("PASS: each slope and intercept center is linear and comaximal")


# The slope correction alone is the affine modification C*V=m(s).  Since
# m(s) is linear with nonzero leading coefficient, this hypersurface is an
# affine space: solve uniquely for s in terms of C*V.
C, V = sp.symbols("C V")
for _, _, slope_shift, _ in data:
    slope_poly = sp.Poly(slope_shift, s)
    leading = slope_poly.nth(1)
    constant = slope_poly.nth(0)
    recovered_s = reduce_a((C * V - constant) / leading, C, V)
    assert reduce_a(
        slope_shift.subs(s, recovered_s) - C * V,
        C,
        V,
    ) == 0

print("PASS: the slope-only modification C*V=m(s) is affine space")


# The full incidence transition also requires
#
#   C^2*U = n(s) + delta*B*C.
#
# Because (m,n)=(1), the two modification equations force C to be a unit.
# This is an exact ring identity, not only a set-theoretic observation.
B, U = sp.symbols("B U")
field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
function_polynomial_ring = number_field.poly_ring(s)

for _, translation, slope_shift, intercept_shift in data:
    slope_over_field = sp.Poly(
        slope_shift.subs(A, alpha),
        s,
        domain=number_field,
    )
    intercept_over_field = sp.Poly(
        intercept_shift.subs(A, alpha),
        s,
        domain=number_field,
    )
    bezout_s, bezout_t, gcd = sp.gcdex(
        slope_over_field,
        intercept_over_field,
    )
    assert gcd.degree() == 0
    bezout_s = bezout_s.quo(gcd)
    bezout_t = bezout_t.quo(gcd)
    assert (
        bezout_s * slope_over_field
        + bezout_t * intercept_over_field
    ).as_expr() == 1

    u = bezout_s.as_expr()
    v = bezout_t.as_expr()
    delta = translation.subs(A, alpha)
    m = slope_shift.subs(A, alpha)
    n = intercept_shift.subs(A, alpha)

    first_relation = C * V - m
    second_relation = C**2 * U - n - delta * B * C
    inverse_C = u * V + v * (C * U - delta * B)
    certificate = sp.expand(
        C * inverse_C
        - 1
        - u * first_relation
        - v * second_relation
    )
    assert sp.Poly(
        certificate,
        s,
        B,
        C,
        U,
        V,
        domain=number_field,
    ).is_zero

print("PASS: the full two-center modification has an explicit inverse for C")
print("PASS: canonical weighted gluing deletes rather than fills C=0")


# Give slope and intercept independent boundary variables C and D.  The
# forced modification then avoids inverting C, but linear elimination puts
# it in one uniform normal form:
#
#   A^2 x {D^2*U-C*R=beta},  beta != 0.
#
# Its threefold factor has class L^3-L, not L^3.
D, R = sp.symbols("D R")
normal_forms = []
for _, translation, slope_shift, intercept_shift in data:
    slope_poly = sp.Poly(slope_shift, s)
    intercept_poly = sp.Poly(intercept_shift, s)
    m1, m0 = slope_poly.nth(1), slope_poly.nth(0)
    n1, n0 = intercept_poly.nth(1), intercept_poly.nth(0)
    alpha_ratio = reduce_a(n1 / m1)
    beta = reduce_a(n0 - alpha_ratio * m0)
    assert beta != 0
    recovered_s = reduce_a((C * V - m0) / m1, C, V)
    assert reduce_a(
        intercept_shift.subs(s, recovered_s)
        - alpha_ratio * C * V
        - beta,
        C,
        V,
    ) == 0
    delta = translation
    new_coordinate = alpha_ratio * V + delta * B
    assert reduce_a(delta) != 0
    assert reduce_a(
        (
            D**2 * U
            - intercept_shift.subs(s, recovered_s)
            - delta * B * C
        )
        - (D**2 * U - C * new_coordinate - beta),
        B,
        C,
        D,
        U,
        V,
    ) == 0
    normal_forms.append(beta)

# Smoothness of D^2*U-C*R=beta for beta a unit.
beta_symbol = sp.symbols("beta_symbol", nonzero=True)
normal_equation = D**2 * U - C * R - beta_symbol
smooth_ideal = sp.groebner(
    [
        normal_equation,
        sp.diff(normal_equation, D),
        sp.diff(normal_equation, U),
        sp.diff(normal_equation, C),
        sp.diff(normal_equation, R),
    ],
    D,
    U,
    C,
    R,
    domain=sp.QQ.frac_field(beta_symbol),
)
assert len(smooth_ideal.polys) == 1
assert smooth_ideal.polys[0].as_expr() == 1

# D!=0 contributes (L-1)L^2.  At D=0, C*R=-beta and U is free,
# contributing (L-1)L.  Hence [X]=L^3-L and two free coordinates give
# L^5-L^3.
L = sp.symbols("L")
class_threefold = sp.expand((L - 1) * L**2 + (L - 1) * L)
assert class_threefold == L**3 - L
assert sp.expand(L**2 * class_threefold) == L**5 - L**3
assert class_threefold != L**3

print("PASS: independent boundary variables reduce to D^2*U-C*R=beta")
print("PASS: the split-boundary target is smooth with class L^5-L^3, not L^5")

# The natural orientation-forgetting involution D -> -D has invariant
# coordinate X=D^2.  Its quotient is X*U-C*R=beta, an SL2-type threefold,
# not affine three-space.  Thus the successful cubic factor-swap descent
# has no direct analogue here.
X = sp.symbols("X")
quotient_equation = X * U - C * R - beta_symbol
quotient_smooth_ideal = sp.groebner(
    [
        quotient_equation,
        sp.diff(quotient_equation, X),
        sp.diff(quotient_equation, U),
        sp.diff(quotient_equation, C),
        sp.diff(quotient_equation, R),
    ],
    X,
    U,
    C,
    R,
    domain=sp.QQ.frac_field(beta_symbol),
)
assert len(quotient_smooth_ideal.polys) == 1
assert quotient_smooth_ideal.polys[0].as_expr() == 1
quotient_class = sp.expand((L - 1) * L**2 + (L - 1) * L)
assert quotient_class == L**3 - L
assert quotient_class != L**3

print("PASS: the D -> -D quotient is SL2-type, not affine three-space")
print("PASS Davenport weighted glue obstruction")

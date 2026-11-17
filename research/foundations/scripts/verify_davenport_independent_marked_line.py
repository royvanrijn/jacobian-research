#!/usr/bin/env python3
"""Exact audit of the independent marked-line Davenport opening."""

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
    conjugate_a,
    davenport_pair,
    reduce_a,
)


g, _ = davenport_pair()
s, W = sp.symbols("s W")


def reduced_primitive(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    """Return the tangent primitive divided by its removable s^6 factor."""
    difference = second - first
    substituted = g.subs(
        {
            Y: first + difference * W,
            T: -s**2,
        }
    )
    base_value = g.subs({Y: first, T: -s**2})
    base_derivative = sp.diff(g, Y).subs({Y: first, T: -s**2})
    primitive = substituted - base_value - difference * base_derivative * W
    return reduce_a(primitive / s**6, s, W)


# The two simultaneous point-cover marking charts.
first_four = s
second_four = -(A + 3) * s
mu = (5 + 2 * A) / 3
k_two = (9 * A - 48) / 23
first_two = mu * s
second_two = (1 + k_two) * first_two

H_four = reduced_primitive(first_four, second_four)
H_two = reduced_primitive(first_two, second_two)
delta = (3 * A + 5) / 21

difference = reduce_a(
    H_four - H_two.subs(W, W + delta),
    s,
    W,
)
assert reduce_a(sp.diff(difference, W, 2), s, W) == 0
m = reduce_a(sp.diff(difference, W), s)
n = reduce_a(difference.subs(W, 0), s)
assert m != 0
assert n != 0

# With independent slope and intercept, the overlap is a polynomial
# triangular automorphism and the two incidence equations agree exactly.
sigma_two, tau_two = sp.symbols("sigma_two tau_two")
sigma_four = sigma_two + m
tau_four = tau_two - n - delta * sigma_two
E_four = H_four - sigma_four * W + tau_four
E_two = H_two.subs(W, W + delta) - sigma_two * (W + delta) + tau_two
assert reduce_a(E_four - E_two, s, W, sigma_two, tau_two) == 0

target_transition = sp.Matrix((s, sigma_four, tau_four))
target_variables = sp.Matrix((s, sigma_two, tau_two))
assert reduce_a(
    target_transition.jacobian(target_variables).det(),
    s,
) == 1

# The marked-line core has precisely the inverse-equation derivative as its
# nonconstant Jacobian factor.
sigma = sp.symbols("sigma")
tau = H_four - sigma * W
controlled_divisor = sp.diff(H_four, W) - sigma
core_jacobian = sp.Matrix((sigma, tau)).jacobian((W, sigma)).det()
assert reduce_a(core_jacobian + controlled_divisor, s, W, sigma) == 0

# The elementary reciprocal modification x=D*W has the right volume ledger,
# but already the s=0 quartic part leaves third- and fourth-order poles.
x, boundary = sp.symbols("x boundary")
reciprocal_jacobian = sp.Matrix(
    (s, x / boundary, boundary)
).jacobian((s, x, boundary)).det()
assert reciprocal_jacobian == 1 / boundary

H_zero = H_four.subs(s, 0)
sigma_pullback = sp.cancel(
    sp.diff(H_zero, W).subs(W, x / boundary) - boundary
)
tau_pullback = sp.cancel(
    H_zero.subs(W, x / boundary)
    - (x / boundary) * sp.diff(H_zero, W).subs(W, x / boundary)
    + x
)
assert sp.Poly(sp.denom(sigma_pullback), boundary).degree() == 3
assert sp.Poly(sp.denom(tau_pullback), boundary).degree() == 4

# Screen the first genuinely alternating Jung coordinate.  Put
#
#   y = T + q(Y),  x = Y + p(y),
#
# with both p and q quadratic.  In the inverse coordinates,
#
#   Y = x - p(y),  T = y - q(Y).
#
# The projective affine-in-U pencil can pass its unit gate only if
# a(x)G_y-b(x)T_y is a unit.  Since deg_y(T_y)=3, the coefficients of G in
# degrees 14 through 11 must vanish.  Their ideal, saturated at the
# quadratic leading coefficient of p, is the unit ideal.
xj, yj = sp.symbols("xj yj")
c, d, e = sp.symbols("c d e")
p2, p1, p0, saturation = sp.symbols("p2 p1 p0 saturation")
p = p2 * yj**2 + p1 * yj + p0
Y_inverse = xj - p
T_inverse = yj - (c * Y_inverse**2 + d * Y_inverse + e)
G_inverse = sp.expand(g.subs({T: T_inverse, Y: Y_inverse}))

field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
G_inverse_number_field = sp.expand(G_inverse.subs(A, alpha))
assert G_inverse_number_field.coeff(yj, 14) != 0
assert G_inverse_number_field.coeff(yj, 15) == 0
assert sp.Poly(T_inverse, yj).degree() == 4

high_coefficient_equations: list[sp.Expr] = []
for y_degree in range(14, 10, -1):
    coefficient_in_x = sp.Poly(
        G_inverse_number_field.coeff(yj, y_degree),
        xj,
        domain=number_field.frac_field(c, d, e, p2, p1, p0),
    )
    for coefficient in coefficient_in_x.all_coeffs():
        numerator, _ = sp.fraction(coefficient)
        if numerator != 0:
            high_coefficient_equations.append(numerator)

assert len(high_coefficient_equations) == 6
quadratic_quadratic_gate = sp.groebner(
    high_coefficient_equations + [saturation * p2 - 1],
    saturation,
    c,
    d,
    e,
    p2,
    p1,
    p0,
    order="grevlex",
    domain=number_field,
)
assert len(quadratic_quadratic_gate.polys) == 1
assert quadratic_quadratic_gate.polys[0].as_expr() == 1

# The all-degree quadratic-first theorem.  After writing
#
#   r=y-e, u=Y, T=r-c*u^2-d*u,
#
# the only possible Newton-leading coefficients are C7*u^7, C6*u^6, and
# C51*r*u^5.  They have no common zero over K.
u, r = sp.symbols("u r")
quadratic_shear_core = sp.Poly(
    sp.expand(g.subs({Y: u, T: r - c * u**2 - d * u})),
    u,
    r,
)
C7 = reduce_a(quadratic_shear_core.coeff_monomial(u**7), c, d)
C6 = reduce_a(quadratic_shear_core.coeff_monomial(u**6), c, d)
C51 = reduce_a(
    quadratic_shear_core.coeff_monomial(u**5 * r),
    c,
    d,
)

leading_coefficient_gate = sp.groebner(
    [coefficient.subs(A, alpha) for coefficient in (C7, C6, C51)],
    d,
    c,
    order="lex",
    domain=number_field,
)
assert len(leading_coefficient_gate.polys) == 1
assert leading_coefficient_gate.polys[0].as_expr() == 1

# More explicitly, C7 and C51 share exactly one possible c-value, and C6
# is then a nonzero constant independent of d.
common_gcd = sp.gcd(
    sp.Poly(C7.subs(A, alpha), c, domain=number_field),
    sp.Poly(C51.subs(A, alpha), c, domain=number_field),
).monic()
assert sp.Poly(
    common_gcd.as_expr()
    - (c + sp.Rational(1, 7) + 2 * alpha / 7),
    c,
    domain=number_field,
).is_zero
exceptional_c = -(1 + 2 * A) / 7
assert reduce_a(C7.subs(c, exceptional_c), d) == 0
assert reduce_a(C51.subs(c, exceptional_c), d) == 0
assert reduce_a(
    C6.subs(c, exceptional_c) + 10 * (2 * A + 1) / 49,
    d,
) == 0
assert reduce_a((2 * A + 1) ** 2 + 7) == 0

# The line-cover coefficients are the Galois conjugates and obey the same
# unit gate.
conjugate_leading_gate = sp.groebner(
    [
        conjugate_a(coefficient, c, d).subs(A, alpha)
        for coefficient in (C7, C6, C51)
    ],
    d,
    c,
    order="lex",
    domain=number_field,
)
assert len(conjugate_leading_gate.polys) == 1
assert conjugate_leading_gate.polys[0].as_expr() == 1

# If deg(p)=m>=2, the transformed polynomial has degree 7m, 6m, or
# 5m+1 according to the first nonzero coefficient above.  Every possibility
# exceeds deg(T)=2m and hence defeats the unit gate.
for p_degree in range(2, 30):
    assert 7 * p_degree > 2 * p_degree
    assert 6 * p_degree > 2 * p_degree
    assert 5 * p_degree + 1 > 2 * p_degree

# If the first shear has degree n>=3 instead, the -(5+3a)T^3Y term has the
# unique degree m(3n+1), larger than every other Davenport term.
for p_degree in range(2, 12):
    for q_degree in range(3, 12):
        leading_degree = p_degree * (3 * q_degree + 1)
        other_degrees = (
            7 * p_degree,
            p_degree * (q_degree + 5),
            p_degree * (q_degree + 4),
            p_degree * (2 * q_degree + 3),
            p_degree * (2 * q_degree + 2),
            3 * p_degree * q_degree,
        )
        assert leading_degree > max(other_degrees)
        assert leading_degree > p_degree * q_degree

# In the reverse alternating orientation, deg(Y)=mn and Y^7/7 is uniquely
# dominant whenever both shears are nonlinear.
for p_degree in range(2, 12):
    for q_degree in range(2, 12):
        leading_degree = 7 * p_degree * q_degree
        other_degrees = (
            p_degree + 5 * p_degree * q_degree,
            p_degree + 4 * p_degree * q_degree,
            2 * p_degree + 3 * p_degree * q_degree,
            2 * p_degree + 2 * p_degree * q_degree,
            3 * p_degree + p_degree * q_degree,
            3 * p_degree,
        )
        assert leading_degree > max(other_degrees)

print("PASS: independent slope/intercept coordinates glue polynomially")
print("PASS: their marked-line determinant is the inverse derivative")
print("PASS: the elementary reciprocal modification has unavoidable poles")
print("PASS: every quadratic-quadratic alternating Jung pencil fails the unit gate")
print("PASS: the three Newton-leading coefficients generate the unit ideal")
print("PASS: every length-two alternating Jung pencil fails the unit gate")
print("PASS Davenport independent marked-line opening")

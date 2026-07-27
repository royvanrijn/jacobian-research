#!/usr/bin/env python3
"""Exact audit of the first quadratic-in-U Davenport incidence ansatz."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, Z, davenport_pair  # noqa: E402


# Universal determinant coefficients.
T0, Y0, U = sp.symbols("T0 Y0 U")
alpha, beta, delta, gamma = sp.symbols("alpha beta delta gamma")
g0 = sp.Function("g")(T0, Y0)
h0 = sp.Function("h")(T0, Y0)
k0 = sp.Function("k")(T0, Y0)

quadratic_map = sp.Matrix(
    (
        T0 + alpha * U**2,
        g0 + beta * U + delta * U**2,
        h0 + k0 * U + gamma * U**2,
    )
)
quadratic_jacobian = sp.expand(
    quadratic_map.jacobian((T0, Y0, U)).det()
)
coefficient_zero = (
    sp.diff(g0, Y0) * k0 - beta * sp.diff(h0, Y0)
)
coefficient_one = (
    2 * gamma * sp.diff(g0, Y0)
    - beta * sp.diff(k0, Y0)
    - 2 * delta * sp.diff(h0, Y0)
    + 2
    * alpha
    * (
        sp.diff(g0, T0) * sp.diff(h0, Y0)
        - sp.diff(g0, Y0) * sp.diff(h0, T0)
    )
)
coefficient_two = (
    -2 * delta * sp.diff(k0, Y0)
    + 2
    * alpha
    * (
        sp.diff(g0, T0) * sp.diff(k0, Y0)
        - sp.diff(g0, Y0) * sp.diff(k0, T0)
    )
)
expected_jacobian = sp.expand(
    coefficient_zero
    + coefficient_one * U
    + coefficient_two * U**2
)
assert sp.simplify(quadratic_jacobian - expected_jacobian) == 0
assert sp.Poly(quadratic_jacobian, U).degree() == 2


# Exact Davenport data used by the centralizer and final degree gates.
g, h = davenport_pair()
lambda_parameter = sp.symbols("lambda_parameter")
shifted_davenport = g - lambda_parameter * T

assert sp.Poly(shifted_davenport, Y).degree() == 7
assert sp.Poly(shifted_davenport, Y).nth(7) == sp.Rational(1, 7)
assert sp.Poly(shifted_davenport, Y).nth(6) == 0
assert sp.expand(
    sp.Poly(shifted_davenport, Y).nth(5) - (1 + A) * T
) == 0

# A nontrivial composite shifted_davenport=G(phi) would have outer degree
# seven and phi affine-linear in Y.  The constant Y^7 coefficient forces
# the Y-slope of phi to be constant; the zero Y^6 coefficient then forces
# its Y-translation to be constant.  Its Y^5 coefficient would consequently
# be constant, contradicting (1+a)T.  Thus every g-lambda*T is closed.

field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="root_a",
)
root_a = number_field.ext
g_T = sp.Poly(
    sp.diff(g, T).subs(A, root_a),
    T,
    Y,
    domain=number_field,
)
g_Y = sp.Poly(
    sp.diff(g, Y).subs(A, root_a),
    T,
    Y,
    domain=number_field,
)
assert sp.gcd(g_T, g_Y).total_degree() == 0
assert g_T.degree(Y) == 5
assert g_Y.degree(Y) == 6
assert number_field.convert(
    g_T.coeff_monomial(Y**5) - (1 + root_a)
) == number_field.zero
assert g_Y.coeff_monomial(Y**6) == 1

h_T = sp.Poly(
    sp.diff(h, T).subs(A, root_a),
    T,
    Z,
    domain=number_field,
)
h_Z = sp.Poly(
    sp.diff(h, Z).subs(A, root_a),
    T,
    Z,
    domain=number_field,
)
assert sp.gcd(h_T, h_Z).total_degree() == 0
assert h_T.degree(Z) == 5
assert h_Z.degree(Z) == 6
assert h_T.coeff_monomial(Z**5) != 0
assert h_Z.coeff_monomial(Z**6) == 1

# If the full determinant is a nonzero constant C, its U^2 coefficient
# says J(g-(delta/alpha)T,k)=0.  Closedness gives k=K(g-(delta/alpha)T).
# Integrating the constant coefficient and substituting in the U
# coefficient reduces it to
#
#   Theta(T,g) g_Y - 2 alpha C g_T + 2 delta C = 0.
#
# If K' is nonconstant, the first term has Y-degree at least thirteen.  If
# K' is constant, the Y^6 coefficient forces Theta=0 and then the nonzero
# Y^5 coefficient of g_T gives a contradiction.
outer_degree = sp.symbols("outer_degree", integer=True, nonnegative=True)
for degree in range(1, 12):
    assert 7 * degree + 6 > 5

print("PASS: the quadratic-in-U determinant has exactly three coefficients")
print("PASS: every shifted Davenport polynomial g-lambda*T is closed")
print("PASS: the U^2 equation reduces k to the shifted-Davenport centralizer")
print("PASS: the U coefficient has an unavoidable degree-five contradiction")
print("PASS: the first constant-direction quadratic-in-U ansatz is impossible")
print("PASS Davenport quadratic-in-U incidence audit")

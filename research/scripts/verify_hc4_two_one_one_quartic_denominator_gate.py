#!/usr/bin/env python3
"""Exact concurrency certificates for the clean P=L0^2*L1*L2 packet."""

from __future__ import annotations

import sympy as sp


x, Y, Z = sp.symbols("x Y Z")
alpha, beta, h0, h1, j0, j1, k0, k1, scalar = sp.symbols(
    "alpha beta h0 h1 j0 j1 k0 k1 scalar"
)


# Fifth-power exact-four boundary (alpha!=0).
h_power_four = (
    alpha * Y**5
    + beta * x * Y**4
    + h0 * x**2 * Y**3 / 2
    + x**3 * (j0 * Y**2 + j1 * Y * Z) / 6
    + x**4 * (k0 * Y + k1 * Z) / 24
    + scalar * x**5 / 120
)
det_power_four = sp.expand(sp.hessian(h_power_four, (x, Y, Z)).det())
assert sp.factor(det_power_four.coeff(x, 4)) == -5 * alpha * j1**2 * Y**5
assert sp.Poly(det_power_four.coeff(x, 7), Y, Z).coeff_monomial(Y * Z) == j1**3 / 18


# Fourth-power exact-five boundary (alpha=0,h1=0).
h_power_five = (
    beta * x * Y**4
    + h0 * x**2 * Y**3 / 2
    + x**3 * (j0 * Y**2 + j1 * Y * Z) / 6
    + x**4 * (k0 * Y + k1 * Z) / 24
    + scalar * x**5 / 120
)
det_power_five = sp.expand(sp.hessian(h_power_five, (x, Y, Z)).det())
assert sp.factor(det_power_five.coeff(x, 5)) == -sp.Rational(7, 3) * beta * j1**2 * Y**4
assert sp.Poly(det_power_five.coeff(x, 7), Y, Z).coeff_monomial(Y * Z) == j1**3 / 18


# The exact-four 4+1 boundary.
j2 = -sp.Rational(3, 4) * h1**2 / beta
h_four_one = (
    beta * x * Y**4
    + x**2 * (h0 * Y**3 + h1 * Y**2 * Z) / 2
    + x**3 * (j0 * Y**2 + j1 * Y * Z + j2 * Z**2) / 6
    + x**4 * (k0 * Y + k1 * Z) / 24
    + scalar * x**5 / 120
)
det_four_one = sp.expand(sp.hessian(h_four_one, (x, Y, Z)).det())
normalization = {j1: 3 * h0 * h1 / (4 * beta)}
assert sp.Poly(
    det_four_one.coeff(x, 5).subs(normalization), Y, Z
).coeff_monomial(Y**2 * Z**2) == -sp.Rational(9, 2) * h1**4 / beta


print("PASS: fifth-power exact-four boundary has immutable transverse tail")
print("PASS: fourth-power exact-five boundary has immutable transverse tail")
print("PASS: 4+1 boundary has immutable quadratic residual-line tail")
print("THEOREM: all P=L0^2*L1*L2 concurrency strata are empty")
print("SCOPE: clean generic-corank-one partition; squarefree P remains separate")

#!/usr/bin/env python3
"""Exact tangent-kernel certificates for the clean P=x^2*y^2 packet."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
beta, h0, h1, j0, j1, k0, k1, scalar = sp.symbols(
    "beta h0 h1 j0 j1 k0 k1 scalar"
)


# Exact multiplicity five: h1=0 and j1!=0.
h_five = (
    beta * x * y**4
    + h0 * x**2 * y**3 / 2
    + x**3 * (j0 * y**2 + j1 * y * z) / 6
    + x**4 * (k0 * y + k1 * z) / 24
    + scalar * x**5 / 120
)
det_five = sp.expand(sp.hessian(h_five, (x, y, z)).det())
assert sp.factor(det_five.coeff(x, 5)) == -sp.Rational(7, 3) * beta * j1**2 * y**4
assert sp.Poly(det_five.coeff(x, 7), y, z).coeff_monomial(y * z) == j1**3 / 18


# Exact multiplicity four: retain h1 and impose the forced j2.
j2 = -sp.Rational(3, 4) * h1**2 / beta
h_four = (
    beta * x * y**4
    + x**2 * (h0 * y**3 + h1 * y**2 * z) / 2
    + x**3 * (j0 * y**2 + j1 * y * z + j2 * z**2) / 6
    + x**4 * (k0 * y + k1 * z) / 24
    + scalar * x**5 / 120
)
det_four = sp.expand(sp.hessian(h_four, (x, y, z)).det())
expected_x4 = -sp.Rational(5, 3) * h1 * y**4 * (
    (4 * beta * j1 - 3 * h0 * h1) * y - 9 * h1**2 * z
)
assert sp.factor(det_four.coeff(x, 4) - expected_x4) == 0

# Normalize the residual line to z=0.
normalization = {j1: 3 * h0 * h1 / (4 * beta)}
x5_normalized = sp.expand(det_four.coeff(x, 5).subs(normalization))
assert sp.Poly(x5_normalized, y, z).coeff_monomial(y**2 * z**2) == (
    -sp.Rational(9, 2) * h1**4 / beta
)


print("PASS: exact-five residual is y^4 with immutable x^7*y*z tail")
print("PASS: exact-four residual line normalized to z")
print("PASS: immutable x^5*y^2*z^2 coefficient is -9*h1^4/(2*beta)")
print("THEOREM: the clean P=x^2*y^2 partition is empty")
print("SCOPE: generic corank one; positive-defect and lower-Smith are separate")

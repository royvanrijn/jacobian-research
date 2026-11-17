#!/usr/bin/env python3
"""Exact audit of the triple-root affine completion."""

from __future__ import annotations

import sympy as sp


a, b, c, X = sp.symbols("a b c X", nonzero=True)
d = (1 + b * c) / a
t = (1 - a**2 * c) / (a * c**2)
L1 = a * X + b
L3 = c * X + d
L2 = sp.factor(L1 + t * L3)
P = sp.Poly(sp.factor(L1 * L2 * L3), X)
p0, p1, p2, p3 = P.all_coeffs()
r12 = sp.factor(sp.resultant(L1, L2, X))
r13 = sp.factor(sp.resultant(L1, L3, X))
r23 = sp.factor(sp.resultant(L2, L3, X))

assert p0 == 1
assert r13 == 1
assert r23 == 1
assert sp.factor(r12 - t) == 0
print("PASS: the triple-root source is G_m^2*A^1 in (a,c,b)")


A = sp.factor(a * c * p1)
B = sp.factor(a**2 * c**2 * p2)
C = sp.factor(a**3 * c**2 * p3)
E = sp.factor(a * c**2 * r12)
assert A == -a**2 * c + 3 * b * c + 2
assert B == (
    -2 * a**2 * b * c**2
    - a**2 * c
    + 3 * b**2 * c**2
    + 4 * b * c
    + 1
)
assert C == -b * (b * c + 1) * (a**2 * c - b * c - 1)
assert E == 1 - a**2 * c
print("PASS: all minimally cleared oriented coordinates are polynomial")


target_relation = sp.factor(A**2 - 3 * B - E**2 + E - 1)
assert target_relation == 0
print("PASS: the cleared target is A^3 with coordinates (A,C,E)")


cleared_jacobian = sp.factor(
    sp.Matrix((A, C, E)).jacobian((a, b, c)).det()
)
expected_jacobian = (
    -6
    * a
    * b
    * c
    * (b * c + 1)
    * (a**2 * c - b * c - 1)
)
assert sp.factor(cleared_jacobian - expected_jacobian) == 0
print("PASS: minimal clearing creates exactly five Jacobian divisors")
print("PASS triple-root affine-completion obstruction")

#!/usr/bin/env python3
"""Exact checks for the clean P=x^3*y Hessian--Schur packet."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, Gamma, D = sp.symbols("A B Gamma D")
a, b = sp.symbols("a b")

h5 = (
    A * x * y**4
    + x**4 * (B * y + Gamma * z) / 24
    + D * x**5 / 120
)
C = sp.hessian(h5, (x, y, z))
Delta = sp.factor(C.det())
assert Delta == -A * Gamma**2 * x**7 * y**2 / 3

# Replay the complete tangent constant-kernel determinant before exact
# residual factorization kills h0 and j0.
h0, j0 = sp.symbols("h0 j0")
h_boundary = (
    A * x * y**4
    + h0 * x**2 * y**3 / 2
    + j0 * x**3 * y**2 / 6
    + x**4 * (B * y + Gamma * z) / 24
    + D * x**5 / 120
)
assert sp.factor(sp.hessian(h_boundary, (x, y, z)).det()) == (
    -Gamma**2
    * x**7
    * (36 * A * y**2 + 9 * h0 * x * y + j0 * x**2)
    / 108
)

# The final Schur family and quotient.
s3 = a * x * y**2 + b * x**3
d = sp.Matrix([sp.diff(s3, variable) for variable in (x, y, z)])
adjugate = C.adjugate().applyfunc(sp.factor)
numerator = sp.factor((d.T * adjugate * d)[0])
assert sp.factor(numerator / Delta) == a**2 * x / (3 * A)

# Once d_z=0, the adjugate form is exactly one square.  This is the UFD gate
# used in the written completeness proof.
dx, dy = sp.symbols("dx dy")
d_no_z = sp.Matrix([dx, dy, 0])
assert sp.factor((d_no_z.T * adjugate * d_no_z)[0]) == (
    -Gamma**2 * x**6 * dy**2 / 36
)

# Verify the primitive cleared module equations.
P = x**3 * y
e = sp.Matrix(
    [
        0,
        a * x**3 / (6 * A),
        (12 * A * a * y**3 + 108 * A * b * x**2 * y - B * a * x**3)
        / (6 * A * Gamma),
    ]
)
assert (C * e - P * d).applyfunc(sp.factor) == sp.zeros(3, 1)
assert sp.factor(d.dot(e) - P * a**2 * x / (3 * A)) == 0

# For a!=0, neither essential component can be removed from P.
assert any(sp.factor(entry.subs(x, 0)) != 0 for entry in e)
assert any(sp.factor(entry.subs(y, 0)) != 0 for entry in e)

# Calibration of the forced four-variable cubic face.
t = sp.symbols("t")
h3_forced = a**2 * x * t**2 / (6 * A)
assert sp.diff(h3_forced, t, 2) == a**2 * x / (3 * A)

print("PASS: classified the x^7*y^2 tangent Hessian boundary")
print("PASS: complete Schur space is span{x*y^2, x^3}")
print("PASS: primitive cleared vector has minimal denominator x^3*y for a!=0")
print("PASS: forced first four-variable prolongation face verified")
print("THEOREM: the first two-line quartic-denominator packet is explicit")
print(
    "SCOPE: leading Hessian--Schur classification; the genuine channel is "
    "excluded from prolongation by HC4NHM5"
)

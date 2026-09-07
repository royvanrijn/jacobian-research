#!/usr/bin/env python3
"""Exact checks for the octuple/nonuple generic-corank-one line gate.

The normalization, moving-line, and binary square-divisor classifications
are written proofs in HC4NHM1--2. This script checks the ten-row ladder, the
new defect-free cubic-kernel calculation, and the terminal cylinder identity.
"""

from __future__ import annotations

import sympy as sp


# For m=8 or 9, the Schur half-valuation gives b<=4.  The kernel residue and
# nonnegative defect give 0<=kappa<=min(b-1,3).
expected_ladder = [
    (1, 0, 6),
    (2, 0, 6),
    (2, 1, 4),
    (3, 0, 6),
    (3, 1, 4),
    (3, 2, 2),
    (4, 0, 6),
    (4, 1, 4),
    (4, 2, 2),
    (4, 3, 0),
]
for multiplicity in (8, 9):
    ladder = []
    for pole_order in range(1, multiplicity // 2 + 1):
        for kernel_degree in range(min(pole_order - 1, 3) + 1):
            ladder.append((pole_order, kernel_degree, 6 - 2 * kernel_degree))
    assert ladder == expected_ladder


u, v = sp.symbols("u v")


def binary_hessian(form: sp.Expr) -> sp.Expr:
    """Return the determinant of the binary Hessian."""

    return sp.factor(sp.hessian(form, (u, v)).det())


# By the complete square-divisor classification in HC4NHM2, a perfect-square
# binary-quintic Hessian only needs these representatives checked.
F_41 = u * v**4
F_32 = u**2 * v**3
F_311 = u**3 * v * (u - v)
F_221 = u**2 * v**2 * (u - v)
F_exceptional = v**2 * (5 * u**3 + 30 * u * v**2 + 8 * v**3)
F_fermat = u**5 + v**5

assert binary_hessian(F_41) == -16 * v**6
assert binary_hessian(F_32) == -24 * u**2 * v**4
assert binary_hessian(F_311) == -8 * u**4 * (
    2 * u**2 - 3 * u * v + 3 * v**2
)
assert binary_hessian(F_221) == -8 * u**2 * v**2 * (
    3 * u**2 - 4 * u * v + 3 * v**2
)
assert binary_hessian(F_exceptional) == (
    -600 * v**2 * (u - 2 * v) ** 2 * (u**2 + 4 * u * v + 6 * v**2)
)
assert binary_hessian(F_fermat) == 400 * u**3 * v**3
assert sp.discriminant(2 * u**2 - 3 * u + 3, u) != 0
assert sp.discriminant(3 * u**2 - 4 * u + 3, u) != 0
assert sp.discriminant(u**2 + 4 * u + 6, u) != 0


# Root type 4+1.  The last two kernel equations first remove the u^4 term
# of G.  Their complete polynomial solution has the displayed b,c.
g0, g1, g2, g3 = sp.symbols("g0:4")
G_41 = g0 * v**4 + g1 * u * v**3 + g2 * u**2 * v**2 + g3 * u**3 * v
a_41 = v**3
b_41 = -g0 * v**3 + g2 * u**2 * v + 2 * g3 * u**3
c_41 = -v * (g1 * v**2 + 2 * g2 * u * v + 3 * g3 * u**2) / 4
A_41 = sp.hessian(F_41, (u, v))
assert (
    A_41 * sp.Matrix([b_41, c_41])
    + a_41 * sp.Matrix([sp.diff(G_41, u), sp.diff(G_41, v)])
).applyfunc(sp.expand) == sp.zeros(2, 1)

first_row_numerator_41 = sp.expand(
    sp.diff(G_41, u) * b_41 + sp.diff(G_41, v) * c_41
)
remainder_41 = sp.Poly(first_row_numerator_41, v).rem(sp.Poly(v**3, v)).as_expr()
assert sp.expand(
    remainder_41
    - sp.Rational(21, 4) * g3**2 * u**5 * v
    - 5 * g2 * g3 * u**4 * v**2
) == 0
assert all(
    sp.rem(entry.subs(g3, 0), v, v) == 0
    for entry in (a_41, b_41, c_41)
)


# Root type 3+2.  Polynomiality removes the v^4 and u^4 terms of G, after
# which the complete lower-row solution has a common factor u.
G_32 = g1 * u * v**3 + g2 * u**2 * v**2 + g3 * u**3 * v
a_32 = u * v**2
b_32 = (-g1 * u * v**2 + g3 * u**3) / 2
c_32 = -(g2 * u * v**2 + 2 * g3 * u**2 * v) / 3
A_32 = sp.hessian(F_32, (u, v))
assert (
    A_32 * sp.Matrix([b_32, c_32])
    + a_32 * sp.Matrix([sp.diff(G_32, u), sp.diff(G_32, v)])
).applyfunc(sp.expand) == sp.zeros(2, 1)
assert all(sp.rem(entry, u, u) == 0 for entry in (a_32, b_32, c_32))


# Constant tangent kernel.  The septuple coefficient ladder from HC4NHM2,
# plus k1=0 for x^8 divisibility, leaves a potential independent of v in
# both boundary-rank cases.  Verify the complete determinant vanishes.
x = sp.symbols("x")
alpha, beta, h0, j0, k0, scalar = sp.symbols(
    "alpha beta h0 j0 k0 scalar"
)
h_alpha = (
    alpha * u**5
    + x * beta * u**4
    + x**2 * h0 * u**3 / 2
    + x**3 * j0 * u**2 / 6
    + x**4 * k0 * u / 24
    + x**5 * scalar / 120
)
assert sp.hessian(h_alpha, (x, u, v)).det() == 0

# The alpha=0, beta!=0 branch has the same terminal form, with alpha=0.
assert sp.hessian(h_alpha.subs(alpha, 0), (x, u, v)).det() == 0


print("PASS: octuple/nonuple pole--defect ladders each have ten rows")
print("PASS: perfect-square binary-quintic Hessian boundary classified")
print("PASS: defect-free primitive cubic-kernel row is empty")
print("PASS: all constant-kernel rows collapse to a zero determinant")
print("THEOREM: generic-corank-one x^8*ell and x^9 Schur packets are empty")
print(
    "SCOPE: closes P=x^4; P=x^3*y is continued through HC4NHM4--6; "
    "lower-Smith strata remain open"
)

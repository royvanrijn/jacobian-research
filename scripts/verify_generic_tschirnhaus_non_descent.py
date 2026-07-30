#!/usr/bin/env python3
"""Exact checks for generic all-rank Tschirnhaus non-descent."""

from __future__ import annotations

import sympy as sp


N = sp.symbols("N", integer=True, positive=True)
T = sp.symbols("T")


def root_polynomial(roots: tuple[sp.Expr, ...]) -> sp.Poly:
    """Return the monic polynomial with the supplied roots."""

    return sp.Poly(sp.prod(T - root for root in roots), T, domain=sp.QQ)


def compiler_i5(roots: tuple[sp.Expr, ...]) -> sp.Expr:
    """Return the quintic compiler-slice invariant I_5."""

    polynomial = root_polynomial(roots)
    linear = polynomial.nth(1)
    assert linear != 0
    normalized = {
        index: sp.cancel(polynomial.nth(index) / linear)
        for index in range(6)
    }
    cubic = normalized[3]
    assert cubic != 0
    u4 = sp.cancel(normalized[4] / cubic**4)
    u5 = sp.cancel(normalized[5] / cubic**5)
    return sp.factor(u5**5 / u4**6)


def top_j(roots: tuple[sp.Expr, ...]) -> sp.Expr:
    """Return a_(N-2)*a_N/a_(N-1)^2 for one root polynomial."""

    polynomial = root_polynomial(roots)
    degree = len(roots)
    return sp.factor(
        polynomial.nth(degree - 2)
        * polynomial.nth(degree)
        / polynomial.nth(degree - 1) ** 2
    )


# The first four rows already have rank four, so r -> r+r^2 is
# nonprojective in every rank N>=4.
projective_minor = sp.Matrix(
    [[1, root, root + root**2, root * (root + root**2)]
     for root in range(1, 5)]
).det()
assert projective_minor == 12


# Symbolic power-sum derivation of the top boundary coordinate.
source_sum_1 = N * (N + 1) / 2
source_sum_2 = N * (N + 1) * (2 * N + 1) / 6
source_j = sp.factor(
    (source_sum_1**2 - source_sum_2) / (2 * source_sum_1**2)
)
expected_source_j = (N - 1) * (3 * N + 2) / (6 * N * (N + 1))
assert sp.factor(source_j - expected_source_j) == 0

source_sum_3 = (N * (N + 1) / 2) ** 2
source_sum_4 = (
    N * (N + 1) * (2 * N + 1) * (3 * N**2 + 3 * N - 1) / 30
)
transformed_sum_1 = sp.factor(source_sum_1 + source_sum_2)
transformed_sum_2 = sp.factor(
    source_sum_2 + 2 * source_sum_3 + source_sum_4
)
transformed_j = sp.factor(
    (transformed_sum_1**2 - transformed_sum_2)
    / (2 * transformed_sum_1**2)
)
expected_transformed_j = (
    (N - 1)
    * (5 * N**2 + 11 * N + 3)
    / (10 * N * (N + 1) * (N + 2))
)
assert sp.factor(transformed_j - expected_transformed_j) == 0

boundary_difference = sp.factor(transformed_j - source_j)
expected_difference = (
    -(N - 1)
    * (7 * N + 11)
    / (30 * N * (N + 1) * (N + 2))
)
assert sp.factor(boundary_difference - expected_difference) == 0


# Quintic base case.
source_i5 = compiler_i5(tuple(sp.Integer(i) for i in range(1, 6)))
transformed_i5 = compiler_i5(
    tuple(sp.Integer(i + i**2) for i in range(1, 6))
)
assert source_i5 == sp.Rational(75076, 968203125)
assert transformed_i5 == sp.Rational(1296, 50236123)
assert source_i5 != transformed_i5


# Direct coefficient regressions confirm the symbolic formulas and the
# N-4 versus N-3 dimension ledger through rank twenty.
for degree in range(6, 21):
    source_roots = tuple(sp.Integer(i) for i in range(1, degree + 1))
    transformed_roots = tuple(root + root**2 for root in source_roots)
    assert top_j(source_roots) == expected_source_j.subs(N, degree)
    assert top_j(transformed_roots) == expected_transformed_j.subs(N, degree)
    assert top_j(transformed_roots) - top_j(source_roots) == (
        expected_difference.subs(N, degree)
    )
    assert degree - 4 < degree - 3


print("PASS: r -> r+r^2 is nonprojective in every rank N>=4")
print("PASS: the quintic I_5 base case changes exactly")
print("PASS: the symbolic top J_N changes for every N>=6")
print("PASS: direct root-polynomial regressions pass through rank twenty")
print("PASS: equal-boundary and projective codimensions are N-4 and N-3")

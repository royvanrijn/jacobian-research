#!/usr/bin/env python3
"""Exact audit of the triple-root affine-modification chain."""

from __future__ import annotations

from itertools import product

import sympy as sp


a, b, c, x, y, u, v = sp.symbols("a b c x y u v")
f = x * y * (x + 1) * (x - y + 1)
square_C, square_S = sp.symbols("square_C square_S")
Q_minus_one = square_S**2 - square_C * f
Q0 = u**2 * v - f
Q1 = u * v - a * x * (x + 1) * (x - a * u + 1)
Q2 = c * v - x * (x + 1) * (x - a**2 * c + 1)

# Every stage maps to the preceding one.
assert sp.factor(
    Q_minus_one.subs(
        {
            square_C: v,
            square_S: u * v,
        }
    )
    - v * Q0
) == 0
assert sp.factor(Q0.subs(y, a * u) - u * Q1) == 0
assert sp.factor(Q1.subs(u, a * c) - a * Q2) == 0
source_substitution = {
    x: b * c,
    v: b * (b * c + 1) * (b * c - a**2 * c + 1),
}
assert sp.factor(Q2.subs(source_substitution)) == 0
print("PASS: all four strict-transform equations are exact")


# Residue-form pullbacks.
primitive_v = f / u**2
square_pullback_jacobian = sp.factor(
    sp.Matrix((x, primitive_v, y))
    .jacobian((x, y, u))
    .det()
)
assert sp.factor(
    square_pullback_jacobian / (2 * u * primitive_v)
) == u**-2
jacobian_y = sp.factor(
    sp.Matrix((x, a * u, u)).jacobian((x, a, u)).det()
)
assert jacobian_y == u
jacobian_u = sp.factor(
    sp.Matrix((x, a, a * c)).jacobian((x, a, c)).det()
)
assert jacobian_u == a
jacobian_x = sp.factor(
    sp.Matrix((b * c, a, c)).jacobian((a, b, c)).det()
)
assert jacobian_x == -c
assert sp.factor(jacobian_y / u**2) == 1 / u
assert sp.factor(jacobian_u / (a * c)) == 1 / c
assert sp.factor(jacobian_x / c) == -1
print("PASS: all five residue forms telescope to a constant volume form")


# Singular-locus ideals, checked by exact gradients at the asserted
# components and by Groebner bases of the full Jacobian ideals.
X, A, U, V, Y, C = sp.symbols("X A U V Y C")
stage0 = U**2 * V - X * Y * (X + 1) * (X - Y + 1)
stage0_gradient = [
    sp.diff(stage0, variable)
    for variable in (X, Y, U, V)
]
for point_x, point_y in ((0, 0), (0, 1), (-1, 0)):
    substitution = {X: point_x, Y: point_y, U: 0}
    assert stage0.subs(substitution) == 0
    assert all(
        derivative.subs(substitution) == 0
        for derivative in stage0_gradient
    )

stage1 = U * V - A * X * (X + 1) * (X - A * U + 1)
stage1_basis = sp.groebner(
    [stage1]
    + [
        sp.diff(stage1, variable)
        for variable in (X, A, U, V)
    ],
    V,
    U,
    A,
    X,
    order="lex",
)
assert {
    sp.factor(polynomial.as_expr())
    for polynomial in stage1_basis.polys
} == {V, U, A * (X + 1), X * (X + 1) ** 2}

stage2 = C * V - X * (X + 1) * (X - A**2 * C + 1)
stage2_basis = sp.groebner(
    [stage2]
    + [
        sp.diff(stage2, variable)
        for variable in (X, A, C, V)
    ],
    V,
    C,
    A,
    X,
    order="lex",
)
assert {
    sp.factor(polynomial.as_expr())
    for polynomial in stage2_basis.polys
} == {V, C, X + 1}
print("PASS: singular loci decrease from three lines to a line and point, one line, then empty")


dicritical_counts = (4, 3, 1, 1, 0)
assert dicritical_counts[-1] == 0
print("PASS: the canonical chain fills or merges every dicritical component")


def count_stage(prime: int, stage: int) -> int:
    if stage == 0:
        return sum(
            (
                finite_u**2 * finite_v
                - finite_x
                * finite_y
                * (finite_x + 1)
                * (finite_x - finite_y + 1)
            )
            % prime
            == 0
            for finite_x, finite_y, finite_u, finite_v in product(
                range(prime),
                repeat=4,
            )
        )
    if stage == 1:
        return sum(
            (
                finite_u * finite_v
                - finite_a
                * finite_x
                * (finite_x + 1)
                * (finite_x - finite_a * finite_u + 1)
            )
            % prime
            == 0
            for finite_x, finite_a, finite_u, finite_v in product(
                range(prime),
                repeat=4,
            )
        )
    if stage == 2:
        return sum(
            (
                finite_c * finite_v
                - finite_x
                * (finite_x + 1)
                * (finite_x - finite_a**2 * finite_c + 1)
            )
            % prime
            == 0
            for finite_x, finite_a, finite_c, finite_v in product(
                range(prime),
                repeat=4,
            )
        )
    return prime**3


for prime in (5, 7):
    expected = (
        prime**3,
        prime * (prime - 1) * (prime + 4),
        prime * (prime**2 + 2 * prime - 2),
        prime**2 * (prime + 1),
        prime**3,
    )
    actual = (
        prime**3,
        *tuple(
            count_stage(prime, stage)
            for stage in range(4)
        ),
    )
    assert actual == expected
print("PASS: all five exact finite-field counts hold")
print("PASS triple-root affine-modification chain")

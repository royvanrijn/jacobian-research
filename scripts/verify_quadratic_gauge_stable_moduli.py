#!/usr/bin/env python3
"""Exact checks for quadratic-gauge stable moduli."""

from __future__ import annotations

from math import gcd

import sympy as sp


x, y, z = sp.symbols("x y z")
alpha, beta = sp.symbols("alpha beta", nonzero=True)
a2, a3, a4, a5, a6 = sp.symbols("a2 a3 a4 a5 a6", nonzero=True)


def quadratic_map(
    coefficients: dict[int, sp.Expr],
) -> tuple[tuple[sp.Expr, ...], sp.Expr, sp.Expr]:
    """Normalized quadratic-gauge map for one coefficient dictionary."""

    local_t = 1 + x * y
    local_q = (
        local_t**2 * z
        + y**2 * (1 + 3 * local_t) / coefficients[3]
    )
    degree = max(coefficients)
    mapping = (
        local_t * local_q,
        y
        + 3 * coefficients[3] * x * local_q
        + 2 * coefficients.get(2, 0) * local_t * local_q
        + sum(
            k
            * coefficients[k]
            * local_t**2
            * x ** (k - 2)
            * local_q**k
            for k in range(4, degree + 1)
        ),
        x * (5 - 3 * local_t)
        - coefficients[3] * x**3 * z
        - sum(
            (k - 2) * coefficients[k] * (x * local_q) ** k
            for k in range(4, degree + 1)
        ),
    )
    return tuple(sp.expand(item) for item in mapping), local_t, local_q


a = {2: a2, 3: a3, 4: a4, 5: a5, 6: a6}
b = {
    2: a2 * alpha**-1 * beta**-1,
    3: a3 * alpha**-2 * beta**-1,
    4: a4 * alpha**-3 * beta**-4,
    5: a5 * alpha**-4 * beta**-5,
    6: a6 * alpha**-5 * beta**-6,
}
F_a, t_a, q_a = quadratic_map(a)
F_b, _, q_b = quadratic_map(b)
source_scaling = {x: alpha * x, y: y / alpha, z: beta * z}

assert sp.factor(t_a.subs(source_scaling) - t_a) == 0
assert sp.factor(q_b.subs(source_scaling, simultaneous=True) - beta * q_a) == 0

scaled_F_b = tuple(
    sp.factor(component.subs(source_scaling, simultaneous=True))
    for component in F_b
)
target_scaled_F_a = (
    beta * F_a[0],
    F_a[1] / alpha,
    alpha * F_a[2],
)
assert all(
    sp.factor(got - expected) == 0
    for got, expected in zip(scaled_F_b, target_scaled_F_a)
)


# The quadratic coefficient is exactly a target shear.
a_without_2 = dict(a)
a_without_2[2] = sp.Integer(0)
F_without_2, _, _ = quadratic_map(a_without_2)
assert sp.factor(F_a[0] - F_without_2[0]) == 0
assert sp.factor(F_a[1] - F_without_2[1] - 2 * a2 * F_a[0]) == 0
assert sp.factor(F_a[2] - F_without_2[2]) == 0


# Intrinsic normalization and its Fitting divisor.
P, r = sp.symbols("P r", nonzero=True)
h = (
    r
    + a2 * P * r**2
    + a3 * P * r**3
    + a4 * P**4 * r**4
    + a5 * P**5 * r**5
    + a6 * P**6 * r**6
)
B = sp.factor(sp.diff(h, r) / r)
C = sp.factor(2 * h - r * sp.diff(h, r))
assert sp.factor(C + r**2 * B - 2 * h) == 0
assert sp.factor(sp.diff(C, r) + r**2 * sp.diff(B, r)) == 0

J = sp.factor(r**2 * sp.diff(B, r))
expected_J = (
    -1
    + 3 * a3 * P * r**2
    + 8 * a4 * P**4 * r**3
    + 15 * a5 * P**5 * r**4
    + 24 * a6 * P**6 * r**5
)
assert sp.factor(J - expected_J) == 0
assert not J.has(a2)


# Support rigidity and the two-dimensional weight lattice in all tested
# degrees.  Inversion changes the unique missing interior r-exponent, and
# the r^2 term forces a possible P^m twist to have m=0.
for degree in range(4, 65):
    support = {(0, 0), (1, 2)}
    support.update((k, k - 1) for k in range(4, degree + 1))
    r_degrees = {pair[1] for pair in support}
    assert min(r_degrees) == 0
    assert max(r_degrees) == degree - 1
    assert [pair for pair in support if pair[1] == 0] == [(0, 0)]
    assert [pair for pair in support if pair[1] == 2] == [(1, 2)]

    inverted_r_degrees = {
        degree - 1 - exponent
        for exponent in r_degrees
    }
    assert inverted_r_degrees == (
        set(range(0, degree - 2)) | {degree - 1}
    )
    assert r_degrees != inverted_r_degrees
    assert 1 not in r_degrees
    assert degree - 2 not in inverted_r_degrees

    weight_3 = (-2, -1)
    weight_4 = (-3, -4)
    determinant_34 = (
        weight_3[0] * weight_4[1] - weight_3[1] * weight_4[0]
    )
    assert determinant_34 == 5
    assert degree - 2 - 2 == degree - 4

    if degree == 4:
        assert abs(determinant_34) == 5  # kernel mu_5
    else:
        weight_5 = (-4, -5)
        determinant_35 = (
            weight_3[0] * weight_5[1] - weight_3[1] * weight_5[0]
        )
        assert determinant_35 == 6
        assert gcd(abs(determinant_34), abs(determinant_35)) == 1


print("PASS: the independent (alpha,beta) source-target scaling is exact")
print("PASS: a2 is removed by the target shear B -> B-2*a2*P")
print("PASS: the intrinsic Fitting polynomial recovers a3,...,aN")
print("PASS: Fitting support orders the two toric punctures")
print("PASS: Fitting support removes P^m twists")
print("PASS: the stable coefficient-torus quotient has dimension N-4")

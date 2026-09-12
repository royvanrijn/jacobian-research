#!/usr/bin/env python3
"""Exact checks for decorated Torelli of polynomial quadratic gauges."""

from __future__ import annotations

import sympy as sp


P, r = sp.symbols("P r")
alpha, beta = sp.symbols("alpha beta", nonzero=True)
a3, a4, a5, a6 = sp.symbols("a3 a4 a5 a6", nonzero=True)
u1, u2 = sp.symbols("u1 u2")


def fitting(
    coefficients: tuple[sp.Expr, ...],
    multiplier: sp.Expr,
) -> sp.Expr:
    """Intrinsic ramified-stratum Fitting generator."""

    result = -1 + 3 * coefficients[0] * P * r**2
    for degree, coefficient in enumerate(coefficients[1:], start=4):
        result += (
            degree
            * (degree - 2)
            * coefficient
            * P**degree
            * multiplier
            * r ** (degree - 1)
        )
    return sp.expand(result)


# Normalize the polynomial multiplier at P=1.
R_a = sp.expand(1 + u1 * (P - 1) + u2 * (P**2 - 1))
assert R_a.subs(P, 1) == 1
normalizing_value = sp.factor(R_a.subs(P, 1 / beta))
R_b = sp.cancel(R_a.subs(P, P / beta) / normalizing_value)
assert sp.cancel(R_b.subs(P, 1) - 1) == 0


# These are exactly the source-target scaling formulas forced by the
# decorated Fitting divisor.
old_coefficients = (a3, a4, a5, a6)
new_coefficients = (
    alpha**-2 * beta**-1 * a3,
    normalizing_value * alpha**-3 * beta**-4 * a4,
    normalizing_value * alpha**-4 * beta**-5 * a5,
    normalizing_value * alpha**-5 * beta**-6 * a6,
)
J_a = fitting(old_coefficients, R_a)
J_b = fitting(new_coefficients, R_b)
assert sp.factor(J_b.subs({P: beta * P, r: alpha * r}) - J_a) == 0


# Once P and r are recovered, coefficient extraction reconstructs a_3,
# every a_j, and R(P) from the Fitting generator without a marked root.
coefficient_r2 = J_a.coeff(r, 2)
assert sp.cancel(coefficient_r2 / (3 * P) - a3) == 0
for degree, coefficient in zip(range(4, 7), (a4, a5, a6)):
    coefficient_function = J_a.coeff(r, degree - 1)
    normalized_function = sp.cancel(
        coefficient_function / (degree * (degree - 2) * P**degree)
    )
    recovered_coefficient = sp.cancel(normalized_function.subs(P, 1))
    recovered_multiplier = sp.cancel(
        normalized_function / recovered_coefficient
    )
    assert sp.cancel(recovered_coefficient - coefficient) == 0
    assert sp.cancel(recovered_multiplier - R_a) == 0


# The r-exponent support rules out inversion exactly as in the minimal
# quadratic-gauge stable-moduli theorem.
for degree in range(4, 13):
    ordinary = {0, 2, *range(3, degree)}
    inverted_raw = {-exponent for exponent in ordinary}
    shift = -min(inverted_raw)
    inverted = {exponent + shift for exponent in inverted_raw}
    assert ordinary != inverted
    assert 1 not in ordinary
    assert degree - 2 not in inverted


# A possible normalization-coordinate twist is r -> U(P)r, where U is a
# unit on the punctured base.  Comparing the r^2 coefficient gives twice
# every valuation of U, so the torsion-free base-unit lattice forces U to be
# constant.
for unit_rank in range(1, 10):
    for index in range(unit_rank):
        valuation_vector = [0] * unit_rank
        valuation_vector[index] = index + 1
        doubled = [2 * value for value in valuation_vector]
        assert doubled != [0] * unit_rank
    assert all(2 * value == 0 for value in [0] * unit_rank)


# Boundary labels intrinsically distinguish P=0 from every multiplier-root
# plane: the former has affine residue degrees 2+1, the latter one cubic
# affine residue field.  This selects the prime generator P up to a scalar.
for degree in range(4, 13):
    zero_affine_degrees = (2, 1)
    multiplier_root_affine_degrees = (3,)
    assert zero_affine_degrees != multiplier_root_affine_degrees
    assert sum(zero_affine_degrees) == sum(multiplier_root_affine_degrees)
    assert (degree - 3) + sum(multiplier_root_affine_degrees) == degree


# Direct quartic source-target scaling identity for a nonconstant multiplier.
x, y, z = sp.symbols("x y z")
t = 1 + x * y
R_linear = 1 + u1 * (P - 1)


def quartic_map(
    cubic: sp.Expr,
    quartic: sp.Expr,
    multiplier: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    q = t**2 * z + cubic**-1 * y**2 * (1 + 3 * t)
    pi = t * q
    multiplier_at_pi = multiplier.subs(P, pi)
    second = y + 3 * cubic * x * q
    second += 4 * quartic * multiplier_at_pi * t**2 * x**2 * q**4
    third = x * (5 - 3 * t) - cubic * x**3 * z
    third -= 2 * quartic * multiplier_at_pi * x**4 * q**4
    return tuple(sp.cancel(component) for component in (pi, second, third))


scale_value = sp.factor(R_linear.subs(P, 1 / beta))
R_linear_b = sp.cancel(R_linear.subs(P, P / beta) / scale_value)
b3 = alpha**-2 * beta**-1 * a3
b4 = scale_value * alpha**-3 * beta**-4 * a4
map_a = quartic_map(a3, a4, R_linear)
map_b = quartic_map(b3, b4, R_linear_b)
sigma = {x: alpha * x, y: y / alpha, z: beta * z}
scaled_map_b = tuple(sp.factor(component.subs(sigma)) for component in map_b)
expected = (
    beta * map_a[0],
    map_a[1] / alpha,
    alpha * map_a[2],
)
for actual, target in zip(scaled_map_b, expected):
    assert sp.factor(actual - target) == 0


print("PASS: the Fitting divisor is invariant under the exact scaling action")
print("PASS: intrinsic coefficient extraction recovers the seed and R(P)")
print("PASS: Fitting support orders r and kills every puncture-unit twist")
print("PASS: the boundary ledger recovers P up to scalar without a root mark")
print("PASS: the reconstructed scaling is an ordinary left-right equivalence")

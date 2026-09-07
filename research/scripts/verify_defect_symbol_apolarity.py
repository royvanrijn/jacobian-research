#!/usr/bin/env python3
"""Exact certificates for the defect-symbol/apolarity package."""

from __future__ import annotations

from itertools import product

import sympy as sp


def complete_intersection_hilbert(degrees: tuple[int, ...]) -> tuple[int, ...]:
    """Coefficients of product_i (1 + T + ... + T^(d_i-1))."""

    coefficients = [1]
    for degree in degrees:
        updated = [0] * (len(coefficients) + degree - 1)
        for left_index, left in enumerate(coefficients):
            for right_index in range(degree):
                updated[left_index + right_index] += left
        coefficients = updated
    return tuple(coefficients)


def monomial_inverse_dimension(generators: tuple[tuple[int, ...], ...]) -> int:
    """Dimension of the derivative closure of monomial DP generators."""

    derivatives: set[tuple[int, ...]] = set()
    for generator in generators:
        ranges = (range(exponent + 1) for exponent in generator)
        derivatives.update(product(*ranges))
    return len(derivatives)


def killed_by_monomial_ideal(
    generators: tuple[tuple[int, ...], ...],
    ideal_generators: tuple[tuple[int, ...], ...],
) -> bool:
    """Check that every monomial operator kills every DP generator."""

    return all(
        any(operator[i] > potential[i] for i in range(len(operator)))
        for operator in ideal_generators
        for potential in generators
    )


# A pair of binary quadrics and its cross-product inverse system.
a, b, c, d, e, f = sp.symbols("a b c d e f")
dual_a = b * f - c * e
dual_b = c * d - a * f
dual_c = a * e - b * d

assert sp.expand(a * dual_a + b * dual_b + c * dual_c) == 0
assert sp.expand(d * dual_a + e * dual_b + f * dual_c) == 0

x = sp.symbols("x")
q = a * x**2 + b * x + c
r = d * x**2 + e * x + f
resultant = sp.factor(sp.resultant(q, r, x))
dual_discriminant = sp.factor(dual_b**2 - dual_a * dual_c)
assert sp.factor(resultant - dual_discriminant) == 0


# Specialize to the degree-42 quadratic Kuranishi pair after removing w1.
e1, e2, t = sp.symbols("e1 e2 t")
q42 = (
    -sp.Rational(9, 2) * e2,
    9 * e2**2 + sp.Rational(3, 2) * e1,
    -sp.Rational(9, 2) * e2**3
    + sp.Rational(15, 4) * e1 * e2
    + sp.Rational(3, 4) * t,
)
r42 = (
    -sp.Rational(3, 4),
    sp.Rational(3, 2) * e2,
    -sp.Rational(3, 4) * e2**2 + sp.Rational(3, 4) * e1,
)
dual42 = (
    sp.factor(q42[1] * r42[2] - q42[2] * r42[1]),
    sp.factor(q42[2] * r42[0] - q42[0] * r42[2]),
    sp.factor(q42[0] * r42[1] - q42[1] * r42[0]),
)
expected_dual42 = (
    -sp.Rational(9, 8) * (-e1**2 + e2 * t),
    -sp.Rational(9, 16) * (-e1 * e2 + t),
    sp.Rational(9, 8) * e1,
)
assert all(
    sp.factor(actual - expected) == 0
    for actual, expected in zip(dual42, expected_dual42, strict=True)
)
delta = (t + e1 * e2) ** 2 - 4 * e1**3
assert sp.factor(
    dual42[1] ** 2
    - dual42[0] * dual42[2]
    - sp.Rational(81, 256) * delta
) == 0


# The terminal binary cubics on the w1=0 branch, after removing w2.
u, v = sp.symbols("u v")
p3 = (
    (sp.Rational(15, 4) * e2**2 + sp.Rational(5, 8) * e1) * u**3
    + (
        -sp.Rational(45, 4) * e2**3
        + sp.Rational(75, 8) * e1 * e2
        + sp.Rational(15, 8) * t
    )
    * u**2
    * v
    + (
        sp.Rational(45, 4) * e2**4
        - sp.Rational(165, 8) * e1 * e2**2
        - sp.Rational(15, 8) * e1**2
        + sp.Rational(75, 8) * e2 * t
    )
    * u
    * v**2
    + (
        -sp.Rational(15, 4) * e2**5
        + sp.Rational(85, 8) * e1 * e2**3
        - sp.Rational(5, 2) * e1**2 * e2
        - sp.Rational(55, 8) * e2**2 * t
        - sp.Rational(5, 4) * e1 * t
    )
    * v**3
)
q3 = (
    sp.Rational(5, 8) * e2 * u**3
    + (-sp.Rational(15, 8) * e2**2 + sp.Rational(15, 8) * e1)
    * u**2
    * v
    + (
        sp.Rational(15, 8) * e2**3
        - sp.Rational(15, 4) * e1 * e2
        + sp.Rational(15, 8) * t
    )
    * u
    * v**2
    + (
        -sp.Rational(5, 8) * e2**4
        + sp.Rational(15, 8) * e1 * e2**2
        - sp.Rational(5, 8) * e1**2
        - sp.Rational(5, 4) * e2 * t
    )
    * v**3
)


def binary_cubic_coefficients(form: sp.Expr) -> list[sp.Expr]:
    expanded = sp.expand(form)
    return [
        expanded.coeff(u, 3 - index).coeff(v, index) for index in range(4)
    ]


p_coefficients = binary_cubic_coefficients(p3)
q_coefficients = binary_cubic_coefficients(q3)
catalecticant = sp.Matrix(
    [
        p_coefficients + [0],
        [0] + p_coefficients,
        q_coefficients + [0],
        [0] + q_coefficients,
    ]
)
dual_quartic = []
for column in range(5):
    retained = [index for index in range(5) if index != column]
    dual_quartic.append(
        sp.factor((-1) ** column * catalecticant[:, retained].det())
    )
assert all(
    sp.factor(entry) == 0
    for entry in catalecticant * sp.Matrix(dual_quartic)
)

factor_a = 4 * e1**3 - e1**2 * e2**2 + e2**3 * t - 6 * e1 * e2 * t
factor_b = (
    e1**4 * e2**4
    - 2 * e1**2 * e2**5 * t
    + e2**6 * t**2
    - 8 * e1**5 * e2**2
    + 20 * e1**3 * e2**3 * t
    - 12 * e1 * e2**4 * t**2
    + 16 * e1**6
    - 48 * e1**4 * e2 * t
    + 27 * e1**2 * e2**2 * t**2
    + 9 * e2**3 * t**3
    + 36 * e1**3 * t**2
    - 54 * e1 * e2 * t**3
    + 27 * t**4
)
f0, f1, f2, f3, f4 = dual_quartic
quartic_i = sp.factor(f0 * f4 - 4 * f1 * f3 + 3 * f2**2)
quartic_j = sp.factor(
    f0 * f2 * f4
    + 2 * f1 * f2 * f3
    - f0 * f3**2
    - f1**2 * f4
    - f2**3
)
expected_i = (
    sp.Rational(3 * 5**8, 2**24)
    * factor_a
    * (factor_a + 3 * t**2)
    * factor_b
)
expected_j = -sp.Rational(5**12, 2**36) * factor_a**2 * factor_b**2
expected_discriminant = (
    sp.Rational(3**6 * 5**24, 2**72)
    * t**6
    * factor_a**3
    * factor_b**3
)
assert sp.factor(quartic_i - expected_i) == 0
assert sp.factor(quartic_j - expected_j) == 0
assert sp.factor(quartic_i**3 - 27 * quartic_j**2 - expected_discriminant) == 0
assert sp.factor((factor_a + 3 * t**2) ** 3 - factor_a * factor_b - 27 * t**6) == 0


# Complete-intersection Hilbert functions and finite-jet closure exponents.
assert complete_intersection_hilbert((2, 2)) == (1, 2, 1)
assert complete_intersection_hilbert((3, 3)) == (1, 2, 3, 2, 1)
assert complete_intersection_hilbert((5,)) == (1, 1, 1, 1, 1)
assert complete_intersection_hilbert((4, 2)) == (1, 2, 2, 2, 1)


# Monomial Gorenstein models.  Containment in the annihilator together with
# equality of inverse-system and quotient dimensions proves equality.
monomial_models = (
    # (inverse generators, ideal generators, quotient length)
    (((1,),), ((2,),), 2),
    (((1, 1),), ((2, 0), (0, 2)), 4),
    (((4,),), ((5,),), 5),
    (((3, 1),), ((4, 0), (0, 2)), 8),
)
for inverse_generators, ideal_generators, quotient_length in monomial_models:
    assert killed_by_monomial_ideal(inverse_generators, ideal_generators)
    assert monomial_inverse_dimension(inverse_generators) == quotient_length

# Three independent primitive mergers give the Boolean Hilbert vector
# (1,3,3,1), with squarefree-monomial inverse system.
boolean_three = ((1, 1, 1),)
assert killed_by_monomial_ideal(
    boolean_three,
    ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
)
assert monomial_inverse_dimension(boolean_three) == 8


# Two nongorenstein boundary models require two inverse-system generators.
common_tangent_inverse = ((1, 0), (0, 2))
common_tangent_ideal = ((2, 0), (1, 1), (0, 3))
assert killed_by_monomial_ideal(common_tangent_inverse, common_tangent_ideal)
assert monomial_inverse_dimension(common_tangent_inverse) == 4

merged_failure_inverse = ((2, 0), (0, 1))
merged_failure_ideal = ((3, 0), (1, 1), (0, 2))
assert killed_by_monomial_ideal(merged_failure_inverse, merged_failure_ideal)
assert monomial_inverse_dimension(merged_failure_inverse) == 4


print("PASS: a binary quadratic pencil is dual to its coefficient cross product")
print("PASS: the dual discriminant equals the binary quadratic resultant")
print("PASS: the degree-42 dual discriminant is (81/256)*Delta")
print("PASS: the degree-42 dual quartic has factored I, J, and discriminant")
print("PASS: complete-intersection Hilbert vectors give closure powers 3 and 5")
print("PASS: the primitive and degree-30 inverse-system models have exact lengths")
print("PASS: both nongorenstein boundary models require two dual generators")

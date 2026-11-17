#!/usr/bin/env python3
"""Exact audit of all linear three-factor coefficient slices."""

from __future__ import annotations

from itertools import product
from math import gcd

import sympy as sp


a, b, c, d, t, X = sp.symbols("a b c d t X")
L1 = a * X + b
L3 = c * X + d
L2 = L1 + t * L3
coefficients = sp.Poly(sp.expand(L1 * L2 * L3), X).all_coeffs()
determinant_relation = a * d - b * c - 1


def reduce_on_sl2(expression: sp.Expr) -> sp.Poly:
    return sp.rem(
        sp.Poly(sp.expand(expression), b),
        sp.Poly(determinant_relation, b),
    )


# Representatives of the triple, double, and squarefree dual-cubic strata.
p0, p1, _, p3 = coefficients
triple_B0 = a**2 * c
triple_B1 = a * c**2
assert reduce_on_sl2(p0 - (triple_B0 + t * triple_B1)).is_zero

double_B0 = a * (3 * a * d - 2)
double_B1 = c * (3 * a * d - 1)
assert reduce_on_sl2(p1 - (double_B0 + t * double_B1)).is_zero

squarefree_B0 = a**2 * c + b**2 * d
squarefree_B1 = a * c**2 + b * d**2
assert reduce_on_sl2(
    p0 + p3 - (squarefree_B0 + t * squarefree_B1)
).is_zero
print("PASS: triple, double, and squarefree slice formulas hold on SL_2")


L = sp.symbols("L")
class_sl2 = L**3 - L
class_triple = sp.expand(class_sl2 - 2 * L * (L - 1))
class_double = sp.expand(
    class_sl2
    - ((L - 1) * L + (L - 1) ** 2)
    + L * (L + (L - 1))
)
S = sp.symbols("S")
class_squarefree = L**3 - 2 * L - S
assert class_triple == L**3 - 2 * L**2 + L
assert class_double == L**3 + L - 1
assert class_squarefree.subs(S, 2) == L**3 - 2 * L - 2
print("PASS: the three geometric classes are exact")


def source_count(
    prime: int,
    coefficient_vector: tuple[int, int, int, int],
) -> int:
    count = 0
    for a0, b0, c0, d0 in product(range(prime), repeat=4):
        if (a0 * d0 - b0 * c0) % prime != 1:
            continue
        base = (
            a0 * a0 * c0,
            a0 * a0 * d0 + 2 * a0 * b0 * c0,
            2 * a0 * b0 * d0 + b0 * b0 * c0,
            b0 * b0 * d0,
        )
        slope = (
            a0 * c0 * c0,
            2 * a0 * c0 * d0 + b0 * c0 * c0,
            a0 * d0 * d0 + 2 * b0 * c0 * d0,
            b0 * d0 * d0,
        )
        base_value = sum(
            coefficient_vector[index] * base[index]
            for index in range(4)
        ) % prime
        slope_value = sum(
            coefficient_vector[index] * slope[index]
            for index in range(4)
        ) % prime
        count += (
            1
            if slope_value
            else prime if base_value == 1 else 0
        )
    return count


for prime in (5, 7, 11):
    assert source_count(prime, (1, 0, 0, 0)) == (
        prime**3 - 2 * prime**2 + prime
    )
    assert source_count(prime, (0, 1, 0, 0)) == (
        prime**3 + prime - 1
    )
    artin_points = sum(
        1
        for value in range(prime)
        if (value**2 - value + 1) % prime == 0
    )
    assert source_count(prime, (1, 0, 0, 1)) == (
        prime**3 - 2 * prime - artin_points
    )
print("PASS: finite-field counts realize all three motivic classes")


q = sp.symbols("q")
hodge_classes = (
    q**3 - 2 * q**2 + q,
    q**3 + q - 1,
    q**3 - 2 * q - 2,
)
assert all(sp.expand(value - q**3) != 0 for value in hodge_classes)
for suspension_dimension in range(8):
    assert all(
        sp.expand(q**suspension_dimension * value - q ** (
            suspension_dimension + 3
        ))
        != 0
        for value in hodge_classes
    )
print("PASS: no linear slice is affine or stably affine")


def canonical_vector(
    vector: tuple[int, int, int, int],
) -> tuple[int, int, int, int] | None:
    common_divisor = 0
    for entry in vector:
        common_divisor = gcd(common_divisor, abs(entry))
    if not common_divisor:
        return None
    primitive = tuple(entry // common_divisor for entry in vector)
    first_nonzero = next(entry for entry in primitive if entry)
    if first_nonzero < 0:
        primitive = tuple(-entry for entry in primitive)
    return primitive


# Bounded regression which led to the orbit classification.
height = 5
coefficient_vectors = {
    canonical_vector(vector)
    for vector in product(range(-height, height + 1), repeat=4)
    if any(vector)
}
assert None not in coefficient_vectors
assert len(coefficient_vectors) == 6928

affine_count_candidates = []
for vector in coefficient_vectors:
    assert vector is not None
    counts = [
        source_count(prime, vector)
        for prime in (5, 7, 11)
    ]
    if counts == [5**3, 7**3, 11**3]:
        affine_count_candidates.append(vector)
assert not affine_count_candidates
print("PASS: no height-five slice has affine counts over F_5,F_7,F_11")
print("PASS linear hyperplane Cox classification")

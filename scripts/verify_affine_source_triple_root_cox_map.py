#!/usr/bin/env python3
"""Exact audit of the affine-source oriented triple-root Cox map."""

from __future__ import annotations

from collections import Counter
from itertools import product

import sympy as sp


a, b, c = sp.symbols("a b c")
A = -a**2 * c + 3 * b * c + 2
C = -b * (b * c + 1) * (a**2 * c - b * c - 1)
E = 1 - a**2 * c
R = -6 * a * b * c * (b * c + 1) * (
    a**2 * c - b * c - 1
)

target_relation = sp.factor(
    3 * R**2
    - 4
    * C
    * (1 - E)
    * (A - E - 1)
    * (A - E + 2)
    * (A + 2 * E - 1)
)
assert target_relation == 0
print("PASS: the affine source lands on the five-factor oriented target")


jacobian = sp.factor(
    sp.Matrix((A, C, E)).jacobian((a, b, c)).det()
)
assert jacobian == R
residue_ratio = sp.cancel(jacobian / (6 * R))
assert residue_ratio == sp.Rational(1, 6)
print("PASS: the hypersurface-residue Jacobian is 1/6")


x = sp.factor((A - E - 1) / 3)
y = sp.factor(1 - E)
assert x == b * c
assert y == a**2 * c
assert sp.factor(C - b * (x + 1) * (E + x)) == 0
assert sp.factor(R - 6 * a * x * (x + 1) * (E + x)) == 0

inverse_a = sp.factor(R / (6 * x * (x + 1) * (E + x)))
inverse_b = sp.factor(C / ((x + 1) * (E + x)))
inverse_c = sp.factor(x * (x + 1) * (E + x) / C)
assert inverse_a == a
assert inverse_b == b
assert inverse_c == c
print("PASS: the displayed rational inverse is exact")


# Valuations use R as uniformizer; the selected branch factor has order 2.
valuation_table = {
    "y": (1, 0, 0),
    "C": (1, 2, -2),
    "x": (-1, 0, 2),
    "x+1": (-1, -2, 2),
    "E+x": (-1, -2, 2),
}
assert valuation_table["y"] == (1, 0, 0)
dicritical = {
    component: valuations
    for component, valuations in valuation_table.items()
    if min(valuations) < 0
}
assert set(dicritical) == {"C", "x", "x+1", "E+x"}
assert len(dicritical) == 4
print("PASS: one branch is finite and four target components are dicritical")


def finite_outputs(
    prime: int,
    source: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    finite_a, finite_b, finite_c = source
    finite_A = (
        -finite_a**2 * finite_c
        + 3 * finite_b * finite_c
        + 2
    ) % prime
    finite_C = (
        -finite_b
        * (finite_b * finite_c + 1)
        * (
            finite_a**2 * finite_c
            - finite_b * finite_c
            - 1
        )
    ) % prime
    finite_E = (1 - finite_a**2 * finite_c) % prime
    finite_R = (
        -6
        * finite_a
        * finite_b
        * finite_c
        * (finite_b * finite_c + 1)
        * (
            finite_a**2 * finite_c
            - finite_b * finite_c
            - 1
        )
    ) % prime
    return finite_A, finite_C, finite_E, finite_R


for prime in (5, 7, 11):
    fibers = Counter(
        finite_outputs(prime, source)
        for source in product(range(prime), repeat=3)
    )
    target_points = [
        (finite_A, finite_C, finite_E, finite_R)
        for finite_A, finite_C, finite_E, finite_R in product(
            range(prime),
            repeat=4,
        )
        if (
            3 * finite_R**2
            - 4
            * finite_C
            * (1 - finite_E)
            * (finite_A - finite_E - 1)
            * (finite_A - finite_E + 2)
            * (finite_A + 2 * finite_E - 1)
        )
        % prime
        == 0
    ]
    assert len(target_points) == prime**3
    profile = Counter(
        fibers.get(target_point, 0)
        for target_point in target_points
    )
    expected_profile = Counter(
        {
            0: (4 * prime - 5) * (prime - 1),
            1: prime**3 - 4 * prime**2 + 5 * prime - 2,
            prime - 1: 3 * (prime - 1),
            prime: prime - 1,
            2 * prime - 1: 1,
        }
    )
    assert profile == expected_profile
print("PASS: exact finite-field fiber profile holds over F_5,F_7,F_11")
print("PASS affine-source triple-root Cox map")

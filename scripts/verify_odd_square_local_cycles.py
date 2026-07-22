#!/usr/bin/env python3
"""Exact arithmetic checks for ODD_SQUARE_LOCAL_CYCLES.md."""
from __future__ import annotations

from fractions import Fraction
import math

import sympy as sp


def valuation(integer: int, prime: int) -> int:
    value = 0
    while integer % prime == 0:
        integer //= prime
        value += 1
    return value


def lower_hull(values: list[int]) -> tuple[tuple[int, int], ...]:
    hull: list[tuple[int, int]] = []
    for x_value, y_value in enumerate(values):
        point = (x_value, y_value)
        while len(hull) >= 2:
            x_0, y_0 = hull[-2]
            x_1, y_1 = hull[-1]
            old_slope = Fraction(y_1 - y_0, x_1 - x_0)
            new_slope = Fraction(y_value - y_1, x_value - x_1)
            if old_slope < new_slope:
                break
            hull.pop()
        hull.append(point)
    return tuple(hull)


def least_prime_divisor(integer: int) -> int:
    return min(sp.factorint(integer))


def q_cycle_criterion(a: int, prime: int) -> bool:
    degree = 2 * a * a - 1
    outer_degree = degree + 2
    assert outer_degree % prime == 0
    assert valuation(outer_degree, prime) == 1
    cycle = prime - 2
    return (
        outer_degree // prime > 1
        and math.gcd(degree, cycle) == 1
        and cycle > degree // least_prime_divisor(degree)
    )


X = sp.symbols("x")

# The Ishida/trinomial identity and the square discriminant formula.
for a in range(3, 24, 2):
    degree = 2 * a * a - 1
    outer_degree = degree + 2
    f_polynomial = sum(
        (index + 1) * X ** (degree - index) for index in range(degree + 1)
    )
    assert sp.expand((X - 1) ** 2 * f_polynomial) == (
        X**outer_degree - outer_degree * X + outer_degree - 1
    )

    discriminant = (
        2
        * outer_degree ** (degree - 1)
        * (degree + 1) ** (degree - 2)
    )
    discriminant_root = math.isqrt(discriminant)
    assert discriminant_root * discriminant_root == discriminant

# The q-adic edge (0,1)--(q-2,0) and the exact block criterion for the
# seven witnesses covered by Corollary 3.
Q_WITNESSES = {
    5: 17,
    7: 11,
    13: 113,
    15: 41,
    17: 193,
    19: 241,
    23: 353,
}

for a, prime in Q_WITNESSES.items():
    degree = 2 * a * a - 1
    outer_degree = degree + 2
    assert sp.isprime(prime)
    assert valuation(outer_degree, prime) == 1
    assert all(
        valuation(math.comb(outer_degree, index), prime) == 1
        for index in range(2, prime)
    )
    assert valuation(math.comb(outer_degree, prime), prime) == 0
    local_values = [
        valuation(math.comb(outer_degree, index + 2), prime)
        for index in range(prime - 1)
    ]
    assert lower_hull(local_values) == ((0, 1), (prime - 2, 0))
    assert q_cycle_criterion(a, prime)

    cofactor = outer_degree // prime
    assert math.gcd(degree, prime - 2) == math.gcd(cofactor - 1, prime - 2)
    print(
        f"PASS m={degree}: q={prime} isolates a {prime - 2}-cycle; "
        "the block criterion proves primitivity"
    )

# In the cofactor-three rows the block inequalities are automatic.
for a in (5, 13, 17, 19, 23, 31, 49):
    outer_degree = 2 * a * a + 1
    prime = outer_degree // 3
    assert outer_degree == 3 * prime and sp.isprime(prime)
    degree = outer_degree - 2
    assert math.gcd(degree, prime - 2) == 1
    assert degree // least_prime_divisor(degree) < prime - 2

# The p-adic polygons (14)--(15), retained as the next local-separation
# target rather than asserted to produce pure cycles.
for a in range(3, 34, 2):
    outer_degree = 2 * a * a + 1
    for prime in sp.factorint(a):
        height = valuation(outer_degree - 1, prime)
        last_power = prime**height

        for index in range(2, last_power + 1):
            expected = height - valuation(index * (index - 1), prime)
            assert valuation(math.comb(outer_degree, index), prime) == expected

        quotient_values = [
            valuation(math.comb(outer_degree, index + 2), prime)
            for index in range(last_power - 1)
        ]
        expected_quotient_hull = ((0, height),) + tuple(
            (prime**power - 2, height - power)
            for power in range(1, height + 1)
        )
        assert lower_hull(quotient_values) == expected_quotient_hull

        # The linear coefficient at a nontrivial residue root vanishes.
        # A sentinel above the hull represents its infinite valuation.
        other_cluster_values = [height, height + 1] + [
            valuation(math.comb(outer_degree, index), prime)
            for index in range(2, last_power + 1)
        ]
        expected_other_hull = ((0, height),) + tuple(
            (prime**power, height - power)
            for power in range(1, height + 1)
        )
        assert lower_hull(other_cluster_values) == expected_other_hull

print("PASS: odd-square local-cycle arithmetic and block criteria")


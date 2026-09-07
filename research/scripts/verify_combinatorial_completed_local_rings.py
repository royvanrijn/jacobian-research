#!/usr/bin/env python3
"""Audit the numerical local-ring tests for nested collisions."""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.collision_local_rings import (
    monomial_blowup,
    nested_node_local_algebra,
    node_local_algebra,
    phase_quotient_representatives,
)

def matrix_rank(rows: list[list[int]]) -> int:
    """Exact rational row rank for the small phase-character matrices."""

    matrix = [[Fraction(entry) for entry in row] for row in rows]
    rank = 0
    width = len(matrix[0]) if matrix else 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [entry / pivot_value for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


def matrix_rank_mod(rows: list[list[int]], prime: int) -> int:
    """Row rank over a prime field."""

    matrix = [[entry % prime for entry in row] for row in rows]
    rank = 0
    width = len(matrix[0]) if matrix else 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [
            entry * inverse % prime for entry in matrix[rank]
        ]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


def weak_compositions(total: int, parts: int):
    """Yield exponent tuples of fixed total degree."""

    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first, *tail)


def square_profile_conductor_order(number_of_square_nodes: int) -> int:
    """Conductor order for s_i^2=q after fixing one global sign.

    On the normalization branches, a degree-d monomial is t^d times a
    character of the phase quotient.  Since s_1 is diagonal t, the first
    degree in which those characters span every branch remains full in all
    higher degrees and is exactly the conductor order.
    """

    branches = tuple(
        (1, *signs)
        for signs in product((-1, 1), repeat=number_of_square_nodes - 1)
    )
    for degree in range(number_of_square_nodes + 1):
        evaluation_rows = []
        for exponents in weak_compositions(degree, number_of_square_nodes):
            evaluation_rows.append(
                [
                    math.prod(
                        sign**exponent
                        for sign, exponent in zip(branch, exponents)
                    )
                    for branch in branches
                ]
            )
        if matrix_rank(evaluation_rows) == len(branches):
            return degree
    raise AssertionError("phase characters did not span the normalization")


def is_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def prime_divisors(number: int) -> tuple[int, ...]:
    result = []
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            result.append(divisor)
            while number % divisor == 0:
                number //= divisor
        divisor += 1
    if number > 1:
        result.append(number)
    return tuple(result)


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = prime_divisors(order)
    return next(
        candidate
        for candidate in range(1, prime)
        if all(
            pow(candidate, order // factor, prime) != 1
            for factor in factors
        )
    )


def splitting_prime(common_index: int) -> int:
    """A small prime containing all common-index roots of unity."""

    multiplier = 1
    while True:
        candidate = multiplier * common_index + 1
        if is_prime(candidate):
            return candidate
        multiplier += 1


def weighted_exponents(
    total: int,
    weights: tuple[int, ...],
):
    """Yield nonnegative exponent vectors of the requested weighted degree."""

    if len(weights) == 1:
        if total % weights[0] == 0:
            yield (total // weights[0],)
        return
    first_weight = weights[0]
    for first in range(total // first_weight + 1):
        for tail in weighted_exponents(
            total - first * first_weight,
            weights[1:],
        ):
            yield (first, *tail)


def phase_character_conductor_order(
    profile: tuple[int, ...],
) -> int:
    """Recover the conductor independently from normalization characters."""

    package = node_local_algebra(profile)
    prime = splitting_prime(package.common_index)
    generator = primitive_root(prime)
    roots = tuple(
        pow(generator, (prime - 1) // index, prime)
        for index in profile
    )
    phases = package.phase_representatives

    def coefficient_rank(order: int) -> int:
        rows = []
        for exponents in weighted_exponents(order, package.source_orders):
            rows.append(
                [
                    math.prod(
                        pow(root, phase * exponent, prime)
                        for root, phase, exponent in zip(
                            roots,
                            branch,
                            exponents,
                        )
                    )
                    % prime
                    for branch in phases
                ]
            )
        return matrix_rank_mod(rows, prime)

    # Multiplication by q raises order by L without changing a phase
    # character.  A full block of L consecutive orders is therefore enough.
    for candidate in range(package.conductor_order + 1):
        if all(
            coefficient_rank(order) == package.branch_count
            for order in range(
                candidate,
                candidate + package.common_index,
            )
        ):
            return candidate
    raise AssertionError("no phase-character conductor was found")


degree_five = monomial_blowup(3, 2)
assert degree_five.primitive_ray == (2, 3)
assert degree_five.chart_indices == (3, 2)

# Ramifying x -> X^2 preserves the limiting cover decorations but changes
# one chart from a quotient singularity of index two to a regular chart.
base_changed = monomial_blowup(6, 2)
assert base_changed.primitive_ray == (1, 3)
assert base_changed.chart_indices == (3, 1)

pair_maxwell = node_local_algebra((1, 1, 2, 2))
assert pair_maxwell.common_index == 2
assert pair_maxwell.branch_count == 2
assert pair_maxwell.raw_fitting_order == 2
assert pair_maxwell.normalized_different_order == 1
assert pair_maxwell.conductor_order == 1

triple_maxwell = node_local_algebra((2, 2, 2))
assert triple_maxwell.common_index == 2
assert triple_maxwell.branch_count == 4
assert triple_maxwell.raw_fitting_order == 3
assert triple_maxwell.normalized_different_order == 1
assert triple_maxwell.conductor_order == 2

caustic = node_local_algebra((1, 1, 1, 3))
assert caustic.common_index == 3
assert caustic.branch_count == 1
assert caustic.raw_fitting_order == 2
assert caustic.normalized_different_order == 2
assert caustic.conductor_order == 0

mixed = nested_node_local_algebra(
    ((1, 1, 2, 2), (1, 1, 1, 3))
)
assert mixed.branch_count == 2
assert mixed.conductor_exponent_vector == (1, 0)
assert mixed.raw_fitting_exponent_vector == (2, 2)
assert mixed.normalized_different_exponent_vector == (1, 2)

# Conductor formula on a selected normalization branch.
assert node_local_algebra((3, 2)).conductor_order == 2  # cusp <2,3>
assert square_profile_conductor_order(2) == 1  # two square branches
assert square_profile_conductor_order(3) == 2  # four square branches

# Exhaustively check the branch count and closed conductor formula.  Identity
# indices are deliberately included: they must cancel from the formula.
profile_count = 0
for width in range(1, 5):
    for profile in product(range(1, 7), repeat=width):
        package = node_local_algebra(profile)
        assert len(phase_quotient_representatives(profile)) == (
            math.prod(profile) // math.lcm(*profile)
        )
        assert package.conductor_order == (
            (width - 1) * package.common_index
            - sum(package.source_orders)
            + 1
        )
        profile_count += 1

phase_conductor_count = 0
for width in range(1, 4):
    for profile in product(range(1, 6), repeat=width):
        assert phase_character_conductor_order(profile) == (
            node_local_algebra(profile).conductor_order
        )
        phase_conductor_count += 1

print("PASS degree-five contact-monoid counterexample: (3,2) versus (6,2)")
print("PASS Maxwell profiles: (1,1,2,2) and (2,2,2)")
print("PASS caustic profile: (1,1,1,3)")
print("PASS mixed Maxwell-caustic tensor factor")
print(
    "PASS closed conductor/Fitting/different identity:",
    profile_count,
    "profiles",
)
print(
    "PASS independent phase-character conductor:",
    phase_conductor_count,
    "profiles",
)
print("COMBINATORIAL_COMPLETED_LOCAL_RINGS_PASS")

#!/usr/bin/env python3
"""Exact finite replays for binary translation-tangent rigidity.

The proof is in
``extended-geometry/BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md``.
This dependency-free checker verifies bounded instances of its row kernel,
the nonprimitive rank jump, the cyclotomic linearized ranks away from the
chosen integer minor, and the abstract factorially weighted module
obstruction.  The bounded checks are not the proof of the unbounded theorem.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, factorial, gcd


def choose(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return comb(n, k)


def tangent_row(degree: int, slope: int, moment: int) -> list[int]:
    return [
        comb(degree, channel)
        * choose(
            degree * (moment - 1),
            slope * moment - channel,
        )
        for channel in range(degree + 1)
    ]


def rational_rank(rows: list[list[int]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    height = len(matrix)
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(pivot_row, height)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [
            value / pivot_value for value in matrix[pivot_row]
        ]
        for row in range(height):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                value - multiplier * pivot_entry
                for value, pivot_entry in zip(
                    matrix[row],
                    matrix[pivot_row],
                    strict=True,
                )
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            multiplier = work[row][column] / pivot_value
            for entry in range(column + 1, size):
                work[row][entry] -= multiplier * work[column][entry]
    assert answer.denominator == 1
    return answer.numerator


def modular_rank(rows: list[list[int]], prime: int) -> int:
    matrix = [
        [value % prime for value in row]
        for row in rows
    ]
    if not matrix:
        return 0
    height = len(matrix)
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(pivot_row, height)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [
            value * inverse % prime for value in matrix[pivot_row]
        ]
        for row in range(height):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                (value - multiplier * pivot_entry) % prime
                for value, pivot_entry in zip(
                    matrix[row],
                    matrix[pivot_row],
                    strict=True,
                )
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def primes_through(bound: int) -> list[int]:
    answer = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in answer if prime * prime <= candidate):
            answer.append(candidate)
    return answer


def verify_tangent_kernel(max_degree: int, prime_bound: int) -> tuple[int, int]:
    primitive_cases = 0
    nonprimitive_cases = 0
    for degree in range(2, max_degree + 1):
        for slope in range(1, degree):
            rows = [
                tangent_row(degree, slope, moment)
                for moment in range(1, 3 * degree + 1)
            ]
            tangent = [
                channel - slope for channel in range(degree + 1)
            ]
            assert all(
                sum(value * direction for value, direction in zip(
                    row,
                    tangent,
                    strict=True,
                )) == 0
                for row in rows
            )

            common_divisor = gcd(degree, slope)
            expected_rank = degree + 1 - common_divisor
            assert rational_rank(rows) == expected_rank
            if common_divisor != 1:
                nonprimitive_cases += 1
                continue

            primitive_cases += 1
            first_rows = rows[:degree]
            # Delete channel zero.  Its tangent coordinate is -slope, so a
            # nonzero minor here is compatible with the proved kernel.
            minor = [row[1:] for row in first_rows]
            delta = determinant(minor)
            assert delta
            for prime in primes_through(prime_bound):
                if delta % prime == 0:
                    continue
                assert modular_rank(first_rows, prime) == degree
    return primitive_cases, nonprimitive_cases


def constant_term_dilated(moment: int, dilation: int) -> int:
    if moment % 2:
        return 0
    # CT((z^d+z^-d)^N) is independent of the nonzero dilation d.
    assert dilation > 0
    return comb(moment, moment // 2)


def shifted_coefficient(moment: int, dilation: int) -> int:
    # CT(z^-1*(z^d+z^-d)^N) is the coefficient of z^1.
    numerator = moment * dilation - 1
    denominator = 2 * dilation
    if numerator % denominator:
        return 0
    selected_positive = numerator // denominator
    if not 0 <= selected_positive <= moment:
        return 0
    return comb(moment, selected_positive)


def verify_module_obstruction(depth: int) -> None:
    for moment in range(1, depth + 1):
        radial_factorial = factorial(moment) ** 2
        first = radial_factorial * constant_term_dilated(moment, 1)
        second = -radial_factorial * constant_term_dilated(moment, 2)
        assert first + second == 0
        if moment % 2 == 0:
            assert first and second

    for moment in range(1, depth + 1, 2):
        difference = (
            shifted_coefficient(moment, 1)
            - shifted_coefficient(moment, 2)
        )
        assert difference == comb(moment, (moment - 1) // 2)
        assert difference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=12)
    parser.add_argument("--prime-bound", type=int, default=97)
    parser.add_argument("--module-depth", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.max_degree < 2:
        raise SystemExit("max-degree must be at least two")
    primitive, nonprimitive = verify_tangent_kernel(
        arguments.max_degree,
        arguments.prime_bound,
    )
    verify_module_obstruction(arguments.module_depth)
    print(
        "PASS primitive translated-binomial tangent kernels: "
        f"{primitive} cases through degree {arguments.max_degree}"
    )
    print(
        "PASS nonprimitive power/subsequence rank jumps: "
        f"{nonprimitive} cases"
    )
    print(
        "PASS finite-field rank reductions away from the displayed "
        f"minors through prime {arguments.prime_bound}"
    )
    print(
        "PASS factorially weighted free-module inheritance obstruction "
        f"through order {arguments.module_depth}"
    )
    print(
        "STATUS: exact regressions for the proved tangent/large-prime "
        "theorems; affine-carry promotion remains open"
    )


if __name__ == "__main__":
    main()

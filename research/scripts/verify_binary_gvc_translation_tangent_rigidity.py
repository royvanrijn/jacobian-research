#!/usr/bin/env python3
"""Exact finite replays for binary translation-tangent rigidity.

The proof is in
``extended-geometry/BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md``.
This dependency-free checker verifies bounded instances of its row kernel,
the nonprimitive rank jump, the cyclotomic linearized ranks away from the
chosen integer minor used at every prime-power digit, the universal blind
rectangular tangent module and a nonflat member, the exact quadratic Hessian
kernel on bounded rectangles, and the abstract factorially weighted module
obstruction.  The bounded checks are not the proofs of the unbounded results.
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


def verify_blind_rectangle_module(
    max_degree: int,
    depth: int,
) -> tuple[int, int]:
    cases = 0
    blind_generators = 0
    for x_degree in range(2, max_degree + 1):
        for y_degree in range(2, max_degree + 1):
            for x_slope in range(1, x_degree):
                for y_slope in range(1, y_degree):
                    cases += 1
                    blind_generators += x_degree + y_degree + 2
                    generators = []
                    for y_basis in range(y_degree + 1):
                        generators.append([
                            (x_channel - x_slope)
                            if y_channel == y_basis
                            else 0
                            for x_channel in range(x_degree + 1)
                            for y_channel in range(y_degree + 1)
                        ])
                    for x_basis in range(x_degree + 1):
                        generators.append([
                            (y_channel - y_slope)
                            if x_channel == x_basis
                            else 0
                            for x_channel in range(x_degree + 1)
                            for y_channel in range(y_degree + 1)
                        ])
                    assert rational_rank(generators) == (
                        x_degree + y_degree + 1
                    )
                    for moment in range(1, depth + 1):
                        x_row = [
                            comb(x_degree, x_channel)
                            * choose(
                                x_degree * (moment - 1),
                                x_slope * moment - x_channel,
                            )
                            for x_channel in range(x_degree + 1)
                        ]
                        y_row = [
                            comb(y_degree, y_channel)
                            * choose(
                                y_degree * (moment - 1),
                                y_slope * moment - y_channel,
                            )
                            for y_channel in range(y_degree + 1)
                        ]
                        assert sum(
                            (x_channel - x_slope) * x_row[x_channel]
                            for x_channel in range(x_degree + 1)
                        ) == 0
                        assert sum(
                            (y_channel - y_slope) * y_row[y_channel]
                            for y_channel in range(y_degree + 1)
                        ) == 0

                        # These coordinate generators span
                        # (j-r)*v_k + u_j*(k-s).
                        for y_basis in range(y_degree + 1):
                            assert sum(
                                (x_channel - x_slope)
                                * x_row[x_channel]
                                * y_row[y_basis]
                                for x_channel in range(x_degree + 1)
                            ) == 0
                        for x_basis in range(x_degree + 1):
                            assert sum(
                                x_row[x_basis]
                                * (y_channel - y_slope)
                                * y_row[y_channel]
                                for y_channel in range(y_degree + 1)
                            ) == 0

                    # The mixed second difference of
                    # (j-r)*(k-s) is one, while every flat affine tangent
                    # has mixed second difference zero.
                    assert (
                        (1 - x_slope) * (1 - y_slope)
                        - (1 - x_slope) * (0 - y_slope)
                        - (0 - x_slope) * (1 - y_slope)
                        + (0 - x_slope) * (0 - y_slope)
                    ) == 1
    return cases, blind_generators


def rectangle_hessian(
    x_degree: int,
    y_degree: int,
    x_slope: int,
    y_slope: int,
    moment: int,
) -> list[list[int]]:
    channels = [
        (x_channel, y_channel)
        for x_channel in range(x_degree + 1)
        for y_channel in range(y_degree + 1)
    ]
    coefficients = [
        comb(x_degree, x_channel) * comb(y_degree, y_channel)
        for x_channel, y_channel in channels
    ]
    size = len(channels)
    matrix = [[0] * size for _ in range(size)]
    for left, (x_channel, y_channel) in enumerate(channels):
        matrix[left][left] += (
            moment
            * coefficients[left]
            * choose(
                x_degree * (moment - 1),
                x_slope * moment - x_channel,
            )
            * choose(
                y_degree * (moment - 1),
                y_slope * moment - y_channel,
            )
        )
    if moment < 2:
        return matrix
    for left, (x_left, y_left) in enumerate(channels):
        for right, (x_right, y_right) in enumerate(channels):
            matrix[left][right] += (
                moment
                * (moment - 1)
                * coefficients[left]
                * coefficients[right]
                * choose(
                    x_degree * (moment - 2),
                    x_slope * moment - x_left - x_right,
                )
                * choose(
                    y_degree * (moment - 2),
                    y_slope * moment - y_left - y_right,
                )
            )
    return matrix


def verify_quadratic_rectangle_rigidity(
    max_degree: int,
    depth: int,
) -> int:
    cases = 0
    for x_degree in range(2, max_degree + 1):
        for y_degree in range(2, max_degree + 1):
            size = (x_degree + 1) * (y_degree + 1)
            for x_slope in range(1, x_degree):
                for y_slope in range(1, y_degree):
                    cases += 1
                    gram = [[0] * size for _ in range(size)]
                    x_flat = []
                    y_flat = []
                    for x_channel in range(x_degree + 1):
                        for y_channel in range(y_degree + 1):
                            x_flat.append(x_channel - x_slope)
                            y_flat.append(y_channel - y_slope)
                    for moment in range(1, depth + 1):
                        hessian = rectangle_hessian(
                            x_degree,
                            y_degree,
                            x_slope,
                            y_slope,
                            moment,
                        )
                        assert all(
                            sum(
                                row[column] * x_flat[column]
                                for column in range(size)
                            ) == 0
                            for row in hessian
                        )
                        assert all(
                            sum(
                                row[column] * y_flat[column]
                                for column in range(size)
                            ) == 0
                            for row in hessian
                        )
                        gram = [
                            [
                                old + new
                                for old, new in zip(
                                    old_row,
                                    new_row,
                                    strict=True,
                                )
                            ]
                            for old_row, new_row in zip(
                                gram,
                                hessian,
                                strict=True,
                            )
                        ]
                    assert rational_rank(gram) == size - 2
    return cases


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
    parser.add_argument("--rectangle-degree", type=int, default=6)
    parser.add_argument("--rectangle-depth", type=int, default=8)
    parser.add_argument("--quadratic-degree", type=int, default=5)
    parser.add_argument("--quadratic-depth", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.max_degree < 2:
        raise SystemExit("max-degree must be at least two")
    if arguments.rectangle_degree < 2 or arguments.quadratic_degree < 2:
        raise SystemExit("rectangle degrees must be at least two")
    if arguments.quadratic_depth < arguments.quadratic_degree:
        raise SystemExit(
            "quadratic-depth must be at least quadratic-degree"
        )
    primitive, nonprimitive = verify_tangent_kernel(
        arguments.max_degree,
        arguments.prime_bound,
    )
    rectangle_cases, blind_generators = verify_blind_rectangle_module(
        arguments.rectangle_degree,
        arguments.rectangle_depth,
    )
    quadratic_cases = verify_quadratic_rectangle_rigidity(
        arguments.quadratic_degree,
        arguments.quadratic_depth,
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
        "PASS universal blind rectangular tangent module: "
        f"{blind_generators} spanning generators on {rectangle_cases} "
        "slope cases through bidegree "
        f"({arguments.rectangle_degree},{arguments.rectangle_degree}), "
        "including a nonflat bilinear member"
    )
    print(
        "PASS quadratic rectangular rigidity: common Hessian kernel is "
        f"the two flat torus directions in {quadratic_cases} slope cases "
        f"through bidegree "
        f"({arguments.quadratic_degree},{arguments.quadratic_degree})"
    )
    print(
        "PASS factorially weighted free-module inheritance obstruction "
        f"through order {arguments.module_depth}"
    )
    print(
        "STATUS: exact regressions for the proved tangent, quadratic, and "
        "large-prime-power theorems; common-quotient inheritance remains "
        "unproved in the parked route and is bypassed by Hall-envelope "
        "separation"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact infrastructure for the bidegree-(4,4) SIC rank frontier.

This checker does not claim to close ranks two through four.  It verifies
the determinantal parametrization, pure and mixed coefficient formulas,
the fixed-flag nilpotence screen, and the known rank-five endpoint.
"""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_frontier.json"
)
Matrix = list[list[Fraction]]
Polynomial = dict[tuple[int, int], Fraction]


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Fraction(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix, strict=True)]


def determinant(matrix: Matrix) -> Fraction:
    if not matrix:
        return Fraction(1)
    total = Fraction(0)
    for order in permutations(range(len(matrix))):
        inversions = sum(
            order[i] > order[j]
            for i in range(len(order))
            for j in range(i + 1, len(order))
        )
        term = Fraction(-1 if inversions % 2 else 1)
        for i, j in enumerate(order):
            term *= matrix[i][j]
        total += term
    return total


def rank(matrix: Matrix) -> int:
    for size in range(min(len(matrix), len(matrix[0])), 0, -1):
        for rows in combinations(range(len(matrix)), size):
            for columns in combinations(range(len(matrix[0])), size):
                minor = [[matrix[i][j] for j in columns] for i in rows]
                if determinant(minor):
                    return size
    return 0


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (i, j), left_coefficient in left.items():
        for (a, b), right_coefficient in right.items():
            exponent = (i + a, j + b)
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {exponent: value for exponent, value in result.items() if value}


def power(polynomial: Polynomial, order: int) -> Polynomial:
    result = {(0, 0): Fraction(1)}
    for _ in range(order):
        result = multiply(result, polynomial)
    return result


def coefficient_polynomial(matrix: Matrix) -> Polynomial:
    return {
        (i, j): matrix[i][j]
        for i in range(5)
        for j in range(5)
        if matrix[i][j]
    }


def diagonal_contraction(polynomial: Polynomial, degree: int) -> Fraction:
    return sum(
        (
            Fraction(factorial(i) * factorial(degree - i))
            * polynomial.get((i, i), Fraction(0))
            for i in range(degree + 1)
        ),
        Fraction(0),
    )


def pure_moment(matrix: Matrix, order: int) -> Fraction:
    return diagonal_contraction(
        power(coefficient_polynomial(matrix), order),
        4 * order,
    )


def mixed_moment(
    matrix: Matrix,
    multiplier: tuple[int, int],
    order: int,
) -> Fraction:
    """Contract xi1^a xi2^(1-a) z1^b z2^(1-b) F^m."""
    return diagonal_contraction(
        multiply(
            {multiplier: Fraction(1)},
            power(coefficient_polynomial(matrix), order),
        ),
        4 * order + 1,
    )


def mixed_shift_formula(
    matrix: Matrix,
    multiplier: tuple[int, int],
    order: int,
) -> Fraction:
    polynomial = power(coefficient_polynomial(matrix), order)
    a, b = multiplier
    degree = 4 * order + 1
    total = Fraction(0)
    for source_row in range(4 * order + 1):
        source_column = source_row + a - b
        if 0 <= source_column <= 4 * order:
            final_index = source_row + a
            total += (
                factorial(final_index)
                * factorial(degree - final_index)
                * polynomial.get((source_row, source_column), Fraction(0))
            )
    return total


def factor_matrix(rank_value: int) -> Matrix:
    """An exact rank-r point presented globally as U*V^T."""
    u = [
        [Fraction((row + 1) ** column) for column in range(rank_value)]
        for row in range(5)
    ]
    v = [
        [Fraction((row + 2) ** column) for column in range(rank_value)]
        for row in range(5)
    ]
    result = matrix_product(u, transpose(v))
    assert rank(result) == rank_value
    return result


def trace(matrix: Matrix) -> Fraction:
    return sum(
        (matrix[index][index] for index in range(len(matrix))),
        Fraction(0),
    )


def odd_double_factorial(order: int) -> int:
    result = 1
    for value in range(1, order + 1, 2):
        result *= value
    return result


def main() -> None:
    dimensions = {str(value): value * (10 - value) for value in range(1, 5)}
    assert dimensions == {"1": 9, "2": 16, "3": 21, "4": 24}

    formula_orders = 4
    for rank_value in range(1, 5):
        matrix = factor_matrix(rank_value)
        for order in range(1, formula_orders + 1):
            polynomial = power(coefficient_polynomial(matrix), order)
            assert pure_moment(matrix, order) == diagonal_contraction(
                polynomial,
                4 * order,
            )
            for multiplier in ((0, 0), (0, 1), (1, 0), (1, 1)):
                assert mixed_moment(
                    matrix,
                    multiplier,
                    order,
                ) == mixed_shift_formula(matrix, multiplier, order)

    one_sided: Matrix = [
        [Fraction(int(row > column)) for column in range(5)]
        for row in range(5)
    ]
    weights: Matrix = [
        [
            Fraction(factorial(i) * factorial(4 - i))
            if row == i
            else Fraction(0)
            for i in range(5)
        ]
        for row in range(5)
    ]
    endomorphism = matrix_product(one_sided, weights)
    endomorphism_power: Matrix = [
        [Fraction(int(row == column)) for column in range(5)]
        for row in range(5)
    ]
    traces: list[int] = []
    for _ in range(5):
        endomorphism_power = matrix_product(endomorphism_power, endomorphism)
        traces.append(int(trace(endomorphism_power)))
    assert traces == [0, 0, 0, 0, 0]

    witness: Matrix = [
        [Fraction(-1), Fraction(2), 0, 0, 0],
        [Fraction(-3, 2), Fraction(2), Fraction(6), 0, 0],
        [Fraction(-1, 2), Fraction(3, 2), Fraction(6), Fraction(6), 0],
        [0, Fraction(1), Fraction(3, 2), Fraction(2), Fraction(2)],
        [0, 0, Fraction(-1, 2), Fraction(-3, 2), Fraction(-1)],
    ]
    assert determinant(witness) == 48
    assert rank(witness) == 5
    for order in range(1, 5):
        assert pure_moment(witness, order) == 0
        assert mixed_moment(witness, (1, 0), order) == Fraction(
            factorial(4 * order + 2) * factorial(order),
            odd_double_factorial(2 * order + 1),
        )

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-frontier-v1",
        "field": "characteristic zero",
        "principal_invariant": "rank of the 5-by-5 coefficient matrix C",
        "certified_frontier": {
            "lower_bound": 2,
            "upper_bound": 5,
            "exact_value": None,
            "status": "open interval",
            "lower_bound_reason": (
                "the split-symbol theorem followed by fixed dual "
                "differentiation excludes every rank-one SIC point"
            ),
            "upper_bound_reason": "the exact witness has det(C)=48",
        },
        "best_open_target": {
            "coefficient_rank": 2,
            "bidegree": [4, 4],
            "ansatz": "F=A_1(xi)P_1(z)+A_2(xi)P_2(z)",
            "status": "open",
        },
        "rank_strata": {
            "1": {"dimension": 9, "status": "excluded exactly"},
            "2": {"dimension": 16, "status": "open"},
            "3": {"dimension": 21, "status": "open"},
            "4": {"dimension": 24, "status": "open"},
            "5": {
                "dimension": 25,
                "status": "counterexample exists",
                "coefficient_matrix_determinant": 48,
            },
        },
        "formula_replay_orders": formula_orders,
        "one_sided_endomorphism_trace_powers_1_through_5": traces,
        "written_source": (
            "extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS bidegree-(4,4) SIC rank frontier: certified interval [2,5]")
    print("PASS rank charts and pure/mixed formulas through order 4")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

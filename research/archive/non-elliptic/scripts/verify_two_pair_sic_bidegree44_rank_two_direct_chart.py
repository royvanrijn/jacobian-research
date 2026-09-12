#!/usr/bin/env python3
"""Verify the direct GL2-quotient chart for rank-two bidegree-(4,4) SIC.

This dependency-free checker has two purposes.

First, it fixes the internal factor gauge before any moment calculation.
On the open where the first two rows of U are invertible, the factorization

    C = U V^T

has the unique representative U=[I_2;A], V^T=B.  Thus the chart has the
expected sixteen coordinates, not the redundant twenty factor entries.
The checker verifies the gauge normalization, the first-moment pivot, and
the pure and mixed relative-period formulas at exact rational points.

Second, it supplies an all-order exact-rank-two control family on a fixed
one-sided flag.  Its period has strictly negative Laurent u-valuation, so
all pure moments vanish and every bidegree-(e,e) mixed sequence vanishes
for m>e.  This is a safe fixed-flag family, not a semistable component or a
global rank-two classification.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_direct_chart.json"
)
DEGREE = 4
Matrix = list[list[Fraction]]
Polynomial = dict[tuple[int, int], Fraction]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix, strict=True)]


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(
                (
                    left[row][inner] * right[inner][column]
                    for inner in range(len(right))
                ),
                Fraction(0),
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def determinant_two(matrix: Matrix) -> Fraction:
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError("expected a two-by-two matrix")
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def inverse_two(matrix: Matrix) -> Matrix:
    determinant = Fraction(determinant_two(matrix))
    if not determinant:
        raise ZeroDivisionError("singular two-by-two pivot")
    return [
        [Fraction(matrix[1][1]) / determinant, -Fraction(matrix[0][1]) / determinant],
        [-Fraction(matrix[1][0]) / determinant, Fraction(matrix[0][0]) / determinant],
    ]


def multiply_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for (left_x, left_y), left_value in left.items():
        for (right_x, right_y), right_value in right.items():
            exponent = (left_x + right_x, left_y + right_y)
            answer[exponent] = (
                answer.get(exponent, Fraction(0))
                + left_value * right_value
            )
    return {exponent: value for exponent, value in answer.items() if value}


def polynomial_power(polynomial: Polynomial, order: int) -> Polynomial:
    answer: Polynomial = {(0, 0): Fraction(1)}
    for _ in range(order):
        answer = multiply_polynomials(answer, polynomial)
    return answer


def coefficient_polynomial(matrix: Matrix) -> Polynomial:
    return {
        (row, column): matrix[row][column]
        for row in range(DEGREE + 1)
        for column in range(DEGREE + 1)
        if matrix[row][column]
    }


def diagonal_contraction(polynomial: Polynomial, degree: int) -> Fraction:
    return sum(
        (
            Fraction(factorial(index) * factorial(degree - index))
            * polynomial.get((index, index), Fraction(0))
            for index in range(degree + 1)
        ),
        Fraction(0),
    )


def mixed_moment(
    matrix: Matrix,
    order: int,
    multiplier_degree: int = 0,
    dual_index: int = 0,
    coordinate_index: int = 0,
) -> Fraction:
    """Return E_2(M*F^order) for a monomial multiplier M.

    Here

      M=xi1^dual_index xi2^(e-dual_index)
        z1^coordinate_index z2^(e-coordinate_index).
    """

    if not (
        0 <= dual_index <= multiplier_degree
        and 0 <= coordinate_index <= multiplier_degree
    ):
        raise ValueError("multiplier indices must lie between zero and e")
    power = polynomial_power(coefficient_polynomial(matrix), order)
    multiplied = multiply_polynomials(
        {(dual_index, coordinate_index): Fraction(1)},
        power,
    )
    return diagonal_contraction(
        multiplied,
        DEGREE * order + multiplier_degree,
    )


def beta_transform(matrix: Matrix) -> Polynomial:
    """Return P(u,t)=Phi(1,u,t,(1-t)/u) as a Laurent polynomial."""

    answer: Polynomial = {}
    for row in range(DEGREE + 1):
        for column in range(DEGREE + 1):
            value = matrix[row][column]
            if not value:
                continue
            for tail in range(DEGREE - column + 1):
                exponent = (column - row, column + tail)
                answer[exponent] = answer.get(exponent, Fraction(0)) + (
                    value
                    * comb(DEGREE - column, tail)
                    * (-1) ** tail
                )
    return {exponent: value for exponent, value in answer.items() if value}


def multiplier_transform(
    multiplier_degree: int,
    dual_index: int,
    coordinate_index: int,
) -> Polynomial:
    answer: Polynomial = {}
    for tail in range(multiplier_degree - coordinate_index + 1):
        exponent = (
            coordinate_index - dual_index,
            coordinate_index + tail,
        )
        answer[exponent] = Fraction(
            comb(multiplier_degree - coordinate_index, tail)
            * (-1) ** tail
        )
    return answer


def relative_period(
    matrix: Matrix,
    order: int,
    multiplier_degree: int = 0,
    dual_index: int = 0,
    coordinate_index: int = 0,
) -> Fraction:
    integrand = multiply_polynomials(
        multiplier_transform(
            multiplier_degree,
            dual_index,
            coordinate_index,
        ),
        polynomial_power(beta_transform(matrix), order),
    )
    return sum(
        (
            value / Fraction(t_exponent + 1)
            for (u_exponent, t_exponent), value in integrand.items()
            if u_exponent == 0
        ),
        Fraction(0),
    )


def normalize_factorization(
    u_matrix: Matrix,
    v_matrix: Matrix,
    pivot_rows: tuple[int, int] = (0, 1),
) -> tuple[Matrix, Matrix]:
    pivot = [list(u_matrix[row]) for row in pivot_rows]
    gauge = inverse_two(pivot)
    normalized_u = matrix_product(u_matrix, gauge)
    normalized_v = matrix_product(v_matrix, transpose(inverse_two(gauge)))
    return normalized_u, transpose(normalized_v)


def direct_chart_matrix(a_matrix: Matrix, b_matrix: Matrix) -> Matrix:
    if len(a_matrix) != 3 or any(len(row) != 2 for row in a_matrix):
        raise ValueError("A must be three-by-two")
    if len(b_matrix) != 2 or any(len(row) != 5 for row in b_matrix):
        raise ValueError("B must be two-by-five")
    u_matrix: Matrix = [
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
        *a_matrix,
    ]
    return matrix_product(u_matrix, b_matrix)


def first_moment_chart_formula(a_matrix: Matrix, b_matrix: Matrix) -> Fraction:
    return (
        24 * b_matrix[0][0]
        + 6 * b_matrix[1][1]
        + 4
        * (
            a_matrix[0][0] * b_matrix[0][2]
            + a_matrix[0][1] * b_matrix[1][2]
        )
        + 6
        * (
            a_matrix[1][0] * b_matrix[0][3]
            + a_matrix[1][1] * b_matrix[1][3]
        )
        + 24
        * (
            a_matrix[2][0] * b_matrix[0][4]
            + a_matrix[2][1] * b_matrix[1][4]
        )
    )


def fixed_flag_matrix(parameters: tuple[int, ...]) -> Matrix:
    if len(parameters) != 7:
        raise ValueError("expected seven fixed-flag parameters")
    a40, a41, b20, b21, b30, b31, b32 = map(Fraction, parameters)
    u_matrix: Matrix = [
        [0, 0],
        [0, 0],
        [1, 0],
        [0, 1],
        [a40, a41],
    ]
    b_matrix: Matrix = [
        [b20, b21, 0, 0, 0],
        [b30, b31, b32, 0, 0],
    ]
    if not b21 * b32:
        raise ValueError("the localized rank-two minor b21*b32 must be nonzero")
    return matrix_product(u_matrix, b_matrix)


def matrix_as_json(matrix: Matrix) -> list[list[int | str]]:
    return [
        [
            int(value)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
            for value in row
        ]
        for row in matrix
    ]


def main() -> None:
    # An exact factor point with an invertible top U block.
    u_matrix: Matrix = [
        [2, 1],
        [1, 1],
        [3, -1],
        [2, 4],
        [-2, 5],
    ]
    v_matrix: Matrix = [
        [1, 2],
        [3, -1],
        [2, 5],
        [-2, 1],
        [4, 3],
    ]
    original = matrix_product(u_matrix, transpose(v_matrix))
    normalized_u, b_matrix = normalize_factorization(u_matrix, v_matrix)
    assert normalized_u[:2] == [[1, 0], [0, 1]]
    assert matrix_product(normalized_u, b_matrix) == original
    assert b_matrix == [original[0], original[1]]

    # Internal gauge changes disappear after normalization.
    internal_gauge: Matrix = [[1, 2], [3, 5]]
    gauged_u = matrix_product(u_matrix, internal_gauge)
    gauged_v = matrix_product(
        v_matrix,
        transpose(inverse_two(internal_gauge)),
    )
    second_u, second_b = normalize_factorization(gauged_u, gauged_v)
    assert second_u == normalized_u
    assert second_b == b_matrix

    a_matrix = normalized_u[2:]
    chart_matrix = direct_chart_matrix(a_matrix, b_matrix)
    assert chart_matrix == original
    chart_minor = determinant_two(
        [[b_matrix[row][column] for column in (0, 1)] for row in (0, 1)]
    )
    assert chart_minor
    assert first_moment_chart_formula(a_matrix, b_matrix) == mixed_moment(
        chart_matrix,
        1,
    )

    # Verify the relative pure and mixed period identities exactly.
    period_checks = 0
    for order in range(5):
        for multiplier_degree in range(3):
            for dual_index in range(multiplier_degree + 1):
                for coordinate_index in range(multiplier_degree + 1):
                    raw = mixed_moment(
                        chart_matrix,
                        order,
                        multiplier_degree,
                        dual_index,
                        coordinate_index,
                    )
                    normalized = raw / factorial(
                        DEGREE * order + multiplier_degree + 1
                    )
                    assert normalized == relative_period(
                        chart_matrix,
                        order,
                        multiplier_degree,
                        dual_index,
                        coordinate_index,
                    )
                    period_checks += 1

    # A seven-parameter exact-rank-two family on the fixed one-sided flag.
    flag_parameters = (2, -1, 1, 2, 3, -1, 1)
    flag_matrix = fixed_flag_matrix(flag_parameters)
    assert determinant_two(
        [[flag_matrix[row][column] for column in (1, 2)] for row in (2, 3)]
    )
    flag_support = sorted(coefficient_polynomial(flag_matrix))
    assert flag_support
    assert all(row > column for row, column in flag_support)
    flag_beta_support = sorted(beta_transform(flag_matrix))
    assert flag_beta_support
    assert max(u_exponent for u_exponent, _ in flag_beta_support) <= -1

    pure_values = [mixed_moment(flag_matrix, order) for order in range(1, 9)]
    assert pure_values == [0] * 8

    mixed_sequences: dict[str, list[int]] = {}
    nonzero_low_mixed = []
    for multiplier_degree in (1, 2):
        for dual_index in range(multiplier_degree + 1):
            for coordinate_index in range(multiplier_degree + 1):
                values = [
                    mixed_moment(
                        flag_matrix,
                        order,
                        multiplier_degree,
                        dual_index,
                        coordinate_index,
                    )
                    for order in range(1, multiplier_degree + 4)
                ]
                key = f"e{multiplier_degree}_a{dual_index}_b{coordinate_index}"
                mixed_sequences[key] = [int(value) for value in values]
                assert all(
                    not value
                    for order, value in enumerate(values, start=1)
                    if order > multiplier_degree
                )
                if any(values):
                    nonzero_low_mixed.append(key)
    assert nonzero_low_mixed

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-direct-chart-v1",
        "field": "characteristic zero",
        "direct_chart": {
            "pivot_rows": [0, 1],
            "normal_form": "U=[I_2;A], C=U*B",
            "coordinates": {
                "A": 6,
                "B": 10,
                "total": 16,
            },
            "internal_gauge_removed": True,
            "rank_two_localization": "a nonzero 2-by-2 column minor of B",
            "sample_minor_columns_0_1": str(chart_minor),
            "mu1_formula": (
                "24*b00+6*b11+4*(a20*b02+a21*b12)"
                "+6*(a30*b03+a31*b13)+24*(a40*b04+a41*b14)"
            ),
            "mu1_constant_pivot": "24*b00",
        },
        "relative_period": {
            "P": "Phi_C(1,u,t,(1-t)/u)",
            "pure": "mu_m/(4m+1)!=CT_u integral_0^1 P^m dt",
            "mixed": (
                "E_2(M_eab*F^m)/(4m+e+1)!="
                "CT_u integral_0^1 u^(b-a)t^b(1-t)^(e-b)P^m dt"
            ),
            "exact_checks": period_checks,
            "orders": [0, 4],
            "multiplier_degrees": [0, 2],
        },
        "fixed_flag_rank_two_control_family": {
            "parameters": ["a40", "a41", "b20", "b21", "b30", "b31", "b32"],
            "localization": "b21*b32 != 0",
            "dimension": 7,
            "support_condition": "every nonzero c_ij has i>j",
            "relative_valuation": "max u-degree(P)<=-1",
            "pure_recurrence": "nu_(m+1)=0 for every m>=0",
            "initial_vanishing": "mu_1=0",
            "mixed_cutoff": "for every bidegree-(e,e) multiplier, eta_m=0 when m>e",
            "sample_matrix": matrix_as_json(flag_matrix),
            "sample_pure_moments_1_through_8": [int(value) for value in pure_values],
            "sample_low_degree_mixed_sequences": mixed_sequences,
            "nonzero_but_terminating_multipliers": nonzero_low_mixed,
            "persistent_mixed_sequence_found": False,
            "conclusion": "exact-rank-two fixed-flag family is SIC-safe",
        },
        "scope": (
            "direct quotient infrastructure and an exact safe control family; "
            "not an exact semistable component, a global rank-two exclusion, "
            "or a counterexample"
        ),
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS direct GL2 gauge quotient has sixteen coordinates")
    print("PASS exact pure and mixed relative-period identities")
    print("PASS fixed-flag exact-rank-two family has pure recurrence S*nu=0")
    print("PASS all bidegree-(e,e) mixed sequences vanish for m>e")


if __name__ == "__main__":
    main()

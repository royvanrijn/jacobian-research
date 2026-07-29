#!/usr/bin/env python3
"""Exact prerequisite audit for the bidegree-(4,4) rank-two problem.

This checker does not claim an all-order rank-two witness or obstruction.
It verifies the exact rank-two period/generating-function setup and records
that the only displayed exact rank-two point in the current frontier note is
a Jacobian transversality point, not a truncated moment survivor.
Recurrence derivation is parked until an explicit exact-rank-two moment
survivor is available.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_all_order_audit.json"
)
Polynomial = dict[tuple[int, int], Fraction]
Matrix = list[list[Fraction]]


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (i, j), left_value in left.items():
        for (a, b), right_value in right.items():
            exponent = (i + a, j + b)
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_value * right_value
            )
    return {exponent: value for exponent, value in result.items() if value}


def power(polynomial: Polynomial, order: int) -> Polynomial:
    result = {(0, 0): Fraction(1)}
    for _ in range(order):
        result = multiply(result, polynomial)
    return result


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


def matrix_rank(matrix: Matrix) -> int:
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                index
                for index in range(row, len(work))
                if work[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_value = work[row][column]
        work[row] = [value / pivot_value for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            scale = work[index][column]
            work[index] = [
                left - scale * right
                for left, right in zip(
                    work[index],
                    work[row],
                    strict=True,
                )
            ]
        row += 1
    return row


def substituted_polynomial(matrix: Matrix) -> Polynomial:
    """Phi(1,u,t,(1-t)/u), with exponents (u,t)."""
    result: Polynomial = {}
    for i in range(5):
        for j in range(5):
            for extra in range(5 - j):
                exponent = (j - i, j + extra)
                result[exponent] = (
                    result.get(exponent, Fraction(0))
                    + matrix[i][j]
                    * (-1) ** extra
                    * comb(4 - j, extra)
                )
    return {exponent: value for exponent, value in result.items() if value}


def factor_substituted_polynomial(u: Matrix, w: Matrix) -> Polynomial:
    result: Polynomial = {}
    for inner in range(2):
        dual = {
            (4 - i, 0): u[i][inner]
            for i in range(5)
            if u[i][inner]
        }
        coordinate: Polynomial = {}
        for j in range(5):
            for extra in range(5 - j):
                exponent = (j - 4, j + extra)
                coordinate[exponent] = (
                    coordinate.get(exponent, Fraction(0))
                    + w[inner][j]
                    * (-1) ** extra
                    * comb(4 - j, extra)
                )
        product = multiply(dual, coordinate)
        for exponent, value in product.items():
            result[exponent] = result.get(exponent, Fraction(0)) + value
    return {exponent: value for exponent, value in result.items() if value}


def constant_term_integral(polynomial: Polynomial) -> Fraction:
    return sum(
        (
            coefficient / Fraction(t_degree + 1)
            for (u_degree, t_degree), coefficient in polynomial.items()
            if u_degree == 0
        ),
        Fraction(0),
    )


def pure_moment(matrix: Matrix, order: int) -> Fraction:
    coefficient_polynomial = {
        (i, j): matrix[i][j]
        for i in range(5)
        for j in range(5)
        if matrix[i][j]
    }
    polynomial_power = power(coefficient_polynomial, order)
    degree = 4 * order
    return sum(
        (
            factorial(index)
            * factorial(degree - index)
            * polynomial_power.get((index, index), Fraction(0))
            for index in range(degree + 1)
        ),
        Fraction(0),
    )


def convex_hull(points: set[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(points)

    def cross(
        origin: tuple[int, int],
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> int:
        return (
            (left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0])
        )

    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def twice_polygon_area(vertices: list[tuple[int, int]]) -> int:
    return abs(
        sum(
            left[0] * right[1] - left[1] * right[0]
            for left, right in zip(
                vertices,
                vertices[1:] + vertices[:1],
                strict=True,
            )
        )
    )


def main() -> None:
    u: Matrix = [
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
        [Fraction(14), Fraction(17)],
        [Fraction(18), Fraction(4)],
        [Fraction(6), Fraction(13)],
    ]
    w: Matrix = [
        [
            Fraction(8),
            Fraction(10),
            Fraction(1),
            Fraction(8),
            Fraction(4),
        ],
        [
            Fraction(19),
            Fraction(1),
            Fraction(4),
            Fraction(6),
            Fraction(17),
        ],
    ]
    matrix = matrix_product(u, w)
    assert matrix_rank(matrix) == 2
    assert matrix == [
        [8, 10, 1, 8, 4],
        [19, 1, 4, 6, 17],
        [435, 157, 82, 214, 345],
        [220, 184, 34, 168, 140],
        [295, 73, 58, 126, 245],
    ]

    substituted = substituted_polynomial(matrix)
    assert substituted == factor_substituted_polynomial(u, w)
    checked_orders = 4
    moments: list[int] = []
    for order in range(1, checked_orders + 1):
        moment = pure_moment(matrix, order)
        period = constant_term_integral(power(substituted, order))
        assert moment == factorial(4 * order + 1) * period
        moments.append(int(moment))
    assert moments == [
        7414,
        3675739680,
        12167497410877440,
        148010006143680629760000,
    ]

    generic_support = {
        (j - i, j + extra)
        for i in range(5)
        for j in range(5)
        for extra in range(5 - j)
    }
    hull = convex_hull(generic_support)
    assert hull == [(-4, 0), (0, 0), (4, 4), (-4, 4)]
    twice_area = twice_polygon_area(hull)
    assert twice_area == 48

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-prerequisite-audit-v2",
        "field": "characteristic zero",
        "factor_chart": "C=U*W with U in Mat(5,2), W in Mat(2,5)",
        "displayed_rank_two_point": {
            "source_equation": (
                "extended-geometry/"
                "TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md (6.1)"
            ),
            "rank": matrix_rank(matrix),
            "pure_moments_1_through_4": moments,
            "role": "Jacobian transversality point, not a moment survivor",
        },
        "period_identity": {
            "P": "Phi_C(1,u,t,(1-t)/u)",
            "normalized_moment": (
                "mu_m/(4m+1)! = CT_u integral_0^1 P(u,t)^m dt"
            ),
            "ordinary_generating_function": (
                "sum_(m>=0) mu_m/(4m+1)! s^m "
                "= CT_u integral_0^1 1/(1-s*P(u,t)) dt"
            ),
            "checked_orders": checked_orders,
        },
        "generic_rank_two_laurent_support": {
            "newton_polygon_vertices": [list(vertex) for vertex in hull],
            "euclidean_area": Fraction(twice_area, 2).numerator,
            "normalized_volume": twice_area,
            "warning": (
                "the normalized volume is not by itself a recurrence-order "
                "or initial-moment bound at the singular expansion s=0"
            ),
        },
        "explicit_survivor_status": (
            "none recorded: Hilbert theory proves only an existential "
            "semistable point on the rank-at-most-two fiber"
        ),
        "rank_one_boundary_status": (
            "one squarefree Rabinowitsch membership remains open"
        ),
        "recurrence_status": (
            "parked until the rank-one boundary is closed, an explicit "
            "closed semistable point or component is extracted, and its "
            "coefficient rank is proved to be exactly two"
        ),
        "tail_evaluation_status": (
            "mu_14 is not evaluated without an explicit exact-rank-two "
            "specialization"
        ),
        "mandatory_gate_sequence": [
            "close the remaining rank-one membership",
            "extract an explicit closed semistable point or component",
            "prove exact coefficient rank two",
            "derive the specialized recurrence and evaluate mu_14",
        ],
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank-two factor and beta/constant-term period identities")
    print("PASS displayed exact rank-two chart point has mu_1=7414, not zero")
    print("PASS generic rank-two Newton polygon has normalized volume 48")
    print("PASS recurrence derivation is parked behind the exact-point gates")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

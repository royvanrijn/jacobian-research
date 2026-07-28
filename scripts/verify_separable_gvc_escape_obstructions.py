#!/usr/bin/env python3
"""Exact finite checks for the separable GVC escape obstructions.

The all-order theorems are proved in SPLIT_SYMBOL_GVC_THEOREM.md and
SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md.  This dependency-free script checks
the finite linear-algebra claims, translated polarization identity, and a
degree-raising factor-unit example.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


Matrix = list[list[Fraction]]
Polynomial = dict[tuple[int, int], Fraction]
SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "separable_gvc_escape_obstructions.json"
)


def rank(matrix: Matrix) -> int:
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_value = work[row][column]
        work[row] = [entry / pivot_value for entry in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            factor = work[index][column]
            work[index] = [
                work[index][j] - factor * work[row][j]
                for j in range(len(work[0]))
            ]
        row += 1
        if row == len(work):
            break
    return row


def determinant(matrix: Matrix) -> Fraction:
    work = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            index
            for index in range(column, len(work))
            if work[index][column]
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for index in range(column + 1, len(work)):
            factor = work[index][column] / pivot_value
            for j in range(column, len(work)):
                work[index][j] -= factor * work[column][j]
    return result


def outer(left: list[int], right: list[int]) -> Matrix:
    return [
        [Fraction(left_entry * right_entry) for right_entry in right]
        for left_entry in left
    ]


def coefficient(power: int, exponent: int) -> int:
    return int(power == exponent)


def polynomial_add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, value in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + value
    return {exponent: value for exponent, value in result.items() if value}


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (left_x, left_y), left_value in left.items():
        for (right_x, right_y), right_value in right.items():
            exponent = (left_x + right_x, left_y + right_y)
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_value * right_value
            )
    return {exponent: value for exponent, value in result.items() if value}


def polynomial_power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result: Polynomial = {(0, 0): Fraction(1)}
    base = polynomial
    value = exponent
    while value:
        if value % 2:
            result = polynomial_multiply(result, base)
        base = polynomial_multiply(base, base)
        value //= 2
    return result


def derivative(polynomial: Polynomial, x_order: int, y_order: int) -> Polynomial:
    result: Polynomial = {}
    for (x_power, y_power), value in polynomial.items():
        if x_power < x_order or y_power < y_order:
            continue
        multiplier = (
            factorial(x_power)
            // factorial(x_power - x_order)
            * factorial(y_power)
            // factorial(y_power - y_order)
        )
        result[(x_power - x_order, y_power - y_order)] = value * multiplier
    return result


def translated_diagonal_coefficient(
    polynomial: Polynomial,
    order: int,
) -> Polynomial:
    """Return [t1^order*t2^order] R(x+t1,y+t2)."""
    result: Polynomial = {}
    for (x_power, y_power), value in polynomial.items():
        if x_power < order or y_power < order:
            continue
        result[(x_power - order, y_power - order)] = (
            value * comb(x_power, order) * comb(y_power, order)
        )
    return result


def factor_unit_operator(polynomial: Polynomial) -> Polynomial:
    """Apply dx*dy*(1+dx+dy^2)."""
    inner = polynomial_add(
        polynomial,
        derivative(polynomial, 1, 0),
        derivative(polynomial, 0, 2),
    )
    return derivative(inner, 1, 1)


def iterate_operator(polynomial: Polynomial, order: int) -> Polynomial:
    result = polynomial
    for _ in range(order):
        result = factor_unit_operator(result)
    return result


def main() -> None:
    witness = [
        [Fraction(-1), Fraction(2), 0, 0, 0],
        [Fraction(-3, 2), Fraction(2), Fraction(6), 0, 0],
        [Fraction(-1, 2), Fraction(3, 2), Fraction(6), Fraction(6), 0],
        [0, Fraction(1), Fraction(3, 2), Fraction(2), Fraction(2)],
        [0, 0, Fraction(-1, 2), Fraction(-3, 2), Fraction(-1)],
    ]
    assert determinant(witness) == 48
    assert rank(witness) == 5

    rank_one_channels = [
        outer([1, 0, 2, -1, 3], [2, -1, 0, 4, 1]),
        outer([0, 1, 1, 2, -2], [1, 3, -1, 0, 2]),
        outer([2, -1, 0, 1, 1], [0, 2, 3, -1, 1]),
        outer([1, 1, -1, 0, 2], [3, 0, 1, 2, -1]),
    ]
    assert all(rank(channel) == 1 for channel in rank_one_channels)
    four_channel_sum = [
        [
            sum(channel[i][j] for channel in rank_one_channels)
            for j in range(5)
        ]
        for i in range(5)
    ]
    assert rank(four_channel_sum) <= 4

    # [u^j] is not multiplicative.  For j=4, u*u^3 is detected although
    # neither factor is; this is the first coefficient width capable of
    # carrying five rank-one summands.
    assert coefficient(1 + 3, 4) == 1
    assert coefficient(1, 4) * coefficient(3, 4) == 0
    assert 4 + 1 == 5

    # The translated coefficient identity retains all degrees of P.
    high_degree = {
        (3, 0): Fraction(1),
        (2, 1): Fraction(2),
        (0, 4): Fraction(-1),
        (0, 0): Fraction(3),
    }
    for order in range(1, 5):
        power = polynomial_power(high_degree, order)
        differentiated = derivative(power, order, order)
        translated = translated_diagonal_coefficient(power, order)
        assert differentiated == {
            exponent: value * factorial(order) ** 2
            for exponent, value in translated.items()
        }

    # Lambda=dx*dy*(1+dx+dy^2) is nonhomogeneous factor-unit, while
    # P=x^3 has degree greater than the order of the split factor dx*dy.
    p = {(3, 0): Fraction(1)}
    q = {(1, 0): Fraction(1), (0, 2): Fraction(1)}
    factor_unit_mixed_nonzero_orders: list[int] = []
    for order in range(1, 7):
        power = polynomial_power(p, order)
        assert iterate_operator(power, order) == {}
        mixed = iterate_operator(polynomial_multiply(q, power), order)
        if mixed:
            factor_unit_mixed_nonzero_orders.append(order)
    assert factor_unit_mixed_nonzero_orders == [1, 2]

    artifact = {
        "format": "separable-gvc-escape-obstructions-v2",
        "field": "characteristic zero",
        "two_pair_witness": {
            "bidegree": [4, 4],
            "coefficient_matrix_determinant": 48,
            "coefficient_matrix_rank": 5,
        },
        "separated_conversion": {
            "one_channel_rank_bound": 1,
            "minimum_additive_channels": 5,
            "one_variable_coefficient_order_lower_bound": 4,
            "coefficient_extraction_multiplicative": False,
        },
        "all_order_theorem": {
            "statement": (
                "For a binary constant-coefficient operator with lowest "
                "positive order r, deg(P)<=r satisfies GVC."
            ),
            "split_symbol_strengthening": (
                "A homogeneous split-symbol operator satisfies GVC for "
                "arbitrary polynomial degree."
            ),
            "factor_unit_extension": (
                "Lambda=Lambda_0*Gamma satisfies GVC for arbitrary P when "
                "Lambda_0 is homogeneous split and Gamma(0)!=0."
            ),
            "single_coefficient_converse": (
                "Fixed linear translation and one diagonal coefficient "
                "represent exactly products of powers of linear symbols."
            ),
            "associated_graded_obstruction": (
                "For arbitrary nonhomogeneous binary Lambda, every fixed "
                "number of leading mixed homogeneous layers vanishes "
                "eventually."
            ),
            "quadratic_heat_class": (
                "For Lambda with linear and quadratic pieces and deg(P)<=2, "
                "the first two pure equations imply GVC."
            ),
            "rank_one_quadratic_heat_class": (
                "If the quadratic symbol has rank at most one, the first "
                "two pure equations imply GVC for arbitrary P."
            ),
            "binary_drift_diffusion_class": (
                "Every binary operator with nonzero linear part and no "
                "piece above order two satisfies GVC after its first two "
                "pure equations, for arbitrary P."
            ),
            "separated_drift_class": (
                "Every Lambda=dx+h(dy), including formal h acting locally "
                "finitely on polynomials, satisfies GVC for arbitrary P "
                "after its first two pure equations."
            ),
            "formal_drift_straightening": (
                "Every binary Lambda with nonzero linear part factors "
                "formally as U(dx,dy)*(dx+q(dy)), with U a locally finite "
                "differential automorphism; hence it satisfies GVC for "
                "arbitrary P after its first two pure equations."
            ),
            "arbitrary_order_quadratic_polynomial_class": (
                "Every binary Lambda with nonzero linear part satisfies "
                "GVC for deg(P)<=2 after its first two pure equations."
            ),
            "arbitrary_order_cubic_polynomial_class": (
                "Every binary Lambda with nonzero linear part satisfies "
                "GVC for deg(P)<=3 after its first two pure equations; "
                "complete order-four and order-five terms do not cancel "
                "the decisive second-moment branches."
            ),
            "quartic_polynomial_regression": (
                "For deg(P)<=4, the complete operator 7-jet gives the "
                "successive square branches 2304*C^4*p4^2, "
                "15552*G^2*p4^2, and 39168*L^2*p4^2."
            ),
            "quadratic_leading_cubic_polynomial_class": (
                "If the lowest positive operator order is two and "
                "deg(P)<=3, the first three pure equations imply GVC. "
                "The double-line cancellation branch is killed by "
                "-4608*C^3*p_xy2^3; every survivor has a strict weighted "
                "degree separator."
            ),
            "proof": (
                "extended-geometry/"
                "SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md"
            ),
            "external_input": (
                "Duistermaat--van der Kallen constant-term theorem"
            ),
            "bounded_computation": False,
        },
        "finite_regressions": {
            "translated_polarization_orders": [1, 2, 3, 4],
            "translated_polynomial_degree": 4,
            "split_operator_order": 2,
            "factor_unit_operator": "dx*dy*(1+dx+dy^2)",
            "factor_unit_polynomial": "x^3",
            "factor_unit_mixed_nonzero_orders_through_6": (
                factor_unit_mixed_nonzero_orders
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS separable escape: witness coefficient rank is five")
    print("PASS separable escape: four rank-one channels have rank at most four")
    print("PASS separable escape: coefficient extraction is not multiplicative")
    print("PASS split symbol: translated polarization retains degree > order")
    print("PASS factor unit: degree-raising mixed values vanish after order two")
    print("SCOPE: all-order conclusions are the written Newton arguments")
    print(f"PASS separable escape: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

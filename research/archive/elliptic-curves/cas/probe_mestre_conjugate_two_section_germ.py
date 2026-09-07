#!/usr/bin/env python3
"""Local continuation of a conjugate-slope two-section Mestre germ.

The post-max-200 ``p=7`` scout finds the normalized root seed
``(0,7,127,128,225,233)`` with lines

    x_1=233/113-97*T/113,  x_2=233/113+97*T/113.

Its labelled incidence Jacobian has rank six.  The two tangent directions
preserve the conjugate form ``x_1=a-b*T, x_2=a+b*T``.  This checker solves the
first six recursive equations in a bivariate formal neighbourhood with
``(a,b)`` as parameters and audits the remaining residual.  It uses only the
recursive coefficient representation; no expanded residual is formed.

The result is local evidence for a two-dimensional germ, not a reconstructed
global component, a Mordell--Weil independence statement, or a height/Shioda
calculation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from probe_mestre_two_section_local_continuation import (
    EQUATION_NAMES,
    Field,
    FormalBivariate,
    FormalSeries,
    Jet,
    VARIABLES,
    residuals,
    residuals_from_jets,
    row_reduce,
    solve_square_over_q,
)
from screen_mestre_two_section_transverse_seeds import normalized_moduli


Q = Fraction
ROOTS = (0, 7, 127, 128, 225, 233)
SEED = normalized_moduli(ROOTS) + (
    Q(233, 113), -Q(97, 113), Q(233, 113), Q(97, 113)
)
PIVOT_ROWS = tuple(range(6))
PIVOT_COLUMNS = tuple(range(6))
FREE_COLUMNS = (6, 7)
GOOD_TANGENT_PRIMES = (17, 23, 29)
ROOT_COORDINATES = tuple(Q(value, 7) for value in (127, 128, 225, 233))
COMPARISON_ROOTS = (0, 21, 151, 169, 200, 239)
COMPARISON_SEED = normalized_moduli(COMPARISON_ROOTS) + (
    Q(239, 109), -Q(31, 109), Q(239, 109), Q(31, 109)
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def pivot_matrix() -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    values = residuals(SEED, Field())
    if any(value.value for value in values):
        raise AssertionError("the conjugate-slope seed no longer solves the recursive system")
    jacobian = [list(value.gradient) for value in values]
    rank, pivots = row_reduce(jacobian, Field())
    if rank != 6 or pivots != list(PIVOT_COLUMNS):
        raise AssertionError("the conjugate-slope tangent rank changed")
    matrix = [[jacobian[row][column] for column in PIVOT_COLUMNS] for row in PIVOT_ROWS]
    if row_reduce(matrix, Field())[0] != 6:
        raise AssertionError("the chosen six-by-six implicit minor became singular")
    return matrix, jacobian


def coordinates_from_roots(values: list[object]) -> tuple[object, ...]:
    """Convert four moving roots and two affine lines to recursive inputs."""

    roots = values[:4]
    zero = type(roots[0]).constant(roots[0].field, 0)
    c1 = zero
    for root in roots:
        c1 = c1 - root
    c2 = zero
    for first in range(4):
        for second in range(first + 1, 4):
            c2 = c2 + roots[first] * roots[second]
    c3 = zero
    for first in range(4):
        for second in range(first + 1, 4):
            for third in range(second + 1, 4):
                c3 = c3 - roots[first] * roots[second] * roots[third]
    c4 = roots[0] * roots[1] * roots[2] * roots[3]
    return c1, c2, c3, c4, *values[4:]


def root_slice_matrix() -> tuple[list[list[Fraction]], list[int]]:
    """Use r3,r6 as local free coordinates and solve for the other six."""

    field = Field()
    values = [
        Jet.variable(field, value, index)
        for index, value in enumerate((*ROOT_COORDINATES, *SEED[4:]))
    ]
    jacobian = [
        list(value.gradient) for value in residuals_from_jets(coordinates_from_roots(values))
    ]
    columns = [1, 2, 4, 5, 6, 7]
    matrix = [[jacobian[row][column] for column in columns] for row in PIVOT_ROWS]
    if row_reduce(matrix, field)[0] != 6:
        raise AssertionError("the root-motion implicit minor became singular")
    return matrix, columns


def root_motion_slice(order: int) -> dict[str, object]:
    """Lift the chart r6=233/7, r3=127/7+t in exact root coordinates."""

    matrix, columns = root_slice_matrix()
    field = Field(series_order=order)
    values = [
        [Q(0)] * (order + 1)
        for _ in range(8)
    ]
    for index, value in enumerate((*ROOT_COORDINATES, *SEED[4:])):
        values[index][0] = value
    values[0][1] = Q(1)
    remaining = []
    for degree in range(1, order + 1):
        current = [FormalSeries.seed(field, row) for row in values]
        residual_values = residuals_from_jets(coordinates_from_roots(current))
        correction = solve_square_over_q(
            matrix, [-residual_values[row].coefficients[degree] for row in PIVOT_ROWS]
        )
        for column, value in zip(columns, correction):
            values[column][degree] = value
        checked = residuals_from_jets(
            coordinates_from_roots([FormalSeries.seed(field, row) for row in values])
        )
        if any(checked[row].coefficients[degree] for row in PIVOT_ROWS):
            raise AssertionError("the root-motion implicit solve failed")
        remaining.append(checked[6].coefficients[degree])
    local_parameter = FormalSeries.seed(field, [Q(0), Q(1), *([Q(0)] * (order - 1))])
    expected_intercept = Q(233) / (Q(113) - 7 * local_parameter)
    first_intercept = FormalSeries.seed(field, values[4])
    second_intercept = FormalSeries.seed(field, values[6])
    first_slope = FormalSeries.seed(field, values[5])
    second_slope = FormalSeries.seed(field, values[7])
    if first_intercept != expected_intercept or second_intercept != expected_intercept:
        raise AssertionError("the root-motion affine intercept recognition failed")
    if first_slope != -second_slope:
        raise AssertionError("the root-motion slopes ceased to be conjugate")
    return {
        "order": order,
        "parameters": "r3=127/7+t, r6=233/7",
        "all_E0_2_coefficients_through_order_vanish": all(value == 0 for value in remaining),
        "recognized_common_intercept": "233/(113-7t)=233/(240-7r3)",
        "conjugate_slope_identity_through_order": order,
        "remaining_E0_2_coefficients": [rational_text(value) for value in remaining],
    }


def comparison_seed_audit(order: int = 3) -> dict[str, object]:
    """Check that a distant conjugate-slope seed has the same local germ."""

    values = residuals(COMPARISON_SEED, Field())
    if any(value.value for value in values):
        raise AssertionError("the comparison conjugate-slope seed left the residual locus")
    jacobian = [list(value.gradient) for value in values]
    rank, pivots = row_reduce(jacobian, Field())
    if rank != 6 or pivots != list(PIVOT_COLUMNS):
        raise AssertionError("the comparison seed tangent rank changed")
    matrix = [[jacobian[row][column] for column in PIVOT_COLUMNS] for row in PIVOT_ROWS]
    field = Field(series_order=order)
    coefficients = [{(0, 0): value} for value in COMPARISON_SEED]
    coefficients[6][(1, 0)] = Q(1)
    coefficients[7][(0, 1)] = Q(1)
    remaining = []
    for degree in range(1, order + 1):
        for first_degree in range(degree + 1):
            monomial = (first_degree, degree - first_degree)
            current = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            correction = solve_square_over_q(
                matrix, [-current[row].coefficient(monomial) for row in PIVOT_ROWS]
            )
            for column, value in zip(PIVOT_COLUMNS, correction):
                coefficients[column][monomial] = value
            checked = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            if any(checked[row].coefficient(monomial) for row in PIVOT_ROWS):
                raise AssertionError("the comparison bivariate solve failed")
            remaining.append(checked[6].coefficient(monomial))
    modular_ranks = {}
    for prime in (11, 17, 19, 23, 29):
        values = residuals(COMPARISON_SEED, Field(prime))
        prime_rank, _ = row_reduce([value.gradient for value in values], Field(prime))
        if any(value.value for value in values) or prime_rank != 6:
            raise AssertionError("the comparison small-prime tangent rank changed")
        modular_ranks[str(prime)] = prime_rank
    return {
        "integer_roots": list(COMPARISON_ROOTS),
        "seed_sections": [
            [rational_text(value) for value in COMPARISON_SEED[4:6]],
            [rational_text(value) for value in COMPARISON_SEED[6:8]],
        ],
        "exact_tangent_rank_over_Q": rank,
        "modular_tangent_ranks": modular_ranks,
        "bivariate_order": order,
        "all_E0_2_coefficients_through_total_order_vanish": all(
            value == 0 for value in remaining
        ),
    }


def tangent_directions() -> list[list[Fraction]]:
    matrix, jacobian = pivot_matrix()
    directions = []
    for free in FREE_COLUMNS:
        correction = solve_square_over_q(
            matrix, [-jacobian[row][free] for row in PIVOT_ROWS]
        )
        direction = [Q(0)] * len(VARIABLES)
        for column, value in zip(PIVOT_COLUMNS, correction):
            direction[column] = value
        direction[free] = Q(1)
        if any(
            sum(value * derivative for value, derivative in zip(direction, row))
            for row in jacobian
        ):
            raise AssertionError("a displayed germ direction left the tangent kernel")
        directions.append(direction)
    # The free coordinates are the second line's intercept/slope.  The first
    # line follows respectively with the same intercept motion and opposite
    # slope motion, exactly the conjugate-slope tangent plane.
    if directions[0][4] != 1 or directions[0][5] != 0:
        raise AssertionError("the intercept tangent ceased to be conjugate")
    if directions[1][4] != 0 or directions[1][5] != -1:
        raise AssertionError("the slope tangent ceased to be conjugate")
    return directions


def bivariate_germ(order: int) -> dict[str, object]:
    """Lift the six-row implicit system in Q[[a-a0,b-b0]]."""

    matrix, _ = pivot_matrix()
    field = Field(series_order=order)
    coefficients = [{(0, 0): value} for value in SEED]
    coefficients[6][(1, 0)] = Q(1)
    coefficients[7][(0, 1)] = Q(1)
    remaining = {}
    for degree in range(1, order + 1):
        for a_degree in range(degree + 1):
            monomial = (a_degree, degree - a_degree)
            values = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            correction = solve_square_over_q(
                matrix, [-values[row].coefficient(monomial) for row in PIVOT_ROWS]
            )
            for column, value in zip(PIVOT_COLUMNS, correction):
                coefficients[column][monomial] = value
            checked = residuals_from_jets(
                [FormalBivariate.seed(field, row) for row in coefficients]
            )
            if any(checked[row].coefficient(monomial) for row in PIVOT_ROWS):
                raise AssertionError("the six-row bivariate implicit solve failed")
            remaining[monomial] = checked[6].coefficient(monomial)
    return {
        "order": order,
        "parameters": "a=x02-233/113, b=x12-97/113",
        "all_E0_2_coefficients_through_total_order_vanish": all(
            value == 0 for value in remaining.values()
        ),
        "nonzero_E0_2_coefficients": {
            f"a^{a_degree}b^{b_degree}": rational_text(value)
            for (a_degree, b_degree), value in sorted(remaining.items())
            if value
        },
    }


def intercept_slice(order: int) -> dict[str, object]:
    """Lift the slice varying the common affine intercept and fixing slope."""

    matrix, _ = pivot_matrix()
    field = Field(series_order=order)
    coefficients = [[Q(0)] * (order + 1) for _ in SEED]
    for index, value in enumerate(SEED):
        coefficients[index][0] = value
    coefficients[6][1] = Q(1)
    remaining = []
    for degree in range(1, order + 1):
        values = residuals_from_jets(
            [FormalSeries.seed(field, row) for row in coefficients]
        )
        correction = solve_square_over_q(
            matrix, [-values[row].coefficients[degree] for row in PIVOT_ROWS]
        )
        for column, value in zip(PIVOT_COLUMNS, correction):
            coefficients[column][degree] = value
        checked = residuals_from_jets(
            [FormalSeries.seed(field, row) for row in coefficients]
        )
        if any(checked[row].coefficients[degree] for row in PIVOT_ROWS):
            raise AssertionError("the intercept slice left the implicit system")
        remaining.append(checked[6].coefficients[degree])
    if coefficients[4][1] != 1 or any(coefficients[4][degree] for degree in range(2, order + 1)):
        raise AssertionError("the first intercept no longer follows the common slice")
    if any(coefficients[5][degree] for degree in range(1, order + 1)):
        raise AssertionError("the first slope no longer stays fixed on the intercept slice")
    return {
        "order": order,
        "parameter": "x02=233/113+t, x12=97/113",
        "all_E0_2_coefficients_through_order_vanish": all(value == 0 for value in remaining),
        "E0_2_coefficients": [rational_text(value) for value in remaining],
        "conjugate_line_slice": {
            "x01": [rational_text(value) for value in coefficients[4]],
            "x11": [rational_text(value) for value in coefficients[5]],
            "x02": [rational_text(value) for value in coefficients[6]],
            "x12": [rational_text(value) for value in coefficients[7]],
        },
    }


def replay(*, bivariate_order: int = 4, slice_order: int = 12) -> dict[str, object]:
    directions = tangent_directions()
    modular_ranks = {}
    for prime in GOOD_TANGENT_PRIMES:
        values = residuals(SEED, Field(prime))
        rank, _ = row_reduce([value.gradient for value in values], Field(prime))
        if any(value.value for value in values) or rank != 6:
            raise AssertionError("a declared good-prime tangent replay changed")
        modular_ranks[str(prime)] = rank
    bivariate = bivariate_germ(bivariate_order)
    slice_data = intercept_slice(slice_order)
    root_data = root_motion_slice(slice_order)
    comparison = comparison_seed_audit()
    if not (
        bivariate["all_E0_2_coefficients_through_total_order_vanish"]
        and slice_data["all_E0_2_coefficients_through_order_vanish"]
        and root_data["all_E0_2_coefficients_through_order_vanish"]
    ):
        raise AssertionError("the conjugate-slope local continuation developed an obstruction")
    return {
        "status": "conjugate-slope two-section local germ replayed",
        "integer_roots": list(ROOTS),
        "normalized_moduli": [rational_text(value) for value in SEED[:4]],
        "seed_sections": [[rational_text(value) for value in SEED[4:6]], [rational_text(value) for value in SEED[6:8]]],
        "exact_tangent_rank_over_Q": 6,
        "modular_tangent_ranks": modular_ranks,
        "free_tangent_coordinates": [VARIABLES[index] for index in FREE_COLUMNS],
        "conjugate_tangent_directions": [[rational_text(value) for value in direction] for direction in directions],
        "bivariate_germ": bivariate,
        "intercept_slice": slice_data,
        "root_motion_slice": root_data,
        "comparison_seed": comparison,
        "conclusion": "the displayed conjugate-slope seed has an unobstructed bivariate formal germ through the declared orders",
        "not_established": [
            "a rational parametrization or global two-section component identity",
            "a Mordell--Weil independence statement, saturation, heights, intersections, or a Shioda Gram matrix",
            "generic rank at least 14",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bivariate-order", type=int, default=4)
    parser.add_argument("--slice-order", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        replay(bivariate_order=args.bivariate_order, slice_order=args.slice_order),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

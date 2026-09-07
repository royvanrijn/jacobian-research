#!/usr/bin/env python3
r"""Build the affine-orbit-normal quadratic Kuranishi slice at F_4.

This is a finite-field research compiler, not a characteristic-zero primary
decomposition certificate.  It performs the following operations over a
selected good prime:

1. construct the full degree-12 linearized Jacobian matrix at the explicit
   integer-root weighted map F_4;
2. retain a relation-tracked 58-vector tangent basis;
3. quotient the 22-dimensional affine left-right orbit;
4. use the normalized quartic seed tangent as the first of 36 normal
   coordinates;
5. compute and row-reduce the full quadratic Kuranishi ideal;
6. optionally emit a Singular program for Hilbert and primary analysis.

The companion exact-QQ checker certifies the ranks used to interpret this
modular exploration.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import w, x, y, z  # noqa: E402
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    explicit_seed,
    mapping_from_primitive,
    monomials_through,
)


VARIABLES = (x, y, z)


def coefficient_mod(value: sp.Expr, prime: int) -> int:
    value = sp.Rational(value)
    return (
        int(value.p)
        * pow(int(value.q) % prime, -1, prime)
        % prime
    )


def add_term(
    polynomial: dict[tuple[int, int, int], int],
    exponent: tuple[int, int, int],
    coefficient: int,
    prime: int,
) -> None:
    value = (polynomial.get(exponent, 0) + coefficient) % prime
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def add_scalar(
    vector: dict[int, int], index: int, coefficient: int, prime: int
) -> None:
    value = (vector.get(index, 0) + coefficient) % prime
    if value:
        vector[index] = value
    else:
        vector.pop(index, None)


def multiply(
    first: dict[tuple[int, int, int], int],
    second: dict[tuple[int, int, int], int],
    prime: int,
) -> dict[tuple[int, int, int], int]:
    output: dict[tuple[int, int, int], int] = {}
    for exponent_first, coefficient_first in first.items():
        for exponent_second, coefficient_second in second.items():
            add_term(
                output,
                tuple(
                    exponent_first[index] + exponent_second[index]
                    for index in range(3)
                ),
                coefficient_first * coefficient_second,
                prime,
            )
    return output


def derivative(
    polynomial: dict[tuple[int, int, int], int],
    variable: int,
    prime: int,
) -> dict[tuple[int, int, int], int]:
    output: dict[tuple[int, int, int], int] = {}
    for exponent, coefficient in polynomial.items():
        if exponent[variable] == 0:
            continue
        new_exponent = list(exponent)
        new_exponent[variable] -= 1
        output[tuple(new_exponent)] = (
            coefficient * exponent[variable] % prime
        )
    return output


def sparse_polynomial(
    expression: sp.Expr, prime: int
) -> dict[tuple[int, int, int], int]:
    return {
        exponent: coefficient_mod(coefficient, prime)
        for exponent, coefficient in sp.Poly(
            expression, *VARIABLES
        ).terms()
        if coefficient_mod(coefficient, prime)
    }


def leading_exponent(
    polynomial: dict[tuple[int, int, int], int],
) -> tuple[int, int, int]:
    return max(
        polynomial, key=lambda exponent: (sum(exponent), exponent)
    )


def reduce_polynomial(
    polynomial: dict[tuple[int, int, int], int],
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    prime: int,
) -> dict[tuple[int, int, int], int]:
    working = dict(polynomial)
    remainder: dict[tuple[int, int, int], int] = {}
    while working:
        lead = leading_exponent(working)
        coefficient = working.pop(lead)
        pivot = pivots.get(lead)
        if pivot is None:
            remainder[lead] = coefficient
            continue
        for exponent, value in pivot.items():
            if exponent == lead:
                continue
            add_term(
                working, exponent, -coefficient * value, prime
            )
    return remainder


def row_reduce_vectors(
    vectors: list[dict[int, int]], prime: int
) -> tuple[list[dict[int, int]], list[int]]:
    """Return normalized echelon rows and their original indices."""
    pivots: dict[int, dict[int, int]] = {}
    selected: list[dict[int, int]] = []
    selected_indices: list[int] = []
    for original_index, original in enumerate(vectors):
        vector = dict(original)
        while vector:
            lead = max(vector)
            coefficient = vector[lead]
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(coefficient, -1, prime)
                normalized = {
                    index: value * inverse % prime
                    for index, value in vector.items()
                }
                pivots[lead] = normalized
                selected.append(normalized)
                selected_indices.append(original_index)
                break
            for index, value in pivot.items():
                add_scalar(
                    vector, index, -coefficient * value, prime
                )
    return selected, selected_indices


def solve_square(
    columns: list[dict[int, int]],
    rows: list[int],
    right_hand_side: dict[int, int],
    prime: int,
) -> list[int]:
    """Solve a square column system restricted to the selected rows."""
    size = len(columns)
    augmented = [
        [
            columns[column].get(rows[row], 0)
            for column in range(size)
        ]
        + [right_hand_side.get(rows[row], 0)]
        for row in range(size)
    ]
    for column in range(size):
        pivot_row = next(
            row
            for row in range(column, size)
            if augmented[row][column] % prime
        )
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[column],
        )
        inverse = pow(augmented[column][column] % prime, -1, prime)
        augmented[column] = [
            value * inverse % prime
            for value in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            coefficient = augmented[row][column] % prime
            if coefficient:
                augmented[row] = [
                    (
                        augmented[row][index]
                        - coefficient * augmented[column][index]
                    )
                    % prime
                    for index in range(size + 1)
                ]
    return [augmented[index][-1] % prime for index in range(size)]


def tangent_echelon(
    adjugate: sp.Matrix, degree: int, prime: int
) -> tuple[
    list[tuple[int, int, int]],
    list[int],
    list[dict[int, int]],
    dict[tuple[int, int, int], dict[tuple[int, int, int], int]],
    dict[tuple[int, int, int], dict[int, int]],
    list[list[dict[tuple[int, int, int], int]]],
]:
    monomials = monomials_through(degree)
    adjugate_terms = [
        [
            sparse_polynomial(adjugate[row, column], prime)
            for column in range(3)
        ]
        for row in range(3)
    ]
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ] = {}
    pivot_combinations: dict[
        tuple[int, int, int], dict[int, int]
    ] = {}
    relations: list[dict[int, int]] = []
    free_columns: list[int] = []
    column_index = 0

    for component in range(3):
        for exponent in monomials:
            column: dict[tuple[int, int, int], int] = {}
            for variable in range(3):
                if exponent[variable] == 0:
                    continue
                base = list(exponent)
                base[variable] -= 1
                for adjugate_exponent, coefficient in (
                    adjugate_terms[variable][component].items()
                ):
                    add_term(
                        column,
                        tuple(
                            base[index] + adjugate_exponent[index]
                            for index in range(3)
                        ),
                        exponent[variable] * coefficient,
                        prime,
                    )
            combination = {column_index: 1}
            while column:
                lead = leading_exponent(column)
                coefficient = column[lead]
                pivot = pivots.get(lead)
                if pivot is None:
                    inverse = pow(coefficient, -1, prime)
                    pivots[lead] = {
                        monomial: value * inverse % prime
                        for monomial, value in column.items()
                    }
                    pivot_combinations[lead] = {
                        index: value * inverse % prime
                        for index, value in combination.items()
                    }
                    break
                for monomial, value in pivot.items():
                    add_term(
                        column,
                        monomial,
                        -coefficient * value,
                        prime,
                    )
                for index, value in pivot_combinations[lead].items():
                    add_scalar(
                        combination,
                        index,
                        -coefficient * value,
                        prime,
                    )
            if not column:
                relations.append(combination)
                free_columns.append(column_index)
            column_index += 1

    return (
        monomials,
        free_columns,
        relations,
        pivots,
        pivot_combinations,
        adjugate_terms,
    )


def relation_matrices(
    relations: list[dict[int, int]],
    monomials: list[tuple[int, int, int]],
    adjugate_terms: list[
        list[dict[tuple[int, int, int], int]]
    ],
    prime: int,
) -> tuple[
    list[list[list[dict[tuple[int, int, int], int]]]],
    list[list[list[dict[tuple[int, int, int], int]]]],
]:
    matrices = []
    jacobians = []
    monomial_count = len(monomials)
    for relation in relations:
        direction = [{} for _ in range(3)]
        for index, coefficient in relation.items():
            component = index // monomial_count
            add_term(
                direction[component],
                monomials[index % monomial_count],
                coefficient,
                prime,
            )
        direction_jacobian = [
            [
                derivative(direction[row], column, prime)
                for column in range(3)
            ]
            for row in range(3)
        ]
        matrix = [[{} for _ in range(3)] for _ in range(3)]
        for row in range(3):
            for column in range(3):
                entry: dict[tuple[int, int, int], int] = {}
                for middle in range(3):
                    product = multiply(
                        adjugate_terms[row][middle],
                        direction_jacobian[middle][column],
                        prime,
                    )
                    for exponent, coefficient in product.items():
                        add_term(
                            entry, exponent, coefficient, prime
                        )
                matrix[row][column] = entry
        matrices.append(matrix)
        jacobians.append(direction_jacobian)
    return matrices, jacobians


def free_coordinates(
    direction: sp.Matrix,
    free_columns: list[int],
    monomials: list[tuple[int, int, int]],
    prime: int,
) -> dict[int, int]:
    coordinate: dict[int, int] = {}
    monomial_count = len(monomials)
    polynomials = [
        sp.Poly(sp.expand(direction[component]), *VARIABLES)
        for component in range(3)
    ]
    for tangent_index, column_index in enumerate(free_columns):
        component = column_index // monomial_count
        exponent = monomials[column_index % monomial_count]
        coefficient = coefficient_mod(
            polynomials[component].coeff_monomial(exponent), prime
        )
        if coefficient:
            coordinate[tangent_index] = coefficient
    return coordinate


def affine_directions(mapping: sp.Matrix, jacobian: sp.Matrix) -> list[sp.Matrix]:
    source_vector = sp.Matrix(VARIABLES)
    directions: list[sp.Matrix] = []
    for index in range(3):
        translation = sp.zeros(3, 1)
        translation[index] = 1
        directions.append(translation)
    for index in range(3):
        translation = sp.zeros(3, 1)
        translation[index] = 1
        directions.append(jacobian * translation)
    for side in ("target", "source"):
        for row in range(3):
            for column in range(3):
                if row == column:
                    continue
                matrix = sp.zeros(3)
                matrix[row, column] = 1
                directions.append(
                    matrix * mapping
                    if side == "target"
                    else jacobian * (matrix * source_vector)
                )
    diagonal_data = (
        (sp.diag(1, 0, -1), sp.zeros(3)),
        (sp.diag(0, 1, -1), sp.zeros(3)),
        (sp.diag(-1, 0, 0), sp.diag(1, 0, 0)),
        (sp.diag(-1, 0, 0), sp.diag(0, 1, 0)),
        (sp.diag(-1, 0, 0), sp.diag(0, 0, 1)),
    )
    for target_matrix, source_matrix in diagonal_data:
        directions.append(
            target_matrix * mapping
            + jacobian * (source_matrix * source_vector)
        )
    return directions


def independent_columns(
    columns: list[dict[int, int]], prime: int
) -> tuple[list[dict[int, int]], list[int]]:
    pivots: dict[int, dict[int, int]] = {}
    selected: list[dict[int, int]] = []
    selected_indices: list[int] = []
    for original_index, original in enumerate(columns):
        vector = dict(original)
        while vector:
            lead = max(vector)
            coefficient = vector[lead]
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(coefficient, -1, prime)
                normalized = {
                    index: value * inverse % prime
                    for index, value in vector.items()
                }
                pivots[lead] = normalized
                selected.append(original)
                selected_indices.append(original_index)
                break
            for index, value in pivot.items():
                add_scalar(
                    vector, index, -coefficient * value, prime
                )
    return selected, selected_indices


def combine_matrices(
    coefficients: dict[int, int],
    matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    prime: int,
) -> list[list[dict[tuple[int, int, int], int]]]:
    output = [[{} for _ in range(3)] for _ in range(3)]
    for matrix_index, scalar in coefficients.items():
        for row in range(3):
            for column in range(3):
                for exponent, coefficient in (
                    matrices[matrix_index][row][column].items()
                ):
                    add_term(
                        output[row][column],
                        exponent,
                        scalar * coefficient,
                        prime,
                    )
    return output


def trace_product(
    first: list[list[dict[tuple[int, int, int], int]]],
    second: list[list[dict[tuple[int, int, int], int]]],
    prime: int,
    scalar: int = 1,
) -> dict[tuple[int, int, int], int]:
    output: dict[tuple[int, int, int], int] = {}
    for row in range(3):
        for column in range(3):
            product = multiply(
                first[row][column],
                second[column][row],
                prime,
            )
            for exponent, coefficient in product.items():
                add_term(
                    output,
                    exponent,
                    scalar * coefficient,
                    prime,
                )
    return output


def determinant_three(
    matrix: list[list[dict[tuple[int, int, int], int]]],
    prime: int,
) -> dict[tuple[int, int, int], int]:
    permutations = (
        ((0, 1, 2), 1),
        ((0, 2, 1), -1),
        ((1, 0, 2), -1),
        ((1, 2, 0), 1),
        ((2, 0, 1), 1),
        ((2, 1, 0), -1),
    )
    output: dict[tuple[int, int, int], int] = {}
    for permutation, sign in permutations:
        product = multiply(
            multiply(
                matrix[0][permutation[0]],
                matrix[1][permutation[1]],
                prime,
            ),
            matrix[2][permutation[2]],
            prime,
        )
        for exponent, coefficient in product.items():
            add_term(
                output, exponent, sign * coefficient, prime
            )
    return output


def solve_linearized_image(
    polynomial: dict[tuple[int, int, int], int],
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    pivot_combinations: dict[
        tuple[int, int, int], dict[int, int]
    ],
    prime: int,
) -> tuple[dict[int, int], dict[tuple[int, int, int], int]]:
    """Write polynomial=L(solution)+remainder in the fixed echelon split."""
    working = dict(polynomial)
    solution: dict[int, int] = {}
    remainder: dict[tuple[int, int, int], int] = {}
    while working:
        lead = leading_exponent(working)
        coefficient = working.pop(lead)
        pivot = pivots.get(lead)
        if pivot is None:
            remainder[lead] = coefficient
            continue
        for exponent, value in pivot.items():
            if exponent != lead:
                add_term(
                    working,
                    exponent,
                    -coefficient * value,
                    prime,
                )
        for index, value in pivot_combinations[lead].items():
            add_scalar(
                solution, index, coefficient * value, prime
            )
    return solution, remainder


def coefficient_direction_jacobian(
    coefficients: dict[int, int],
    monomials: list[tuple[int, int, int]],
    prime: int,
) -> list[list[dict[tuple[int, int, int], int]]]:
    direction = [{} for _ in range(3)]
    monomial_count = len(monomials)
    for index, coefficient in coefficients.items():
        add_term(
            direction[index // monomial_count],
            monomials[index % monomial_count],
            coefficient,
            prime,
        )
    return [
        [
            derivative(direction[row], column, prime)
            for column in range(3)
        ]
        for row in range(3)
    ]


def left_multiply_adjugate(
    adjugate_terms: list[
        list[dict[tuple[int, int, int], int]]
    ],
    matrix: list[list[dict[tuple[int, int, int], int]]],
    prime: int,
) -> list[list[dict[tuple[int, int, int], int]]]:
    output = [[{} for _ in range(3)] for _ in range(3)]
    for row in range(3):
        for column in range(3):
            for middle in range(3):
                product = multiply(
                    adjugate_terms[row][middle],
                    matrix[middle][column],
                    prime,
                )
                for exponent, coefficient in product.items():
                    add_term(
                        output[row][column],
                        exponent,
                        coefficient,
                        prime,
                    )
    return output


def vector_span_membership(
    vectors: list[dict[tuple[int, int, int], int]],
    target: dict[tuple[int, int, int], int],
    prime: int,
) -> tuple[int, bool]:
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ] = {}
    for original in vectors:
        vector = dict(original)
        while vector:
            lead = leading_exponent(vector)
            coefficient = vector[lead]
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(coefficient, -1, prime)
                pivots[lead] = {
                    exponent: value * inverse % prime
                    for exponent, value in vector.items()
                }
                break
            for exponent, value in pivot.items():
                add_term(
                    vector,
                    exponent,
                    -coefficient * value,
                    prime,
                )
    reduced = dict(target)
    while reduced:
        lead = leading_exponent(reduced)
        coefficient = reduced[lead]
        pivot = pivots.get(lead)
        if pivot is None:
            return len(pivots), False
        for exponent, value in pivot.items():
            add_term(
                reduced,
                exponent,
                -coefficient * value,
                prime,
            )
    return len(pivots), True


def vector_span_solution(
    vectors: list[dict[tuple[int, int, int], int]],
    target: dict[tuple[int, int, int], int],
    prime: int,
) -> dict[int, int] | None:
    pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ] = {}
    pivot_combinations: dict[
        tuple[int, int, int], dict[int, int]
    ] = {}
    for vector_index, original in enumerate(vectors):
        vector = dict(original)
        combination = {vector_index: 1}
        while vector:
            lead = leading_exponent(vector)
            coefficient = vector[lead]
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(coefficient, -1, prime)
                pivots[lead] = {
                    exponent: value * inverse % prime
                    for exponent, value in vector.items()
                }
                pivot_combinations[lead] = {
                    index: value * inverse % prime
                    for index, value in combination.items()
                }
                break
            for exponent, value in pivot.items():
                add_term(
                    vector,
                    exponent,
                    -coefficient * value,
                    prime,
                )
            for index, value in pivot_combinations[lead].items():
                add_scalar(
                    combination,
                    index,
                    -coefficient * value,
                    prime,
                )
    working = dict(target)
    solution: dict[int, int] = {}
    while working:
        lead = leading_exponent(working)
        coefficient = working[lead]
        pivot = pivots.get(lead)
        if pivot is None:
            return None
        for exponent, value in pivot.items():
            add_term(
                working,
                exponent,
                -coefficient * value,
                prime,
            )
        for index, value in pivot_combinations[lead].items():
            add_scalar(
                solution, index, coefficient * value, prime
            )
    return solution


def add_matrix(
    first: list[list[dict[tuple[int, int, int], int]]],
    second: list[list[dict[tuple[int, int, int], int]]],
    prime: int,
) -> list[list[dict[tuple[int, int, int], int]]]:
    output = [
        [dict(first[row][column]) for column in range(3)]
        for row in range(3)
    ]
    for row in range(3):
        for column in range(3):
            for exponent, coefficient in second[row][column].items():
                add_term(
                    output[row][column],
                    exponent,
                    coefficient,
                    prime,
                )
    return output


def determinant_series_coefficient(
    series: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    order: int,
    prime: int,
) -> dict[tuple[int, int, int], int]:
    permutations = (
        ((0, 1, 2), 1),
        ((0, 2, 1), -1),
        ((1, 0, 2), -1),
        ((1, 2, 0), 1),
        ((2, 0, 1), 1),
        ((2, 1, 0), -1),
    )
    output: dict[tuple[int, int, int], int] = {}
    maximum_order = len(series) - 1
    for permutation, sign in permutations:
        for first_order in range(
            min(order, maximum_order) + 1
        ):
            for second_order in range(
                min(order - first_order, maximum_order) + 1
            ):
                third_order = order - first_order - second_order
                if not 0 <= third_order <= maximum_order:
                    continue
                product = multiply(
                    multiply(
                        series[first_order][0][permutation[0]],
                        series[second_order][1][permutation[1]],
                        prime,
                    ),
                    series[third_order][2][permutation[2]],
                    prime,
                )
                for exponent, coefficient in product.items():
                    add_term(
                        output,
                        exponent,
                        sign * coefficient,
                        prime,
                    )
    return output


def jet_lift_axis(
    name: str,
    first_jacobian: list[
        list[dict[tuple[int, int, int], int]]
    ],
    first_matrix: list[
        list[dict[tuple[int, int, int], int]]
    ],
    jacobian_terms: list[
        list[dict[tuple[int, int, int], int]]
    ],
    tangent_jacobians: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    tangent_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    monomials: list[tuple[int, int, int]],
    linear_pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    pivot_combinations: dict[
        tuple[int, int, int], dict[int, int]
    ],
    prime: int,
    maximum_order: int,
) -> dict[str, object]:
    gauge_columns = [
        reduce_polynomial(
            trace_product(
                first_matrix,
                tangent_matrix,
                prime,
                scalar=-1,
            ),
            linear_pivots,
            prime,
        )
        for tangent_matrix in tangent_matrices
    ]
    series = [jacobian_terms, first_jacobian]
    adjustment_counts: list[int] = []
    correction_sizes: list[int] = []
    for order in range(2, maximum_order + 1):
        obstruction = determinant_series_coefficient(
            series, order, prime
        )
        solution, remainder = solve_linearized_image(
            obstruction,
            linear_pivots,
            pivot_combinations,
            prime,
        )
        if remainder:
            return {
                "axis": name,
                "lifts_through_order": order - 1,
                "first_failed_order": order,
                "failure_remainder_terms": len(remainder),
                "tangent_adjustment_counts": adjustment_counts,
                "correction_coefficient_counts": correction_sizes,
            }
        correction_coefficients = {
            index: -coefficient % prime
            for index, coefficient in solution.items()
        }
        correction_sizes.append(len(correction_coefficients))
        correction_jacobian = coefficient_direction_jacobian(
            correction_coefficients, monomials, prime
        )
        series.append(correction_jacobian)
        if order == maximum_order:
            break
        next_obstruction = determinant_series_coefficient(
            series, order + 1, prime
        )
        next_remainder = reduce_polynomial(
            next_obstruction, linear_pivots, prime
        )
        target = {
            exponent: -coefficient % prime
            for exponent, coefficient in next_remainder.items()
        }
        tangent_solution = vector_span_solution(
            gauge_columns, target, prime
        )
        if tangent_solution is None:
            return {
                "axis": name,
                "lifts_through_order": order,
                "first_failed_order": order + 1,
                "failure_remainder_terms": len(next_remainder),
                "tangent_adjustment_counts": adjustment_counts,
                "correction_coefficient_counts": correction_sizes,
            }
        adjustment_counts.append(len(tangent_solution))
        tangent_correction = combine_matrices(
            tangent_solution, tangent_jacobians, prime
        )
        series[-1] = add_matrix(
            series[-1], tangent_correction, prime
        )
        adjusted_next = determinant_series_coefficient(
            series, order + 1, prime
        )
        assert not reduce_polynomial(
            adjusted_next, linear_pivots, prime
        ), (
            f"failed to apply tangent lookahead for {name} "
            f"at order {order + 1}"
        )
    return {
        "axis": name,
        "lifts_through_order": maximum_order,
        "first_failed_order": None,
        "failure_remainder_terms": 0,
        "tangent_adjustment_counts": adjustment_counts,
        "correction_coefficient_counts": correction_sizes,
    }


def cubic_axis_screen(
    variable_names: list[str],
    normal_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    normal_jacobians: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    tangent_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    monomials: list[tuple[int, int, int]],
    linear_pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    pivot_combinations: dict[
        tuple[int, int, int], dict[int, int]
    ],
    adjugate_terms: list[
        list[dict[tuple[int, int, int], int]]
    ],
    prime: int,
) -> list[dict[str, object]]:
    output = []
    inverse_two = pow(2, -1, prime)
    for index, name in enumerate(variable_names):
        matrix = normal_matrices[index]
        quadratic = trace_product(
            matrix, matrix, prime, scalar=-inverse_two
        )
        order_two_solution, quadratic_remainder = solve_linearized_image(
            quadratic,
            linear_pivots,
            pivot_combinations,
            prime,
        )
        if quadratic_remainder:
            continue
        correction_coefficients = {
            coefficient_index: -coefficient % prime
            for coefficient_index, coefficient in order_two_solution.items()
        }
        correction_jacobian = coefficient_direction_jacobian(
            correction_coefficients, monomials, prime
        )
        correction_matrix = left_multiply_adjugate(
            adjugate_terms, correction_jacobian, prime
        )
        cubic = trace_product(
            matrix, correction_matrix, prime, scalar=-1
        )
        determinant = determinant_three(
            normal_jacobians[index], prime
        )
        for exponent, coefficient in determinant.items():
            add_term(cubic, exponent, coefficient, prime)
        cubic_remainder = reduce_polynomial(
            cubic, linear_pivots, prime
        )
        gauge_columns = [
            reduce_polynomial(
                trace_product(
                    matrix,
                    tangent_matrix,
                    prime,
                    scalar=-1,
                ),
                linear_pivots,
                prime,
            )
            for tangent_matrix in tangent_matrices
        ]
        gauge_rank, liftable = vector_span_membership(
            gauge_columns, cubic_remainder, prime
        )
        output.append(
            {
                "axis": name,
                "quadratic_lift_coefficient_count": len(
                    correction_coefficients
                ),
                "cubic_remainder_term_count": len(cubic_remainder),
                "cubic_gauge_rank": gauge_rank,
                "lifts_through_order_three": liftable,
            }
        )
    return output


def quadratic_remainders(
    matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    linear_pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    prime: int,
) -> tuple[
    list[tuple[int, int]],
    list[dict[tuple[int, int, int], int]],
]:
    pairs: list[tuple[int, int]] = []
    remainders: list[dict[tuple[int, int, int], int]] = []
    inverse_two = pow(2, -1, prime)
    for first in range(len(matrices)):
        for second in range(first, len(matrices)):
            polynomial: dict[tuple[int, int, int], int] = {}
            for row in range(3):
                for column in range(3):
                    product = multiply(
                        matrices[first][row][column],
                        matrices[second][column][row],
                        prime,
                    )
                    for exponent, coefficient in product.items():
                        add_term(
                            polynomial,
                            exponent,
                            -coefficient,
                            prime,
                        )
            if first == second:
                polynomial = {
                    exponent: coefficient * inverse_two % prime
                    for exponent, coefficient in polynomial.items()
                }
            pairs.append((first, second))
            remainders.append(
                reduce_polynomial(polynomial, linear_pivots, prime)
            )
    return pairs, remainders


def independent_quadrics(
    pairs: list[tuple[int, int]],
    remainders: list[dict[tuple[int, int, int], int]],
    prime: int,
) -> list[dict[int, int]]:
    equations_by_cokernel_monomial: dict[
        tuple[int, int, int], dict[int, int]
    ] = {}
    for pair_index, remainder in enumerate(remainders):
        for exponent, coefficient in remainder.items():
            add_scalar(
                equations_by_cokernel_monomial.setdefault(exponent, {}),
                pair_index,
                coefficient,
                prime,
            )
    equations = [
        equations_by_cokernel_monomial[exponent]
        for exponent in sorted(
            equations_by_cokernel_monomial,
            key=lambda value: (sum(value), value),
        )
    ]
    independent, _indices = row_reduce_vectors(equations, prime)
    return independent


def singular_polynomial(
    equation: dict[int, int],
    pairs: list[tuple[int, int]],
    variable_names: list[str],
    prime: int,
) -> str:
    terms = []
    for pair_index, coefficient in sorted(equation.items()):
        first, second = pairs[pair_index]
        symmetric = coefficient if coefficient <= prime // 2 else coefficient - prime
        monomial = variable_names[first]
        if first == second:
            monomial += "^2"
        else:
            monomial += f"*{variable_names[second]}"
        terms.append(f"({symmetric})*{monomial}")
    return "+".join(terms) if terms else "0"


def emit_singular(
    path: Path,
    prime: int,
    variable_names: list[str],
    pairs: list[tuple[int, int]],
    equations: list[dict[int, int]],
) -> None:
    polynomial_lines = [
        singular_polynomial(
            equation, pairs, variable_names, prime
        )
        for equation in equations
    ]
    text = "\n".join(
        [
            f"ring r={prime},({','.join(variable_names)}),dp;",
            "ideal I=",
            ",\n".join(polynomial_lines) + ";",
            'print("KURANISHI_IDEAL_GENERATORS");',
            "print(size(I));",
            "ideal G=slimgb(I);",
            'print("KURANISHI_STANDARD_BASIS");',
            "print(size(G));",
            'print("KURANISHI_DIMENSION");',
            "print(dim(G));",
            'print("KURANISHI_DEGREE");',
            "print(deg(G));",
            'print("KURANISHI_HILBERT");',
            "print(hilb(G));",
            'print("KURANISHI_MINASS_START");',
            'LIB "primdec.lib";',
            "list L=minAssGTZ(I);",
            "print(size(L));",
            "for (int k=1; k<=size(L); k++)",
            "{",
            '  print("MINASS_COMPONENT");',
            "  print(k);",
            "  print(dim(std(L[k])));",
            "  print(deg(std(L[k])));",
            "  print(L[k]);",
            "}",
            'print("KURANISHI_MINASS_END");',
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def emit_macaulay2(
    path: Path,
    prime: int,
    variable_names: list[str],
    pairs: list[tuple[int, int]],
    equations: list[dict[int, int]],
) -> None:
    polynomial_lines = [
        singular_polynomial(
            equation, pairs, variable_names, prime
        ).replace("^", "^")
        for equation in equations
    ]
    text = "\n".join(
        [
            "needsPackage \"PrimaryDecomposition\";",
            (
                f"R=GF({prime})["
                + ",".join(variable_names)
                + ",MonomialOrder=>GRevLex];"
            ),
            "I=ideal(",
            ",\n".join(polynomial_lines) + ");",
            'print "KURANISHI_IDEAL_GENERATORS";',
            "print numgens source gens I;",
            "elapsed G=gb I;",
            'print "KURANISHI_DIMENSION";',
            "print dim I;",
            'print "KURANISHI_DEGREE";',
            "print degree I;",
            'print "KURANISHI_HILBERT";',
            "print hilbertSeries I;",
            'print "KURANISHI_MINASS_START";',
            "elapsed L=minimalPrimes I;",
            "print #L;",
            "scan(L,Q -> (print dim Q; print degree Q; print Q));",
            'print "KURANISHI_MINASS_END";',
            "exit 0;",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--singular-output", type=Path)
    parser.add_argument("--macaulay2-output", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument(
        "--jet-order",
        type=int,
        default=8,
        help="largest bounded-box jet order for coordinate-axis lifting",
    )
    args = parser.parse_args()
    prime = args.prime
    if not sp.isprime(prime) or prime in (2, 3, 5, 17):
        parser.error("choose a good odd prime away from 3, 5, and 17")
    if args.jet_order < 2:
        parser.error("--jet-order must be at least two")

    mapping = mapping_from_primitive(explicit_seed(4))
    jacobian = mapping.jacobian(VARIABLES)
    jacobian_terms = [
        [
            sparse_polynomial(jacobian[row, column], prime)
            for column in range(3)
        ]
        for row in range(3)
    ]
    (
        monomials,
        free_columns,
        relations,
        linear_pivots,
        pivot_combinations,
        adjugate_terms,
    ) = tangent_echelon(jacobian.adjugate(), 12, prime)
    assert len(relations) == 58
    relation_matrix_entries, relation_jacobians = relation_matrices(
        relations, monomials, adjugate_terms, prime
    )

    affine_columns = [
        free_coordinates(
            direction, free_columns, monomials, prime
        )
        for direction in affine_directions(mapping, jacobian)
    ]
    independent_affine, independent_affine_indices = independent_columns(
        affine_columns, prime
    )
    assert len(independent_affine) == 22
    affine_rows, orbit_pivot_rows = row_reduce_vectors(
        [
            {
                column_index: column.get(row_index, 0)
                for column_index, column in enumerate(independent_affine)
                if column.get(row_index, 0)
            }
            for row_index in range(58)
        ],
        prime,
    )
    del affine_rows
    assert len(orbit_pivot_rows) == 22
    normal_indices = [
        index for index in range(58) if index not in orbit_pivot_rows
    ]
    assert len(normal_indices) == 36

    seed_parameter = sp.symbols("seed_parameter")
    seed_family = mapping_from_primitive(
        explicit_seed(4)
        + seed_parameter * w**2 * (w - 1) ** 2
    )
    seed_direction = sp.Matrix(
        [
            sp.diff(component, seed_parameter).subs(seed_parameter, 0)
            for component in seed_family
        ]
    )
    seed_coordinate = free_coordinates(
        seed_direction, free_columns, monomials, prime
    )
    affine_solution = solve_square(
        independent_affine,
        orbit_pivot_rows,
        seed_coordinate,
        prime,
    )
    normal_seed = dict(seed_coordinate)
    for column, scalar in zip(
        independent_affine, affine_solution, strict=True
    ):
        for index, coefficient in column.items():
            add_scalar(
                normal_seed, index, -scalar * coefficient, prime
            )
    assert all(
        normal_seed.get(index, 0) == 0
        for index in orbit_pivot_rows
    )
    removed_normal_index = next(
        index for index in normal_indices if normal_seed.get(index, 0)
    )

    normal_basis = [normal_seed] + [
        {index: 1}
        for index in normal_indices
        if index != removed_normal_index
    ]
    assert len(normal_basis) == 36
    normal_matrices = [
        combine_matrices(
            coordinate, relation_matrix_entries, prime
        )
        for coordinate in normal_basis
    ]
    normal_jacobians = [
        combine_matrices(
            coordinate, relation_jacobians, prime
        )
        for coordinate in normal_basis
    ]
    pairs, remainders = quadratic_remainders(
        normal_matrices, linear_pivots, prime
    )
    equations = independent_quadrics(pairs, remainders, prime)
    assert len(equations) == 53, (
        "unexpected normal-slice quadratic rank "
        f"{len(equations)} modulo {prime}"
    )

    variable_names = ["seed"] + [
        f"n{index}" for index in range(1, 36)
    ]
    seed_pair = pairs.index((0, 0))
    assert all(
        equation.get(seed_pair, 0) == 0 for equation in equations
    )
    cubic_screen = cubic_axis_screen(
        variable_names,
        normal_matrices,
        normal_jacobians,
        relation_matrix_entries,
        monomials,
        linear_pivots,
        pivot_combinations,
        adjugate_terms,
        prime,
    )
    unobstructed_axis_indices = [
        index
        for index in range(len(variable_names))
        if not remainders[pairs.index((index, index))]
    ]
    jet_screen = [
        jet_lift_axis(
            variable_names[index],
            normal_jacobians[index],
            normal_matrices[index],
            jacobian_terms,
            relation_jacobians,
            relation_matrix_entries,
            monomials,
            linear_pivots,
            pivot_combinations,
            prime,
            args.jet_order,
        )
        for index in unobstructed_axis_indices
    ]

    summary = {
        "prime": prime,
        "ambient_tangent_dimension": 58,
        "affine_lr_parameter_count": 23,
        "affine_lr_orbit_dimension": 22,
        "affine_stabilizer_dimension": 1,
        "normal_slice_dimension": 36,
        "visible_seed_coordinate": variable_names[0],
        "quadratic_pair_count": len(pairs),
        "nonzero_quadratic_pair_count": sum(
            bool(remainder) for remainder in remainders
        ),
        "quadratically_unobstructed_coordinate_axes": [
            variable_names[index]
            for index in range(len(variable_names))
            if not remainders[pairs.index((index, index))]
        ],
        "nonzero_seed_cross_pair_count": sum(
            bool(remainders[pairs.index((0, index))])
            for index in range(1, len(variable_names))
        ),
        "quadratic_kuranishi_rank": len(equations),
        "cokernel_monomial_count": len(
            {
                exponent
                for remainder in remainders
                for exponent in remainder
            }
        ),
        "independent_affine_parameter_indices": independent_affine_indices,
        "orbit_pivot_tangent_indices": orbit_pivot_rows,
        "normal_tangent_indices": normal_indices,
        "removed_normal_index_for_seed": removed_normal_index,
        "normal_seed_support": {
            str(index): coefficient
            for index, coefficient in sorted(normal_seed.items())
        },
        "cubic_coordinate_axis_screen": cubic_screen,
        "bounded_jet_coordinate_axis_screen": jet_screen,
        "bounded_jet_maximum_order": args.jet_order,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.singular_output is not None:
        emit_singular(
            args.singular_output,
            prime,
            variable_names,
            pairs,
            equations,
        )
    if args.macaulay2_output is not None:
        emit_macaulay2(
            args.macaulay2_output,
            prime,
            variable_names,
            pairs,
            equations,
        )
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()

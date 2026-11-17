#!/usr/bin/env python3
r"""Exact QQ certificate for the full-box quadratic Kuranishi map at F_4.

The coefficient box is X(3,12).  The checker constructs a relation-tracked
58-vector basis of ker(L_F), computes all 1711 polarized quadratic
determinant classes, reduces them by a canonical full echelon normal form in
coker(L_F), and certifies that their span has rank 53 over QQ.

It also certifies the dimensions of three explicit reduced families through
the point:

* affine left-right orbit: 22 (one-dimensional weighted stabilizer);
* affine orbit plus the four target shears adding C^2,C^3 to A or B: 26;
* the preceding family plus the normalized quartic seed direction: 27.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import w, x, y, z  # noqa: E402
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    explicit_seed,
    mapping_from_primitive,
    monomials_through,
    rational,
)


VARIABLES = (x, y, z)
Exponent = tuple[int, int, int]
Polynomial = dict[Exponent, Fraction]


def add_term(
    polynomial: Polynomial, exponent: Exponent, coefficient: Fraction
) -> None:
    value = polynomial.get(exponent, Fraction()) + coefficient
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def multiply(first: Polynomial, second: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for exponent_first, coefficient_first in first.items():
        for exponent_second, coefficient_second in second.items():
            add_term(
                output,
                tuple(
                    exponent_first[index] + exponent_second[index]
                    for index in range(3)
                ),
                coefficient_first * coefficient_second,
            )
    return output


def derivative(polynomial: Polynomial, variable: int) -> Polynomial:
    output: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        if exponent[variable] == 0:
            continue
        new_exponent = list(exponent)
        new_exponent[variable] -= 1
        output[tuple(new_exponent)] = (
            coefficient * exponent[variable]
        )
    return output


def sparse_polynomial(expression: sp.Expr) -> Polynomial:
    return {
        exponent: rational(coefficient)
        for exponent, coefficient in sp.Poly(
            expression, *VARIABLES
        ).terms()
    }


def leading_exponent(polynomial: Polynomial) -> Exponent:
    return max(
        polynomial, key=lambda exponent: (sum(exponent), exponent)
    )


def tangent_echelon(
    adjugate: sp.Matrix, degree: int
) -> tuple[
    list[Exponent],
    list[int],
    list[dict[int, Fraction]],
    dict[Exponent, Polynomial],
]:
    monomials = monomials_through(degree)
    adjugate_terms = [
        [
            sparse_polynomial(adjugate[row, column])
            for column in range(3)
        ]
        for row in range(3)
    ]
    pivots: dict[Exponent, Polynomial] = {}
    pivot_combinations: dict[Exponent, dict[int, Fraction]] = {}
    relations: list[dict[int, Fraction]] = []
    free_columns: list[int] = []
    column_index = 0

    for component in range(3):
        for exponent in monomials:
            column: Polynomial = {}
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
                    )
            combination = {column_index: Fraction(1)}
            while column:
                lead = leading_exponent(column)
                coefficient = column[lead]
                pivot = pivots.get(lead)
                if pivot is None:
                    pivots[lead] = {
                        monomial: value / coefficient
                        for monomial, value in column.items()
                    }
                    pivot_combinations[lead] = {
                        index: value / coefficient
                        for index, value in combination.items()
                    }
                    break
                for monomial, value in pivot.items():
                    add_term(
                        column, monomial, -coefficient * value
                    )
                for index, value in pivot_combinations[lead].items():
                    new_value = (
                        combination.get(index, Fraction())
                        - coefficient * value
                    )
                    if new_value:
                        combination[index] = new_value
                    else:
                        combination.pop(index, None)
            if not column:
                relations.append(combination)
                free_columns.append(column_index)
            column_index += 1

    return monomials, free_columns, relations, pivots


def full_remainder(
    polynomial: Polynomial, pivots: dict[Exponent, Polynomial]
) -> Polynomial:
    working = dict(polynomial)
    remainder: Polynomial = {}
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
                    working, exponent, -coefficient * value
                )
    return remainder


def relation_matrices(
    relations: list[dict[int, Fraction]],
    monomials: list[Exponent],
    adjugate: sp.Matrix,
) -> list[list[list[Polynomial]]]:
    adjugate_terms = [
        [
            sparse_polynomial(adjugate[row, column])
            for column in range(3)
        ]
        for row in range(3)
    ]
    matrices = []
    monomial_count = len(monomials)
    for relation in relations:
        direction = [{} for _ in range(3)]
        for index, coefficient in relation.items():
            add_term(
                direction[index // monomial_count],
                monomials[index % monomial_count],
                coefficient,
            )
        direction_jacobian = [
            [
                derivative(direction[row], column)
                for column in range(3)
            ]
            for row in range(3)
        ]
        matrix = [[{} for _ in range(3)] for _ in range(3)]
        for row in range(3):
            for column in range(3):
                for middle in range(3):
                    product = multiply(
                        adjugate_terms[row][middle],
                        direction_jacobian[middle][column],
                    )
                    for exponent, coefficient in product.items():
                        add_term(
                            matrix[row][column],
                            exponent,
                            coefficient,
                        )
        matrices.append(matrix)
    return matrices


def quadratic_class(
    first: list[list[Polynomial]],
    second: list[list[Polynomial]],
    diagonal: bool,
) -> Polynomial:
    output: Polynomial = {}
    for row in range(3):
        for column in range(3):
            product = multiply(
                first[row][column], second[column][row]
            )
            for exponent, coefficient in product.items():
                add_term(output, exponent, -coefficient)
    if diagonal:
        output = {
            exponent: coefficient / 2
            for exponent, coefficient in output.items()
        }
    return output


def quadratic_rank(
    matrices: list[list[list[Polynomial]]],
    linear_pivots: dict[Exponent, Polynomial],
) -> tuple[int, int, int, int]:
    quadratic_pivots: dict[Exponent, Polynomial] = {}
    nonzero_pairs = 0
    cokernel_coordinates: set[Exponent] = set()
    maximum_bit_size = 0
    for first in range(len(matrices)):
        for second in range(first, len(matrices)):
            column = full_remainder(
                quadratic_class(
                    matrices[first],
                    matrices[second],
                    first == second,
                ),
                linear_pivots,
            )
            nonzero_pairs += bool(column)
            cokernel_coordinates.update(column)
            while column:
                lead = leading_exponent(column)
                coefficient = column[lead]
                pivot = quadratic_pivots.get(lead)
                if pivot is None:
                    normalized = {
                        exponent: value / coefficient
                        for exponent, value in column.items()
                    }
                    quadratic_pivots[lead] = normalized
                    maximum_bit_size = max(
                        maximum_bit_size,
                        max(
                            max(
                                abs(value.numerator).bit_length(),
                                value.denominator.bit_length(),
                            )
                            for value in normalized.values()
                        ),
                    )
                    break
                for exponent, value in pivot.items():
                    add_term(
                        column,
                        exponent,
                        -coefficient * value,
                    )
    return (
        len(quadratic_pivots),
        nonzero_pairs,
        len(cokernel_coordinates),
        maximum_bit_size,
    )


def free_coordinates(
    direction: sp.Matrix,
    free_columns: list[int],
    monomials: list[Exponent],
) -> sp.Matrix:
    monomial_count = len(monomials)
    polynomials = [
        sp.Poly(sp.expand(direction[component]), *VARIABLES)
        for component in range(3)
    ]
    return sp.Matrix(
        [
            polynomials[column_index // monomial_count].coeff_monomial(
                monomials[column_index % monomial_count]
            )
            for column_index in free_columns
        ]
    )


def affine_directions(
    mapping: sp.Matrix, jacobian: sp.Matrix
) -> list[sp.Matrix]:
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
    for target_side in (True, False):
        for row in range(3):
            for column in range(3):
                if row == column:
                    continue
                matrix = sp.zeros(3)
                matrix[row, column] = 1
                directions.append(
                    matrix * mapping
                    if target_side
                    else jacobian * (matrix * source_vector)
                )
    for target_matrix, source_matrix in (
        (sp.diag(1, 0, -1), sp.zeros(3)),
        (sp.diag(0, 1, -1), sp.zeros(3)),
        (sp.diag(-1, 0, 0), sp.diag(1, 0, 0)),
        (sp.diag(-1, 0, 0), sp.diag(0, 1, 0)),
        (sp.diag(-1, 0, 0), sp.diag(0, 0, 1)),
    ):
        directions.append(
            target_matrix * mapping
            + jacobian * (source_matrix * source_vector)
        )
    return directions


def reduced_family_ranks(
    mapping: sp.Matrix,
    jacobian: sp.Matrix,
    free_columns: list[int],
    monomials: list[Exponent],
) -> tuple[int, int, int]:
    affine = sp.Matrix.hstack(
        *[
            free_coordinates(
                direction, free_columns, monomials
            )
            for direction in affine_directions(mapping, jacobian)
        ]
    )
    shears = []
    for target_component in (0, 1):
        for power in (2, 3):
            direction = sp.zeros(3, 1)
            direction[target_component] = mapping[2] ** power
            shears.append(
                free_coordinates(
                    direction, free_columns, monomials
                )
            )
    shear_matrix = sp.Matrix.hstack(*shears)
    parameter = sp.symbols("seed_parameter")
    family = mapping_from_primitive(
        explicit_seed(4)
        + parameter * w**2 * (w - 1) ** 2
    )
    seed_direction = sp.Matrix(
        [
            sp.diff(component, parameter).subs(parameter, 0)
            for component in family
        ]
    )
    seed_column = free_coordinates(
        seed_direction, free_columns, monomials
    )
    return (
        affine.rank(),
        affine.row_join(shear_matrix).rank(),
        affine.row_join(shear_matrix).row_join(seed_column).rank(),
    )


def main() -> None:
    mapping = mapping_from_primitive(explicit_seed(4))
    jacobian = mapping.jacobian(VARIABLES)
    assert sp.factor(jacobian.det()) == 1
    monomials, free_columns, relations, linear_pivots = tangent_echelon(
        jacobian.adjugate(), 12
    )
    assert len(linear_pivots) == 1307
    assert len(relations) == 58
    matrices = relation_matrices(
        relations, monomials, jacobian.adjugate()
    )
    rank, nonzero_pairs, coordinates, maximum_bits = quadratic_rank(
        matrices, linear_pivots
    )
    assert rank == 53
    assert nonzero_pairs == 727
    assert coordinates == 338
    assert maximum_bits == 214

    affine_rank, shear_rank, seed_rank = reduced_family_ranks(
        mapping, jacobian, free_columns, monomials
    )
    assert (affine_rank, shear_rank, seed_rank) == (22, 26, 27)

    print("PASS: dim_QQ ker(L_F4)=58 and rank_QQ(L_F4)=1307")
    print("PASS: the full quadratic Kuranishi map has exact QQ-rank 53")
    print(
        "PASS: 727 nonzero polarized pairs use 338 canonical "
        "cokernel coordinates"
    )
    print(
        "PASS: affine, target-shear, and seed family ranks are "
        "22, 26, and 27"
    )


if __name__ == "__main__":
    main()

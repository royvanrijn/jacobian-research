#!/usr/bin/env python3
"""Exact generic-family tangent profiles for F_4, F_5, and F_6.

For inverse degree N, put d=5N-8 and kmax=floor(d/4).  Starting from a
generic rational normalized seed, add C^k to each of the first two target
coordinates for 2<=k<=kmax.  Together with the affine left-right orbit these
give an explicit reduced family of tangent rank

    22 + 2(kmax-1) + (N-3).

The checker certifies that rank and computes the full X(3,d) tangent
dimension at one rational point of each family.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import w, x, y, z  # noqa: E402
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    exact_sparse_rank,
    explicit_seed,
    mapping_from_primitive,
    monomials_through,
)


VARIABLES = (x, y, z)
EXPECTED = {
    4: {
        "coefficient_degree": 12,
        "family_rank": 27,
        "tangent_dimension": 49,
    },
    5: {
        "coefficient_degree": 17,
        "family_rank": 30,
        "tangent_dimension": 80,
    },
    6: {
        "coefficient_degree": 22,
        "family_rank": 33,
        "tangent_dimension": 109,
    },
}


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


def coefficient_column_mod(
    direction: sp.Matrix,
    monomials: list[tuple[int, int, int]],
    prime: int,
) -> dict[int, int]:
    monomial_indices = {
        exponent: index for index, exponent in enumerate(monomials)
    }
    polynomials = [
        sp.Poly(sp.expand(direction[component]), *VARIABLES)
        for component in range(3)
    ]
    column: dict[int, int] = {}
    for component, polynomial in enumerate(polynomials):
        for exponent, coefficient in polynomial.terms():
            coefficient = sp.Rational(coefficient)
            value = (
                int(coefficient.p)
                * pow(int(coefficient.q) % prime, -1, prime)
                % prime
            )
            if value:
                column[
                    component * len(monomials)
                    + monomial_indices[exponent]
                ] = value
    return column


def modular_column_rank(
    columns: list[dict[int, int]], prime: int
) -> int:
    pivots: dict[int, dict[int, int]] = {}
    rank = 0
    for original in columns:
        column = dict(original)
        while column:
            lead = max(column)
            coefficient = column[lead]
            pivot = pivots.get(lead)
            if pivot is None:
                inverse = pow(coefficient, -1, prime)
                pivots[lead] = {
                    row: value * inverse % prime
                    for row, value in column.items()
                }
                rank += 1
                break
            for row, value in pivot.items():
                new_value = (
                    column.get(row, 0) - coefficient * value
                ) % prime
                if new_value:
                    column[row] = new_value
                else:
                    column.pop(row, None)
    return rank


def audit_degree(degree: int) -> None:
    primitive = explicit_seed(degree) + sum(
        (index + 1) * w ** (index + 2) * (w - 1) ** 2
        for index in range(degree - 3)
    )
    seed_mapping = sp.Matrix(mapping_from_primitive(primitive))
    coefficient_degree = 5 * degree - 8
    kmax = coefficient_degree // 4
    first, second, third = seed_mapping
    first_shear_coefficients = {
        index: index - 1 for index in range(2, kmax + 1)
    }
    second_shear_coefficients = {
        index: index + kmax for index in range(2, kmax + 1)
    }
    mapping = sp.Matrix(
        [
            first
            + sum(
                first_shear_coefficients[index] * third**index
                for index in range(2, kmax + 1)
            ),
            second
            + sum(
                second_shear_coefficients[index] * third**index
                for index in range(2, kmax + 1)
            ),
            third,
        ]
    ).applyfunc(sp.expand)
    jacobian = mapping.jacobian(VARIABLES)
    assert sp.factor(jacobian.det()) == 1
    assert max(
        sp.Poly(component, *VARIABLES).total_degree()
        for component in mapping
    ) == coefficient_degree

    parameter_directions: list[sp.Matrix] = []
    for target_component in (0, 1):
        for power in range(2, kmax + 1):
            direction = sp.zeros(3, 1)
            direction[target_component] = third**power
            parameter_directions.append(direction)
    seed_parameter = sp.symbols("seed_parameter")
    for index in range(degree - 3):
        varied_mapping = sp.Matrix(
            mapping_from_primitive(
                primitive
                + seed_parameter
                * w ** (index + 2)
                * (w - 1) ** 2
            )
        )
        raw_direction = varied_mapping.diff(
            seed_parameter
        ).subs(seed_parameter, 0)
        raw_first, raw_second, raw_third = raw_direction
        parameter_directions.append(
            sp.Matrix(
                [
                    raw_first
                    + sum(
                        coefficient
                        * power
                        * third ** (power - 1)
                        * raw_third
                        for power, coefficient in (
                            first_shear_coefficients.items()
                        )
                    ),
                    raw_second
                    + sum(
                        coefficient
                        * power
                        * third ** (power - 1)
                        * raw_third
                        for power, coefficient in (
                            second_shear_coefficients.items()
                        )
                    ),
                    raw_third,
                ]
            ).applyfunc(sp.expand)
        )
    family_directions = (
        affine_directions(mapping, jacobian) + parameter_directions
    )
    monomials = monomials_through(coefficient_degree)
    prime = 1_000_003
    family_rank = modular_column_rank(
        [
            coefficient_column_mod(
                direction, monomials, prime
            )
            for direction in family_directions
        ],
        prime,
    )
    expected_family_rank = (
        22 + 2 * (kmax - 1) + (degree - 3)
    )
    # There is one explicit kernel direction among the 23 affine parameters:
    # the weighted source vector field minus its target weight field.  Hence
    # the rational rank is at most parameter_count-1=expected_family_rank.
    # The modular rank is a lower bound for the rational rank.
    assert len(family_directions) == expected_family_rank + 1
    assert family_rank == expected_family_rank

    tangent_rank, _rows, _bits, _pivots = exact_sparse_rank(
        jacobian.adjugate(), coefficient_degree
    )
    coefficient_variables = 3 * len(monomials)
    tangent_dimension = coefficient_variables - tangent_rank
    expected = EXPECTED[degree]
    assert coefficient_degree == expected["coefficient_degree"]
    assert family_rank == expected["family_rank"]
    assert tangent_dimension == expected["tangent_dimension"]

    print(
        f"N={degree}: d={coefficient_degree}, kmax={kmax}, "
        f"reduced-family rank={family_rank}, "
        f"full tangent dimension={tangent_dimension}, "
        f"transverse tangent excess={tangent_dimension - family_rank}"
    )


def main() -> None:
    for degree in (4, 5, 6):
        audit_degree(degree)
    print(
        "PASS: generic reduced-family ranks are 27, 30, and 33"
    )
    print(
        "PASS: generic full tangent dimensions are 49, 80, and 109"
    )


if __name__ == "__main__":
    main()

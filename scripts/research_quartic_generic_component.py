#!/usr/bin/env python3
"""Normal quadratic slice at a generic point of the quartic reduced family.

The selected rational point has seed parameter 1 and target-shear parameters
(1,2,3,4).  The exact family consists of the affine left-right orbit, the
four target shears adding C^2,C^3 to the first two outputs, and the normalized
quartic seed parameter.  Its tangent rank is 27.

This research compiler works modulo a good prime.  It quotients those 27
directions from the full 49-dimensional coefficient tangent space and emits
the resulting 22-variable quadratic Kuranishi ideal.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import w
from research_quartic_coefficient_kuranishi import (
    VARIABLES,
    affine_directions,
    combine_matrices,
    emit_macaulay2,
    emit_singular,
    free_coordinates,
    independent_columns,
    independent_quadrics,
    quadratic_remainders,
    relation_matrices,
    row_reduce_vectors,
    sparse_polynomial,
    tangent_echelon,
)
from verify_all_degree_coefficient_tangents import (
    explicit_seed,
    mapping_from_primitive,
)


def generic_map_and_family_directions() -> tuple[
    sp.Matrix, list[sp.Matrix]
]:
    seed_parameter = sp.symbols("seed_parameter")
    seed_family = mapping_from_primitive(
        explicit_seed(4)
        + seed_parameter * w**2 * (w - 1) ** 2
    )
    first, second, third = seed_family
    shear_parameters = sp.symbols("a2 a3 b2 b3")
    a2, a3, b2, b3 = shear_parameters
    family = sp.Matrix(
        [
            first + a2 * third**2 + a3 * third**3,
            second + b2 * third**2 + b3 * third**3,
            third,
        ]
    )
    point_substitution = {
        seed_parameter: 1,
        a2: 1,
        a3: 2,
        b2: 3,
        b3: 4,
    }
    mapping = family.subs(point_substitution).applyfunc(sp.expand)
    directions = [
        family.diff(parameter).subs(point_substitution).applyfunc(sp.expand)
        for parameter in (*shear_parameters, seed_parameter)
    ]
    return mapping, directions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--singular-output", type=Path)
    parser.add_argument("--macaulay2-output", type=Path)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()
    prime = args.prime
    if not sp.isprime(prime) or prime in (2, 3, 5, 17):
        parser.error("choose a good odd prime away from 3, 5, and 17")

    mapping, nonlinear_family_directions = (
        generic_map_and_family_directions()
    )
    jacobian = mapping.jacobian(VARIABLES)
    assert sp.factor(jacobian.det()) == 1
    assert tuple(
        sp.Poly(component, *VARIABLES).total_degree()
        for component in mapping
    ) == (12, 12, 4)

    (
        monomials,
        free_columns,
        relations,
        linear_pivots,
        _pivot_combinations,
        adjugate_terms,
    ) = tangent_echelon(jacobian.adjugate(), 12, prime)
    assert len(linear_pivots) == 1316
    assert len(relations) == 49
    tangent_matrices, _tangent_jacobians = relation_matrices(
        relations, monomials, adjugate_terms, prime
    )

    family_directions = (
        affine_directions(mapping, jacobian)
        + nonlinear_family_directions
    )
    family_columns = [
        free_coordinates(
            direction, free_columns, monomials, prime
        )
        for direction in family_directions
    ]
    independent_family, independent_parameter_indices = (
        independent_columns(family_columns, prime)
    )
    assert len(independent_family) == 27
    _rows, family_pivot_rows = row_reduce_vectors(
        [
            {
                column_index: column.get(row_index, 0)
                for column_index, column in enumerate(independent_family)
                if column.get(row_index, 0)
            }
            for row_index in range(49)
        ],
        prime,
    )
    assert len(family_pivot_rows) == 27
    normal_indices = [
        index for index in range(49) if index not in family_pivot_rows
    ]
    assert len(normal_indices) == 22
    normal_matrices = [
        combine_matrices(
            {index: 1}, tangent_matrices, prime
        )
        for index in normal_indices
    ]
    pairs, remainders = quadratic_remainders(
        normal_matrices, linear_pivots, prime
    )
    equations = independent_quadrics(pairs, remainders, prime)
    assert len(equations) == 22
    variable_names = [f"n{index}" for index in range(22)]
    square_pair_indices = {
        pairs.index((index, index)): index for index in range(22)
    }
    square_rows = [
        {
            square_pair_indices[pair_index]: coefficient
            for pair_index, coefficient in equation.items()
            if pair_index in square_pair_indices
        }
        for equation in equations
    ]
    square_echelon, _square_indices = row_reduce_vectors(
        square_rows, prime
    )
    square_coefficient_rank = len(square_echelon)

    summary = {
        "prime": prime,
        "selected_seed_parameter": 1,
        "selected_target_shear_parameters": [1, 2, 3, 4],
        "full_tangent_dimension": 49,
        "reduced_family_tangent_dimension": 27,
        "normal_tangent_dimension": 22,
        "quadratic_kuranishi_rank": 22,
        "quadratic_generator_count": len(equations),
        "quadratic_pair_count": len(pairs),
        "square_coefficient_rank": square_coefficient_rank,
        "nonzero_quadratic_pair_count": sum(
            bool(remainder) for remainder in remainders
        ),
        "cokernel_monomial_count": len(
            {
                exponent
                for remainder in remainders
                for exponent in remainder
            }
        ),
        "independent_family_parameter_indices": (
            independent_parameter_indices
        ),
        "family_pivot_tangent_indices": family_pivot_rows,
        "normal_tangent_indices": normal_indices,
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

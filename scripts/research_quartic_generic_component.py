#!/usr/bin/env python3
"""Normal quadratic slice at a generic point of the quartic reduced family.

The selected rational point has seed parameter 1 and target-shear parameters
(1,2,3,4).  The exact family consists of the affine left-right orbit, the
four target shears adding C^2,C^3 to the first two outputs, and the normalized
quartic seed parameter.  Its tangent rank is 27.

This research compiler works modulo a good prime.  It quotients those 27
directions from the full 49-dimensional coefficient tangent space and emits
the resulting 22-variable quadratic Kuranishi ideal.

Passing --singular-order3-output also compiles the complete canonical cubic
homogeneous layer.  That slower output is a research input, not a completed
local primary decomposition.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import w
from research_quartic_coefficient_kuranishi import (
    VARIABLES,
    add_scalar,
    add_term,
    coefficient_direction_jacobian,
    affine_directions,
    combine_matrices,
    cubic_axis_screen,
    emit_macaulay2,
    emit_singular,
    free_coordinates,
    independent_columns,
    independent_quadrics,
    jet_lift_axis,
    left_multiply_adjugate,
    multiply,
    quadratic_remainders,
    reduce_polynomial,
    relation_matrices,
    row_reduce_vectors,
    singular_polynomial,
    solve_linearized_image,
    sparse_polynomial,
    tangent_echelon,
    trace_product,
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


def quadratic_correction_matrices(
    normal_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    pairs: list[tuple[int, int]],
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
) -> list[list[list[dict[tuple[int, int, int], int]]]]:
    """Canonical order-two corrections, retaining the cokernel equations."""
    corrections = []
    inverse_two = pow(2, -1, prime)
    for first, second in pairs:
        polynomial = trace_product(
            normal_matrices[first],
            normal_matrices[second],
            prime,
            scalar=-1,
        )
        if first == second:
            polynomial = {
                exponent: coefficient * inverse_two % prime
                for exponent, coefficient in polynomial.items()
            }
        solution, _remainder = solve_linearized_image(
            polynomial,
            linear_pivots,
            pivot_combinations,
            prime,
        )
        correction_jacobian = coefficient_direction_jacobian(
            {
                index: -coefficient % prime
                for index, coefficient in solution.items()
            },
            monomials,
            prime,
        )
        corrections.append(
            left_multiply_adjugate(
                adjugate_terms, correction_jacobian, prime
            )
        )
    return corrections


def mixed_determinant_coefficient(
    triple: tuple[int, int, int],
    normal_jacobians: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    prime: int,
) -> dict[tuple[int, int, int], int]:
    """Coefficient of n_i*n_j*n_k in det(sum n_i D G_i)."""
    column_permutations = (
        ((0, 1, 2), 1),
        ((0, 2, 1), -1),
        ((1, 0, 2), -1),
        ((1, 2, 0), 1),
        ((2, 0, 1), 1),
        ((2, 1, 0), -1),
    )
    output: dict[tuple[int, int, int], int] = {}
    for row_assignment in set(itertools.permutations(triple)):
        for column_assignment, sign in column_permutations:
            product = multiply(
                multiply(
                    normal_jacobians[row_assignment[0]][0][
                        column_assignment[0]
                    ],
                    normal_jacobians[row_assignment[1]][1][
                        column_assignment[1]
                    ],
                    prime,
                ),
                normal_jacobians[row_assignment[2]][2][
                    column_assignment[2]
                ],
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


def cubic_kuranishi_remainders(
    normal_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    normal_jacobians: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    pairs: list[tuple[int, int]],
    correction_matrices: list[
        list[list[dict[tuple[int, int, int], int]]]
    ],
    linear_pivots: dict[
        tuple[int, int, int], dict[tuple[int, int, int], int]
    ],
    prime: int,
) -> tuple[
    list[tuple[int, int, int]],
    list[dict[tuple[int, int, int], int]],
]:
    """Canonical order-three cokernel equations on the normal variables."""
    pair_indices = {
        pair: index for index, pair in enumerate(pairs)
    }
    triples = list(
        itertools.combinations_with_replacement(
            range(len(normal_matrices)), 3
        )
    )
    remainders = []
    for triple in triples:
        polynomial = mixed_determinant_coefficient(
            triple, normal_jacobians, prime
        )
        for first in range(len(normal_matrices)):
            remainder_pair = list(triple)
            try:
                remainder_pair.remove(first)
            except ValueError:
                continue
            pair = tuple(sorted(remainder_pair))
            pair_index = pair_indices[pair]
            cross_term = trace_product(
                normal_matrices[first],
                correction_matrices[pair_index],
                prime,
                scalar=-1,
            )
            for exponent, coefficient in cross_term.items():
                add_term(
                    polynomial, exponent, coefficient, prime
                )
        remainders.append(
            reduce_polynomial(
                polynomial, linear_pivots, prime
            )
        )
    return triples, remainders


def independent_homogeneous_equations(
    remainders: list[dict[tuple[int, int, int], int]],
    prime: int,
) -> list[dict[int, int]]:
    equations_by_cokernel_monomial: dict[
        tuple[int, int, int], dict[int, int]
    ] = {}
    for monomial_index, remainder in enumerate(remainders):
        for exponent, coefficient in remainder.items():
            add_scalar(
                equations_by_cokernel_monomial.setdefault(
                    exponent, {}
                ),
                monomial_index,
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


def singular_homogeneous_polynomial(
    equation: dict[int, int],
    monomials: list[tuple[int, ...]],
    variable_names: list[str],
    prime: int,
) -> str:
    terms = []
    for monomial_index, coefficient in sorted(equation.items()):
        symmetric = (
            coefficient
            if coefficient <= prime // 2
            else coefficient - prime
        )
        factors = [
            variable_names[index] for index in monomials[monomial_index]
        ]
        terms.append(f"({symmetric})*" + "*".join(factors))
    return "+".join(terms) if terms else "0"


def emit_singular_through_cubic(
    path: Path,
    prime: int,
    variable_names: list[str],
    pairs: list[tuple[int, int]],
    quadratic_equations: list[dict[int, int]],
    triples: list[tuple[int, int, int]],
    cubic_equations: list[dict[int, int]],
) -> None:
    generators = [
        singular_polynomial(
            equation, pairs, variable_names, prime
        )
        for equation in quadratic_equations
    ] + [
        singular_homogeneous_polynomial(
            equation, triples, variable_names, prime
        )
        for equation in cubic_equations
    ]
    text = "\n".join(
        [
            f"ring r={prime},({','.join(variable_names)}),dp;",
            "ideal I=",
            ",\n".join(generators) + ";",
            'print("KURANISHI_ORDER23_GENERATORS");',
            "print(size(I));",
            "ideal G=slimgb(I);",
            'print("KURANISHI_ORDER23_STANDARD_BASIS");',
            "print(size(G));",
            'print("KURANISHI_ORDER23_DIMENSION");',
            "print(dim(G));",
            'print("KURANISHI_ORDER23_DEGREE");',
            "print(deg(G));",
            'print("KURANISHI_ORDER23_HILBERT");',
            "print(hilb(G));",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--singular-output", type=Path)
    parser.add_argument("--singular-order3-output", type=Path)
    parser.add_argument("--macaulay2-output", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument(
        "--greedy-jet-order",
        type=int,
        default=6,
        help=(
            "largest bounded-box order in the noncanonical coordinate-axis "
            "lookahead"
        ),
    )
    args = parser.parse_args()
    prime = args.prime
    if not sp.isprime(prime) or prime in (2, 3, 5, 17):
        parser.error("choose a good odd prime away from 3, 5, and 17")
    if args.greedy_jet_order < 2:
        parser.error("--greedy-jet-order must be at least two")

    mapping, nonlinear_family_directions = (
        generic_map_and_family_directions()
    )
    jacobian = mapping.jacobian(VARIABLES)
    jacobian_terms = [
        [
            sparse_polynomial(jacobian[row, column], prime)
            for column in range(3)
        ]
        for row in range(3)
    ]
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
        pivot_combinations,
        adjugate_terms,
    ) = tangent_echelon(jacobian.adjugate(), 12, prime)
    assert len(linear_pivots) == 1316
    assert len(relations) == 49
    tangent_matrices, tangent_jacobians = relation_matrices(
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
    normal_jacobians = [
        combine_matrices(
            {index: 1}, tangent_jacobians, prime
        )
        for index in normal_indices
    ]
    pairs, remainders = quadratic_remainders(
        normal_matrices, linear_pivots, prime
    )
    equations = independent_quadrics(pairs, remainders, prime)
    assert len(equations) == 22
    triples: list[tuple[int, int, int]] = []
    cubic_remainders: list[
        dict[tuple[int, int, int], int]
    ] = []
    cubic_equations: list[dict[int, int]] = []
    if args.singular_order3_output is not None:
        correction_matrices = quadratic_correction_matrices(
            normal_matrices,
            pairs,
            monomials,
            linear_pivots,
            pivot_combinations,
            adjugate_terms,
            prime,
        )
        triples, cubic_remainders = cubic_kuranishi_remainders(
            normal_matrices,
            normal_jacobians,
            pairs,
            correction_matrices,
            linear_pivots,
            prime,
        )
        cubic_equations = independent_homogeneous_equations(
            cubic_remainders, prime
        )
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
    cubic_screen = cubic_axis_screen(
        variable_names,
        normal_matrices,
        normal_jacobians,
        tangent_matrices,
        monomials,
        linear_pivots,
        pivot_combinations,
        adjugate_terms,
        prime,
    )
    if triples:
        for entry in cubic_screen:
            axis_index = variable_names.index(str(entry["axis"]))
            cubic_index = triples.index(
                (axis_index, axis_index, axis_index)
            )
            assert len(cubic_remainders[cubic_index]) == entry[
                "cubic_remainder_term_count"
            ]
    unobstructed_axis_indices = [
        index
        for index in range(len(variable_names))
        if not remainders[pairs.index((index, index))]
    ]
    greedy_jet_screen = [
        jet_lift_axis(
            variable_names[index],
            normal_jacobians[index],
            normal_matrices[index],
            jacobian_terms,
            tangent_jacobians,
            tangent_matrices,
            monomials,
            linear_pivots,
            pivot_combinations,
            prime,
            args.greedy_jet_order,
        )
        for index in unobstructed_axis_indices
    ]

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
        "quadratically_unobstructed_coordinate_axes": [
            entry["axis"] for entry in cubic_screen
        ],
        "coordinate_axis_cubic_screen": cubic_screen,
        "greedy_coordinate_axis_jet_screen": greedy_jet_screen,
        "greedy_coordinate_axis_maximum_order": (
            args.greedy_jet_order
        ),
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
    if triples:
        summary.update(
            {
                "cubic_monomial_count": len(triples),
                "nonzero_cubic_monomial_count": sum(
                    bool(remainder)
                    for remainder in cubic_remainders
                ),
                "cubic_equation_rank": len(
                    cubic_equations
                ),
                "cubic_cokernel_monomial_count": len(
                    {
                        exponent
                        for remainder in cubic_remainders
                        for exponent in remainder
                    }
                ),
            }
        )
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.singular_output is not None:
        emit_singular(
            args.singular_output,
            prime,
            variable_names,
            pairs,
            equations,
        )
    if args.singular_order3_output is not None:
        emit_singular_through_cubic(
            args.singular_order3_output,
            prime,
            variable_names,
            pairs,
            equations,
            triples,
            cubic_equations,
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

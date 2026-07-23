#!/usr/bin/env python3
"""First exact regressions for the Hessian--Ritt deformation complex.

This checker separates three levels which the deformation programme must not
conflate:

* the tree-local composition differential;
* the linear cotangent complex at the reduced augmentation;
* the completed Artin algebra, which retains higher nilpotence.

It also checks the ordinary cellular cohomology of the filled degree-thirty
braid as a combinatorial baseline.  That last calculation is deliberately
not treated as a substitute for the coefficient-decorated complex.
"""

from __future__ import annotations

import sys
from collections import Counter
from itertools import product
from math import prod
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.ritt_complex import (  # noqa: E402
    MoveType,
    compose_factors,
    composition_differential,
    degree_thirty_braid_decorations,
    symmetric_braid_complex,
)


X = sp.symbols("x")


def audit_composition_differential() -> None:
    """Compare the tree formula with literal dual-number differentiation."""

    factor_degrees = (2, 3, 5)
    parameter_blocks = [
        sp.symbols(f"a{position}_1:{degree}")
        for position, degree in enumerate(factor_degrees)
    ]
    variation_blocks = [
        sp.symbols(f"da{position}_1:{degree}")
        for position, degree in enumerate(factor_degrees)
    ]
    factors = tuple(
        X**degree
        + sum(
            parameter * X**power
            for power, parameter in enumerate(parameters, 1)
        )
        for degree, parameters in zip(factor_degrees, parameter_blocks)
    )
    variations = tuple(
        sum(
            parameter * X**power
            for power, parameter in enumerate(parameters, 1)
        )
        for parameters in variation_blocks
    )

    tree_derivative = composition_differential(factors, variations, X)
    epsilon = sp.symbols("epsilon")
    dual_number_composition = compose_factors(
        tuple(
            factor + epsilon * variation
            for factor, variation in zip(factors, variations)
        ),
        X,
    )
    literal_derivative = sp.diff(dual_number_composition, epsilon).subs(
        epsilon, 0
    )
    assert sp.expand(tree_derivative - literal_derivative) == 0

    variation_parameters = tuple(
        parameter
        for block in variation_blocks
        for parameter in block
    )
    columns = [
        [
            sp.Poly(
                sp.diff(tree_derivative, parameter), X
            ).nth(coefficient_degree)
            for coefficient_degree in range(1, prod(factor_degrees))
        ]
        for parameter in variation_parameters
    ]
    tangent_matrix = sp.Matrix.hstack(
        *(sp.Matrix(column) for column in columns)
    )
    specialization = {
        parameter: value
        for block in parameter_blocks
        for parameter, value in zip(
            block,
            range(1, len(block) + 1),
        )
    }
    assert tangent_matrix.subs(specialization).rank() == sum(
        degree - 1 for degree in factor_degrees
    )


def cellular_coboundaries(complex_) -> tuple[sp.Matrix, sp.Matrix]:
    """Return ``C^0 -> C^1 -> C^2`` for the unlabelled filled braid."""

    vertices = tuple(sorted(complex_.words))
    vertex_index = {word: index for index, word in enumerate(vertices)}
    edges = tuple(complex_.edges)
    edge_index = {
        frozenset(edge.endpoints): index for index, edge in enumerate(edges)
    }

    delta_zero = sp.zeros(len(edges), len(vertices))
    for row, edge in enumerate(edges):
        left, right = edge.endpoints
        delta_zero[row, vertex_index[left]] = -1
        delta_zero[row, vertex_index[right]] = 1

    delta_one = sp.zeros(len(complex_.two_cells), len(edges))
    for row, cell in enumerate(complex_.two_cells):
        first, second = cell.paths
        oriented_boundary = list(zip(first, first[1:])) + list(
            zip(reversed(second), reversed(second[:-1]))
        )
        for left, right in oriented_boundary:
            column = edge_index[frozenset((left, right))]
            stored_left, stored_right = edges[column].endpoints
            delta_one[row, column] += (
                1 if (left, right) == (stored_left, stored_right) else -1
            )
    return delta_zero, delta_one


def audit_filled_braid_baseline() -> None:
    """Check contractibility before coefficient decorations are attached."""

    complex_ = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    delta_zero, delta_one = cellular_coboundaries(complex_)
    assert delta_one * delta_zero == sp.zeros(1, 6)
    ranks = (delta_zero.rank(), delta_one.rank())
    assert ranks == (5, 1)
    cohomology_dimensions = (
        6 - ranks[0],
        6 - ranks[0] - ranks[1],
        1 - ranks[1],
    )
    assert cohomology_dimensions == (1, 0, 0)


def hilbert_vector(exponents: tuple[int, ...]) -> tuple[int, ...]:
    """Count the standard monomials of a monomial complete intersection."""

    degrees = Counter(
        sum(monomial)
        for monomial in product(
            *(range(exponent) for exponent in exponents)
        )
    )
    return tuple(degrees[degree] for degree in range(max(degrees) + 1))


def audit_transverse_cotangent_complexes() -> None:
    """Recover linear cotangent homology and distinguish higher thickenings."""

    decorations = degree_thirty_braid_decorations()
    linear_signatures = []
    completed_signatures = []
    for decoration in decorations:
        exponents = decoration.transverse_slice.exponents
        variables = sp.symbols(f"u0:{len(exponents)}")
        relations = [
            variable**exponent
            for variable, exponent in zip(variables, exponents)
        ]
        jacobian = sp.Matrix(
            [
                [sp.diff(relation, variable) for variable in variables]
                for relation in relations
            ]
        )
        augmentation = {variable: 0 for variable in variables}
        point_differential = jacobian.subs(augmentation)
        assert point_differential == sp.zeros(
            len(relations), len(variables)
        )

        h_zero = len(variables) - point_differential.rank()
        h_one = len(relations) - point_differential.rank()
        linear_signatures.append((h_zero, h_one))

        vector = hilbert_vector(exponents)
        assert vector == decoration.transverse_slice.hilbert_vector
        assert sum(vector) == decoration.transverse_slice.length
        maximal_ideal_nilpotence = len(vector)
        completed_signatures.append(
            (vector, maximal_ideal_nilpotence)
        )

    assert tuple(linear_signatures) == ((1, 1), (2, 2), (2, 2))
    assert tuple(completed_signatures) == (
        ((1, 1, 1, 1, 1), 5),
        ((1, 2, 1), 3),
        ((1, 2, 2, 2, 1), 5),
    )

    # The last two sectors have identical point-cotangent complexes but
    # different completed algebras.  This is the exact reason the linear
    # complex alone cannot encode the braid thickening.
    assert linear_signatures[1] == linear_signatures[2]
    assert completed_signatures[1] != completed_signatures[2]

    # The nilradical of the path scheme and the maximal ideal of its
    # conductor slice are different filtrations; their indices need not
    # agree in the first and third sectors.
    assert tuple(
        decoration.nilpotence_index for decoration in decorations
    ) == (4, 3, 4)
    assert tuple(
        signature[1] for signature in completed_signatures
    ) == (5, 3, 5)


def main() -> None:
    audit_composition_differential()
    audit_filled_braid_baseline()
    audit_transverse_cotangent_complexes()
    print("PASS: the tree-local composition differential equals dual-number differentiation")
    print("PASS: the normalized 2 o 3 o 5 factor tangent map has rank seven")
    print("PASS: the unlabelled filled braid has cohomology dimensions (1,0,0)")
    print("PASS: the three sector cotangent homology ranks are (1,1),(2,2),(2,2)")
    print("PASS: completed Hilbert data distinguish sectors with identical linear complexes")
    print("PASS: path-nilradical and conductor-slice filtrations remain distinct")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Bounded sparse-cubic nonlinear shear search in the two exceptional charts.

The complete linear calculation leaves only singular families in charts 1110
and 1111, represented by the quadratic potential ell*u_4^2/2.  This search
adds one, two, or three cubic monomials:

    V = ell*u_4^2/2 + c_1*m_1 + c_2*m_2,

For supports of size one or two the cubic coefficients lie in
{-2,-1,1,2}; for size three they lie in {-1,1}.  In every case ell lies in
[-2,2].  All 20 degree-three monomials are allowed.  A modular multi-point
determinant filter is followed by exact symbolic verification of every
survivor.
"""

from __future__ import annotations

import argparse
from itertools import combinations, product
import runpy

import sympy as sp


graph_namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
utility_namespace = runpy.run_path("scripts/search_hc4_lagrangian_shears.py")

h_variables = graph_namespace["h_variables"]
Q = list(graph_namespace["position_coordinates_h"])
M = list(graph_namespace["momentum_coordinates_h"])
rational_mod = utility_namespace["rational_mod"]
determinant_mod = utility_namespace["determinant_mod"]

u = sp.symbols("u0:4")
cubic_monomials = [
    (exponents, sp.prod(u[index] ** exponents[index] for index in range(4)))
    for exponents in product(range(4), repeat=4)
    if sum(exponents) == 3
]
cubic_hessians = [
    sp.hessian(monomial, u) for _, monomial in cubic_monomials
]


def matrix_mod(
    matrix: sp.Matrix, substitution: dict[sp.Symbol, int], prime: int
) -> list[list[int]]:
    return [
        [
            rational_mod(matrix[row, column].subs(substitution), prime)
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ]


def multiply_mod(
    left: list[list[int]], right: list[list[int]], prime: int
) -> list[list[int]]:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(4))
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1_000_003)
    args = parser.parse_args()

    sample_points = (
        (0, 0, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
        (1, 1, 1, 1),
        (-1, 2, -2, 1),
        (2, -1, 1, -2),
        (3, 1, -1, 2),
    )
    coefficient_values = (-2, -1, 1, 2)
    unit_coefficient_values = (-1, 1)
    ell_values = range(-2, 3)
    support_choices = [
        support
        for support_size in (1, 2, 3)
        for support in combinations(range(len(cubic_monomials)), support_size)
    ]
    candidates_per_chart = sum(
        (
            len(coefficient_values)
            if len(support) <= 2
            else len(unit_coefficient_values)
        )
        ** len(support)
        * len(tuple(ell_values))
        for support in support_choices
    )

    modular_survivors = []
    exact_survivors = []

    for mask in ((1, 1, 1, 0), (1, 1, 1, 1)):
        q0 = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
        m0 = sp.Matrix(
            [-Q[index] if mask[index] else M[index] for index in range(4)]
        )
        jq0 = q0.jacobian(h_variables)
        jm0 = m0.jacobian(h_variables)
        evaluated_points = []

        for point in sample_points:
            substitution = dict(zip(h_variables, point, strict=True))
            q0_value = [
                q0[index].subs(substitution) for index in range(4)
            ]
            u_substitution = dict(zip(u, q0_value, strict=True))
            jq_value = matrix_mod(jq0, substitution, args.prime)
            jm_value = matrix_mod(jm0, substitution, args.prime)

            base_hessian = [
                [1 if row == column == 3 else 0 for column in range(4)]
                for row in range(4)
            ]
            base_contribution = multiply_mod(
                base_hessian, jq_value, args.prime
            )
            cubic_contributions = []
            for hessian in cubic_hessians:
                hessian_value = [
                    [
                        rational_mod(
                            hessian[row, column].subs(u_substitution),
                            args.prime,
                        )
                        for column in range(4)
                    ]
                    for row in range(4)
                ]
                cubic_contributions.append(
                    multiply_mod(hessian_value, jq_value, args.prime)
                )
            evaluated_points.append(
                (jm_value, base_contribution, cubic_contributions)
            )

        chart_survivors = []
        for support in support_choices:
            support_coefficient_values = (
                coefficient_values
                if len(support) <= 2
                else unit_coefficient_values
            )
            for coefficients in product(
                support_coefficient_values, repeat=len(support)
            ):
                for ell_value in ell_values:
                    expected_determinant = None
                    passes = True
                    for (
                        jm_value,
                        base_contribution,
                        cubic_contributions,
                    ) in evaluated_points:
                        candidate_matrix = [
                            [
                                (
                                    jm_value[row][column]
                                    + ell_value * base_contribution[row][column]
                                    + sum(
                                        coefficient
                                        * cubic_contributions[monomial_index][row][
                                            column
                                        ]
                                        for monomial_index, coefficient in zip(
                                            support, coefficients, strict=True
                                        )
                                    )
                                )
                                % args.prime
                                for column in range(4)
                            ]
                            for row in range(4)
                        ]
                        determinant = determinant_mod(
                            candidate_matrix, args.prime
                        )
                        if expected_determinant is None:
                            expected_determinant = determinant
                            if determinant == 0:
                                passes = False
                                break
                        elif determinant != expected_determinant:
                            passes = False
                            break
                    if passes:
                        chart_survivors.append(
                            (support, coefficients, ell_value)
                        )

        modular_survivors.extend((mask, *candidate) for candidate in chart_survivors)

        for support, coefficients, ell_value in chart_survivors:
            potential = sp.Rational(ell_value, 2) * u[3] ** 2
            potential += sum(
                coefficient * cubic_monomials[monomial_index][1]
                for monomial_index, coefficient in zip(
                    support, coefficients, strict=True
                )
            )
            gradient = sp.Matrix(
                [sp.diff(potential, variable) for variable in u]
            ).subs(dict(zip(u, q0, strict=True)), simultaneous=True)
            determinant = sp.factor(
                (m0 + gradient)
                .jacobian(h_variables)
                .det(method="berkowitz")
            )
            if determinant != 0 and not determinant.free_symbols:
                exact_survivors.append(
                    (mask, support, coefficients, ell_value, determinant)
                )

        label = "".join(map(str, mask))
        print(
            f"chart {label}: {candidates_per_chart} candidates, "
            f"{len(chart_survivors)} modular survivors"
        )

    print(f"total modular survivors: {len(modular_survivors)}")
    print(f"exact constant-nonzero survivors: {len(exact_survivors)}")
    for survivor in exact_survivors:
        print(survivor)


if __name__ == "__main__":
    main()

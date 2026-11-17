#!/usr/bin/env python3
"""Deterministic dense random search for nonlinear graph shears in all charts.

In each of the 16 Lagrangian charts, sample potentials V=V_2+V_3 with all
ten quadratic Hessian coefficients and all twenty homogeneous cubic
coefficients independently chosen from [-2,2].  The random seed is fixed by
default.  Candidates are filtered at six graph points modulo a large prime;
every survivor is reconstructed and verified symbolically over Q.

This is exploratory evidence, not an exhaustive coefficient classification.
"""

from __future__ import annotations

import argparse
from itertools import product
import random
import runpy

import sympy as sp


graph_namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
utility_namespace = runpy.run_path("scripts/search_hc4_lagrangian_shears.py")

h_variables = graph_namespace["h_variables"]
Q = list(graph_namespace["position_coordinates_h"])
M = list(graph_namespace["momentum_coordinates_h"])
rational_mod = utility_namespace["rational_mod"]
determinant_mod = utility_namespace["determinant_mod"]

prime = 1_000_003
u = sp.symbols("u0:4")
quadratic_entries = [
    (row, column) for row in range(4) for column in range(row, 4)
]
quadratic_hessians = []
quadratic_potentials = []
for row, column in quadratic_entries:
    hessian = sp.zeros(4)
    hessian[row, column] = 1
    hessian[column, row] = 1
    quadratic_hessians.append(hessian)
    quadratic_potentials.append(
        u[row] ** 2 / 2 if row == column else u[row] * u[column]
    )

cubic_monomials = [
    sp.prod(u[index] ** exponents[index] for index in range(4))
    for exponents in product(range(4), repeat=4)
    if sum(exponents) == 3
]
cubic_hessians = [sp.hessian(monomial, u) for monomial in cubic_monomials]


def rational_matrix_mod(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [
            rational_mod(matrix[row, column], prime)
            for column in range(matrix.cols)
        ]
        for row in range(matrix.rows)
    ]


def multiply_mod(
    left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(4))
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]


def linear_combination(
    base: list[list[int]],
    contributions: list[list[list[int]]],
    coefficients: tuple[int, ...],
) -> list[list[int]]:
    return [
        [
            (
                base[row][column]
                + sum(
                    coefficient * contribution[row][column]
                    for coefficient, contribution in zip(
                        coefficients, contributions, strict=True
                    )
                    if coefficient
                )
            )
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-chart", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=260722198)
    args = parser.parse_args()
    generator = random.Random(args.seed)

    sample_points = (
        (0, 0, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (1, 1, 1, 1),
        (-1, 2, -2, 1),
    )
    all_modular_survivors = []
    exact_survivors = []

    for mask in product((0, 1), repeat=4):
        q0 = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
        m0 = sp.Matrix(
            [-Q[index] if mask[index] else M[index] for index in range(4)]
        )
        jq0 = q0.jacobian(h_variables)
        jm0 = m0.jacobian(h_variables)
        point_data = []
        for point in sample_points:
            substitution = dict(zip(h_variables, point, strict=True))
            q0_value = q0.subs(substitution)
            u_substitution = dict(zip(u, list(q0_value), strict=True))
            jq_value = rational_matrix_mod(jq0.subs(substitution))
            jm_value = rational_matrix_mod(jm0.subs(substitution))
            quadratic_contributions = [
                multiply_mod(rational_matrix_mod(hessian), jq_value)
                for hessian in quadratic_hessians
            ]
            cubic_contributions = [
                multiply_mod(
                    rational_matrix_mod(hessian.subs(u_substitution)),
                    jq_value,
                )
                for hessian in cubic_hessians
            ]
            point_data.append(
                (jm_value, quadratic_contributions, cubic_contributions)
            )

        chart_survivors = []
        for _ in range(args.samples_per_chart):
            quadratic_coefficients = tuple(
                generator.randint(-2, 2) for _ in quadratic_hessians
            )
            cubic_coefficients = tuple(
                generator.randint(-2, 2) for _ in cubic_hessians
            )
            if not any(cubic_coefficients):
                continue

            expected = None
            passes = True
            all_coefficients = quadratic_coefficients + cubic_coefficients
            for (
                jm_value,
                quadratic_contributions,
                cubic_contributions,
            ) in point_data:
                matrix = linear_combination(
                    jm_value,
                    quadratic_contributions + cubic_contributions,
                    all_coefficients,
                )
                determinant = determinant_mod(matrix, prime)
                if expected is None:
                    expected = determinant
                    if determinant == 0:
                        passes = False
                        break
                elif determinant != expected:
                    passes = False
                    break
            if passes:
                chart_survivors.append(
                    (quadratic_coefficients, cubic_coefficients)
                )

        all_modular_survivors.extend(
            (mask, *survivor) for survivor in chart_survivors
        )

        for quadratic_coefficients, cubic_coefficients in chart_survivors:
            potential = sum(
                coefficient * monomial
                for coefficient, monomial in zip(
                    quadratic_coefficients,
                    quadratic_potentials,
                    strict=True,
                )
            )
            potential += sum(
                coefficient * monomial
                for coefficient, monomial in zip(
                    cubic_coefficients, cubic_monomials, strict=True
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
                    (
                        mask,
                        quadratic_coefficients,
                        cubic_coefficients,
                        determinant,
                    )
                )

        print(
            f"chart {''.join(map(str, mask))}: "
            f"{args.samples_per_chart} samples, "
            f"{len(chart_survivors)} modular survivors"
        )

    print(
        f"total samples: {16 * args.samples_per_chart}, "
        f"seed: {args.seed}"
    )
    print(f"total modular survivors: {len(all_modular_survivors)}")
    print(f"exact constant-nonzero survivors: {len(exact_survivors)}")
    for survivor in exact_survivors:
        print(survivor)


if __name__ == "__main__":
    main()

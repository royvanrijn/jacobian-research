#!/usr/bin/env python3
"""Search genuinely nonlinear sparse shears in the exceptional graph charts.

For charts 1110 and 1111, search potentials

    V(u) = V_2(u) + V_3(u),

where the symmetric Hessian of V_2 has at most two nonzero entries and V_3
has one or two cubic monomials.  Every nonzero coefficient is +1 or -1.
Unlike a search around the singular residual linear family, this permits a
nondegenerate determinant at the origin.

A modular nine-point filter is followed by exact symbolic verification of
every survivor.
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
    if row == column:
        potential = u[row] ** 2 / 2
    else:
        potential = u[row] * u[column]
    quadratic_hessians.append(hessian)
    quadratic_potentials.append(potential)

cubic_monomials = [
    (exponents, sp.prod(u[index] ** exponents[index] for index in range(4)))
    for exponents in product(range(4), repeat=4)
    if sum(exponents) == 3
]
cubic_hessians = [
    sp.hessian(monomial, u) for _, monomial in cubic_monomials
]


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


def add_contributions(
    base: list[list[int]],
    indexed_contributions: list[list[list[int]]],
    support: tuple[int, ...],
    coefficients: tuple[int, ...],
) -> list[list[int]]:
    return [
        [
            (
                base[row][column]
                + sum(
                    coefficient
                    * indexed_contributions[index][row][column]
                    for index, coefficient in zip(
                        support, coefficients, strict=True
                    )
                )
            )
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quadratic-support", type=int, default=2)
    parser.add_argument("--cubic-support", type=int, default=2)
    args = parser.parse_args()
    if not 0 <= args.quadratic_support <= len(quadratic_entries):
        raise ValueError("quadratic support is out of range")
    if not 1 <= args.cubic_support <= len(cubic_monomials):
        raise ValueError("cubic support is out of range")

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
    signs = (-1, 1)
    quadratic_choices = [
        (support, coefficients)
        for support_size in range(args.quadratic_support + 1)
        for support in combinations(range(len(quadratic_entries)), support_size)
        for coefficients in product(signs, repeat=support_size)
    ]
    cubic_choices = [
        (support, coefficients)
        for support_size in range(1, args.cubic_support + 1)
        for support in combinations(range(len(cubic_monomials)), support_size)
        for coefficients in product(signs, repeat=support_size)
    ]
    candidates_per_chart = len(quadratic_choices) * len(cubic_choices)

    all_modular_survivors = []
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
            q0_value = [q0[index].subs(substitution) for index in range(4)]
            u_substitution = dict(zip(u, q0_value, strict=True))
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
            evaluated_points.append(
                (
                    jm_value,
                    quadratic_contributions,
                    cubic_contributions,
                )
            )

        # Cubic Hessians vanish at the origin, so filter quadratic choices
        # once before taking their product with all cubic choices.
        valid_quadratic_choices = []
        origin_jm, origin_quadratic, _ = evaluated_points[0]
        for support, coefficients in quadratic_choices:
            origin_matrix = add_contributions(
                origin_jm, origin_quadratic, support, coefficients
            )
            origin_determinant = determinant_mod(origin_matrix, prime)
            if origin_determinant != 0:
                valid_quadratic_choices.append(
                    (support, coefficients, origin_determinant)
                )

        chart_survivors = []
        for (
            quadratic_support,
            quadratic_coefficients,
            expected_determinant,
        ) in valid_quadratic_choices:
            for cubic_support, cubic_coefficients in cubic_choices:
                passes = True
                for (
                    jm_value,
                    quadratic_contributions,
                    cubic_contributions,
                ) in evaluated_points[1:]:
                    candidate_matrix = add_contributions(
                        jm_value,
                        quadratic_contributions,
                        quadratic_support,
                        quadratic_coefficients,
                    )
                    candidate_matrix = add_contributions(
                        candidate_matrix,
                        cubic_contributions,
                        cubic_support,
                        cubic_coefficients,
                    )
                    if (
                        determinant_mod(candidate_matrix, prime)
                        != expected_determinant
                    ):
                        passes = False
                        break
                if passes:
                    chart_survivors.append(
                        (
                            quadratic_support,
                            quadratic_coefficients,
                            cubic_support,
                            cubic_coefficients,
                        )
                    )

        all_modular_survivors.extend(
            (mask, *survivor) for survivor in chart_survivors
        )

        for (
            quadratic_support,
            quadratic_coefficients,
            cubic_support,
            cubic_coefficients,
        ) in chart_survivors:
            potential = sum(
                coefficient * quadratic_potentials[index]
                for index, coefficient in zip(
                    quadratic_support,
                    quadratic_coefficients,
                    strict=True,
                )
            )
            potential += sum(
                coefficient * cubic_monomials[index][1]
                for index, coefficient in zip(
                    cubic_support, cubic_coefficients, strict=True
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
                        quadratic_support,
                        quadratic_coefficients,
                        cubic_support,
                        cubic_coefficients,
                        determinant,
                    )
                )

        label = "".join(map(str, mask))
        print(
            f"chart {label}: {candidates_per_chart} candidates "
            f"({len(valid_quadratic_choices)} origin-nondegenerate "
            f"quadratic parts), {len(chart_survivors)} modular survivors"
        )

    print(f"total modular survivors: {len(all_modular_survivors)}")
    print(f"exact constant-nonzero survivors: {len(exact_survivors)}")
    for survivor in exact_survivors:
        print(survivor)


if __name__ == "__main__":
    main()

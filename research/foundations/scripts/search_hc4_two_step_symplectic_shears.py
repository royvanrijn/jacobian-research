#!/usr/bin/env python3
"""Search two-step nonlinear symplectic polarizations of the PC(2) graph.

Starting from ambient Darboux coordinates (Q,M), apply

    M1 = M + grad V(Q),
    Q1 = Q + grad W(M1).

The components of Q1 Poisson-commute, so B=Q1 is a nonlinear Lagrangian
projection.  The searched potentials are

    V = 0 or +/- one cubic monomial,
    W = (1/2) M1^T C M1 + (0 or +/- one cubic monomial),

where the symmetric matrix C has at most two nonzero entries, each +/-1.
The determinant is filtered at nine exact graph points modulo a large prime,
then every survivor is verified symbolically over Q.
"""

from __future__ import annotations

from itertools import combinations, product
import runpy

import sympy as sp


graph_namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
utility_namespace = runpy.run_path("scripts/search_hc4_lagrangian_shears.py")

h_variables = graph_namespace["h_variables"]
Q = sp.Matrix(graph_namespace["position_coordinates_h"])
M = sp.Matrix(graph_namespace["momentum_coordinates_h"])
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
cubic_gradients = [
    sp.Matrix([sp.diff(monomial, variable) for variable in u])
    for monomial in cubic_monomials
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


def rational_vector_mod(vector: sp.Matrix) -> list[int]:
    return [rational_mod(vector[index], prime) for index in range(4)]


def matrix_multiply_mod(
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


def add_scaled_matrix(
    base: list[list[int]],
    contribution: list[list[int]],
    coefficient: int,
) -> list[list[int]]:
    return [
        [
            (base[row][column] + coefficient * contribution[row][column])
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]


def main() -> None:
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
    c_choices = [
        (support, coefficients)
        for support_size in (0, 1, 2)
        for support in combinations(range(len(quadratic_entries)), support_size)
        for coefficients in product(signs, repeat=support_size)
    ]
    cubic_choices = [(None, 0)] + [
        (index, sign)
        for index in range(len(cubic_monomials))
        for sign in signs
    ]

    jq = Q.jacobian(h_variables)
    jm = M.jacobian(h_variables)
    evaluated_points = []
    for point in sample_points:
        substitution = dict(zip(h_variables, point, strict=True))
        q_value = Q.subs(substitution)
        m_value = M.subs(substitution)
        jq_value = rational_matrix_mod(jq.subs(substitution))
        jm_value = rational_matrix_mod(jm.subs(substitution))
        evaluated_points.append(
            (q_value, m_value, jq_value, jm_value)
        )

    # Precompute M1, dM1, the C-basis contributions, and every W3
    # contribution for each V choice and sample point.
    precomputed = {}
    for v_choice in cubic_choices:
        v_index, v_sign = v_choice
        point_data = []
        for q_value, m_value, jq_value, jm_value in evaluated_points:
            if v_index is None:
                gradient_value = sp.zeros(4, 1)
                hessian_value_mod = [[0] * 4 for _ in range(4)]
            else:
                q_substitution = dict(zip(u, list(q_value), strict=True))
                gradient_value = (
                    v_sign
                    * cubic_gradients[v_index].subs(q_substitution)
                )
                hessian_value_mod = rational_matrix_mod(
                    v_sign
                    * cubic_hessians[v_index].subs(q_substitution)
                )

            m1_value = m_value + gradient_value
            v_contribution = matrix_multiply_mod(
                hessian_value_mod, jq_value
            )
            jm1_value = [
                [
                    (jm_value[row][column] + v_contribution[row][column])
                    % prime
                    for column in range(4)
                ]
                for row in range(4)
            ]
            c_contributions = [
                matrix_multiply_mod(
                    rational_matrix_mod(hessian), jm1_value
                )
                for hessian in quadratic_hessians
            ]
            w_contributions = {}
            m1_substitution = dict(zip(u, list(m1_value), strict=True))
            for w_index, w_sign in cubic_choices:
                if w_index is None:
                    w_contributions[(w_index, w_sign)] = [
                        [0] * 4 for _ in range(4)
                    ]
                else:
                    w_hessian_mod = rational_matrix_mod(
                        w_sign
                        * cubic_hessians[w_index].subs(m1_substitution)
                    )
                    w_contributions[(w_index, w_sign)] = matrix_multiply_mod(
                        w_hessian_mod, jm1_value
                    )
            point_data.append(
                (jq_value, c_contributions, w_contributions)
            )
        precomputed[v_choice] = point_data

    # At the origin all cubic Hessians vanish.  Filter C once.
    origin_jq = precomputed[(None, 0)][0][0]
    origin_c_contributions = precomputed[(None, 0)][0][1]
    valid_c_choices = []
    for support, coefficients in c_choices:
        matrix = [row[:] for row in origin_jq]
        for index, coefficient in zip(support, coefficients, strict=True):
            matrix = add_scaled_matrix(
                matrix, origin_c_contributions[index], coefficient
            )
        determinant = determinant_mod(matrix, prime)
        if determinant != 0:
            valid_c_choices.append((support, coefficients, determinant))

    modular_survivors = []
    for c_support, c_coefficients, expected_determinant in valid_c_choices:
        for v_choice in cubic_choices:
            for w_choice in cubic_choices:
                passes = True
                for jq_value, c_contributions, w_contributions in precomputed[
                    v_choice
                ][1:]:
                    matrix = [row[:] for row in jq_value]
                    for index, coefficient in zip(
                        c_support, c_coefficients, strict=True
                    ):
                        matrix = add_scaled_matrix(
                            matrix, c_contributions[index], coefficient
                        )
                    matrix = add_scaled_matrix(
                        matrix, w_contributions[w_choice], 1
                    )
                    if determinant_mod(matrix, prime) != expected_determinant:
                        passes = False
                        break
                if passes:
                    modular_survivors.append(
                        (c_support, c_coefficients, v_choice, w_choice)
                    )

    exact_survivors = []
    for c_support, c_coefficients, v_choice, w_choice in modular_survivors:
        v_index, v_sign = v_choice
        w_index, w_sign = w_choice
        v_potential = (
            0
            if v_index is None
            else v_sign * cubic_monomials[v_index]
        )
        c_potential = sum(
            coefficient * quadratic_potentials[index]
            for index, coefficient in zip(
                c_support, c_coefficients, strict=True
            )
        )
        w_potential = c_potential + (
            0
            if w_index is None
            else w_sign * cubic_monomials[w_index]
        )

        gradient_v = sp.Matrix(
            [sp.diff(v_potential, variable) for variable in u]
        ).subs(dict(zip(u, Q, strict=True)), simultaneous=True)
        m1 = M + gradient_v
        gradient_w = sp.Matrix(
            [sp.diff(w_potential, variable) for variable in u]
        ).subs(dict(zip(u, m1, strict=True)), simultaneous=True)
        b_map = Q + gradient_w
        determinant = sp.factor(
            b_map.jacobian(h_variables).det(method="berkowitz")
        )
        if determinant != 0 and not determinant.free_symbols:
            exact_survivors.append(
                (
                    c_support,
                    c_coefficients,
                    v_choice,
                    w_choice,
                    determinant,
                    b_map,
                    m1,
                )
            )

    candidate_count = (
        len(c_choices) * len(cubic_choices) * len(cubic_choices)
    )
    print(
        f"two-step candidates: {candidate_count} "
        f"({len(valid_c_choices)} origin-nondegenerate C choices)"
    )
    print(f"modular survivors: {len(modular_survivors)}")
    print(f"exact constant-nonzero survivors: {len(exact_survivors)}")
    for survivor in exact_survivors:
        print(survivor[:-2])


if __name__ == "__main__":
    main()

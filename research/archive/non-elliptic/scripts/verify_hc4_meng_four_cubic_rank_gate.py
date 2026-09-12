#!/usr/bin/env python3
"""Finite-field odd-layer gate for four-monomial cubic corrections.

The sparse quartic checker leaves 234 principal parts.  This continuation
enumerates all 4,845 quadruples of the 20 cubic monomials over each quartic.
Determinant degrees seven and one form a linear system in the four cubic
coefficients.  Rank four is inconsistent.  Rank-three systems leave a line:
coordinate-boundary lines are covered by the characteristic-zero
three-cubic theorem, while every genuine four-support line is tested by the
full determinant gcd over F_1000003.

Lower-rank coefficient spaces are counted and split into boundary and
genuine loci for the next nonlinear lift.  This checker is an exact
finite-field computation, not yet a characteristic-zero theorem.
"""

from __future__ import annotations

import contextlib
import io
from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp


PARENT = Path(__file__).with_name(
    "verify_hc4_meng_sparse_quartic_obstruction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    parent = runpy.run_path(str(PARENT))

PRIME: int = parent["PRIME"]
VARIABLE_COUNT: int = parent["VARIABLE_COUNT"]
quartics = parent["all_principal_quartics_mod"]
quartic_exponents = parent["quartic_exponents"]
cubic_exponents = parent["cubic_exponents"]
sample_points = parent["sample_points"]
points = parent["two_cubic_points"]
quartic_monomial_hessians = parent["two_cubic_quartic_hessians_mod"]
cubic_monomial_hessians = parent["two_cubic_cubic_hessians_mod"]
base_hessian = parent["base_hessian_mod"]
determinant_mod = parent["determinant_mod"]
interpolate_degree_four_mod = parent["interpolate_degree_four_mod"]
polynomial_gcd_mod = parent["polynomial_gcd_mod"]
scaling_coefficient_mod = parent["scaling_coefficient_mod"]


def inverse_mod(value: int) -> int:
    return pow(value % PRIME, -1, PRIME)


def rank_and_null_basis(
    columns: tuple[tuple[int, ...], ...],
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    column_count = len(columns)
    rows = [list(row) for row in zip(*columns, strict=True)]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(rows))
                if rows[row][column] % PRIME
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = inverse_mod(rows[pivot_row][column])
        rows[pivot_row] = [
            entry * scale % PRIME for entry in rows[pivot_row]
        ]
        for row in range(len(rows)):
            if row == pivot_row or rows[row][column] % PRIME == 0:
                continue
            scale = rows[row][column] % PRIME
            rows[row] = [
                (
                    rows[row][index]
                    - scale * rows[pivot_row][index]
                )
                % PRIME
                for index in range(column_count)
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == column_count:
            break

    free_columns = [
        column for column in range(column_count) if column not in pivots
    ]
    basis = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -rows[row][free] % PRIME
        basis.append(tuple(vector))
    return len(pivots), tuple(basis)


def add_scaled(
    target: list[list[int]],
    scale: int,
    contribution: tuple[tuple[int, ...], ...],
) -> None:
    for row in range(VARIABLE_COUNT):
        for column in range(VARIABLE_COUNT):
            target[row][column] = (
                target[row][column]
                + scale * contribution[row][column]
            ) % PRIME


degree_point_indices = tuple(range(4, len(sample_points)))
rank_counts = {rank: 0 for rank in range(5)}
rank_three_boundary_lines = 0
rank_three_genuine_lines = 0
rank_three_survivors = []
rank_two_boundary_planes = 0
rank_two_genuine_planes = 0
rank_two_survivors = []
maximum_rank_two_points = 0
rank_one_boundary_spaces = 0
rank_one_genuine_spaces = 0
rank_one_survivors = []
maximum_rank_one_points = 0
lambda_parameter, mu_parameter, nu_parameter = sp.symbols("lambda mu nu")

for quartic_index, (support, coefficients) in enumerate(quartics):
    quartic_hessians = []
    base_matrices = []
    for point_index in range(len(points)):
        quartic_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            add_scaled(
                quartic_hessian,
                coefficient,
                quartic_monomial_hessians[point_index][monomial_index],
            )
        quartic_hessians.append(quartic_hessian)
        base_matrices.append(
            [
                [
                    (
                        base_hessian[row][column]
                        + quartic_hessian[row][column]
                    )
                    % PRIME
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
        )

    signatures = []
    for cubic_index in range(len(cubic_exponents)):
        signature = []
        for point_index in degree_point_indices:
            values = []
            for coefficient in range(5):
                values.append(
                    determinant_mod(
                        [
                            [
                                (
                                    quartic_hessians[point_index][row][column]
                                    + coefficient
                                    * cubic_monomial_hessians[point_index][
                                        cubic_index
                                    ][row][column]
                                )
                                % PRIME
                                for column in range(VARIABLE_COUNT)
                            ]
                            for row in range(VARIABLE_COUNT)
                        ]
                    )
                )
            polynomial = interpolate_degree_four_mod(values)
            signature.append(
                polynomial[1] if len(polynomial) > 1 else 0
            )
        for point_index in degree_point_indices:
            signature.append(
                scaling_coefficient_mod(
                    quartic_hessians[point_index],
                    cubic_monomial_hessians[point_index][cubic_index],
                    1,
                )
            )
        signatures.append(tuple(signature))

    for cubic_quadruple in combinations(range(len(cubic_exponents)), 4):
        rank, null_basis = rank_and_null_basis(
            tuple(signatures[index] for index in cubic_quadruple)
        )
        rank_counts[rank] += 1
        if rank == 4:
            continue
        if rank == 3:
            assert len(null_basis) == 1
            direction = null_basis[0]
            if 0 in direction:
                rank_three_boundary_lines += 1
                continue
            rank_three_genuine_lines += 1
            common_gcd = None
            for point_index in range(len(points)):
                directed_hessian = [
                    [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
                ]
                for coefficient, cubic_index in zip(
                    direction, cubic_quadruple, strict=True
                ):
                    add_scaled(
                        directed_hessian,
                        coefficient,
                        cubic_monomial_hessians[point_index][cubic_index],
                    )
                values = []
                for parameter in range(5):
                    values.append(
                        (
                            determinant_mod(
                                [
                                    [
                                        (
                                            base_matrices[point_index][row][
                                                column
                                            ]
                                            + parameter
                                            * directed_hessian[row][column]
                                        )
                                        % PRIME
                                        for column in range(VARIABLE_COUNT)
                                    ]
                                    for row in range(VARIABLE_COUNT)
                                ]
                            )
                            - 64
                        )
                        % PRIME
                    )
                equation = interpolate_degree_four_mod(values)
                if len(equation) == 1 and equation[0] == 0:
                    continue
                common_gcd = (
                    equation
                    if common_gcd is None
                    else polynomial_gcd_mod(common_gcd, equation)
                )
                if len(common_gcd) == 1:
                    break
            if common_gcd is None or len(common_gcd) > 1:
                rank_three_survivors.append(
                    (quartic_index, cubic_quadruple, direction, common_gcd)
                )
            continue

        if rank == 2:
            assert len(null_basis) == 2
            if any(
                all(vector[index] == 0 for vector in null_basis)
                for index in range(4)
            ):
                rank_two_boundary_planes += 1
            else:
                rank_two_genuine_planes += 1
                direction_hessians = []
                for direction in null_basis:
                    point_hessians = []
                    for point_index in range(len(points)):
                        directed_hessian = [
                            [0] * VARIABLE_COUNT
                            for _ in range(VARIABLE_COUNT)
                        ]
                        for coefficient, cubic_index in zip(
                            direction, cubic_quadruple, strict=True
                        ):
                            add_scaled(
                                directed_hessian,
                                coefficient,
                                cubic_monomial_hessians[point_index][
                                    cubic_index
                                ],
                            )
                        point_hessians.append(directed_hessian)
                    direction_hessians.append(point_hessians)

                equations = []
                unit = False
                for point_index in range(len(points)):
                    candidate = sp.Matrix(base_matrices[point_index])
                    candidate += lambda_parameter * sp.Matrix(
                        direction_hessians[0][point_index]
                    )
                    candidate += mu_parameter * sp.Matrix(
                        direction_hessians[1][point_index]
                    )
                    equation = sp.Poly(
                        sp.expand(
                            candidate.det(method="berkowitz") - 64
                        ),
                        lambda_parameter,
                        mu_parameter,
                        modulus=PRIME,
                    )
                    if not equation.is_zero:
                        equations.append(equation.as_expr())
                    if point_index < 2:
                        continue
                    basis = sp.groebner(
                        equations,
                        lambda_parameter,
                        mu_parameter,
                        modulus=PRIME,
                        order="grevlex",
                    )
                    unit = (
                        len(basis.polys) == 1
                        and sp.expand(basis.polys[0].as_expr()) == 1
                    )
                    if unit:
                        maximum_rank_two_points = max(
                            maximum_rank_two_points, point_index + 1
                        )
                        break
                if not unit:
                    rank_two_survivors.append(
                        (quartic_index, cubic_quadruple)
                    )
            continue

        if rank == 1:
            assert len(null_basis) == 3
            if any(
                all(vector[index] == 0 for vector in null_basis)
                for index in range(4)
            ):
                rank_one_boundary_spaces += 1
            else:
                rank_one_genuine_spaces += 1
                direction_hessians = []
                for direction in null_basis:
                    point_hessians = []
                    for point_index in range(len(points)):
                        directed_hessian = [
                            [0] * VARIABLE_COUNT
                            for _ in range(VARIABLE_COUNT)
                        ]
                        for coefficient, cubic_index in zip(
                            direction, cubic_quadruple, strict=True
                        ):
                            add_scaled(
                                directed_hessian,
                                coefficient,
                                cubic_monomial_hessians[point_index][
                                    cubic_index
                                ],
                            )
                        point_hessians.append(directed_hessian)
                    direction_hessians.append(point_hessians)

                equations = []
                unit = False
                parameters = (
                    lambda_parameter,
                    mu_parameter,
                    nu_parameter,
                )
                for point_index in range(len(points)):
                    candidate = sp.Matrix(base_matrices[point_index])
                    for parameter, direction in zip(
                        parameters, direction_hessians, strict=True
                    ):
                        candidate += parameter * sp.Matrix(
                            direction[point_index]
                        )
                    equation = sp.Poly(
                        sp.expand(
                            candidate.det(method="berkowitz") - 64
                        ),
                        *parameters,
                        modulus=PRIME,
                    )
                    if not equation.is_zero:
                        equations.append(equation.as_expr())
                    if point_index < 2:
                        continue
                    basis = sp.groebner(
                        equations,
                        *parameters,
                        modulus=PRIME,
                        order="grevlex",
                    )
                    unit = (
                        len(basis.polys) == 1
                        and sp.expand(basis.polys[0].as_expr()) == 1
                    )
                    if unit:
                        maximum_rank_one_points = max(
                            maximum_rank_one_points, point_index + 1
                        )
                        break
                if not unit:
                    rank_one_survivors.append(
                        (quartic_index, cubic_quadruple)
                    )


assert rank_counts == {
    0: 5_430,
    1: 79_396,
    2: 353_740,
    3: 504_818,
    4: 190_346,
}
assert not rank_three_survivors
assert not rank_two_survivors
assert not rank_one_survivors

print(
    "PASS: all 1133730 quartic/quadruple pairs have odd-layer ranks "
    "r0=5430, r1=79396, r2=353740, r3=504818, r4=190346"
)
print(
    "PASS: every genuine rank-three four-support line has unit "
    "full-determinant gcd modulo 1000003"
)
print(
    f"DETAIL: rank-three boundary={rank_three_boundary_lines}, "
    f"genuine={rank_three_genuine_lines}"
)
print(
    "PASS: all 6082 genuine rank-two planes have unit full-determinant "
    f"ideals modulo 1000003 using at most {maximum_rank_two_points} points"
)
print(
    f"DETAIL: rank-two boundary={rank_two_boundary_planes}, "
    f"genuine={rank_two_genuine_planes}"
)
print(
    "PASS: all 7956 genuine rank-one three-spaces have unit "
    "full-determinant ideals modulo 1000003 using at most "
    f"{maximum_rank_one_points} points"
)
print(
    f"DETAIL: rank-one boundary={rank_one_boundary_spaces}, "
    f"genuine={rank_one_genuine_spaces}"
)
print(
    f"FRONTIER: the {rank_counts[0]} rank-zero four-parameter spaces "
    "remain"
)

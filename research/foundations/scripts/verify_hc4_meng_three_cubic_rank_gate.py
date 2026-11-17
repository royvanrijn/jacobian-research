#!/usr/bin/env python3
"""Screen three-monomial cubic corrections by the odd determinant layers.

This checker continues ``verify_hc4_meng_sparse_quartic_obstruction.py``.
That certificate leaves 234 quartic principal parts.  For each one and each
of the 1,140 triples of cubic monomials, this script forms the coefficient
matrix of determinant degrees seven and one.  Both layers are linear in the
three cubic coefficients.

Rank three is immediately inconsistent.  Rank two leaves a single cubic
direction and is tested by univariate gcds.  Rank one leaves a plane and is
tested by bivariate determinant ideals.  Rank zero proceeds directly to the
full three-parameter determinant ideal, interpolated in the 35 monomials of
total degree at most four.  Boundary lines and planes of support at most two
are delegated to the preceding certificate.

The calculation is over F_1000003.  It is an exact finite-field screen and a
reproducible search result; it is not, by itself, a characteristic-zero
proof for arbitrary rational cubic coefficients.
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

# The parent is itself a certificate.  Running it supplies its certified
# quartic list and small exact-arithmetic helpers; suppress its PASS lines so
# this checker's output states only the new result.
with contextlib.redirect_stdout(io.StringIO()):
    parent = runpy.run_path(str(PARENT))

PRIME: int = parent["PRIME"]
VARIABLE_COUNT: int = parent["VARIABLE_COUNT"]
all_principal_quartics_mod = parent["all_principal_quartics_mod"]
quartic_exponents = parent["quartic_exponents"]
cubic_exponents = parent["cubic_exponents"]
sample_points = parent["sample_points"]
two_cubic_points = parent["two_cubic_points"]
quartic_monomial_hessians = parent["two_cubic_quartic_hessians_mod"]
cubic_monomial_hessians = parent["two_cubic_cubic_hessians_mod"]
base_hessian_mod = parent["base_hessian_mod"]
determinant_mod = parent["determinant_mod"]
interpolate_degree_four_mod = parent["interpolate_degree_four_mod"]
polynomial_gcd_mod = parent["polynomial_gcd_mod"]
scaling_coefficient_mod = parent["scaling_coefficient_mod"]
interpolation_inverse = parent["interpolation_inverse"]
invert_matrix_mod = parent["invert_matrix_mod"]


def inverse_mod(value: int) -> int:
    return pow(value % PRIME, -1, PRIME)


def matrix_rank_and_null_vector_mod(
    columns: tuple[tuple[int, ...], ...],
) -> tuple[int, tuple[int, ...] | None]:
    """Return column rank and, only in corank one, a null vector."""

    column_count = len(columns)
    rows = [list(row) for row in zip(*columns, strict=True)]
    pivot_columns: list[int] = []
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
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == column_count:
            break

    rank = len(pivot_columns)
    if rank != column_count - 1:
        return rank, None

    free_column = next(
        column
        for column in range(column_count)
        if column not in pivot_columns
    )
    null_vector = [0] * column_count
    null_vector[free_column] = 1
    for row, column in enumerate(pivot_columns):
        null_vector[column] = -rows[row][free_column] % PRIME
    return rank, tuple(null_vector)


def add_scaled_matrix(
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


def rank_one_null_basis(
    columns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]]:
    """Return one nonzero row equation and a basis of its null plane."""

    equation = next(
        tuple(column[row] for column in columns)
        for row in range(len(columns[0]))
        if any(column[row] for column in columns)
    )
    pivot = next(index for index, entry in enumerate(equation) if entry)
    basis: list[tuple[int, ...]] = []
    for free in range(len(columns)):
        if free == pivot:
            continue
        vector = [0] * len(columns)
        vector[free] = 1
        vector[pivot] = (
            -equation[free] * inverse_mod(equation[pivot])
        ) % PRIME
        basis.append(tuple(vector))
    assert len(basis) == 2
    return equation, (basis[0], basis[1])


def bivariate_determinant_equation(
    base_matrix: list[list[int]],
    left_direction: list[list[int]],
    right_direction: list[list[int]],
) -> sp.Expr:
    """Interpolate det(base+lambda*left+mu*right)-64."""

    values = []
    for left_value in range(5):
        row_values = []
        for right_value in range(5):
            candidate = [
                [
                    (
                        base_matrix[row][column]
                        + left_value * left_direction[row][column]
                        + right_value * right_direction[row][column]
                    )
                    % PRIME
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
            row_values.append((determinant_mod(candidate) - 64) % PRIME)
        values.append(row_values)

    coefficients = [
        [
            sum(
                interpolation_inverse[left_degree][left_value]
                * interpolation_inverse[right_degree][right_value]
                * values[left_value][right_value]
                for left_value in range(5)
                for right_value in range(5)
            )
            % PRIME
            for right_degree in range(5)
        ]
        for left_degree in range(5)
    ]
    assert all(
        coefficients[left_degree][right_degree] == 0
        for left_degree in range(5)
        for right_degree in range(5)
        if left_degree + right_degree > 4
    )
    return sp.Add(
        *(
            coefficient
            * left_parameter**left_degree
            * right_parameter**right_degree
            for left_degree, row in enumerate(coefficients)
            for right_degree, coefficient in enumerate(row)
            if coefficient
        )
    )


parameter_monomials = tuple(
    (left_degree, middle_degree, right_degree)
    for total_degree in range(5)
    for left_degree in range(total_degree + 1)
    for middle_degree in range(total_degree - left_degree + 1)
    for right_degree in (total_degree - left_degree - middle_degree,)
)
assert len(parameter_monomials) == 35
parameter_evaluation_points = parameter_monomials
trivariate_vandermonde = [
    [
        (
            pow(point[0], exponents[0], PRIME)
            * pow(point[1], exponents[1], PRIME)
            * pow(point[2], exponents[2], PRIME)
        )
        % PRIME
        for exponents in parameter_monomials
    ]
    for point in parameter_evaluation_points
]
trivariate_interpolation_inverse = invert_matrix_mod(
    trivariate_vandermonde
)


def trivariate_determinant_equation(
    base_matrix: list[list[int]],
    directions: tuple[
        list[list[int]], list[list[int]], list[list[int]]
    ],
) -> sp.Expr:
    """Interpolate the total-degree-four determinant in three parameters."""

    values = []
    for parameter_point in parameter_evaluation_points:
        candidate = [
            [
                (
                    base_matrix[row][column]
                    + sum(
                        parameter_point[index]
                        * directions[index][row][column]
                        for index in range(3)
                    )
                )
                % PRIME
                for column in range(VARIABLE_COUNT)
            ]
            for row in range(VARIABLE_COUNT)
        ]
        values.append((determinant_mod(candidate) - 64) % PRIME)
    coefficients = [
        sum(
            trivariate_interpolation_inverse[row][column] * values[column]
            for column in range(35)
        )
        % PRIME
        for row in range(35)
    ]
    return sp.Add(
        *(
            coefficient
            * first_parameter**exponents[0]
            * second_parameter**exponents[1]
            * third_parameter**exponents[2]
            for coefficient, exponents in zip(
                coefficients, parameter_monomials, strict=True
            )
            if coefficient
        )
    )


degree_seven_point_indices = tuple(range(4, len(sample_points)))
rank_counts = {0: 0, 1: 0, 2: 0, 3: 0}
rank_two_boundary_lines = 0
rank_two_genuine_lines = 0
rank_two_full_determinant_survivors: list[
    tuple[int, tuple[int, int, int], tuple[int, int, int], list[int] | None]
] = []
rank_one_boundary_planes = 0
rank_one_genuine_planes = 0
rank_one_full_determinant_survivors: list[
    tuple[int, tuple[int, int, int], tuple[int, ...]]
] = []
left_parameter, right_parameter = sp.symbols("lambda mu")
first_parameter, second_parameter, third_parameter = sp.symbols(
    "lambda mu nu"
)
rank_zero_full_determinant_survivors: list[
    tuple[int, tuple[int, int, int]]
] = []

for quartic_index, (support, coefficients) in enumerate(
    all_principal_quartics_mod
):
    quartic_hessians: list[list[list[int]]] = []
    base_matrices: list[list[list[int]]] = []
    for point_index in range(len(two_cubic_points)):
        quartic_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            add_scaled_matrix(
                quartic_hessian,
                coefficient,
                quartic_monomial_hessians[point_index][monomial_index],
            )
        quartic_hessians.append(quartic_hessian)
        base_matrices.append(
            [
                [
                    (
                        base_hessian_mod[row][column]
                        + quartic_hessian[row][column]
                    )
                    % PRIME
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
        )

    odd_layer_signatures: list[tuple[int, ...]] = []
    for cubic_index in range(len(cubic_exponents)):
        signature: list[int] = []
        for point_index in degree_seven_point_indices:
            determinant_values = []
            for cubic_coefficient in range(5):
                candidate = [
                    [
                        (
                            quartic_hessians[point_index][row][column]
                            + cubic_coefficient
                            * cubic_monomial_hessians[point_index][
                                cubic_index
                            ][row][column]
                        )
                        % PRIME
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                determinant_values.append(determinant_mod(candidate))
            interpolated = interpolate_degree_four_mod(determinant_values)
            signature.append(
                interpolated[1] if len(interpolated) > 1 else 0
            )
        for point_index in degree_seven_point_indices:
            signature.append(
                scaling_coefficient_mod(
                    quartic_hessians[point_index],
                    cubic_monomial_hessians[point_index][cubic_index],
                    1,
                )
            )
        odd_layer_signatures.append(tuple(signature))

    for cubic_triple in combinations(range(len(cubic_exponents)), 3):
        rank, null_vector = matrix_rank_and_null_vector_mod(
            tuple(
                odd_layer_signatures[cubic_index]
                for cubic_index in cubic_triple
            )
        )
        rank_counts[rank] += 1
        if rank == 0:
            equations = []
            for point_index in range(len(two_cubic_points)):
                directions = tuple(
                    [
                        list(row)
                        for row in cubic_monomial_hessians[point_index][
                            cubic_index
                        ]
                    ]
                    for cubic_index in cubic_triple
                )
                equations.append(
                    trivariate_determinant_equation(
                        base_matrices[point_index], directions
                    )
                )
            equations = [equation for equation in equations if equation != 0]
            basis = sp.groebner(
                equations,
                first_parameter,
                second_parameter,
                third_parameter,
                modulus=PRIME,
                order="grevlex",
            )
            if {
                sp.expand(polynomial.as_expr())
                for polynomial in basis.polys
            } != {sp.Integer(1)}:
                rank_zero_full_determinant_survivors.append(
                    (quartic_index, cubic_triple)
                )
            continue
        if rank == 1:
            equation, null_basis = rank_one_null_basis(
                tuple(
                    odd_layer_signatures[cubic_index]
                    for cubic_index in cubic_triple
                )
            )
            if sum(entry != 0 for entry in equation) == 1:
                rank_one_boundary_planes += 1
                continue

            rank_one_genuine_planes += 1
            cubic_directions = []
            for direction in null_basis:
                directed_hessians = []
                for point_index in range(len(two_cubic_points)):
                    directed_hessian = [
                        [0] * VARIABLE_COUNT
                        for _ in range(VARIABLE_COUNT)
                    ]
                    for coefficient, cubic_index in zip(
                        direction, cubic_triple, strict=True
                    ):
                        add_scaled_matrix(
                            directed_hessian,
                            coefficient,
                            cubic_monomial_hessians[point_index][
                                cubic_index
                            ],
                        )
                    directed_hessians.append(directed_hessian)
                cubic_directions.append(directed_hessians)

            equations = [
                bivariate_determinant_equation(
                    base_matrices[point_index],
                    cubic_directions[0][point_index],
                    cubic_directions[1][point_index],
                )
                for point_index in range(len(two_cubic_points))
            ]
            equations = [equation for equation in equations if equation != 0]
            basis = sp.groebner(
                equations,
                left_parameter,
                right_parameter,
                modulus=PRIME,
                order="grevlex",
            )
            if {
                sp.expand(polynomial.as_expr())
                for polynomial in basis.polys
            } != {sp.Integer(1)}:
                rank_one_full_determinant_survivors.append(
                    (quartic_index, cubic_triple, equation)
                )
            continue
        if rank != 2:
            continue
        assert null_vector is not None
        if 0 in null_vector:
            rank_two_boundary_lines += 1
            continue

        rank_two_genuine_lines += 1
        common_gcd: list[int] | None = None
        for point_index in range(len(two_cubic_points)):
            directed_cubic_hessian = [
                [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
            ]
            for coefficient, cubic_index in zip(
                null_vector, cubic_triple, strict=True
            ):
                add_scaled_matrix(
                    directed_cubic_hessian,
                    coefficient,
                    cubic_monomial_hessians[point_index][cubic_index],
                )

            determinant_values = []
            for line_parameter in range(5):
                candidate = [
                    [
                        (
                            base_matrices[point_index][row][column]
                            + line_parameter
                            * directed_cubic_hessian[row][column]
                        )
                        % PRIME
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                determinant_values.append(
                    (determinant_mod(candidate) - 64) % PRIME
                )
            equation = interpolate_degree_four_mod(determinant_values)
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
            rank_two_full_determinant_survivors.append(
                (
                    quartic_index,
                    cubic_triple,
                    null_vector,
                    common_gcd,
                )
            )


assert sum(rank_counts.values()) == 234 * 1_140
assert rank_counts == {
    0: 5_480,
    1: 53_364,
    2: 130_508,
    3: 77_408,
}
assert (
    rank_two_boundary_lines + rank_two_genuine_lines
    == rank_counts[2]
)
assert not rank_two_full_determinant_survivors
assert rank_one_boundary_planes == 50_412
assert rank_one_genuine_planes == 2_952
assert not rank_one_full_determinant_survivors
assert not rank_zero_full_determinant_survivors

print(
    "PASS: all 266760 quartic/triple pairs have odd-layer ranks "
    "r0=5480, r1=53364, r2=130508, r3=77408"
)
print(
    "PASS: every rank-two line is either a support-<=2 boundary line or "
    "has unit full-determinant gcd modulo 1000003"
)
print(
    f"DETAIL: rank-two boundary lines={rank_two_boundary_lines}, "
    f"genuine three-support lines={rank_two_genuine_lines}"
)
print(
    "PASS: rank-one planes split into 50412 support-<=2 boundaries and "
    "2952 genuine three-support planes; every genuine plane has unit "
    "full-determinant ideal modulo 1000003"
)
print(
    "PASS: all 5480 rank-zero triples have unit three-parameter "
    "full-determinant ideals modulo 1000003"
)
print(
    "SCOPE: this checker closes cubic support <=3 only over the certificate "
    "field; the companion QQ checker supplies the characteristic-zero lift"
)

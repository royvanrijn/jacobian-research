#!/usr/bin/env python3
"""Promote the mixed quartic--sextic plane ideals to QQ.

Exactly four sparse Meng quartic principal parts have zero immutable
determinant-degree-two signature.  A zero-gradient sextic supported on four
collinear exponent vectors has a two-dimensional coefficient kernel; there
are 519 such supports.

For every one of the 519 supports and every one of the four quartics, this
checker constructs the full determinant equations at three exact spatial
points as bivariate polynomials over QQ.  The determinant of a 4x4 matrix
affine in the two sextic parameters is expanded directly through the 24
permutations.  Singular then verifies that every three-equation ideal is
the unit ideal over QQ.

Requires Singular on PATH.  Use --limit for a non-certifying benchmark.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd, lcm
import shutil
import subprocess

import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument(
    "--limit",
    type=int,
    default=None,
    help="benchmark only the first N plane supports",
)
parser.add_argument("--batch-size", type=int, default=100)
arguments = parser.parse_args()

assert shutil.which("Singular") is not None, "Singular is required on PATH"

PRIME = 1_000_003
VARIABLE_COUNT = 4
sextic_exponents = tuple(
    (first, second, third, 6 - first - second - third)
    for first in range(7)
    for second in range(7 - first)
    for third in range(7 - first - second)
)
quartic_exponents = tuple(
    (first, second, third, 4 - first - second - third)
    for first in range(5)
    for second in range(5 - first)
    for third in range(5 - first - second)
)
assert len(sextic_exponents) == 84
assert len(quartic_exponents) == 35

collision_point = (
    Fraction(1),
    Fraction(-3, 2),
    Fraction(6),
    Fraction(81, 8),
)
sample_points = (
    (Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    (Fraction(1), Fraction(2), Fraction(3), Fraction(5)),
    (Fraction(2), Fraction(3), Fraction(5), Fraction(7)),
)
base_hessian = (
    (Fraction(0), Fraction(0), Fraction(0), Fraction(4)),
    (Fraction(0), Fraction(0), Fraction(2), Fraction(0)),
    (Fraction(0), Fraction(2), Fraction(0), Fraction(0)),
    (Fraction(4), Fraction(0), Fraction(0), Fraction(0)),
)
quartic_data = (
    (
        (0, 1, 5, 34),
        (
            Fraction(-5632, 1594323),
            Fraction(512, 177147),
            Fraction(-2048, 177147),
            Fraction(-81, 8),
        ),
    ),
    (
        (0, 32, 33, 34),
        (
            Fraction(-512, 531441),
            Fraction(3),
            Fraction(-12),
            Fraction(-297, 8),
        ),
    ),
    (
        (3, 4, 14, 18),
        (
            Fraction(-1, 54),
            Fraction(29, 576),
            Fraction(8, 9),
            Fraction(-3, 16),
        ),
    ),
    (
        (4, 12, 14, 24),
        (
            Fraction(1, 288),
            Fraction(32, 27),
            Fraction(116, 9),
            Fraction(12),
        ),
    ),
)


def inverse_mod(value: int) -> int:
    return pow(value % PRIME, -1, PRIME)


def zero_kernel_mod(
    support: tuple[int, ...],
) -> list[list[int]]:
    variable_count = len(support)
    matrix = [
        [
            sextic_exponents[monomial][row] % PRIME
            for monomial in support
        ]
        for row in range(VARIABLE_COUNT)
    ]
    pivot_columns = []
    pivot_row = 0
    for column in range(variable_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, VARIABLE_COUNT)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = (
            matrix[pivot],
            matrix[pivot_row],
        )
        scale = inverse_mod(matrix[pivot_row][column])
        matrix[pivot_row] = [
            entry * scale % PRIME for entry in matrix[pivot_row]
        ]
        for row in range(VARIABLE_COUNT):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                (
                    matrix[row][index]
                    - scale * matrix[pivot_row][index]
                )
                % PRIME
                for index in range(variable_count)
            ]
        pivot_columns.append(column)
        pivot_row += 1
    free_columns = [
        column
        for column in range(variable_count)
        if column not in pivot_columns
    ]
    kernel = []
    for free_column in free_columns:
        vector = [0] * variable_count
        vector[free_column] = 1
        for row, column in enumerate(pivot_columns):
            vector[column] = -matrix[row][free_column] % PRIME
        kernel.append(vector)
    return kernel


plane_supports = []
for support in combinations(range(len(sextic_exponents)), 4):
    if len(zero_kernel_mod(support)) == 2:
        plane_supports.append(support)
assert len(plane_supports) == 519
if arguments.limit is not None:
    assert arguments.limit > 0
    plane_supports = plane_supports[: arguments.limit]


def zero_kernel_exact(
    support: tuple[int, ...],
) -> list[list[Fraction]]:
    variable_count = len(support)
    matrix = [
        [
            Fraction(sextic_exponents[monomial][row])
            for monomial in support
        ]
        for row in range(VARIABLE_COUNT)
    ]
    pivot_columns = []
    pivot_row = 0
    for column in range(variable_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, VARIABLE_COUNT)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = (
            matrix[pivot],
            matrix[pivot_row],
        )
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [
            entry / scale for entry in matrix[pivot_row]
        ]
        for row in range(VARIABLE_COUNT):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                matrix[row][index]
                - scale * matrix[pivot_row][index]
                for index in range(variable_count)
            ]
        pivot_columns.append(column)
        pivot_row += 1
    free_columns = [
        column
        for column in range(variable_count)
        if column not in pivot_columns
    ]
    kernel = []
    for free_column in free_columns:
        vector = [Fraction(0)] * variable_count
        vector[free_column] = Fraction(1)
        for row, column in enumerate(pivot_columns):
            vector[column] = -matrix[row][free_column]
        kernel.append(vector)
    return kernel


def monomial_value(
    point: tuple[Fraction, ...],
    exponents: tuple[int, ...],
) -> Fraction:
    value = Fraction(1)
    for coordinate, exponent in zip(point, exponents, strict=True):
        value *= coordinate**exponent
    return value


collision_monomial_inverses = tuple(
    1 / monomial_value(collision_point, exponents)
    for exponents in sextic_exponents
)


def monomial_hessian(
    exponents: tuple[int, ...],
    point: tuple[Fraction, ...],
    scale: Fraction = Fraction(1),
) -> tuple[tuple[Fraction, ...], ...]:
    matrix = []
    for row in range(VARIABLE_COUNT):
        matrix_row = []
        for column in range(VARIABLE_COUNT):
            coefficient = exponents[row] * (
                exponents[column] - int(row == column)
            )
            if coefficient == 0:
                matrix_row.append(Fraction(0))
                continue
            reduced = list(exponents)
            reduced[row] -= 1
            reduced[column] -= 1
            value = Fraction(coefficient) * scale
            for coordinate, exponent in zip(
                point, reduced, strict=True
            ):
                value *= coordinate**exponent
            matrix_row.append(value)
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


sextic_hessian_basis = tuple(
    tuple(
        monomial_hessian(
            exponents,
            point,
            collision_monomial_inverses[index],
        )
        for index, exponents in enumerate(sextic_exponents)
    )
    for point in sample_points
)
quartic_hessian_basis = tuple(
    tuple(
        monomial_hessian(exponents, point)
        for exponents in quartic_exponents
    )
    for point in sample_points
)

quartic_hessians = []
for support, coefficients in quartic_data:
    point_hessians = []
    for point_index in range(len(sample_points)):
        matrix = [
            [Fraction(0)] * VARIABLE_COUNT
            for _ in range(VARIABLE_COUNT)
        ]
        for monomial, coefficient in zip(
            support, coefficients, strict=True
        ):
            contribution = quartic_hessian_basis[point_index][monomial]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    matrix[row][column] += (
                        coefficient * contribution[row][column]
                    )
        point_hessians.append(matrix)
    quartic_hessians.append(tuple(point_hessians))


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


permutation_data = tuple(
    (permutation, permutation_sign(permutation))
    for permutation in permutations(range(VARIABLE_COUNT))
)
Polynomial = dict[tuple[int, int], Fraction]


def polynomial_multiply(
    left: Polynomial, right: Polynomial
) -> Polynomial:
    result: Polynomial = {}
    for (left_u, left_v), left_coefficient in left.items():
        for (right_u, right_v), right_coefficient in right.items():
            exponent = (left_u + right_u, left_v + right_v)
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def affine_determinant_polynomial(
    constant: list[list[Fraction]],
    first: list[list[Fraction]],
    second: list[list[Fraction]],
) -> Polynomial:
    determinant: Polynomial = {}
    for permutation, sign in permutation_data:
        term: Polynomial = {(0, 0): Fraction(sign)}
        for row, column in enumerate(permutation):
            entry = {
                (0, 0): constant[row][column],
                (1, 0): first[row][column],
                (0, 1): second[row][column],
            }
            term = polynomial_multiply(term, entry)
        for exponent, coefficient in term.items():
            determinant[exponent] = (
                determinant.get(exponent, Fraction(0)) + coefficient
            )
    determinant[(0, 0)] = (
        determinant.get((0, 0), Fraction(0)) - 64
    )
    return {
        exponent: coefficient
        for exponent, coefficient in determinant.items()
        if coefficient
    }


def primitive_integer_polynomial(polynomial: Polynomial) -> str:
    denominator = 1
    for coefficient in polynomial.values():
        denominator = lcm(denominator, coefficient.denominator)
    integer_coefficients = {
        exponent: coefficient.numerator
        * (denominator // coefficient.denominator)
        for exponent, coefficient in polynomial.items()
    }
    content = 0
    for coefficient in integer_coefficients.values():
        content = gcd(content, abs(coefficient))
    assert content
    integer_coefficients = {
        exponent: coefficient // content
        for exponent, coefficient in integer_coefficients.items()
    }
    leading_exponent = max(integer_coefficients)
    if integer_coefficients[leading_exponent] < 0:
        integer_coefficients = {
            exponent: -coefficient
            for exponent, coefficient in integer_coefficients.items()
        }
    terms = []
    for (u_exponent, v_exponent), coefficient in sorted(
        integer_coefficients.items(), reverse=True
    ):
        factors = [str(coefficient)]
        if u_exponent:
            factors.append(
                "u" if u_exponent == 1 else f"u^{u_exponent}"
            )
        if v_exponent:
            factors.append(
                "v" if v_exponent == 1 else f"v^{v_exponent}"
            )
        terms.append("*".join(factors))
    return "+".join(terms).replace("+-", "-")


ideals = []
for plane_index, support in enumerate(plane_supports):
    kernel = zero_kernel_exact(support)
    assert len(kernel) == 2
    direction_hessians = []
    for point_index in range(len(sample_points)):
        point_directions = []
        for direction in kernel:
            point_directions.append(
                [
                    [
                        sum(
                            (
                                coefficient
                                * sextic_hessian_basis[point_index][
                                    monomial
                                ][row][column]
                                for coefficient, monomial in zip(
                                    direction, support, strict=True
                                )
                            ),
                            Fraction(0),
                        )
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
            )
        direction_hessians.append(point_directions)

    for quartic_index in range(len(quartic_data)):
        equations = []
        for point_index in range(len(sample_points)):
            constant = [
                [
                    base_hessian[row][column]
                    + quartic_hessians[quartic_index][point_index][
                        row
                    ][column]
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
            equations.append(
                primitive_integer_polynomial(
                    affine_determinant_polynomial(
                        constant,
                        direction_hessians[point_index][0],
                        direction_hessians[point_index][1],
                    )
                )
            )
        ideals.append(
            (
                plane_index,
                quartic_index,
                tuple(equations),
            )
        )


nonunit_ideals = []
for batch_start in range(0, len(ideals), arguments.batch_size):
    batch = ideals[batch_start : batch_start + arguments.batch_size]
    script = "ring rr=0,(u,v),dp; option(redSB);\n"
    labels = []
    for local_index, (plane_index, quartic_index, equations) in enumerate(
        batch
    ):
        label = f"{plane_index}_{quartic_index}"
        labels.append(label)
        script += (
            f"ideal I{local_index}="
            + ",".join(equations)
            + ";\n"
            + f"ideal G{local_index}=std(I{local_index});\n"
            + f"if(size(G{local_index})==1)"
            + "{"
            + f"if(G{local_index}[1]!=1)"
            + "{"
            + f'print(\"NONUNIT {label}\");'
            + "};"
            + "}else{"
            + f'print(\"NONUNIT {label}\");'
            + "};\n"
        )
    script += 'print("BATCH PASS");quit;\n'
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=900,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "?" not in result.stdout, result.stdout
    for line in result.stdout.splitlines():
        if line.startswith("NONUNIT "):
            nonunit_ideals.append(line.removeprefix("NONUNIT "))
    assert "BATCH PASS" in result.stdout

assert not nonunit_ideals, nonunit_ideals
if arguments.limit is None:
    assert len(ideals) == 519 * 4
    print(
        "PASS: all 519 zero-gradient sextic planes are reconstructed "
        "over QQ"
    )
    print(
        "PASS: all 2076 mixed quartic--sextic three-evaluation ideals "
        "are unit ideals over QQ"
    )
    print(
        "SCOPE: characteristic-zero plane-family promotion for the four "
        "low-layer quartics"
    )
else:
    print(
        f"BENCHMARK PASS: first {len(plane_supports)} plane supports, "
        f"{len(ideals)} ideals over QQ"
    )

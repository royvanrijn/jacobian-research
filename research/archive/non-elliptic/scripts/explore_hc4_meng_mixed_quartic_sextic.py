#!/usr/bin/env python3
"""Explore sparse sextics over the four low-layer Meng quartics.

The full cubic-kernel theorem leaves 234 collision-compatible quartic
principal parts.  Exactly four have zero determinant-degree-two signature,
which is necessary before adding a sextic without a cubic.  Their exact
supports and coefficients are pinned below and rechecked symbolically.

For a sextic added to one of these quartics, collision requires

    grad(h6)(p) = 0.

This script exhausts supports of at most four in that homogeneous exponent
kernel.  One-dimensional kernel directions are screened first by
det Hess(h6), then by univariate determinant-evaluation gcds.  A projective
certificate promotes 52,686 lines per quartic directly from F_1000003; the
remaining 68,460 lower-scale-degree lines per quartic are checked by exact
rational gcds.  Two-dimensional kernels are checked here by bivariate
Groebner bases over the certificate field and in the companion checker over
QQ.

With --lines-only this is a characteristic-zero line-family checker.  The
default also replays the finite-field plane calculation.
"""

from __future__ import annotations

import argparse
import contextlib
from fractions import Fraction
import io
from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument(
    "--lines-only",
    action="store_true",
    help="stop after the characteristic-zero line-family certificates",
)
arguments = parser.parse_args()

SEXTIC_CHECKER = Path(__file__).with_name(
    "verify_hc4_meng_sparse_sextic_obstruction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    sextic = runpy.run_path(str(SEXTIC_CHECKER))

PRIME = sextic["PRIME"]
VARIABLE_COUNT = sextic["VARIABLE_COUNT"]
sextic_exponents = sextic["sextic_exponents"]
sample_points = sextic["sample_points"][:5]
sample_hessians = tuple(
    point_hessians
    for point_hessians in sextic["sample_hessians"][:5]
)
exact_sample_hessians = sextic["exact_sample_hessians"][:2]
base_hessian = sextic["base_hessian"]
determinant_mod = sextic["determinant_mod"]
determinant_exact = sextic["determinant_exact"]
interpolate_degree_four_exact = sextic[
    "interpolate_degree_four_exact"
]
polynomial_gcd_exact = sextic["polynomial_gcd_exact"]
inverse_mod = sextic["inverse_mod"]

quartic_exponents = tuple(
    (first, second, third, 4 - first - second - third)
    for first in range(5)
    for second in range(5 - first)
    for third in range(5 - first - second)
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


def rational_mod(value: Fraction) -> int:
    return value.numerator * inverse_mod(value.denominator) % PRIME


def monomial_hessian_mod(
    exponents: tuple[int, ...],
    point: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    matrix = []
    for row in range(VARIABLE_COUNT):
        matrix_row = []
        for column in range(VARIABLE_COUNT):
            coefficient = exponents[row] * (
                exponents[column] - int(row == column)
            )
            if coefficient == 0:
                matrix_row.append(0)
                continue
            reduced = list(exponents)
            reduced[row] -= 1
            reduced[column] -= 1
            value = coefficient
            for coordinate, exponent in zip(
                point, reduced, strict=True
            ):
                value = value * pow(coordinate, exponent, PRIME) % PRIME
            matrix_row.append(value)
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


quartic_monomial_hessians = tuple(
    tuple(
        monomial_hessian_mod(exponents, point)
        for exponents in quartic_exponents
    )
    for point in sample_points
)
quartic_hessians = []
for support, coefficients in quartic_data:
    point_hessians = []
    for point_index in range(len(sample_points)):
        matrix = [[0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)]
        for monomial, coefficient in zip(
            support, coefficients, strict=True
        ):
            contribution = quartic_monomial_hessians[point_index][
                monomial
            ]
            scalar = rational_mod(coefficient)
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    matrix[row][column] = (
                        matrix[row][column]
                        + scalar * contribution[row][column]
                    ) % PRIME
        point_hessians.append(matrix)
    quartic_hessians.append(tuple(point_hessians))


# Recheck the pinned four quartics exactly: collision, quartic principal
# determinant, and immutable determinant-degree-two signature.
x, y, r, s = spatial_variables = sp.symbols("x y r s")
collision_point_exact = (
    sp.Rational(1),
    sp.Rational(-3, 2),
    sp.Rational(6),
    sp.Rational(81, 8),
)
base_hessian_exact = sp.Matrix(base_hessian)
collision_target_exact = -base_hessian_exact * sp.Matrix(
    collision_point_exact
)
for support, coefficients in quartic_data:
    quartic = sp.expand(
        sum(
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(
                    spatial_variables,
                    quartic_exponents[monomial],
                    strict=True,
                )
            )
            for monomial, coefficient in zip(
                support, coefficients, strict=True
            )
        )
    )
    assert sp.Matrix(
        [
            sp.diff(quartic, variable).subs(
                dict(zip(spatial_variables, collision_point_exact))
            )
            for variable in spatial_variables
        ]
    ) == collision_target_exact
    quartic_hessian = sp.hessian(quartic, spatial_variables)
    assert sp.expand(quartic_hessian.det(method="berkowitz")) == 0
    assert sp.expand(
        sp.trace(base_hessian_exact.adjugate() * quartic_hessian)
    ) == 0


def zero_gradient_kernel(
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


def directed_hessian(
    point_index: int,
    support: tuple[int, ...],
    direction: list[int],
) -> list[list[int]]:
    return [
        [
            sum(
                coefficient
                * sample_hessians[point_index][monomial][row][column]
                for coefficient, monomial in zip(
                    direction, support, strict=True
                )
            )
            % PRIME
            for column in range(VARIABLE_COUNT)
        ]
        for row in range(VARIABLE_COUNT)
    ]


def invert_matrix_mod(matrix: list[list[int]]) -> list[list[int]]:
    size = len(matrix)
    augmented = [
        row[:]
        + [int(row_index == column) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row
            for row in range(column, size)
            if augmented[row][column]
        )
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = inverse_mod(augmented[column][column])
        augmented[column] = [
            entry * scale % PRIME for entry in augmented[column]
        ]
        for row in range(size):
            if row == column or augmented[row][column] == 0:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                (
                    augmented[row][index]
                    - scale * augmented[column][index]
                )
                % PRIME
                for index in range(2 * size)
            ]
    return [row[size:] for row in augmented]


interpolation_inverse = invert_matrix_mod(
    [
        [pow(value, degree, PRIME) for degree in range(5)]
        for value in range(5)
    ]
)


def interpolate_parameter(values: list[int]) -> list[int]:
    coefficients = [
        sum(
            interpolation_inverse[row][column] * values[column]
            for column in range(5)
        )
        % PRIME
        for row in range(5)
    ]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_remainder(
    dividend: list[int], divisor: list[int]
) -> list[int]:
    remainder = dividend[:]
    leading_inverse = inverse_mod(divisor[-1])
    while len(remainder) >= len(divisor) and not (
        len(remainder) == 1 and remainder[0] == 0
    ):
        offset = len(remainder) - len(divisor)
        scale = remainder[-1] * leading_inverse % PRIME
        for index, coefficient in enumerate(divisor):
            remainder[index + offset] = (
                remainder[index + offset] - scale * coefficient
            ) % PRIME
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return remainder


def polynomial_gcd(left: list[int], right: list[int]) -> list[int]:
    while not (len(right) == 1 and right[0] == 0):
        left, right = right, polynomial_remainder(left, right)
    scale = inverse_mod(left[-1])
    return [coefficient * scale % PRIME for coefficient in left]


line_support_counts = {3: 0, 4: 0}
principal_line_counts = {3: 0, 4: 0}
line_unit_prefixes = [
    [0] * (len(sample_points) + 1) for _ in quartic_data
]
line_projective_certificates = [0] * len(quartic_data)
line_projective_exceptions = [[] for _ in quartic_data]
plane_families = []

for support_size in (3, 4):
    for support in combinations(
        range(len(sextic_exponents)), support_size
    ):
        kernel = zero_gradient_kernel(support)
        if not kernel:
            continue
        if len(kernel) == 2:
            assert support_size == 4
            plane_families.append((support, kernel))
            continue
        assert len(kernel) == 1
        line_support_counts[support_size] += 1
        direction = kernel[0]
        direction_hessians = []
        for point_index in range(len(sample_points)):
            hessian = directed_hessian(
                point_index, support, direction
            )
            direction_hessians.append(hessian)
            if determinant_mod(hessian):
                break
        else:
            principal_line_counts[support_size] += 1
            for quartic_index in range(len(quartic_data)):
                common_gcd = None
                equations = []
                for point_index in range(len(sample_points)):
                    constant_hessian = [
                        [
                            (
                                base_hessian[row][column]
                                + quartic_hessians[quartic_index][
                                    point_index
                                ][row][column]
                            )
                            % PRIME
                            for column in range(VARIABLE_COUNT)
                        ]
                        for row in range(VARIABLE_COUNT)
                    ]
                    values = []
                    for parameter in range(5):
                        candidate = [
                            [
                                (
                                    constant_hessian[row][column]
                                    + parameter
                                    * direction_hessians[point_index][
                                        row
                                    ][column]
                                )
                                % PRIME
                                for column in range(VARIABLE_COUNT)
                            ]
                            for row in range(VARIABLE_COUNT)
                        ]
                        values.append(
                            (determinant_mod(candidate) - 64) % PRIME
                        )
                    equation = interpolate_parameter(values)
                    equations.append(equation)
                    if len(equation) == 1 and equation[0] == 0:
                        continue
                    common_gcd = (
                        equation
                        if common_gcd is None
                        else polynomial_gcd(common_gcd, equation)
                    )
                    if len(common_gcd) == 1:
                        # An exact candidate has det Hess(h6)=0, so its
                        # determinant polynomial in the sextic scale has
                        # degree at most three.  A nonzero cubic coefficient
                        # rejects the unique projective point at infinity.
                        # Together with the affine unit gcd, proper reduction
                        # over Z_(p) rules out a characteristic-zero root.
                        projectively_complete = any(
                            len(previous) == 4 and previous[3]
                            for previous in equations
                        )
                        if projectively_complete:
                            line_projective_certificates[
                                quartic_index
                            ] += 1
                        else:
                            line_projective_exceptions[
                                quartic_index
                            ].append(
                                (
                                    support,
                                    tuple(direction),
                                    tuple(
                                        tuple(polynomial)
                                        for polynomial in equations
                                    ),
                                )
                            )
                        line_unit_prefixes[quartic_index][
                            point_index + 1
                        ] += 1
                        break
                else:
                    raise AssertionError(
                        (quartic_index, support, common_gcd)
                    )

assert line_support_counts == {3: 976, 4: 205494}
assert principal_line_counts == {3: 768, 4: 120378}
assert len(plane_families) == 519
expected_line_prefixes = [0, 0, 121146, 0, 0, 0]
assert all(
    prefixes == expected_line_prefixes
    for prefixes in line_unit_prefixes
)
assert all(
    certified + len(exceptions) == 121146
    for certified, exceptions in zip(
        line_projective_certificates,
        line_projective_exceptions,
        strict=True,
    )
)


def zero_gradient_kernel_exact(
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


def quartic_monomial_hessian_exact(
    exponents: tuple[int, ...],
    point: tuple[Fraction, ...],
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
            value = Fraction(coefficient)
            for coordinate, exponent in zip(
                point, reduced, strict=True
            ):
                value *= coordinate**exponent
            matrix_row.append(value)
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


exact_points = tuple(
    tuple(Fraction(coordinate) for coordinate in point)
    for point in sample_points[:2]
)
exact_quartic_monomial_hessians = tuple(
    tuple(
        quartic_monomial_hessian_exact(exponents, point)
        for exponents in quartic_exponents
    )
    for point in exact_points
)
exact_quartic_hessians = []
for support, coefficients in quartic_data:
    point_hessians = []
    for point_index in range(len(exact_points)):
        matrix = [
            [Fraction(0)] * VARIABLE_COUNT
            for _ in range(VARIABLE_COUNT)
        ]
        for monomial, coefficient in zip(
            support, coefficients, strict=True
        ):
            contribution = exact_quartic_monomial_hessians[
                point_index
            ][monomial]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    matrix[row][column] += (
                        coefficient * contribution[row][column]
                    )
        point_hessians.append(matrix)
    exact_quartic_hessians.append(tuple(point_hessians))


exact_direction_hessian_cache = {}


def exact_direction_hessians(
    support: tuple[int, ...],
) -> tuple[list[list[Fraction]], ...]:
    cached = exact_direction_hessian_cache.get(support)
    if cached is not None:
        return cached
    kernel = zero_gradient_kernel_exact(support)
    assert len(kernel) == 1
    direction = kernel[0]
    point_hessians = []
    for point_index in range(len(exact_points)):
        point_hessians.append(
            [
                [
                    sum(
                        (
                            coefficient
                            * exact_sample_hessians[point_index][
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
    result = tuple(point_hessians)
    exact_direction_hessian_cache[support] = result
    return result


exact_line_exception_counts = [0] * len(quartic_data)
for quartic_index, exceptions in enumerate(
    line_projective_exceptions
):
    for support, _, _ in exceptions:
        direction_hessians = exact_direction_hessians(support)
        common_gcd = None
        for point_index in range(len(exact_points)):
            constant_hessian = [
                [
                    Fraction(base_hessian[row][column])
                    + exact_quartic_hessians[quartic_index][
                        point_index
                    ][row][column]
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
            values = []
            for parameter in range(5):
                candidate = [
                    [
                        constant_hessian[row][column]
                        + parameter
                        * direction_hessians[point_index][row][column]
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                values.append(determinant_exact(candidate) - 64)
            equation = interpolate_degree_four_exact(values)
            if len(equation) == 1 and equation[0] == 0:
                continue
            common_gcd = (
                equation
                if common_gcd is None
                else polynomial_gcd_exact(common_gcd, equation)
            )
            if len(common_gcd) == 1:
                exact_line_exception_counts[quartic_index] += 1
                break
        else:
            raise AssertionError(
                (quartic_index, support, common_gcd)
            )

assert exact_line_exception_counts == [68460] * 4

if arguments.lines_only:
    print(
        "PASS: affine line gcds are units for all 121146 principal "
        "directions and all four quartics"
    )
    print(
        f"DETAIL: cubic-leading projective certificates "
        f"{line_projective_certificates}"
    )
    print(
        "DETAIL: lower-scale-degree exceptions "
        f"{[len(exceptions) for exceptions in line_projective_exceptions]}"
    )
    print(
        "PASS: every lower-scale-degree exception has unit exact QQ gcd "
        "at the same two evaluations"
    )
    raise SystemExit(0)


plane_parameters = sp.symbols("u v")
plane_unit_prefixes = [
    [0] * (len(sample_points) + 1) for _ in quartic_data
]
for support, kernel in plane_families:
    directed_hessians = [
        [
            directed_hessian(point_index, support, direction)
            for direction in kernel
        ]
        for point_index in range(len(sample_points))
    ]
    for quartic_index in range(len(quartic_data)):
        equations = []
        for point_index in range(len(sample_points)):
            entries = []
            for row in range(VARIABLE_COUNT):
                matrix_row = []
                for column in range(VARIABLE_COUNT):
                    entry = (
                        base_hessian[row][column]
                        + quartic_hessians[quartic_index][point_index][
                            row
                        ][column]
                    )
                    for parameter, hessian in zip(
                        plane_parameters,
                        directed_hessians[point_index],
                        strict=True,
                    ):
                        entry += parameter * hessian[row][column]
                    matrix_row.append(entry)
                entries.append(matrix_row)
            equation = sp.Poly(
                sp.Matrix(entries).det(method="berkowitz") - 64,
                *plane_parameters,
                modulus=PRIME,
            ).as_expr()
            equations.append(equation)
            basis = sp.groebner(
                equations,
                *plane_parameters,
                modulus=PRIME,
            )
            if basis.contains(sp.Integer(1)):
                plane_unit_prefixes[quartic_index][
                    point_index + 1
                ] += 1
                break
        else:
            raise AssertionError(
                (quartic_index, support, tuple(basis.polys))
            )

expected_plane_prefixes = [0, 0, 0, 519, 0, 0]
assert all(
    prefixes == expected_plane_prefixes
    for prefixes in plane_unit_prefixes
)

print(
    "PASS: exactly four quartic principal parts have zero immutable "
    "degree-two signature"
)
print(
    "PASS: zero-gradient sextic support census gives "
    f"{line_support_counts} line supports and 519 planes"
)
print(
    "PASS: principal Hessian cancellation leaves "
    f"{principal_line_counts[3] + principal_line_counts[4]} lines"
)
print(
    "PASS: every principal line is excluded in characteristic zero by "
    "projective reduction or an exact rational gcd"
)
print(
    "PASS: all 519 plane ideals are units by evaluation three for all "
    "four quartics"
)
print(
    "SCOPE: line families are excluded in characteristic zero; the 519 "
    "plane ideals per quartic remain an exact F_1000003 experiment"
)

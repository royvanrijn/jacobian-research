#!/usr/bin/env python3
"""Exclude genuinely mixed quartic--sextic corrections of total support <=4.

Use the normalized Meng potential psi_0=2*y*r+4*x*s and collision point
p=(1,-3/2,6,81/8).  A mixed correction h4+h6 satisfies the collision iff

    sum d_e*e = (-81/2,18,18,-81/2),  d_e=c_e*p^e,

where quartic and sextic exponent columns occur in the same four-row
system.  This checker exhausts every support of total size at most four
containing at least one monomial of each degree.

Isolated collision points are screened modulo 1,000,003.  For each spatial
point, det(H0+z*H4+z^2*H6)-64 is interpolated in z, so the determinant
layers are tested from degree sixteen down to degree two.  Every affine
collision line is then reconstructed over QQ and rejected by exact rational
evaluation gcds.  The 34 collision planes are handled by exact bivariate
Groebner bases over QQ.

This is a characteristic-zero bounded-total-support theorem.  It does not
cover total support at least five or simultaneous cubic corrections.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import sympy as sp


PRIME = 1_000_003
VARIABLE_COUNT = 4
TOTAL_SUPPORT_BOUND = 4


def inverse_mod(value: int) -> int:
    return pow(value % PRIME, -1, PRIME)


def determinant_mod(matrix: list[list[int]]) -> int:
    work = [[entry % PRIME for entry in row] for row in matrix]
    determinant = 1
    for column in range(VARIABLE_COUNT):
        pivot = next(
            (
                row
                for row in range(column, VARIABLE_COUNT)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % PRIME
        pivot_inverse = inverse_mod(pivot_value)
        for row in range(column + 1, VARIABLE_COUNT):
            if work[row][column] == 0:
                continue
            scale = work[row][column] * pivot_inverse % PRIME
            for index in range(column, VARIABLE_COUNT):
                work[row][index] = (
                    work[row][index]
                    - scale * work[column][index]
                ) % PRIME
    return determinant % PRIME


def determinant_exact(
    matrix: list[list[Fraction]],
) -> Fraction:
    work = [row[:] for row in matrix]
    determinant = Fraction(1)
    for column in range(VARIABLE_COUNT):
        pivot = next(
            (
                row
                for row in range(column, VARIABLE_COUNT)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, VARIABLE_COUNT):
            if work[row][column] == 0:
                continue
            scale = work[row][column] / pivot_value
            for index in range(column, VARIABLE_COUNT):
                work[row][index] -= scale * work[column][index]
    return determinant


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


def invert_matrix_exact(
    matrix: list[list[Fraction]],
) -> list[list[Fraction]]:
    size = len(matrix)
    augmented = [
        row[:]
        + [
            Fraction(int(row_index == column))
            for column in range(size)
        ]
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
        scale = augmented[column][column]
        augmented[column] = [
            entry / scale for entry in augmented[column]
        ]
        for row in range(size):
            if row == column or augmented[row][column] == 0:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                augmented[row][index]
                - scale * augmented[column][index]
                for index in range(2 * size)
            ]
    return [row[size:] for row in augmented]


scaling_inverse = invert_matrix_mod(
    [
        [pow(value, degree, PRIME) for degree in range(9)]
        for value in range(9)
    ]
)
parameter_inverse_exact = invert_matrix_exact(
    [
        [Fraction(value**degree) for degree in range(5)]
        for value in range(5)
    ]
)


def scaling_coefficients(values: list[int]) -> list[int]:
    assert len(values) == 9
    return [
        sum(
            scaling_inverse[row][column] * values[column]
            for column in range(9)
        )
        % PRIME
        for row in range(9)
    ]


def parameter_coefficients_exact(
    values: list[Fraction],
) -> list[Fraction]:
    assert len(values) == 5
    coefficients = [
        sum(
            (
                parameter_inverse_exact[row][column] * values[column]
                for column in range(5)
            ),
            Fraction(0),
        )
        for row in range(5)
    ]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_remainder_exact(
    dividend: list[Fraction],
    divisor: list[Fraction],
) -> list[Fraction]:
    remainder = dividend[:]
    while len(remainder) > 1 and remainder[-1] == 0:
        remainder.pop()
    while len(remainder) >= len(divisor) and not (
        len(remainder) == 1 and remainder[0] == 0
    ):
        offset = len(remainder) - len(divisor)
        scale = remainder[-1] / divisor[-1]
        for index, coefficient in enumerate(divisor):
            remainder[index + offset] -= scale * coefficient
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return remainder


def polynomial_gcd_exact(
    left: list[Fraction],
    right: list[Fraction],
) -> list[Fraction]:
    while not (len(right) == 1 and right[0] == 0):
        left, right = right, polynomial_remainder_exact(left, right)
    scale = left[-1]
    return [coefficient / scale for coefficient in left]


quartic_exponents = tuple(
    (first, second, third, 4 - first - second - third)
    for first in range(5)
    for second in range(5 - first)
    for third in range(5 - first - second)
)
sextic_exponents = tuple(
    (first, second, third, 6 - first - second - third)
    for first in range(7)
    for second in range(7 - first)
    for third in range(7 - first - second)
)
all_exponents = quartic_exponents + sextic_exponents
degree_split = len(quartic_exponents)
assert degree_split == 35
assert len(all_exponents) == 119

inverse_two = inverse_mod(2)
collision_target_mod = (
    -81 * inverse_two % PRIME,
    18,
    18,
    -81 * inverse_two % PRIME,
)
collision_target_exact = (
    Fraction(-81, 2),
    Fraction(18),
    Fraction(18),
    Fraction(-81, 2),
)
collision_point_mod = (
    1,
    -3 * inverse_two % PRIME,
    6,
    81 * inverse_mod(8) % PRIME,
)
collision_point_exact = (
    Fraction(1),
    Fraction(-3, 2),
    Fraction(6),
    Fraction(81, 8),
)
base_hessian_mod = (
    (0, 0, 0, 4),
    (0, 0, 2, 0),
    (0, 2, 0, 0),
    (4, 0, 0, 0),
)
base_hessian_exact = tuple(
    tuple(Fraction(entry) for entry in row)
    for row in base_hessian_mod
)
sample_points_mod = (
    (1, 1, 1, 1),
    (1, 2, 3, 5),
    (2, 3, 5, 7),
)
sample_points_exact = tuple(
    tuple(Fraction(entry) for entry in point)
    for point in sample_points_mod
)

# Hadamard bounds preserve every collision rank on reduction.
doubled_target = (-81, 36, 36, -81)
assert 6**4 < PRIME
assert (
    (6**3) ** 2 * sum(entry * entry for entry in doubled_target)
    < PRIME**2
)


def monomial_value_mod(
    point: tuple[int, ...], exponents: tuple[int, ...]
) -> int:
    value = 1
    for coordinate, exponent in zip(point, exponents, strict=True):
        value = value * pow(coordinate, exponent, PRIME) % PRIME
    return value


def monomial_value_exact(
    point: tuple[Fraction, ...],
    exponents: tuple[int, ...],
) -> Fraction:
    value = Fraction(1)
    for coordinate, exponent in zip(point, exponents, strict=True):
        value *= coordinate**exponent
    return value


collision_inverses_mod = tuple(
    inverse_mod(monomial_value_mod(collision_point_mod, exponents))
    for exponents in all_exponents
)
collision_inverses_exact = tuple(
    1 / monomial_value_exact(collision_point_exact, exponents)
    for exponents in all_exponents
)


def monomial_hessian_mod(
    point: tuple[int, ...],
    exponents: tuple[int, ...],
    scale: int,
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
            value = coefficient * scale % PRIME
            for coordinate, exponent in zip(
                point, reduced, strict=True
            ):
                value = value * pow(coordinate, exponent, PRIME) % PRIME
            matrix_row.append(value)
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


def monomial_hessian_exact(
    point: tuple[Fraction, ...],
    exponents: tuple[int, ...],
    scale: Fraction,
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


hessian_basis_mod = tuple(
    tuple(
        monomial_hessian_mod(
            point,
            exponents,
            collision_inverses_mod[index],
        )
        for index, exponents in enumerate(all_exponents)
    )
    for point in sample_points_mod
)
hessian_basis_exact = tuple(
    tuple(
        monomial_hessian_exact(
            point,
            exponents,
            collision_inverses_exact[index],
        )
        for index, exponents in enumerate(all_exponents)
    )
    for point in sample_points_exact
)


def solve_collision_mod(
    support: tuple[int, ...],
) -> tuple[list[int], list[list[int]]] | None:
    variable_count = len(support)
    matrix = [
        [
            all_exponents[monomial][row] % PRIME
            for monomial in support
        ]
        + [collision_target_mod[row]]
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
                for index in range(variable_count + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1
    for row in range(pivot_row, VARIABLE_COUNT):
        if (
            all(matrix[row][column] == 0 for column in range(variable_count))
            and matrix[row][-1]
        ):
            return None
    point = [0] * variable_count
    for row, column in enumerate(pivot_columns):
        point[column] = matrix[row][-1]
    free_columns = [
        column
        for column in range(variable_count)
        if column not in pivot_columns
    ]
    directions = []
    for free_column in free_columns:
        direction = [0] * variable_count
        direction[free_column] = 1
        for row, column in enumerate(pivot_columns):
            direction[column] = -matrix[row][free_column] % PRIME
        directions.append(direction)
    return point, directions


def solve_collision_exact(
    support: tuple[int, ...],
) -> tuple[list[Fraction], list[list[Fraction]]]:
    variable_count = len(support)
    matrix = [
        [
            Fraction(all_exponents[monomial][row])
            for monomial in support
        ]
        + [collision_target_exact[row]]
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
                for index in range(variable_count + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1
    point = [Fraction(0)] * variable_count
    for row, column in enumerate(pivot_columns):
        point[column] = matrix[row][-1]
    free_columns = [
        column
        for column in range(variable_count)
        if column not in pivot_columns
    ]
    directions = []
    for free_column in free_columns:
        direction = [Fraction(0)] * variable_count
        direction[free_column] = Fraction(1)
        for row, column in enumerate(pivot_columns):
            direction[column] = -matrix[row][free_column]
        directions.append(direction)
    return point, directions


def split_hessians_mod(
    point_index: int,
    support: tuple[int, ...],
    coefficients: list[int],
) -> tuple[list[list[int]], list[list[int]]]:
    quartic = [[0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)]
    sextic = [[0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)]
    for coefficient, monomial in zip(
        coefficients, support, strict=True
    ):
        destination = quartic if monomial < degree_split else sextic
        contribution = hessian_basis_mod[point_index][monomial]
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                destination[row][column] = (
                    destination[row][column]
                    + coefficient * contribution[row][column]
                ) % PRIME
    return quartic, sextic


def full_hessian_exact(
    point_index: int,
    support: tuple[int, ...],
    coefficients: list[Fraction],
) -> list[list[Fraction]]:
    matrix = [list(row) for row in base_hessian_exact]
    for coefficient, monomial in zip(
        coefficients, support, strict=True
    ):
        contribution = hessian_basis_exact[point_index][monomial]
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                matrix[row][column] += (
                    coefficient * contribution[row][column]
                )
    return matrix


support_counts = {
    size: {"inconsistent": 0, "isolated": 0, "line": 0, "plane": 0}
    for size in range(2, TOTAL_SUPPORT_BOUND + 1)
}
isolated_layer_rejections = {degree: 0 for degree in range(2, 17, 2)}
isolated_survivors = []
line_supports = []
plane_supports = []

for support_size in range(2, TOTAL_SUPPORT_BOUND + 1):
    for support in combinations(range(len(all_exponents)), support_size):
        if support[-1] < degree_split or support[0] >= degree_split:
            continue
        solution = solve_collision_mod(support)
        if solution is None:
            support_counts[support_size]["inconsistent"] += 1
            continue
        coefficients, directions = solution
        if len(directions) == 1:
            support_counts[support_size]["line"] += 1
            line_supports.append(support)
            continue
        if len(directions) == 2:
            support_counts[support_size]["plane"] += 1
            plane_supports.append(support)
            continue
        assert not directions
        support_counts[support_size]["isolated"] += 1

        rejected = False
        for point_index in range(len(sample_points_mod)):
            quartic_hessian, sextic_hessian = split_hessians_mod(
                point_index, support, coefficients
            )
            values = []
            for scale in range(9):
                candidate = [
                    [
                        (
                            base_hessian_mod[row][column]
                            + scale * quartic_hessian[row][column]
                            + scale
                            * scale
                            * sextic_hessian[row][column]
                        )
                        % PRIME
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                values.append(
                    (determinant_mod(candidate) - 64) % PRIME
                )
            coefficients_by_scale = scaling_coefficients(values)
            first_nonzero = next(
                (
                    degree
                    for degree in range(8, 0, -1)
                    if coefficients_by_scale[degree]
                ),
                None,
            )
            if first_nonzero is not None:
                isolated_layer_rejections[2 * first_nonzero] += 1
                rejected = True
                break
        if not rejected:
            isolated_survivors.append((support, coefficients))


expected_support_counts = {
    2: {"inconsistent": 2930, "isolated": 10, "line": 0, "plane": 0},
    3: {
        "inconsistent": 165514,
        "isolated": 6446,
        "line": 30,
        "plane": 0,
    },
    4: {
        "inconsistent": 695358,
        "isolated": 5219228,
        "line": 44270,
        "plane": 34,
    },
}
assert support_counts == expected_support_counts
assert not isolated_survivors


line_unit_prefixes = [0] * (len(sample_points_exact) + 1)
line_survivors = []
for support in line_supports:
    collision_coefficients, directions = solve_collision_exact(support)
    assert len(directions) == 1
    direction = directions[0]
    common_gcd = None
    for point_index in range(len(sample_points_exact)):
        values = []
        for parameter in range(5):
            coefficients = [
                collision_coefficients[index]
                + parameter * direction[index]
                for index in range(len(support))
            ]
            values.append(
                determinant_exact(
                    full_hessian_exact(
                        point_index, support, coefficients
                    )
                )
                - 64
            )
        equation = parameter_coefficients_exact(values)
        if len(equation) == 1 and equation[0] == 0:
            continue
        common_gcd = (
            equation
            if common_gcd is None
            else polynomial_gcd_exact(common_gcd, equation)
        )
        if len(common_gcd) == 1:
            line_unit_prefixes[point_index + 1] += 1
            break
    else:
        line_survivors.append((support, common_gcd))

assert not line_survivors


plane_parameters = sp.symbols("u v")
plane_unit_prefixes = [0] * (len(sample_points_exact) + 1)
plane_survivors = []
for support in plane_supports:
    collision_coefficients, directions = solve_collision_exact(support)
    assert len(directions) == 2
    equations = []
    for point_index in range(len(sample_points_exact)):
        entries = []
        for row in range(VARIABLE_COUNT):
            matrix_row = []
            for column in range(VARIABLE_COUNT):
                entry = base_hessian_exact[row][column]
                for coefficient_index, monomial in enumerate(support):
                    coefficient = sp.Rational(
                        collision_coefficients[
                            coefficient_index
                        ].numerator,
                        collision_coefficients[
                            coefficient_index
                        ].denominator,
                    )
                    for parameter, direction in zip(
                        plane_parameters, directions, strict=True
                    ):
                        coefficient += parameter * sp.Rational(
                            direction[coefficient_index].numerator,
                            direction[coefficient_index].denominator,
                        )
                    contribution = hessian_basis_exact[point_index][
                        monomial
                    ][row][column]
                    entry += coefficient * sp.Rational(
                        contribution.numerator,
                        contribution.denominator,
                    )
                matrix_row.append(entry)
            entries.append(matrix_row)
        equation = sp.Poly(
            sp.Matrix(entries).det(method="berkowitz") - 64,
            *plane_parameters,
            domain=sp.QQ,
        ).as_expr()
        equations.append(equation)
        basis = sp.groebner(
            equations,
            *plane_parameters,
            domain=sp.QQ,
        )
        if basis.contains(sp.Integer(1)):
            plane_unit_prefixes[point_index + 1] += 1
            break
    else:
        plane_survivors.append((support, tuple(basis.polys)))

assert not plane_survivors

print("PASS: coupled quartic--sextic collision equation reconstructed")
print(f"PASS: exact mixed support census {support_counts}")
print(
    "PASS: every isolated collision point is rejected by descending "
    f"determinant layers {isolated_layer_rejections}"
)
print(
    "PASS: all 44300 collision lines have unit exact QQ gcds; "
    f"prefix counts {line_unit_prefixes}"
)
print(
    "PASS: all 34 collision planes have unit exact QQ ideals; "
    f"prefix counts {plane_unit_prefixes}"
)
print(
    "SCOPE: genuinely mixed quartic--sextic corrections of total "
    "support at most four, with no cubic"
)

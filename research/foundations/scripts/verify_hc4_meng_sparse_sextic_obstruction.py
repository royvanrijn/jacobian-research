#!/usr/bin/env python3
"""Exclude sextic-only Meng corrections with support at most four.

Work in the collision-normalized four-variable chart

    psi_0 = 2*y*r + 4*x*s,
    p = (1, -3/2, 6, 81/8).

For a homogeneous sextic h, collision of the gradients at +/-p is

    grad(h)(p) = -H_0*p.

Writing d_e=c_e*p^e turns this into the exponent-vector system

    sum_e d_e*e = (-81/2, 18, 18, -81/2).

This checker exhausts all supports of at most four among the 84 sextic
monomials.  Isolated points are screened modulo the exact certificate prime
1,000,003.  A Hadamard bound on every coefficient and augmented collision
minor proves that reduction preserves the isolated collision solutions, and
the collision-point denominators are nonzero modulo the prime.

Principal-part cancellation is imposed before the lower determinant layers:

    det Hess(h) = 0

is necessary for det Hess(psi_0+h) to be constant.  Isolated collision
solutions are rejected by a nonzero principal evaluation or, for the
principal survivors, by a nonzero full-determinant evaluation.  All
positive-dimensional collision families are promoted to characteristic
zero.  Exact rational univariate gcds treat the 7,566 lines; the unique
plane is handled by an exact Groebner basis over QQ.  Only two lines pass
the principal layer, at parameters -81/8 and 9/2, and determinant degree
twelve rejects both.

The scope is the sextic-only chart with sextic support at most four.  It
does not combine a sextic with quartic or cubic corrections and does not
address dense sextics.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import sympy as sp


PRIME = 1_000_003
VARIABLE_COUNT = 4
SUPPORT_BOUND = 4
DEGREE = 6


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


exact_interpolation_inverse = invert_matrix_exact(
    [
        [Fraction(value**degree) for degree in range(5)]
        for value in range(5)
    ]
)


def interpolate_degree_four_exact(
    values: list[Fraction],
) -> list[Fraction]:
    coefficients = [
        sum(
            (
                exact_interpolation_inverse[row][column]
                * values[column]
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
    divisor = divisor[:]
    while len(divisor) > 1 and divisor[-1] == 0:
        divisor.pop()
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
    result = [coefficient / scale for coefficient in left]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


sextic_exponents = tuple(
    (first, second, third, DEGREE - first - second - third)
    for first in range(DEGREE + 1)
    for second in range(DEGREE + 1 - first)
    for third in range(DEGREE + 1 - first - second)
)
assert len(sextic_exponents) == 84

inverse_two = inverse_mod(2)
collision_target = (
    -81 * inverse_two % PRIME,
    18,
    18,
    -81 * inverse_two % PRIME,
)
# After multiplying the collision target by two, Hadamard's inequality
# bounds every augmented 4x4 minor by
# 6^3*sqrt(2*81^2+2*36^2), well below the certificate prime.  Coefficient
# minors are bounded by 6^4.  Thus reduction preserves every rank used in
# the collision census.
doubled_target = (-81, 36, 36, -81)
assert DEGREE**4 < PRIME
assert (
    (DEGREE**3) ** 2 * sum(entry * entry for entry in doubled_target)
    < PRIME**2
)
collision_point = (
    1,
    -3 * inverse_two % PRIME,
    6,
    81 * inverse_mod(8) % PRIME,
)
base_hessian = (
    (0, 0, 0, 4),
    (0, 0, 2, 0),
    (0, 2, 0, 0),
    (4, 0, 0, 0),
)
assert determinant_mod([list(row) for row in base_hessian]) == 64

sample_points = (
    (1, 1, 1, 1),
    (1, 2, 3, 5),
    (2, 3, 5, 7),
    (3, 5, 7, 11),
    (5, 7, 11, 13),
    (7, 11, 13, 17),
    (11, 13, 17, 19),
    (13, 17, 19, 23),
    (17, 19, 23, 29),
)
principal_point_count = 5


def monomial_value_mod(
    point: tuple[int, ...], exponents: tuple[int, ...]
) -> int:
    value = 1
    for coordinate, exponent in zip(point, exponents, strict=True):
        value = value * pow(coordinate, exponent, PRIME) % PRIME
    return value


collision_monomial_inverses = tuple(
    inverse_mod(monomial_value_mod(collision_point, exponents))
    for exponents in sextic_exponents
)


def normalized_monomial_hessian_mod(
    point: tuple[int, ...],
    exponents: tuple[int, ...],
    collision_inverse: int,
) -> tuple[tuple[int, ...], ...]:
    value = monomial_value_mod(point, exponents)
    coordinate_inverses = tuple(inverse_mod(entry) for entry in point)
    matrix = []
    for row in range(VARIABLE_COUNT):
        matrix_row = []
        for column in range(VARIABLE_COUNT):
            coefficient = exponents[row] * (
                exponents[column] - int(row == column)
            )
            matrix_row.append(
                coefficient
                * value
                * coordinate_inverses[row]
                * coordinate_inverses[column]
                * collision_inverse
                % PRIME
            )
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


sample_hessians = tuple(
    tuple(
        normalized_monomial_hessian_mod(
            point,
            exponents,
            collision_monomial_inverses[index],
        )
        for index, exponents in enumerate(sextic_exponents)
    )
    for point in sample_points
)

collision_target_exact = (
    Fraction(-81, 2),
    Fraction(18),
    Fraction(18),
    Fraction(-81, 2),
)
collision_point_exact = (
    Fraction(1),
    Fraction(-3, 2),
    Fraction(6),
    Fraction(81, 8),
)


def monomial_value_exact(
    point: tuple[Fraction, ...],
    exponents: tuple[int, ...],
) -> Fraction:
    value = Fraction(1)
    for coordinate, exponent in zip(point, exponents, strict=True):
        value *= coordinate**exponent
    return value


collision_monomial_inverses_exact = tuple(
    1 / monomial_value_exact(collision_point_exact, exponents)
    for exponents in sextic_exponents
)


def normalized_monomial_hessian_exact(
    point: tuple[Fraction, ...],
    exponents: tuple[int, ...],
    collision_inverse: Fraction,
) -> tuple[tuple[Fraction, ...], ...]:
    matrix = []
    for row in range(VARIABLE_COUNT):
        matrix_row = []
        for column in range(VARIABLE_COUNT):
            coefficient = exponents[row] * (
                exponents[column] - int(row == column)
            )
            reduced_exponents = list(exponents)
            reduced_exponents[row] -= 1
            reduced_exponents[column] -= 1
            if coefficient == 0:
                matrix_row.append(Fraction(0))
                continue
            value = Fraction(coefficient) * collision_inverse
            for coordinate, exponent in zip(
                point, reduced_exponents, strict=True
            ):
                value *= coordinate**exponent
            matrix_row.append(value)
        matrix.append(tuple(matrix_row))
    return tuple(matrix)


exact_sample_hessians = tuple(
    tuple(
        normalized_monomial_hessian_exact(
            tuple(Fraction(entry) for entry in point),
            exponents,
            collision_monomial_inverses_exact[index],
        )
        for index, exponents in enumerate(sextic_exponents)
    )
    for point in sample_points[:principal_point_count]
)


def solve_collision_support(
    support: tuple[int, ...],
) -> tuple[list[int], list[list[int]]] | None:
    variable_count = len(support)
    matrix = [
        [
            sextic_exponents[monomial][row] % PRIME
            for monomial in support
        ]
        + [collision_target[row]]
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


def solve_collision_support_exact(
    support: tuple[int, ...],
) -> tuple[list[Fraction], list[list[Fraction]]]:
    variable_count = len(support)
    matrix = [
        [
            Fraction(sextic_exponents[monomial][row])
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


def hessian_from_coefficients(
    point_index: int,
    support: tuple[int, ...],
    coefficients: list[int],
    *,
    include_base: bool,
) -> list[list[int]]:
    matrix = [
        [
            base_hessian[row][column] if include_base else 0
            for column in range(VARIABLE_COUNT)
        ]
        for row in range(VARIABLE_COUNT)
    ]
    for coefficient, monomial in zip(
        coefficients, support, strict=True
    ):
        contribution = sample_hessians[point_index][monomial]
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                matrix[row][column] = (
                    matrix[row][column]
                    + coefficient * contribution[row][column]
                ) % PRIME
    return matrix


def hessian_from_coefficients_exact(
    point_index: int,
    support: tuple[int, ...],
    coefficients: list[Fraction],
    *,
    include_base: bool,
) -> list[list[Fraction]]:
    matrix = [
        [
            (
                Fraction(base_hessian[row][column])
                if include_base
                else Fraction(0)
            )
            for column in range(VARIABLE_COUNT)
        ]
        for row in range(VARIABLE_COUNT)
    ]
    for coefficient, monomial in zip(
        coefficients, support, strict=True
    ):
        contribution = exact_sample_hessians[point_index][monomial]
        for row in range(VARIABLE_COUNT):
            for column in range(VARIABLE_COUNT):
                matrix[row][column] += (
                    coefficient * contribution[row][column]
                )
    return matrix


def line_equation_exact(
    point_index: int,
    support: tuple[int, ...],
    collision_point_coefficients: list[Fraction],
    direction: list[Fraction],
    *,
    include_base: bool,
) -> list[Fraction]:
    values = []
    for parameter in range(5):
        coefficients = [
            (
                collision_point_coefficients[index]
                + parameter * direction[index]
            )
            for index in range(len(support))
        ]
        value = determinant_exact(
            hessian_from_coefficients_exact(
                point_index,
                support,
                coefficients,
                include_base=include_base,
            )
        )
        if include_base:
            value -= 64
        values.append(value)
    return interpolate_degree_four_exact(values)


support_counts = {
    size: {"inconsistent": 0, "isolated": 0, "line": 0, "plane": 0}
    for size in range(1, SUPPORT_BOUND + 1)
}
isolated_principal_survivors = []
line_families = []
plane_families = []

for support_size in range(1, SUPPORT_BOUND + 1):
    for support in combinations(range(len(sextic_exponents)), support_size):
        solution = solve_collision_support(support)
        if solution is None:
            support_counts[support_size]["inconsistent"] += 1
            continue
        collision_coefficients, directions = solution
        if not directions:
            support_counts[support_size]["isolated"] += 1
            for point_index in range(principal_point_count):
                hessian = hessian_from_coefficients(
                    point_index,
                    support,
                    collision_coefficients,
                    include_base=False,
                )
                if determinant_mod(hessian):
                    break
            else:
                isolated_principal_survivors.append(
                    (support, collision_coefficients)
                )
            continue
        if len(directions) == 1:
            support_counts[support_size]["line"] += 1
            line_families.append(
                (support, collision_coefficients, directions[0])
            )
            continue
        assert len(directions) == 2
        support_counts[support_size]["plane"] += 1
        plane_families.append(
            (support, collision_coefficients, directions)
        )


expected_support_counts = {
    1: {"inconsistent": 84, "isolated": 0, "line": 0, "plane": 0},
    2: {"inconsistent": 3480, "isolated": 6, "line": 0, "plane": 0},
    3: {
        "inconsistent": 92936,
        "isolated": 2344,
        "line": 4,
        "plane": 0,
    },
    4: {
        "inconsistent": 198450,
        "isolated": 1723488,
        "line": 7562,
        "plane": 1,
    },
}
assert support_counts == expected_support_counts
assert len(isolated_principal_survivors) == 748
assert len(line_families) == 7566
assert len(plane_families) == 1


isolated_full_survivors = []
isolated_rejection_points = [0] * len(sample_points)
for support, coefficients in isolated_principal_survivors:
    for point_index in range(len(sample_points)):
        full_hessian = hessian_from_coefficients(
            point_index,
            support,
            coefficients,
            include_base=True,
        )
        if (determinant_mod(full_hessian) - 64) % PRIME:
            isolated_rejection_points[point_index] += 1
            break
    else:
        isolated_full_survivors.append((support, coefficients))


line_principal_survivor_supports = []
line_principal_gcds = {}
line_full_survivors = []
line_unit_points = [0] * (principal_point_count + 1)
for support, _, _ in line_families:
    collision_coefficients, directions = solve_collision_support_exact(
        support
    )
    assert len(directions) == 1
    direction = directions[0]
    principal_gcd: list[Fraction] | None = None
    for point_index in range(principal_point_count):
        equation = line_equation_exact(
            point_index,
            support,
            collision_coefficients,
            direction,
            include_base=False,
        )
        if len(equation) == 1 and equation[0] == 0:
            continue
        principal_gcd = (
            equation
            if principal_gcd is None
            else polynomial_gcd_exact(principal_gcd, equation)
        )
        if len(principal_gcd) == 1:
            break
    if principal_gcd is None or len(principal_gcd) > 1:
        line_principal_survivor_supports.append(support)
        line_principal_gcds[support] = principal_gcd
    else:
        continue

    full_gcd = principal_gcd
    for point_index in range(principal_point_count):
        equation = line_equation_exact(
            point_index,
            support,
            collision_coefficients,
            direction,
            include_base=True,
        )
        if len(equation) == 1 and equation[0] == 0:
            continue
        full_gcd = (
            equation
            if full_gcd is None
            else polynomial_gcd_exact(full_gcd, equation)
        )
        if len(full_gcd) == 1:
            line_unit_points[point_index + 1] += 1
            break
    else:
        line_full_survivors.append(
            (support, collision_coefficients, direction, full_gcd)
        )


def exact_line_layer_gcds(
    support: tuple[int, ...],
) -> dict[int, sp.Poly]:
    parameter = sp.symbols("t")
    collision_coefficients, directions = solve_collision_support_exact(
        support
    )
    assert len(directions) == 1
    spatial_variables = sp.symbols("x y r s")
    sextic = sp.Rational(0)
    for index, monomial in enumerate(support):
        scaled_coefficient = sp.Rational(
            collision_coefficients[index].numerator,
            collision_coefficients[index].denominator,
        ) + parameter * sp.Rational(
            directions[0][index].numerator,
            directions[0][index].denominator,
        )
        coefficient = scaled_coefficient / sp.Rational(
            monomial_value_exact(
                collision_point_exact,
                sextic_exponents[monomial],
            )
        )
        sextic += coefficient * sp.prod(
            variable**exponent
            for variable, exponent in zip(
                spatial_variables,
                sextic_exponents[monomial],
                strict=True,
            )
        )
    x_symbol, y_symbol, r_symbol, s_symbol = spatial_variables
    determinant = sp.Poly(
        sp.expand(
            sp.hessian(
                2 * y_symbol * r_symbol
                + 4 * x_symbol * s_symbol
                + sextic,
                spatial_variables,
            ).det(method="berkowitz")
            - 64
        ),
        *spatial_variables,
    )
    layers: dict[int, list[sp.Poly]] = {}
    for spatial_exponents, coefficient in determinant.terms():
        spatial_degree = sum(spatial_exponents)
        layers.setdefault(spatial_degree, []).append(
            sp.Poly(coefficient, parameter, domain=sp.QQ)
        )
    gcds = {}
    for degree, polynomials in layers.items():
        gcd = polynomials[0]
        for polynomial in polynomials[1:]:
            gcd = sp.gcd(gcd, polynomial)
        gcds[degree] = gcd.monic()
    return gcds


def spatial_coefficient_equations(
    expression: sp.Expr,
    spatial_variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:
    polynomial = sp.Poly(
        expression,
        *spatial_variables,
        *parameters,
        domain=sp.QQ,
    )
    grouped: dict[tuple[int, ...], sp.Expr] = {}
    for exponents, coefficient in polynomial.terms():
        spatial_exponents = exponents[: len(spatial_variables)]
        parameter_exponents = exponents[len(spatial_variables) :]
        term = sp.Rational(coefficient)
        for parameter, exponent in zip(
            parameters, parameter_exponents, strict=True
        ):
            term *= parameter**exponent
        grouped[spatial_exponents] = (
            grouped.get(spatial_exponents, 0) + term
        )
    return [
        sp.expand(coefficient)
        for coefficient in grouped.values()
        if sp.expand(coefficient) != 0
    ]


plane_principal_survivors = 0
plane_full_survivors = []
plane_parameters = sp.symbols("u v")
spatial_symbols = sp.symbols("x y r s")
for support, _, _ in plane_families:
    collision_coefficients, directions = solve_collision_support_exact(
        support
    )
    assert len(directions) == 2
    sextic = sp.Rational(0)
    for index, monomial in enumerate(support):
        scaled_coefficient = (
            sp.Rational(collision_coefficients[index].numerator)
            / collision_coefficients[index].denominator
        )
        for parameter, direction in zip(
            plane_parameters, directions, strict=True
        ):
            scaled_coefficient += parameter * (
                sp.Rational(direction[index].numerator)
                / direction[index].denominator
            )
        coefficient = scaled_coefficient / sp.Rational(
            monomial_value_exact(
                collision_point_exact,
                sextic_exponents[monomial],
            )
        )
        sextic += coefficient * sp.prod(
            variable**exponent
            for variable, exponent in zip(
                spatial_symbols,
                sextic_exponents[monomial],
                strict=True,
            )
        )
    principal_determinant = sp.expand(
        sp.hessian(sextic, spatial_symbols).det(method="berkowitz")
    )
    principal_equations = spatial_coefficient_equations(
        principal_determinant,
        spatial_symbols,
        plane_parameters,
    )
    principal_basis = sp.groebner(
        principal_equations,
        *plane_parameters,
        domain=sp.QQ,
    )
    principal_is_unit = principal_basis.contains(sp.Integer(1))
    if not principal_is_unit:
        plane_principal_survivors += 1
    if principal_is_unit:
        continue

    x_symbol, y_symbol, r_symbol, s_symbol = spatial_symbols
    base_potential = (
        2 * y_symbol * r_symbol + 4 * x_symbol * s_symbol
    )
    full_determinant = sp.expand(
        sp.hessian(
            base_potential + sextic,
            spatial_symbols,
        ).det(method="berkowitz")
        - 64
    )
    full_equations = principal_equations + spatial_coefficient_equations(
        full_determinant,
        spatial_symbols,
        plane_parameters,
    )
    full_basis = sp.groebner(
        full_equations,
        *plane_parameters,
        domain=sp.QQ,
    )
    if not full_basis.contains(sp.Integer(1)):
        plane_full_survivors.append(
            (
                support,
                collision_coefficients,
                directions,
                tuple(full_basis.polys),
            )
        )


assert not isolated_full_survivors
assert not line_full_survivors
assert not plane_full_survivors
expected_line_principal_supports = [
    (15, 21, 41, 60),
    (51, 55, 58, 64),
]
assert line_principal_survivor_supports == expected_line_principal_supports
assert line_principal_gcds[expected_line_principal_supports[0]] == [
    Fraction(81, 8),
    Fraction(1),
]
assert line_principal_gcds[expected_line_principal_supports[1]] == [
    Fraction(-9, 2),
    Fraction(1),
]
for support in expected_line_principal_supports:
    layer_gcds = exact_line_layer_gcds(support)
    assert set(layer_gcds) == {4, 8, 12, 16}
    assert layer_gcds[12].is_one
assert plane_families[0][0] == (21, 41, 55, 64)

print("PASS: normalized Meng collision equation reconstructed in degree six")
print(f"PASS: exact support census {support_counts}")
print(
    "PASS: principal Hessian screen leaves "
    f"{len(isolated_principal_survivors)} isolated points, "
    f"{len(line_principal_survivor_supports)} lines, and "
    f"{plane_principal_survivors} planes"
)
print(
    "PASS: every principal survivor is rejected by lower determinant "
    "layers; positive-dimensional families are certified over QQ"
)
print(f"DETAIL: isolated rejection points {isolated_rejection_points}")
print(f"DETAIL: line unit-prefix counts {line_unit_points}")
print(
    "DETAIL: principal line supports "
    f"{line_principal_survivor_supports}"
)
print(
    "PASS: the two principal roots t=-81/8 and t=9/2 are both "
    "rejected at determinant degree twelve"
)
print(
    "PASS: the unique collision plane is the binary cubic in yr and xs; "
    "its principal coefficient ideal is already the unit ideal over QQ"
)
print(
    "SCOPE: no sextic-only correction supported on at most four "
    "monomials retains the Meng collision with constant Hessian determinant"
)

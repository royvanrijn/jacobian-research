#!/usr/bin/env python3
"""Exclude sparse homogeneous quartic corrections in the Meng HC(4) descent.

The nonlinear toric descent checker normalizes the corrected four-variable
potential, after a polynomial base gauge, to

    psi_0 = 2*y*r + 4*x*s

and sends the two Meng--Yang points to the antipodal pair +/-p, where

    p = (1, -3/2, 6, 81/8).

Adding a polynomial h(x,y,r,s) before Schur descent is a vertical Hamiltonian
correction independent of the eliminated coordinate t, so the unit quadratic
pivot and the polynomial critical solution are unchanged.

For a homogeneous quartic h, collision of the two descended gradients is
equivalent to

    grad(h)(p) = -H_0*p,

where H_0=Hess(psi_0).  For a monomial c*w^e, multiply the i-th gradient
equation by p_i and put d=c*p^e.  The collision equations become the tiny
linear system

    sum_e d_e*e = (-81/2, 18, 18, -81/2).

This checker exhausts every support of at most four of the 35 quartic
monomials.  Collision is imposed before determinant work.  Isolated
coefficient solutions are rejected by exact reduction modulo 1,000,003 at
nine points.  The exponent minors and all denominators are smaller than the
prime, so any rational constant-determinant solution would survive the
reduction.  Rank-deficient collision supports are checked separately over
QQ: all are one-parameter families, and exact gcds of their determinant
evaluation equations are units.

The same collision solutions are then screened by the degree-eight
principal part det Hess(h_4).  Only 232 isolated quartics and two exact
members of the one-parameter families survive.  For each of these, the
checker adjoins an arbitrary scalar multiple of each one of the 20 cubic
monomials.  Modular or rational gcds of the resulting determinant evaluation
equations reject all of them.

It then exhausts all 190 supports of two cubic monomials over each of the 234
quartic principal-part survivors.  Determinant degrees seven and one give a
linear rank gate.  Rank-one systems reduce to univariate gcds.  Rank-zero
systems pass through degree-six conic linearization, conic rulings, and a
degree-four/two lift.  Only four bivariate coefficient families reach the
terminal calculation, and their full determinant equations have unit
Groebner basis modulo 1,000,003.

This is a bounded-support theorem, not a classification of quartic HC(4)
potentials.  Dense quartics, cubic supports of size at least three, quadratic
renormalizations, degree at least six, and non-coordinate coisotropic
embeddings remain open.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


PRIME = 1_000_003
VARIABLE_COUNT = 4
SUPPORT_BOUND = 4


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
            scale = work[row][column] * pivot_inverse % PRIME
            for index in range(column, VARIABLE_COUNT):
                work[row][index] = (
                    work[row][index] - scale * work[column][index]
                ) % PRIME
    return determinant % PRIME


def trim_polynomial_mod(coefficients: list[int]) -> list[int]:
    coefficients = [coefficient % PRIME for coefficient in coefficients]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_remainder_mod(
    dividend: list[int], divisor: list[int]
) -> list[int]:
    remainder = trim_polynomial_mod(dividend[:])
    divisor = trim_polynomial_mod(divisor[:])
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
        remainder = trim_polynomial_mod(remainder)
    return remainder


def polynomial_gcd_mod(
    left: list[int], right: list[int]
) -> list[int]:
    left = trim_polynomial_mod(left[:])
    right = trim_polynomial_mod(right[:])
    while not (len(right) == 1 and right[0] == 0):
        left, right = right, polynomial_remainder_mod(left, right)
    leading_inverse = inverse_mod(left[-1])
    return trim_polynomial_mod(
        [coefficient * leading_inverse % PRIME for coefficient in left]
    )


def invert_matrix_mod(matrix: list[list[int]]) -> list[list[int]]:
    size = len(matrix)
    augmented = [
        row[:]
        + [1 if row_index == column else 0 for column in range(size)]
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


# A determinant of a 4x4 matrix affine-linear in one parameter has degree at
# most four.  Interpolate it from the values at 0,1,2,3,4.
interpolation_vandermonde = [
    [pow(value, degree, PRIME) for degree in range(5)]
    for value in range(5)
]
interpolation_inverse = invert_matrix_mod(interpolation_vandermonde)
scaling_vandermonde = [
    [pow(value, degree, PRIME) for degree in range(9)]
    for value in range(9)
]
scaling_interpolation_inverse = invert_matrix_mod(scaling_vandermonde)


def interpolate_degree_four_mod(values: list[int]) -> list[int]:
    assert len(values) == 5
    return trim_polynomial_mod(
        [
            sum(
                interpolation_inverse[row][column] * values[column]
                for column in range(5)
            )
            % PRIME
            for row in range(5)
        ]
    )


def scaling_coefficient_mod(
    quartic_hessian: list[list[int]],
    cubic_hessian: list[list[int]],
    degree: int,
) -> int:
    """Coefficient of rho^degree in det(H0+rho*H3+rho^2*H4)."""

    values = []
    for scale in range(9):
        candidate = [
            [
                (
                    base_hessian_mod[row][column]
                    + scale * cubic_hessian[row][column]
                    + scale * scale * quartic_hessian[row][column]
                )
                % PRIME
                for column in range(VARIABLE_COUNT)
            ]
            for row in range(VARIABLE_COUNT)
        ]
        values.append(determinant_mod(candidate))
    return sum(
        scaling_interpolation_inverse[degree][index] * values[index]
        for index in range(9)
    ) % PRIME


def solve_affine_mod(
    rows: list[list[int]], variable_count: int
) -> tuple[str, list[int] | None, list[list[int]]]:
    """RREF an affine system; return status, a point, and nullspace basis."""

    matrix = [
        [entry % PRIME for entry in row] for row in rows
    ]
    pivot_row = 0
    pivot_columns: list[int] = []
    for column in range(variable_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(matrix))
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
        for row in range(len(matrix)):
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

    for row in range(pivot_row, len(matrix)):
        if (
            all(matrix[row][column] == 0 for column in range(variable_count))
            and matrix[row][variable_count] != 0
        ):
            return "inconsistent", None, []

    point = [0] * variable_count
    for row, column in enumerate(pivot_columns):
        point[column] = matrix[row][variable_count]
    free_columns = [
        column
        for column in range(variable_count)
        if column not in pivot_columns
    ]
    nullspace = []
    for free_column in free_columns:
        vector = [0] * variable_count
        vector[free_column] = 1
        for row, column in enumerate(pivot_columns):
            vector[column] = -matrix[row][free_column] % PRIME
        nullspace.append(vector)
    status = "unique" if not nullspace else "family"
    return status, point, nullspace


def solve_collision_mod(
    columns: tuple[tuple[int, ...], ...],
    target: tuple[int, ...],
) -> tuple[str, list[int] | None]:
    """Solve a four-row system, distinguishing unique and family solutions."""

    column_count = len(columns)
    matrix = [
        [columns[column][row] % PRIME for column in range(column_count)]
        + [target[row] % PRIME]
        for row in range(VARIABLE_COUNT)
    ]
    pivot_row = 0
    pivot_columns: list[int] = []
    for column in range(column_count):
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
                (matrix[row][index] - scale * matrix[pivot_row][index])
                % PRIME
                for index in range(column_count + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1

    for row in range(pivot_row, VARIABLE_COUNT):
        if (
            all(matrix[row][column] == 0 for column in range(column_count))
            and matrix[row][column_count] != 0
        ):
            return "inconsistent", None
    if len(pivot_columns) != column_count:
        return "family", None

    solution = [0] * column_count
    for row, column in enumerate(pivot_columns):
        solution[column] = matrix[row][column_count]
    return "unique", solution


def monomial_value_mod(
    exponents: tuple[int, ...], point: tuple[int, ...]
) -> int:
    value = 1
    for exponent, coordinate in zip(exponents, point, strict=True):
        value = value * pow(coordinate % PRIME, exponent, PRIME) % PRIME
    return value


def monomial_hessian_mod(
    exponents: tuple[int, ...], point: tuple[int, ...]
) -> list[list[int]]:
    hessian = [[0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)]
    for row in range(VARIABLE_COUNT):
        for column in range(VARIABLE_COUNT):
            reduced = list(exponents)
            if row == column:
                if exponents[row] < 2:
                    continue
                coefficient = exponents[row] * (exponents[row] - 1)
                reduced[row] -= 2
            else:
                if exponents[row] == 0 or exponents[column] == 0:
                    continue
                coefficient = exponents[row] * exponents[column]
                reduced[row] -= 1
                reduced[column] -= 1
            hessian[row][column] = (
                coefficient
                * monomial_value_mod(tuple(reduced), point)
                % PRIME
            )
    return hessian


def monomial_value_exact(
    exponents: tuple[int, ...], point: tuple[sp.Rational, ...]
) -> sp.Rational:
    return sp.prod(
        coordinate**exponent
        for exponent, coordinate in zip(exponents, point, strict=True)
    )


def monomial_hessian_exact(
    exponents: tuple[int, ...], point: tuple[sp.Rational, ...]
) -> sp.Matrix:
    entries: list[sp.Rational] = []
    for row in range(VARIABLE_COUNT):
        for column in range(VARIABLE_COUNT):
            reduced = list(exponents)
            if row == column:
                if exponents[row] < 2:
                    entries.append(sp.Rational(0))
                    continue
                coefficient = exponents[row] * (exponents[row] - 1)
                reduced[row] -= 2
            else:
                if exponents[row] == 0 or exponents[column] == 0:
                    entries.append(sp.Rational(0))
                    continue
                coefficient = exponents[row] * exponents[column]
                reduced[row] -= 1
                reduced[column] -= 1
            entries.append(
                coefficient * monomial_value_exact(tuple(reduced), point)
            )
    return sp.Matrix(VARIABLE_COUNT, VARIABLE_COUNT, entries)


quartic_exponents = tuple(
    (a, b, c, 4 - a - b - c)
    for a in range(5)
    for b in range(5 - a)
    for c in range(5 - a - b)
)
assert len(quartic_exponents) == 35

collision_point = (
    sp.Rational(1),
    -sp.Rational(3, 2),
    sp.Rational(6),
    sp.Rational(81, 8),
)
base_hessian = sp.Matrix(
    [
        [0, 0, 0, 4],
        [0, 0, 2, 0],
        [0, 2, 0, 0],
        [4, 0, 0, 0],
    ]
)
assert base_hessian.det() == 64
collision_gradient_target = -base_hessian * sp.Matrix(collision_point)
scaled_collision_target = tuple(
    sp.factor(collision_gradient_target[index] * collision_point[index])
    for index in range(VARIABLE_COUNT)
)
assert scaled_collision_target == (
    -sp.Rational(81, 2),
    sp.Rational(18),
    sp.Rational(18),
    -sp.Rational(81, 2),
)

collision_point_mod = (
    1,
    -3 * inverse_mod(2) % PRIME,
    6,
    81 * inverse_mod(8) % PRIME,
)
scaled_target_mod = tuple(
    int(value.p) * inverse_mod(int(value.q)) % PRIME
    for value in scaled_collision_target
)
collision_monomials_mod = tuple(
    monomial_value_mod(exponents, collision_point_mod)
    for exponents in quartic_exponents
)
assert all(collision_monomials_mod)

sample_points = (
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (1, 1, 1, 1),
    (1, 2, 3, 4),
    (-1, 2, -2, 1),
    (2, -1, 1, -2),
    (2, 3, 5, 7),
)
sample_hessians_mod = tuple(
    tuple(
        monomial_hessian_mod(
            exponents, tuple(coordinate % PRIME for coordinate in point)
        )
        for exponents in quartic_exponents
    )
    for point in sample_points
)
base_hessian_mod = [
    [int(base_hessian[row, column]) % PRIME for column in range(4)]
    for row in range(4)
]


support_counts: dict[int, tuple[int, int]] = {}
unique_collision_count = 0
family_supports: list[tuple[int, ...]] = []
modular_survivors: list[tuple[int, ...]] = []
principal_unique_survivors: list[
    tuple[tuple[int, ...], tuple[int, ...]]
] = []

for support_size in range(1, SUPPORT_BOUND + 1):
    size_unique = 0
    size_families = 0
    for support in combinations(range(len(quartic_exponents)), support_size):
        status, scaled_solution = solve_collision_mod(
            tuple(quartic_exponents[index] for index in support),
            scaled_target_mod,
        )
        if status == "inconsistent":
            continue
        if status == "family":
            family_supports.append(support)
            size_families += 1
            continue
        assert scaled_solution is not None
        coefficients = [
            scaled_solution[position]
            * inverse_mod(collision_monomials_mod[index])
            % PRIME
            for position, index in enumerate(support)
        ]
        size_unique += 1

        # A cubic correction cannot cancel the degree-eight term
        # det Hess(h_4).  Retain only quartics for which that principal part
        # vanishes at the four generic sample points.
        principal_survives = True
        for point_index in range(4, len(sample_points)):
            quartic_hessian = [
                [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
            ]
            for coefficient, monomial_index in zip(
                coefficients, support, strict=True
            ):
                contribution = sample_hessians_mod[point_index][
                    monomial_index
                ]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        quartic_hessian[row][column] = (
                            quartic_hessian[row][column]
                            + coefficient * contribution[row][column]
                        ) % PRIME
            if determinant_mod(quartic_hessian) != 0:
                principal_survives = False
                break
        if principal_survives:
            principal_unique_survivors.append(
                (support, tuple(coefficients))
            )

        survives = True
        for point_index in range(len(sample_points)):
            candidate_hessian = [row[:] for row in base_hessian_mod]
            for coefficient, monomial_index in zip(
                coefficients, support, strict=True
            ):
                contribution = sample_hessians_mod[point_index][
                    monomial_index
                ]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        candidate_hessian[row][column] = (
                            candidate_hessian[row][column]
                            + coefficient * contribution[row][column]
                        ) % PRIME
            if determinant_mod(candidate_hessian) != 64:
                survives = False
                break
        if survives:
            modular_survivors.append(support)
    support_counts[support_size] = (size_unique, size_families)
    unique_collision_count += size_unique

assert support_counts == {
    1: (0, 0),
    2: (3, 0),
    3: (328, 1),
    4: (42_622, 514),
}
assert unique_collision_count == 42_953
assert len(family_supports) == 515
assert not modular_survivors
assert len(principal_unique_survivors) == 232


# Every modular family is verified to be an exact rational one-parameter
# collision family.  Determinant evaluation equations at the same points
# have unit gcd in QQ[tau], so no family contains a constant determinant.
tau = sp.symbols("tau")
exact_sample_points = tuple(
    tuple(sp.Rational(coordinate) for coordinate in point)
    for point in sample_points
)
exact_sample_hessians = tuple(
    tuple(
        monomial_hessian_exact(exponents, point)
        for exponents in quartic_exponents
    )
    for point in exact_sample_points
)
collision_monomials_exact = tuple(
    monomial_value_exact(exponents, collision_point)
    for exponents in quartic_exponents
)
scaled_target_matrix = sp.Matrix(scaled_collision_target)
family_gcd_survivors: list[tuple[tuple[int, ...], sp.Expr]] = []
principal_family_survivors: list[
    tuple[tuple[int, ...], tuple[sp.Expr, ...], sp.Poly]
] = []

for support in family_supports:
    exponent_matrix = sp.Matrix.hstack(
        *(sp.Matrix(quartic_exponents[index]) for index in support)
    )
    solution_set = sp.linsolve((exponent_matrix, scaled_target_matrix))
    solution = next(iter(solution_set))
    free_symbols = sorted(
        set().union(*(entry.free_symbols for entry in solution)), key=str
    )
    assert len(free_symbols) == 1
    scaled_solution = tuple(
        sp.expand(entry.subs(free_symbols[0], tau)) for entry in solution
    )
    coefficients = tuple(
        scaled_solution[position] / collision_monomials_exact[index]
        for position, index in enumerate(support)
    )

    principal_gcd: sp.Poly | None = None
    for point_index in range(4, len(exact_sample_points)):
        quartic_hessian = sp.zeros(VARIABLE_COUNT)
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            quartic_hessian += (
                coefficient
                * exact_sample_hessians[point_index][monomial_index]
            )
        principal_equation = sp.Poly(
            sp.expand(quartic_hessian.det(method="berkowitz")),
            tau,
            domain=sp.QQ,
        )
        if principal_equation.is_zero:
            continue
        principal_gcd = (
            principal_equation
            if principal_gcd is None
            else sp.gcd(principal_gcd, principal_equation)
        )
        if principal_gcd.degree() == 0:
            break

    common_gcd: sp.Poly | None = None
    for point_index in range(len(exact_sample_points)):
        candidate_hessian = base_hessian.copy()
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            contribution = (
                coefficient
                * exact_sample_hessians[point_index][monomial_index]
            )
            candidate_hessian += contribution

        equation = sp.Poly(
            sp.expand(candidate_hessian.det(method="berkowitz") - 64),
            tau,
            domain=sp.QQ,
        )
        if equation.is_zero:
            continue
        common_gcd = (
            equation if common_gcd is None else sp.gcd(common_gcd, equation)
        )
        if common_gcd.degree() == 0:
            break
    assert common_gcd is not None
    if common_gcd.degree() > 0:
        family_gcd_survivors.append((support, common_gcd.as_expr()))
    if principal_gcd is not None and principal_gcd.degree() > 0:
        principal_family_survivors.append(
            (support, coefficients, principal_gcd.monic())
        )

assert not family_gcd_survivors
assert len(principal_family_survivors) == 2
principal_family_summary = [
    (support, sp.factor(gcd.as_expr()))
    for support, _, gcd in principal_family_survivors
]
expected_principal_family_summary = [
    ((6, 11, 20, 29), (8 * tau + 81) / 8),
    ((17, 20, 22, 25), (2 * tau - 9) / 2),
]
assert [
    support for support, _ in principal_family_summary
] == [
    support for support, _ in expected_principal_family_summary
]
assert all(
    sp.expand(actual - expected) == 0
    for (_, actual), (_, expected) in zip(
        principal_family_summary,
        expected_principal_family_summary,
        strict=True,
    )
)


# Mixed cubic--quartic extension with one cubic monomial.  The cubic does not
# affect the antipodal collision because its gradient is even.  Its Hessian
# is linear in the variables, so det Hess(h_4) is still the uncancellable
# degree-eight principal part.  Only the 232 isolated principal survivors
# and the two exact family members can proceed.
cubic_exponents = tuple(
    (a, b, c, 3 - a - b - c)
    for a in range(4)
    for b in range(4 - a)
    for c in range(4 - a - b)
)
assert len(cubic_exponents) == 20

mixed_sample_points = tuple(
    tuple(scale * coordinate for coordinate in (1, 2, 3, 4))
    for scale in range(1, 10)
) + (
    (1, 1, 1, 1),
    (-1, 2, -2, 1),
    (2, -1, 1, -2),
)
mixed_quartic_hessians_mod = tuple(
    tuple(
        monomial_hessian_mod(
            exponents, tuple(coordinate % PRIME for coordinate in point)
        )
        for exponents in quartic_exponents
    )
    for point in mixed_sample_points
)
mixed_cubic_hessians_mod = tuple(
    tuple(
        monomial_hessian_mod(
            exponents, tuple(coordinate % PRIME for coordinate in point)
        )
        for exponents in cubic_exponents
    )
    for point in mixed_sample_points
)

mixed_unique_survivors: list[
    tuple[tuple[int, ...], int, list[int]]
] = []
for support, coefficients in principal_unique_survivors:
    base_matrices: list[list[list[int]]] = []
    for point_index in range(len(mixed_sample_points)):
        matrix = [row[:] for row in base_hessian_mod]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            contribution = mixed_quartic_hessians_mod[point_index][
                monomial_index
            ]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    matrix[row][column] = (
                        matrix[row][column]
                        + coefficient * contribution[row][column]
                    ) % PRIME
        base_matrices.append(matrix)

    for cubic_index in range(len(cubic_exponents)):
        common_gcd_mod: list[int] | None = None
        for point_index in range(len(mixed_sample_points)):
            determinant_values = []
            cubic_hessian = mixed_cubic_hessians_mod[point_index][
                cubic_index
            ]
            for cubic_coefficient in range(5):
                candidate = [
                    [
                        (
                            base_matrices[point_index][row][column]
                            + cubic_coefficient
                            * cubic_hessian[row][column]
                        )
                        % PRIME
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                determinant_values.append(
                    (determinant_mod(candidate) - 64) % PRIME
                )
            equation_mod = interpolate_degree_four_mod(
                determinant_values
            )
            if len(equation_mod) == 1 and equation_mod[0] == 0:
                continue
            common_gcd_mod = (
                equation_mod
                if common_gcd_mod is None
                else polynomial_gcd_mod(common_gcd_mod, equation_mod)
            )
            if len(common_gcd_mod) == 1:
                break
        if common_gcd_mod is not None and len(common_gcd_mod) > 1:
            mixed_unique_survivors.append(
                (support, cubic_index, common_gcd_mod)
            )

assert not mixed_unique_survivors


# The two family principal roots are checked over QQ.  Their quartic Hessian
# determinants vanish identically, but no scalar multiple of a single cubic
# monomial makes the full determinant constant.
lambda_parameter = sp.symbols("lambda")
mixed_exact_points = tuple(
    tuple(sp.Rational(coordinate) for coordinate in point)
    for point in mixed_sample_points
)
mixed_quartic_hessians_exact = tuple(
    tuple(
        monomial_hessian_exact(exponents, point)
        for exponents in quartic_exponents
    )
    for point in mixed_exact_points
)
mixed_cubic_hessians_exact = tuple(
    tuple(
        monomial_hessian_exact(exponents, point)
        for exponents in cubic_exponents
    )
    for point in mixed_exact_points
)
mixed_family_survivors: list[
    tuple[tuple[int, ...], int, sp.Expr]
] = []
formal_variables = sp.symbols("w0:4")

for support, coefficient_family, principal_gcd in principal_family_survivors:
    principal_roots = sp.solve(principal_gcd.as_expr(), tau)
    assert len(principal_roots) == 1
    root = principal_roots[0]
    coefficients = tuple(
        sp.factor(coefficient.subs(tau, root))
        for coefficient in coefficient_family
    )
    quartic = sum(
        coefficient
        * sp.prod(
            formal_variables[index] ** quartic_exponents[monomial_index][
                index
            ]
            for index in range(VARIABLE_COUNT)
        )
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        )
    )
    assert sp.factor(sp.hessian(quartic, formal_variables).det()) == 0

    for cubic_index in range(len(cubic_exponents)):
        common_gcd: sp.Poly | None = None
        for point_index in range(len(mixed_exact_points)):
            candidate_hessian = base_hessian.copy()
            for coefficient, monomial_index in zip(
                coefficients, support, strict=True
            ):
                candidate_hessian += (
                    coefficient
                    * mixed_quartic_hessians_exact[point_index][
                        monomial_index
                    ]
                )
            candidate_hessian += (
                lambda_parameter
                * mixed_cubic_hessians_exact[point_index][cubic_index]
            )
            equation = sp.Poly(
                sp.expand(
                    candidate_hessian.det(method="berkowitz") - 64
                ),
                lambda_parameter,
                domain=sp.QQ,
            )
            if equation.is_zero:
                continue
            common_gcd = (
                equation
                if common_gcd is None
                else sp.gcd(common_gcd, equation)
            )
            if common_gcd.degree() == 0:
                break
        assert common_gcd is not None
        if common_gcd.degree() > 0:
            mixed_family_survivors.append(
                (support, cubic_index, common_gcd.as_expr())
            )

assert not mixed_family_survivors


# First two-cubic layer.  For h_3=lambda*m_a+mu*m_b, the degree-seven
# determinant coefficient is linear in (lambda,mu).  Rank two excludes the
# pair immediately.  Rank one fixes a unique cubic direction and leaves one
# scalar parameter; exact modular gcds of the full determinant evaluations
# then test that line.  Rank-zero pairs require the degree-six conic analysis
# and are deliberately left to the next checker.
two_cubic_points = sample_points + mixed_sample_points
two_cubic_quartic_hessians_mod = tuple(
    tuple(
        monomial_hessian_mod(
            exponents, tuple(coordinate % PRIME for coordinate in point)
        )
        for exponents in quartic_exponents
    )
    for point in two_cubic_points
)
two_cubic_cubic_hessians_mod = tuple(
    tuple(
        monomial_hessian_mod(
            exponents, tuple(coordinate % PRIME for coordinate in point)
        )
        for exponents in cubic_exponents
    )
    for point in two_cubic_points
)


def rational_mod(value: sp.Expr) -> int:
    value = sp.Rational(value)
    return int(value.p) * inverse_mod(int(value.q)) % PRIME


all_principal_quartics_mod = list(principal_unique_survivors)
for support, coefficient_family, principal_gcd in principal_family_survivors:
    roots = sp.solve(principal_gcd.as_expr(), tau)
    assert len(roots) == 1
    all_principal_quartics_mod.append(
        (
            support,
            tuple(
                rational_mod(coefficient.subs(tau, roots[0]))
                for coefficient in coefficient_family
            ),
        )
    )
assert len(all_principal_quartics_mod) == 234

degree_seven_point_indices = tuple(range(4, len(sample_points)))
two_cubic_rank_counts = {0: 0, 1: 0, 2: 0}
two_cubic_rank_zero_pairs: list[tuple[int, int, int]] = []
two_cubic_rank_one_survivors: list[
    tuple[tuple[int, ...], tuple[int, int], list[int]]
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
            contribution = two_cubic_quartic_hessians_mod[point_index][
                monomial_index
            ]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    quartic_hessian[row][column] = (
                        quartic_hessian[row][column]
                        + coefficient * contribution[row][column]
                    ) % PRIME
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

    degree_seven_signatures: list[tuple[int, ...]] = []
    for cubic_index in range(len(cubic_exponents)):
        signature = []
        for point_index in degree_seven_point_indices:
            determinant_values = []
            quartic_hessian = quartic_hessians[point_index]
            cubic_hessian = two_cubic_cubic_hessians_mod[point_index][
                cubic_index
            ]
            for cubic_coefficient in range(5):
                candidate = [
                    [
                        (
                            quartic_hessian[row][column]
                            + cubic_coefficient
                            * cubic_hessian[row][column]
                        )
                        % PRIME
                        for column in range(VARIABLE_COUNT)
                    ]
                    for row in range(VARIABLE_COUNT)
                ]
                determinant_values.append(determinant_mod(candidate))
            interpolated = interpolate_degree_four_mod(
                determinant_values
            )
            signature.append(
                interpolated[1] if len(interpolated) > 1 else 0
            )
        # Degree one is also linear in the cubic coefficients.  Append it to
        # the same signature before taking the two-column rank.
        for point_index in degree_seven_point_indices:
            signature.append(
                scaling_coefficient_mod(
                    quartic_hessians[point_index],
                    two_cubic_cubic_hessians_mod[point_index][
                        cubic_index
                    ],
                    1,
                )
            )
        degree_seven_signatures.append(tuple(signature))

    for left_cubic, right_cubic in combinations(
        range(len(cubic_exponents)), 2
    ):
        left_signature = degree_seven_signatures[left_cubic]
        right_signature = degree_seven_signatures[right_cubic]
        if any(
            (
                left_signature[left_point]
                * right_signature[right_point]
                - left_signature[right_point]
                * right_signature[left_point]
            )
            % PRIME
            for left_point in range(len(left_signature))
            for right_point in range(left_point)
        ):
            two_cubic_rank_counts[2] += 1
            continue
        if not any(left_signature) and not any(right_signature):
            two_cubic_rank_counts[0] += 1
            two_cubic_rank_zero_pairs.append(
                (quartic_index, left_cubic, right_cubic)
            )
            continue

        two_cubic_rank_counts[1] += 1
        pivot = next(
            index
            for index in range(len(left_signature))
            if left_signature[index] or right_signature[index]
        )
        left_direction = right_signature[pivot]
        right_direction = -left_signature[pivot] % PRIME

        common_gcd_mod: list[int] | None = None
        for point_index in range(len(two_cubic_points)):
            directed_cubic_hessian = [
                [
                    (
                        left_direction
                        * two_cubic_cubic_hessians_mod[point_index][
                            left_cubic
                        ][row][column]
                        + right_direction
                        * two_cubic_cubic_hessians_mod[point_index][
                            right_cubic
                        ][row][column]
                    )
                    % PRIME
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
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
            equation_mod = interpolate_degree_four_mod(
                determinant_values
            )
            if len(equation_mod) == 1 and equation_mod[0] == 0:
                continue
            common_gcd_mod = (
                equation_mod
                if common_gcd_mod is None
                else polynomial_gcd_mod(common_gcd_mod, equation_mod)
            )
            if len(common_gcd_mod) == 1:
                break
        if common_gcd_mod is not None and len(common_gcd_mod) > 1:
            two_cubic_rank_one_survivors.append(
                (
                    support,
                    (left_cubic, right_cubic),
                    common_gcd_mod,
                )
            )

assert sum(two_cubic_rank_counts.values()) == 234 * 190
assert not two_cubic_rank_one_survivors


# Rank-zero degree-seven pairs proceed to degree six.  Its dependence on the
# cubic coefficients is affine-linear in
#
#     L=lambda^2, M=lambda*mu, N=mu^2,
#
# subject to M^2=L*N.  Solve the linear equations first.  Unique points off
# the conic are impossible; unique points on it give at most two F_p
# coefficient pairs and are checked against the full determinant.
assert PRIME % 4 == 3


def square_roots_mod(value: int) -> tuple[int, ...]:
    value %= PRIME
    if value == 0:
        return (0,)
    if pow(value, (PRIME - 1) // 2, PRIME) != 1:
        return ()
    root = pow(value, (PRIME + 1) // 4, PRIME)
    return (root, -root % PRIME)


def quadratic_roots_mod(
    quadratic: int, linear: int, constant: int
) -> tuple[int, ...] | None:
    """Return F_p roots; None means the polynomial is identically zero."""

    quadratic %= PRIME
    linear %= PRIME
    constant %= PRIME
    if quadratic == 0:
        if linear == 0:
            return None if constant == 0 else ()
        return (-constant * inverse_mod(linear) % PRIME,)
    discriminant = (
        linear * linear - 4 * quadratic * constant
    ) % PRIME
    roots = square_roots_mod(discriminant)
    denominator_inverse = inverse_mod(2 * quadratic)
    return tuple(
        (-linear + root) * denominator_inverse % PRIME for root in roots
    )


def coefficient_pairs_from_quadratics(
    quadratic_values: list[int],
) -> tuple[tuple[int, int], ...]:
    left_square, product_value, right_square = quadratic_values
    if product_value * product_value % PRIME != (
        left_square * right_square % PRIME
    ):
        return ()
    if left_square:
        return tuple(
            (left, product_value * inverse_mod(left) % PRIME)
            for left in square_roots_mod(left_square)
        )
    if product_value:
        return ()
    return tuple(
        (0, right) for right in square_roots_mod(right_square)
    )


degree_six_point_indices = tuple(range(4, len(sample_points)))
degree_six_counts = {
    "inconsistent": 0,
    "off_conic": 0,
    "finite": 0,
    "family": 0,
}
degree_six_finite_survivors: list[
    tuple[int, int, int, int, int]
] = []
degree_six_families: list[
    tuple[int, int, int, list[int], list[list[int]]]
] = []

for quartic_index, left_cubic, right_cubic in two_cubic_rank_zero_pairs:
    support, coefficients = all_principal_quartics_mod[quartic_index]
    equations: list[list[int]] = []
    for point_index in degree_six_point_indices:
        quartic_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            contribution = two_cubic_quartic_hessians_mod[point_index][
                monomial_index
            ]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    quartic_hessian[row][column] = (
                        quartic_hessian[row][column]
                        + coefficient * contribution[row][column]
                    ) % PRIME
        left_hessian = two_cubic_cubic_hessians_mod[point_index][
            left_cubic
        ]
        right_hessian = two_cubic_cubic_hessians_mod[point_index][
            right_cubic
        ]
        zero_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        constant_term = scaling_coefficient_mod(
            quartic_hessian, zero_hessian, 6
        )
        left_value = scaling_coefficient_mod(
            quartic_hessian, left_hessian, 6
        )
        right_value = scaling_coefficient_mod(
            quartic_hessian, right_hessian, 6
        )
        sum_hessian = [
            [
                (
                    left_hessian[row][column]
                    + right_hessian[row][column]
                )
                % PRIME
                for column in range(VARIABLE_COUNT)
            ]
            for row in range(VARIABLE_COUNT)
        ]
        sum_value = scaling_coefficient_mod(
            quartic_hessian, sum_hessian, 6
        )
        left_square_coefficient = (
            left_value - constant_term
        ) % PRIME
        right_square_coefficient = (
            right_value - constant_term
        ) % PRIME
        product_coefficient = (
            sum_value
            - constant_term
            - left_square_coefficient
            - right_square_coefficient
        ) % PRIME
        equations.append(
            [
                left_square_coefficient,
                product_coefficient,
                right_square_coefficient,
                -constant_term % PRIME,
            ]
        )

    status, point, nullspace = solve_affine_mod(equations, 3)
    if status == "inconsistent":
        degree_six_counts["inconsistent"] += 1
        continue
    assert point is not None
    if status == "family":
        degree_six_counts["family"] += 1
        degree_six_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue

    coefficient_pairs = coefficient_pairs_from_quadratics(point)
    if not coefficient_pairs:
        degree_six_counts["off_conic"] += 1
        continue
    degree_six_counts["finite"] += 1

    for left_coefficient, right_coefficient in coefficient_pairs:
        survives = True
        for point_index in range(len(two_cubic_points)):
            candidate_hessian = [
                row[:] for row in base_hessian_mod
            ]
            for coefficient, monomial_index in zip(
                coefficients, support, strict=True
            ):
                contribution = two_cubic_quartic_hessians_mod[
                    point_index
                ][monomial_index]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        candidate_hessian[row][column] = (
                            candidate_hessian[row][column]
                            + coefficient * contribution[row][column]
                        ) % PRIME
            for cubic_coefficient, cubic_index in (
                (left_coefficient, left_cubic),
                (right_coefficient, right_cubic),
            ):
                contribution = two_cubic_cubic_hessians_mod[point_index][
                    cubic_index
                ]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        candidate_hessian[row][column] = (
                            candidate_hessian[row][column]
                            + cubic_coefficient
                            * contribution[row][column]
                        ) % PRIME
            if determinant_mod(candidate_hessian) != 64:
                survives = False
                break
        if survives:
            degree_six_finite_survivors.append(
                (
                    quartic_index,
                    left_cubic,
                    right_cubic,
                    left_coefficient,
                    right_coefficient,
                )
            )

assert sum(degree_six_counts.values()) == two_cubic_rank_counts[0]
assert not degree_six_finite_survivors


# Intersect every positive-dimensional affine solution in (L,M,N) with the
# conic M^2=L*N.  A nullity-one affine line meets it in at most two F_p
# points unless the whole line lies on the conic.  Test every finite
# intersection against the full determinant.
degree_six_conic_finite_count = 0
degree_six_conic_finite_survivors: list[
    tuple[int, int, int, int, int]
] = []
degree_six_conic_families: list[
    tuple[int, int, int, list[int], list[list[int]]]
] = []

for (
    quartic_index,
    left_cubic,
    right_cubic,
    point,
    nullspace,
) in degree_six_families:
    if len(nullspace) != 1:
        degree_six_conic_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue
    direction = nullspace[0]
    constant = (
        point[1] * point[1] - point[0] * point[2]
    ) % PRIME
    linear = (
        2 * point[1] * direction[1]
        - point[0] * direction[2]
        - direction[0] * point[2]
    ) % PRIME
    quadratic = (
        direction[1] * direction[1]
        - direction[0] * direction[2]
    ) % PRIME
    roots = quadratic_roots_mod(quadratic, linear, constant)
    if roots is None:
        degree_six_conic_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue

    support, coefficients = all_principal_quartics_mod[quartic_index]
    for parameter in roots:
        quadratic_values = [
            (point[index] + parameter * direction[index]) % PRIME
            for index in range(3)
        ]
        coefficient_pairs = coefficient_pairs_from_quadratics(
            quadratic_values
        )
        degree_six_conic_finite_count += len(coefficient_pairs)
        for left_coefficient, right_coefficient in coefficient_pairs:
            survives = True
            for point_index in range(len(two_cubic_points)):
                candidate_hessian = [
                    row[:] for row in base_hessian_mod
                ]
                for coefficient, monomial_index in zip(
                    coefficients, support, strict=True
                ):
                    contribution = two_cubic_quartic_hessians_mod[
                        point_index
                    ][monomial_index]
                    for row in range(VARIABLE_COUNT):
                        for column in range(VARIABLE_COUNT):
                            candidate_hessian[row][column] = (
                                candidate_hessian[row][column]
                                + coefficient
                                * contribution[row][column]
                            ) % PRIME
                for cubic_coefficient, cubic_index in (
                    (left_coefficient, left_cubic),
                    (right_coefficient, right_cubic),
                ):
                    contribution = two_cubic_cubic_hessians_mod[
                        point_index
                    ][cubic_index]
                    for row in range(VARIABLE_COUNT):
                        for column in range(VARIABLE_COUNT):
                            candidate_hessian[row][column] = (
                                candidate_hessian[row][column]
                                + cubic_coefficient
                                * contribution[row][column]
                            ) % PRIME
                if determinant_mod(candidate_hessian) != 64:
                    survives = False
                    break
            if survives:
                degree_six_conic_finite_survivors.append(
                    (
                        quartic_index,
                        left_cubic,
                        right_cubic,
                        left_coefficient,
                        right_coefficient,
                    )
                )

assert not degree_six_conic_finite_survivors


# A line contained in the quadratic cone M^2=L*N is a ruling through the
# vertex.  It therefore corresponds to one fixed cubic direction
# (lambda,mu)=kappa*(alpha,beta).  Recover that direction without extracting
# a square root and apply the full univariate determinant gcd.
degree_six_ruling_survivors: list[
    tuple[int, int, int, list[int]]
] = []
degree_six_nonruling_families: list[
    tuple[int, int, int, list[int], list[list[int]]]
] = []

for (
    quartic_index,
    left_cubic,
    right_cubic,
    point,
    nullspace,
) in degree_six_conic_families:
    if len(nullspace) != 1:
        degree_six_nonruling_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue
    direction = nullspace[0]
    pivot = next(
        (index for index in range(3) if direction[index]), None
    )
    assert pivot is not None
    vertex_parameter = (
        -point[pivot] * inverse_mod(direction[pivot])
    ) % PRIME
    if any(
        (
            point[index]
            + vertex_parameter * direction[index]
        )
        % PRIME
        for index in range(3)
    ):
        degree_six_nonruling_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue
    assert (
        direction[1] * direction[1]
        - direction[0] * direction[2]
    ) % PRIME == 0

    if direction[0]:
        left_direction = 1
        right_direction = (
            direction[1] * inverse_mod(direction[0])
        ) % PRIME
    else:
        assert direction[1] == 0 and direction[2] != 0
        left_direction = 0
        right_direction = 1

    support, coefficients = all_principal_quartics_mod[quartic_index]
    common_gcd_mod: list[int] | None = None
    for point_index in range(len(two_cubic_points)):
        candidate_base = [row[:] for row in base_hessian_mod]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            contribution = two_cubic_quartic_hessians_mod[point_index][
                monomial_index
            ]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    candidate_base[row][column] = (
                        candidate_base[row][column]
                        + coefficient * contribution[row][column]
                    ) % PRIME
        directed_cubic_hessian = [
            [
                (
                    left_direction
                    * two_cubic_cubic_hessians_mod[point_index][
                        left_cubic
                    ][row][column]
                    + right_direction
                    * two_cubic_cubic_hessians_mod[point_index][
                        right_cubic
                    ][row][column]
                )
                % PRIME
                for column in range(VARIABLE_COUNT)
            ]
            for row in range(VARIABLE_COUNT)
        ]
        determinant_values = []
        for line_parameter in range(5):
            candidate = [
                [
                    (
                        candidate_base[row][column]
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
        equation_mod = interpolate_degree_four_mod(
            determinant_values
        )
        if len(equation_mod) == 1 and equation_mod[0] == 0:
            continue
        common_gcd_mod = (
            equation_mod
            if common_gcd_mod is None
            else polynomial_gcd_mod(common_gcd_mod, equation_mod)
        )
        if len(common_gcd_mod) == 1:
            break
    if common_gcd_mod is not None and len(common_gcd_mod) > 1:
        degree_six_ruling_survivors.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                common_gcd_mod,
            )
        )

assert not degree_six_ruling_survivors
nonruling_nullities = [
    len(nullspace)
    for _, _, _, _, nullspace in degree_six_nonruling_families
]
nonruling_nullity_counts = {
    nullity: nonruling_nullities.count(nullity)
    for nullity in sorted(set(nonruling_nullities))
}


# Degree four on the conic.  Modulo M^2=L*N, its coefficient is linear in
#
#   L^2, L*M, L*N, M*N, N^2, L, M, N, 1.
#
# Interpolate this coefficient from nine actual cubic coefficient pairs,
# append the degree-six linear equations, and solve the lifted affine system.
degree_four_parameter_pairs = (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 0),
    (2, 1),
)


def conic_lift_row(left: int, right: int) -> list[int]:
    left_square = left * left % PRIME
    product_value = left * right % PRIME
    right_square = right * right % PRIME
    return [
        left_square * left_square % PRIME,
        left_square * product_value % PRIME,
        left_square * right_square % PRIME,
        product_value * right_square % PRIME,
        right_square * right_square % PRIME,
        left_square,
        product_value,
        right_square,
        1,
    ]


degree_four_interpolation_inverse = invert_matrix_mod(
    [
        conic_lift_row(left, right)
        for left, right in degree_four_parameter_pairs
    ]
)
degree_four_counts = {
    "inconsistent": 0,
    "off_lift": 0,
    "finite": 0,
    "family": 0,
}
degree_four_finite_survivors: list[
    tuple[int, int, int, int, int]
] = []
degree_four_families: list[
    tuple[int, int, int, list[int], list[list[int]]]
] = []

for (
    quartic_index,
    left_cubic,
    right_cubic,
    _,
    _,
) in degree_six_nonruling_families:
    support, coefficients = all_principal_quartics_mod[quartic_index]
    lifted_equations: list[list[int]] = []
    for point_index in degree_six_point_indices:
        quartic_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            contribution = two_cubic_quartic_hessians_mod[point_index][
                monomial_index
            ]
            for row in range(VARIABLE_COUNT):
                for column in range(VARIABLE_COUNT):
                    quartic_hessian[row][column] = (
                        quartic_hessian[row][column]
                        + coefficient * contribution[row][column]
                    ) % PRIME
        left_hessian = two_cubic_cubic_hessians_mod[point_index][
            left_cubic
        ]
        right_hessian = two_cubic_cubic_hessians_mod[point_index][
            right_cubic
        ]

        # Retain the degree-six affine-linear constraints in L,M,N.
        zero_hessian = [
            [0] * VARIABLE_COUNT for _ in range(VARIABLE_COUNT)
        ]
        constant_six = scaling_coefficient_mod(
            quartic_hessian, zero_hessian, 6
        )
        left_six = scaling_coefficient_mod(
            quartic_hessian, left_hessian, 6
        )
        right_six = scaling_coefficient_mod(
            quartic_hessian, right_hessian, 6
        )
        sum_hessian = [
            [
                (
                    left_hessian[row][column]
                    + right_hessian[row][column]
                )
                % PRIME
                for column in range(VARIABLE_COUNT)
            ]
            for row in range(VARIABLE_COUNT)
        ]
        sum_six = scaling_coefficient_mod(
            quartic_hessian, sum_hessian, 6
        )
        left_six_coefficient = (left_six - constant_six) % PRIME
        right_six_coefficient = (
            right_six - constant_six
        ) % PRIME
        product_six_coefficient = (
            sum_six
            - constant_six
            - left_six_coefficient
            - right_six_coefficient
        ) % PRIME
        lifted_equations.append(
            [
                0,
                0,
                0,
                0,
                0,
                left_six_coefficient,
                product_six_coefficient,
                right_six_coefficient,
                -constant_six % PRIME,
            ]
        )

        # Degree two is again affine-linear in L,M,N.
        constant_two = scaling_coefficient_mod(
            quartic_hessian, zero_hessian, 2
        )
        left_two = scaling_coefficient_mod(
            quartic_hessian, left_hessian, 2
        )
        right_two = scaling_coefficient_mod(
            quartic_hessian, right_hessian, 2
        )
        sum_two = scaling_coefficient_mod(
            quartic_hessian, sum_hessian, 2
        )
        left_two_coefficient = (left_two - constant_two) % PRIME
        right_two_coefficient = (
            right_two - constant_two
        ) % PRIME
        product_two_coefficient = (
            sum_two
            - constant_two
            - left_two_coefficient
            - right_two_coefficient
        ) % PRIME
        lifted_equations.append(
            [
                0,
                0,
                0,
                0,
                0,
                left_two_coefficient,
                product_two_coefficient,
                right_two_coefficient,
                -constant_two % PRIME,
            ]
        )

        degree_four_values = []
        for left_value, right_value in degree_four_parameter_pairs:
            cubic_hessian = [
                [
                    (
                        left_value * left_hessian[row][column]
                        + right_value
                        * right_hessian[row][column]
                    )
                    % PRIME
                    for column in range(VARIABLE_COUNT)
                ]
                for row in range(VARIABLE_COUNT)
            ]
            degree_four_values.append(
                scaling_coefficient_mod(
                    quartic_hessian, cubic_hessian, 4
                )
            )
        degree_four_coefficients = [
            sum(
                degree_four_interpolation_inverse[row][column]
                * degree_four_values[column]
                for column in range(9)
            )
            % PRIME
            for row in range(9)
        ]
        lifted_equations.append(
            degree_four_coefficients[:8]
            + [-degree_four_coefficients[8] % PRIME]
        )

    status, point, nullspace = solve_affine_mod(
        lifted_equations, 8
    )
    if status == "inconsistent":
        degree_four_counts["inconsistent"] += 1
        continue
    assert point is not None
    if status == "family":
        degree_four_counts["family"] += 1
        degree_four_families.append(
            (
                quartic_index,
                left_cubic,
                right_cubic,
                point,
                nullspace,
            )
        )
        continue

    (
        left_fourth,
        left_cube_right,
        left_square_right_square,
        left_right_cube,
        right_fourth,
        left_square,
        product_value,
        right_square,
    ) = point
    lift_is_consistent = (
        left_fourth == left_square * left_square % PRIME
        and left_cube_right == left_square * product_value % PRIME
        and left_square_right_square
        == left_square * right_square % PRIME
        and left_right_cube == product_value * right_square % PRIME
        and right_fourth == right_square * right_square % PRIME
        and product_value * product_value % PRIME
        == left_square * right_square % PRIME
    )
    if not lift_is_consistent:
        degree_four_counts["off_lift"] += 1
        continue
    coefficient_pairs = coefficient_pairs_from_quadratics(
        [left_square, product_value, right_square]
    )
    if not coefficient_pairs:
        degree_four_counts["off_lift"] += 1
        continue
    degree_four_counts["finite"] += 1

    for left_coefficient, right_coefficient in coefficient_pairs:
        survives = True
        for point_index in range(len(two_cubic_points)):
            candidate_hessian = [
                row[:] for row in base_hessian_mod
            ]
            for coefficient, monomial_index in zip(
                coefficients, support, strict=True
            ):
                contribution = two_cubic_quartic_hessians_mod[
                    point_index
                ][monomial_index]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        candidate_hessian[row][column] = (
                            candidate_hessian[row][column]
                            + coefficient * contribution[row][column]
                        ) % PRIME
            for cubic_coefficient, cubic_index in (
                (left_coefficient, left_cubic),
                (right_coefficient, right_cubic),
            ):
                contribution = two_cubic_cubic_hessians_mod[point_index][
                    cubic_index
                ]
                for row in range(VARIABLE_COUNT):
                    for column in range(VARIABLE_COUNT):
                        candidate_hessian[row][column] = (
                            candidate_hessian[row][column]
                            + cubic_coefficient
                            * contribution[row][column]
                        ) % PRIME
            if determinant_mod(candidate_hessian) != 64:
                survives = False
                break
        if survives:
            degree_four_finite_survivors.append(
                (
                    quartic_index,
                    left_cubic,
                    right_cubic,
                    left_coefficient,
                    right_coefficient,
                )
            )

assert sum(degree_four_counts.values()) == len(
    degree_six_nonruling_families
)
assert not degree_four_finite_survivors
assert len(degree_four_families) == 4


# Terminal exact modular ideals for the four surviving two-cubic families.
# At this point generic elimination is unnecessary: form the full determinant
# equations in the two actual cubic coefficients and certify the unit ideal.
left_parameter, right_parameter = sp.symbols(
    "left_parameter right_parameter"
)
terminal_two_cubic_summaries: list[
    tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]
] = []

for (
    quartic_index,
    left_cubic,
    right_cubic,
    _,
    _,
) in degree_four_families:
    support, coefficients = all_principal_quartics_mod[quartic_index]
    equations: list[sp.Expr] = []
    for point_index in range(len(two_cubic_points)):
        candidate_hessian = sp.Matrix(base_hessian_mod)
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            candidate_hessian += coefficient * sp.Matrix(
                two_cubic_quartic_hessians_mod[point_index][
                    monomial_index
                ]
            )
        candidate_hessian += left_parameter * sp.Matrix(
            two_cubic_cubic_hessians_mod[point_index][left_cubic]
        )
        candidate_hessian += right_parameter * sp.Matrix(
            two_cubic_cubic_hessians_mod[point_index][right_cubic]
        )
        equations.append(
            sp.Poly(
                sp.expand(
                    candidate_hessian.det(method="berkowitz") - 64
                ),
                left_parameter,
                right_parameter,
                modulus=PRIME,
            ).as_expr()
        )

    terminal_basis = sp.groebner(
        equations,
        left_parameter,
        right_parameter,
        modulus=PRIME,
        order="grevlex",
    )
    terminal_expressions = {
        sp.expand(polynomial.as_expr())
        for polynomial in terminal_basis.polys
    }
    assert terminal_expressions == {sp.Integer(1)}
    terminal_two_cubic_summaries.append(
        (
            tuple(quartic_exponents[index] for index in support),
            cubic_exponents[left_cubic],
            cubic_exponents[right_cubic],
        )
    )

assert len(terminal_two_cubic_summaries) == 4


print("PASS: normalized Meng collision points are the antipodal pair +/-p")
print(
    "PASS: collision-first quartic supports <=4 give 42953 isolated choices "
    "and 515 one-parameter families"
)
print(
    "PASS: exact modular determinant evaluations reject every isolated "
    "collision solution"
)
print(
    "PASS: rational determinant-equation gcds reject all 515 collision "
    "families"
)
print(
    "PASS: the degree-eight principal part leaves 232 isolated quartics and "
    "two family members; none admits a one-monomial cubic correction"
)
print(
    "PASS: for two cubic monomials, every rank-two linear odd-layer system "
    "is inconsistent and every rank-one system fails the full determinant gcd"
)
print(
    "PASS: degree-six conic linearization rejects every isolated finite "
    "rank-zero coefficient pair"
)
print(
    "PASS: every finite conic intersection of a degree-six affine family "
    "fails the full determinant test"
)
print(
    "PASS: every conic family is a ruling through the vertex and its full "
    "univariate determinant gcd is a unit"
)
print(
    "PASS: degree-four conic-lift linearization rejects every isolated "
    "finite higher-nullity coefficient pair"
)
print(
    "PASS: the four remaining two-cubic families have exact modular unit "
    "determinant ideals"
)
print(
    "SCOPE: dense quartics, cubic supports >=3, quadratic renormalizations, "
    "degree >=6, and non-coordinate coisotropic embeddings remain open"
)

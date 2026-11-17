#!/usr/bin/env python3
"""Exact local geometry of the displayed bidegree-(4,4) SIC(2) witness."""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
from fractions import Fraction
from math import factorial
from pathlib import Path

import sympy as sp
from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
SEED_SCRIPT = ROOT / "scripts" / "verify_two_pair_image_mathieu_counterexample.py"
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_local_moduli.json"
)

spec = importlib.util.spec_from_file_location("sic2_seed", SEED_SCRIPT)
assert spec and spec.loader
seed = importlib.util.module_from_spec(spec)
spec.loader.exec_module(seed)

BASIS = [
    (i, 4 - i, j, 4 - j)
    for i in range(5)
    for j in range(5)
]


def rational(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def scalar_contraction(polynomial: seed.Polynomial) -> sp.Rational:
    return rational(seed.contraction(polynomial).get(seed.ZERO, Fraction(0)))


def linear_row(f_power: seed.Polynomial) -> list[sp.Rational]:
    return [
        scalar_contraction(seed.multiply(seed.monomial(exponent), f_power))
        for exponent in BASIS
    ]


def bilinear_matrix(
    f_power: seed.Polynomial,
    f_order: int,
) -> sp.Matrix:
    """Matrix E(e_i e_j F^n), using balance rather than polynomial products."""
    result = sp.zeros(25)
    for left, left_exponent in enumerate(BASIS):
        for right in range(left, 25):
            right_exponent = BASIS[right]
            dual_one = left_exponent[0] + right_exponent[0]
            coordinate_one = left_exponent[2] + right_exponent[2]
            value = Fraction(0)
            for (xi1, _xi2, z1, _z2), coefficient in f_power.items():
                if xi1 + dual_one != z1 + coordinate_one:
                    continue
                total_one = xi1 + dual_one
                total_two = 4 * f_order - xi1 + 8 - dual_one
                value += coefficient * factorial(total_one) * factorial(total_two)
            entry = rational(value)
            result[left, right] = entry
            result[right, left] = entry
    return result


def matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, col]) for col in range(matrix.cols)] for row in range(matrix.rows)]


def primitive_quadratic(matrix: sp.Matrix) -> list[str]:
    coefficients = [
        matrix[i, j] if i == j else 2 * matrix[i, j]
        for i in range(matrix.rows)
        for j in range(i, matrix.cols)
    ]
    denominators = [int(sp.denom(value)) for value in coefficients]
    common = sp.ilcm(*denominators)
    integers = [int(value * common) for value in coefficients]
    divisor = 0
    for value in integers:
        divisor = int(sp.igcd(divisor, value))
    if divisor:
        integers = [value // divisor for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return [str(value) for value in integers]


def symbolic_multiply(
    left: dict[seed.Exponent, sp.Expr],
    right: dict[seed.Exponent, sp.Expr],
) -> dict[seed.Exponent, sp.Expr]:
    result: dict[seed.Exponent, sp.Expr] = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(4)
            )
            result[exponent] = (
                result.get(exponent, sp.Integer(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponent: sp.factor(coefficient)
        for exponent, coefficient in result.items()
        if coefficient
    }


def symbolic_contraction(
    numerator: dict[seed.Exponent, sp.Expr],
    f_power: seed.Polynomial,
) -> sp.Expr:
    result = sp.Integer(0)
    for numerator_exponent, numerator_coefficient in numerator.items():
        for f_exponent, f_coefficient in f_power.items():
            xi1 = numerator_exponent[0] + f_exponent[0]
            xi2 = numerator_exponent[1] + f_exponent[1]
            z1 = numerator_exponent[2] + f_exponent[2]
            z2 = numerator_exponent[3] + f_exponent[3]
            if xi1 != z1 or xi2 != z2:
                continue
            result += (
                numerator_coefficient
                * rational(f_coefficient)
                * factorial(xi1)
                * factorial(xi2)
            )
    return sp.factor(result)


def main() -> None:
    f, q, generators = seed.witness()

    # Exact SL2 stabilizer: b=c=0, a=d, d^2=1.
    xi1, xi2, z1, z2, aa, bb, cc, dd = sp.symbols(
        "xi1 xi2 z1 z2 aa bb cc dd"
    )
    r = xi1 * z1 + xi2 * z2
    z = xi1 * z2
    w = 2 * xi2 * z1
    t = xi1 * z1 - xi2 * z2
    f_symbolic = sp.expand(
        (r + z) * (r**2 * w - sp.Rational(1, 2) * (2 * r + z) * t**2)
    )
    transformed = sp.expand(
        f_symbolic.subs(
            {
                z1: aa * z1 + bb * z2,
                z2: cc * z1 + dd * z2,
                xi1: dd * xi1 - cc * xi2,
                xi2: -bb * xi1 + aa * xi2,
            },
            simultaneous=True,
        )
        - f_symbolic
    )
    stabilizer_equations = sp.Poly(transformed, xi1, xi2, z1, z2).coeffs()
    stabilizer_equations.append(aa * dd - bb * cc - 1)
    stabilizer_basis = sp.groebner(
        stabilizer_equations, aa, bb, cc, dd, order="grevlex"
    )
    reduced_stabilizer = {
        sp.factor(item.as_expr()) for item in stabilizer_basis.polys
    }
    assert reduced_stabilizer == {
        bb,
        cc,
        aa - dd,
        (dd - 1) * (dd + 1),
    }

    # Coefficient-operator invariant I2=tr((DC)^2).
    coefficient_matrix = sp.Matrix(
        [
            [
                rational(
                    f.get((i, 4 - i, j, 4 - j), Fraction(0))
                )
                for j in range(5)
            ]
            for i in range(5)
        ]
    )
    factorial_diagonal = sp.diag(
        *[factorial(i) * factorial(4 - i) for i in range(5)]
    )
    invariant_i2 = sp.trace((factorial_diagonal * coefficient_matrix) ** 2)
    assert invariant_i2 == 1152

    # Exact moment differentials.  The tail calculation below needs F^25.
    powers = [seed.monomial(seed.ZERO)]
    for _ in range(25):
        powers.append(seed.multiply(powers[-1], f))
    rows = [linear_row(power) for power in powers]
    first_rows = sp.Matrix(rows[:12])
    assert first_rows.rank() == 12
    assert sp.Matrix(rows).rank() == 12
    tangent_basis = sp.Matrix.hstack(*first_rows.nullspace())
    assert tangent_basis.shape == (25, 13)

    # Four tail generators, and the exact projective tail identity.
    tail = sp.Matrix(rows[12:16])
    assert tail.rank() == 4
    tail_columns = list(tail.rref()[1])
    assert tail_columns == [5, 10, 15, 20]
    tail_minor_inverse = tail[:, tail_columns].inv()
    tail_coordinates = {}
    for order in range(16, 26):
        coordinates = (
            sp.Matrix([[rows[order][index] for index in tail_columns]])
            * tail_minor_inverse
        )
        ratios = [sp.factor(coordinates[index] / coordinates[0]) for index in range(4)]
        expected = [
            1,
            -sp.Rational(order - 12, 1521520 * (order - 13)),
            sp.Rational(order - 12, 9199596806400 * (order - 14)),
            -sp.Rational(
                order - 12,
                216992729791918080000 * (order - 15),
            ),
        ]
        assert ratios == expected
        tail_coordinates[str(order)] = [str(item) for item in coordinates]

    # Lie orbit has dimension three and radial direction is independent.
    e = sp.Matrix([[0, 1], [0, 0]])
    h = sp.Matrix([[1, 0], [0, -1]])
    lower = sp.Matrix([[0, 0], [1, 0]])
    parameter = sp.symbols("parameter")

    def infinitesimal(generator: sp.Matrix) -> sp.Matrix:
        substitutions = {
            z1: z1 + parameter * (generator[0, 0] * z1 + generator[0, 1] * z2),
            z2: z2 + parameter * (generator[1, 0] * z1 + generator[1, 1] * z2),
            xi1: xi1 - parameter * (generator[0, 0] * xi1 + generator[1, 0] * xi2),
            xi2: xi2 - parameter * (generator[0, 1] * xi1 + generator[1, 1] * xi2),
        }
        derivative = sp.expand(
            sp.diff(
                f_symbolic.subs(substitutions, simultaneous=True),
                parameter,
            ).subs(parameter, 0)
        )
        polynomial = sp.Poly(derivative, xi1, xi2, z1, z2)
        return sp.Matrix(
            [
                polynomial.coeff_monomial(
                    xi1**i * xi2 ** (4 - i) * z1**j * z2 ** (4 - j)
                )
                for i in range(5)
                for j in range(5)
            ]
        )

    orbit_vectors = [infinitesimal(generator) for generator in (e, h, lower)]
    radial_vector = sp.Matrix(
        [coefficient_matrix[i, j] for i in range(5) for j in range(5)]
    )
    orbit_radial = sp.Matrix.hstack(*orbit_vectors, radial_vector)
    assert orbit_radial.rank() == 4
    assert first_rows * orbit_radial == sp.zeros(12, 4)
    orbit_radial_coordinates = tangent_basis.gauss_jordan_solve(
        orbit_radial
    )[0]
    quotient_coordinate_basis = orbit_radial_coordinates.copy()
    for index in range(13):
        candidate = quotient_coordinate_basis.row_join(sp.eye(13)[:, index])
        if candidate.rank() > quotient_coordinate_basis.cols:
            quotient_coordinate_basis = candidate
        if quotient_coordinate_basis.cols == 13:
            break
    assert quotient_coordinate_basis.det() != 0

    # Quadratic lifting obstructions after eliminating K with rows 0,...,11.
    pivot_columns = list(first_rows.rref()[1])
    pivot_inverse = first_rows[:, pivot_columns].inv()
    row_coordinates = []
    for order, row in enumerate(rows):
        coordinates = (
            sp.Matrix([[row[index] for index in pivot_columns]])
            * pivot_inverse
        )
        assert sp.Matrix([row]) == coordinates * first_rows
        row_coordinates.append(coordinates)

    q_matrices = [sp.zeros(13)]
    for order in range(1, 26):
        restricted = tangent_basis.T * bilinear_matrix(powers[order - 1], order - 1) * tangent_basis
        q_matrices.append(sp.Rational(order, 2) * restricted)

    obstruction_matrices = []
    for order in range(12, 26):
        obstruction = q_matrices[order].copy()
        for index in range(12):
            obstruction -= row_coordinates[order][index] * q_matrices[index]
        obstruction_matrices.append(sp.simplify(obstruction))
    obstruction_vectors = sp.Matrix(
        [primitive_quadratic(matrix) for matrix in obstruction_matrices[:7]]
    )
    assert obstruction_vectors.rank() == 7
    all_obstruction_vectors = sp.Matrix(
        [
            [
                matrix[i, j] if i == j else 2 * matrix[i, j]
                for i in range(13)
                for j in range(i, 13)
            ]
            for matrix in obstruction_matrices
        ]
    )
    assert all_obstruction_vectors.rank() == 7
    obstruction_pivots = list(obstruction_vectors.rref()[1])
    obstruction_minor_inverse = obstruction_vectors[:, obstruction_pivots].inv()
    later_obstruction_coordinates = {}
    for offset, vector in enumerate(all_obstruction_vectors[7:, :].tolist()):
        coordinates = (
            sp.Matrix([[vector[index] for index in obstruction_pivots]])
            * obstruction_minor_inverse
        )
        assert sp.Matrix([vector]) == coordinates * obstruction_vectors
        later_obstruction_coordinates[str(offset + 19)] = [
            str(item) for item in coordinates
        ]
    quotient_obstruction_matrices = []
    for matrix in obstruction_matrices[:7]:
        assert (
            orbit_radial_coordinates.T
            * matrix
            * orbit_radial_coordinates
        ) == sp.zeros(4)
        transformed = sp.simplify(
            quotient_coordinate_basis.T
            * matrix
            * quotient_coordinate_basis
        )
        assert transformed[:4, :] == sp.zeros(4, 13)
        assert transformed[:, :4] == sp.zeros(13, 4)
        quotient_obstruction_matrices.append(transformed[4:, 4:])
    quotient_obstruction_vectors = sp.Matrix(
        [
            primitive_quadratic(matrix)
            for matrix in quotient_obstruction_matrices
        ]
    )
    assert quotient_obstruction_vectors.rank() == 7

    # The seven quotient quadrics form a degree-three thickening of one
    # five-plane.  Singular supplies an exact characteristic-zero radical,
    # not a numerical or finite-field decomposition.
    quotient_variables = sp.symbols("v0:9")
    quotient_monomials = [
        quotient_variables[i] * quotient_variables[j]
        for i in range(9)
        for j in range(i, 9)
    ]
    quotient_quadrics = [
        sp.expand(
            sum(
                int(quotient_obstruction_vectors[row, column])
                * quotient_monomials[column]
                for column in range(len(quotient_monomials))
            )
        )
        for row in range(7)
    ]
    radical_linear_forms = [
        quotient_variables[0],
        (
            35 * quotient_variables[1]
            - 8 * quotient_variables[4]
            - 70 * quotient_variables[6]
            + 105 * quotient_variables[7]
            - 105 * quotient_variables[8]
        ),
        (
            28 * quotient_variables[2]
            - 43 * quotient_variables[4]
            + 168 * quotient_variables[5]
            - 336 * quotient_variables[6]
            + 336 * quotient_variables[7]
            - 336 * quotient_variables[8]
        ),
        (
            105 * quotient_variables[3]
            - 251 * quotient_variables[4]
            + 840 * quotient_variables[5]
            - 1260 * quotient_variables[6]
            + 1260 * quotient_variables[7]
            - 1365 * quotient_variables[8]
        ),
    ]
    for quadratic in quotient_quadrics:
        assert sp.expand(
            quadratic.subs(
                {
                    quotient_variables[0]: 0,
                    quotient_variables[1]: (
                        8 * quotient_variables[4]
                        + 70 * quotient_variables[6]
                        - 105 * quotient_variables[7]
                        + 105 * quotient_variables[8]
                    )
                    / 35,
                    quotient_variables[2]: (
                        43 * quotient_variables[4]
                        - 168 * quotient_variables[5]
                        + 336 * quotient_variables[6]
                        - 336 * quotient_variables[7]
                        + 336 * quotient_variables[8]
                    )
                    / 28,
                    quotient_variables[3]: (
                        251 * quotient_variables[4]
                        - 840 * quotient_variables[5]
                        + 1260 * quotient_variables[6]
                        - 1260 * quotient_variables[7]
                        + 1365 * quotient_variables[8]
                    )
                    / 105,
                },
                simultaneous=True,
            )
        ) == 0
    singular_polynomials = [
        str(polynomial).replace("**", "^")
        for polynomial in quotient_quadrics
    ]
    singular_linear_forms = [
        str(polynomial).replace("**", "^")
        for polynomial in radical_linear_forms
    ]
    singular_script = "\n".join(
        [
            'LIB "primdec.lib";',
            "ring r=0,(v0,v1,v2,v3,v4,v5,v6,v7,v8),dp;",
            f"ideal I={','.join(singular_polynomials)};",
            f"ideal L={','.join(singular_linear_forms)};",
            "ideal G=std(I);",
            "ideal R=radical(I);",
            'if (dim(G)==5 && deg(G)==3 && '
            "size(reduce(R,std(L)))==0 && "
            "size(reduce(L,std(R)))==0) "
            '{ print("RADICAL_OK"); }',
        ]
    )
    singular_result = subprocess.run(
        ["/opt/homebrew/bin/Singular", "-q"],
        input=singular_script,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    assert "RADICAL_OK" in singular_result.stdout

    # Parametrize the reduced five-plane by h0,...,h4.  The following
    # exact calculation constructs a quadratic tangent correction that
    # kills every cubic lifting obstruction.
    residual_parameters = sp.symbols("h0:5")
    residual_parametrization = sp.zeros(9, 5)
    for index in range(5):
        residual_parametrization[4 + index, index] = 1
    residual_parametrization[1, :] = sp.Matrix(
        [[sp.Rational(8, 35), 0, 2, -3, 3]]
    )
    residual_parametrization[2, :] = sp.Matrix(
        [[sp.Rational(43, 28), -6, 12, -12, 12]]
    )
    residual_parametrization[3, :] = sp.Matrix(
        [[sp.Rational(251, 105), -8, 12, -12, 13]]
    )
    residual_vectors = (
        tangent_basis
        * quotient_coordinate_basis[:, 4:]
        * residual_parametrization
    )
    residual_h = residual_vectors * sp.Matrix(residual_parameters)

    # The apolar adjoint, composed with diag(1,-1), fixes F.  It explains
    # eight of the ten all-moment quotient-tangent directions and splits
    # the reduced quadratic plane into one even and four odd directions
    # after quotienting by the orbit and radial line.
    apolar = sp.zeros(5)
    for index, value in enumerate(
        (1, -sp.Rational(1, 4), sp.Rational(1, 6), -sp.Rational(1, 4), 1)
    ):
        apolar[index, 4 - index] = value
    parity = sp.diag(1, -1, 1, -1, 1)
    local_involution = sp.zeros(25)
    for basis_index in range(25):
        basis_coefficient = sp.zeros(5)
        basis_coefficient[basis_index // 5, basis_index % 5] = 1
        basis_operator = basis_coefficient.T * factorial_diagonal
        transformed_operator = (
            parity
            * apolar.inv()
            * basis_operator.T
            * apolar
            * parity
        )
        transformed_coefficient = (
            factorial_diagonal.inv() * transformed_operator.T
        )
        local_involution[:, basis_index] = sp.Matrix(
            [
                transformed_coefficient[row, column]
                for row in range(5)
                for column in range(5)
            ]
        )
    assert local_involution**2 == sp.eye(25)
    assert local_involution * radial_vector == radial_vector
    assert first_rows * local_involution == first_rows
    plus_space = sp.Matrix.hstack(
        *(local_involution - sp.eye(25)).nullspace()
    )
    minus_space = sp.Matrix.hstack(
        *(local_involution + sp.eye(25)).nullspace()
    )
    assert (plus_space.cols, minus_space.cols) == (15, 10)
    orbit_plus_rank = ((sp.eye(25) + local_involution) * orbit_radial[:, :3]).rank()
    orbit_minus_rank = ((sp.eye(25) - local_involution) * orbit_radial[:, :3]).rank()
    assert (orbit_plus_rank, orbit_minus_rank) == (1, 2)
    quotient_fiber_eigenspaces = (
        plus_space.cols - (first_rows * plus_space).rank() - orbit_plus_rank,
        minus_space.cols - orbit_minus_rank,
    )
    assert quotient_fiber_eigenspaces == (2, 8)
    assert sp.Matrix.hstack(
        orbit_radial, residual_vectors, local_involution * residual_vectors
    ).rank() == sp.Matrix.hstack(orbit_radial, residual_vectors).rank()
    reduced_plane_eigenspaces = (
        sp.Matrix.hstack(
            orbit_radial,
            (sp.eye(25) + local_involution) * residual_vectors,
        ).rank()
        - orbit_radial.rank(),
        sp.Matrix.hstack(
            orbit_radial,
            (sp.eye(25) - local_involution) * residual_vectors,
        ).rank()
        - orbit_radial.rank(),
    )
    assert reduced_plane_eigenspaces == (1, 4)

    residual_q = [sp.Integer(0)]
    for order in range(1, 12):
        residual_q.append(
            sp.Rational(order, 2)
            * (
                residual_h.T
                * bilinear_matrix(powers[order - 1], order - 1)
                * residual_h
            )[0]
        )
    second_correction = sp.zeros(25, 1)
    second_pivot_values = -pivot_inverse * sp.Matrix(residual_q)
    for index, pivot in enumerate(pivot_columns):
        second_correction[pivot] = second_pivot_values[index]

    residual_dictionaries = [
        {
            BASIS[row]: residual_vectors[row, column]
            for row in range(25)
            if residual_vectors[row, column]
        }
        for column in range(5)
    ]
    residual_triples: dict[tuple[int, int, int], dict[seed.Exponent, sp.Expr]] = {}
    for first in range(5):
        for second in range(first, 5):
            for third in range(second, 5):
                residual_triples[first, second, third] = symbolic_multiply(
                    symbolic_multiply(
                        residual_dictionaries[first],
                        residual_dictionaries[second],
                    ),
                    residual_dictionaries[third],
                )

    residual_cubics = []
    for f_order in range(24):
        cubic = sp.Integer(0)
        for (first, second, third), numerator in residual_triples.items():
            multiplicity = (
                1
                if first == second == third
                else (3 if first == second or second == third else 6)
            )
            cubic += (
                multiplicity
                * symbolic_contraction(numerator, powers[f_order])
                * residual_parameters[first]
                * residual_parameters[second]
                * residual_parameters[third]
            )
        residual_cubics.append(sp.expand(cubic))

    tangent_lift_variables = sp.symbols("s0:13")
    lifted_second_correction = (
        second_correction
        + tangent_basis * sp.Matrix(tangent_lift_variables)
    )
    third_known = []
    for order in range(26):
        if order == 0:
            third_known.append(sp.Integer(0))
            continue
        value = order * (
            residual_h.T
            * bilinear_matrix(powers[order - 1], order - 1)
            * lifted_second_correction
        )[0]
        if order >= 2:
            value += (
                sp.Rational(order * (order - 1), 6)
                * residual_cubics[order - 2]
            )
        third_known.append(sp.expand(value))

    third_obstructions = []
    for order in range(12, 26):
        obstruction = third_known[order]
        for index in range(12):
            obstruction -= (
                row_coordinates[order][index] * third_known[index]
            )
        third_obstructions.append(sp.expand(obstruction))

    third_matrix = sp.Matrix(
        [
            [
                obstruction.coeff(variable)
                for variable in tangent_lift_variables
            ]
            for obstruction in third_obstructions
        ]
    )
    third_constant = sp.Matrix(
        [
            obstruction.subs(
                {variable: 0 for variable in tangent_lift_variables}
            )
            for obstruction in third_obstructions
        ]
    )
    quadratic_parameter_monomials = [
        residual_parameters[first] * residual_parameters[second]
        for first in range(5)
        for second in range(first, 5)
    ]
    correction_coefficients = sp.symbols(
        f"a0:{13 * len(quadratic_parameter_monomials)}"
    )
    correction_ansatz = sp.Matrix(
        [
            sum(
                correction_coefficients[
                    row * len(quadratic_parameter_monomials) + column
                ]
                * monomial
                for column, monomial in enumerate(
                    quadratic_parameter_monomials
                )
            )
            for row in range(13)
        ]
    )
    cubic_parameter_monomials = [
        residual_parameters[first]
        * residual_parameters[second]
        * residual_parameters[third]
        for first in range(5)
        for second in range(first, 5)
        for third in range(second, 5)
    ]
    correction_equations = []
    for polynomial in (
        third_matrix * correction_ansatz + third_constant
    ):
        polynomial_object = sp.Poly(polynomial, *residual_parameters)
        correction_equations.extend(
            polynomial_object.coeff_monomial(monomial)
            for monomial in cubic_parameter_monomials
        )
    correction_matrix, correction_right = sp.linear_eq_to_matrix(
        correction_equations,
        correction_coefficients,
    )
    augmented_correction = correction_matrix.row_join(correction_right)
    flint_augmented = fmpq_mat(
        [
            [
                str(augmented_correction[row, column])
                for column in range(augmented_correction.cols)
            ]
            for row in range(augmented_correction.rows)
        ]
    )
    correction_rref, correction_rank = flint_augmented.rref()
    correction_pivots = []
    for row in range(correction_rank):
        pivot = next(
            column
            for column in range(correction_matrix.cols)
            if correction_rref[row, column]
        )
        correction_pivots.append(pivot)
    correction_solution = [
        sp.Integer(0) for _ in range(correction_matrix.cols)
    ]
    for row, pivot in enumerate(correction_pivots):
        correction_solution[pivot] = sp.Rational(
            str(correction_rref[row, correction_matrix.cols])
        )
    polynomial_third_correction = sp.Matrix(
        [
            sum(
                correction_solution[
                    row * len(quadratic_parameter_monomials) + column
                ]
                * monomial
                for column, monomial in enumerate(
                    quadratic_parameter_monomials
                )
            )
            for row in range(13)
        ]
    )
    assert all(
        sp.expand(polynomial) == 0
        for polynomial in (
            third_matrix * polynomial_third_correction
            + third_constant
        )
    )

    # The exact two-parameter homogeneous family.
    a, b = sp.symbols("a b")
    family = sp.expand(
        sp.Rational(1, 2)
        * (a * r + b * z)
        * (
            2 * w * (a * r + b * z) ** 2
            - 2 * a * b * r**3
            - b**2 * r**2 * z
        )
    )
    assert sp.expand(family.subs({a: 1, b: 1}) - f_symbolic) == 0
    family_matrix = sp.Matrix(
        [
            [
                sp.Poly(family, xi1, xi2, z1, z2).coeff_monomial(
                    xi1**i * xi2 ** (4 - i) * z1**j * z2 ** (4 - j)
                )
                for j in range(5)
            ]
            for i in range(5)
        ]
    )
    assert sp.factor(sp.trace((factorial_diagonal * family_matrix) ** 2)) == (
        1152 * a**4 * b**2
    )

    # Exact sparse specializations replay the all-order formulas at sample
    # parameters; the written coefficient-extraction proof is symbolic in m.
    for a_value, b_value in [(1, 1), (2, 3), (-1, 2)]:
        specialized = {
            exponent: Fraction(
                sp.Poly(family.subs({a: a_value, b: b_value}), xi1, xi2, z1, z2)
                .coeff_monomial(
                    xi1**exponent[0]
                    * xi2**exponent[1]
                    * z1**exponent[2]
                    * z2**exponent[3]
                )
            )
            for exponent in BASIS
        }
        specialized = {key: value for key, value in specialized.items() if value}
        specialized_power = seed.monomial(seed.ZERO)
        for order in range(1, 7):
            specialized_power = seed.multiply(specialized_power, specialized)
            assert seed.contraction(specialized_power) == {}
            expected = Fraction(
                factorial(4 * order + 2)
                * a_value ** (2 * order + 1)
                * b_value ** (order - 1)
                * factorial(order),
                seed.double_factorial_odd(2 * order + 1),
            )
            assert seed.contraction(seed.multiply(q, specialized_power)) == {
                seed.ZERO: expected
            }

    artifact = {
        "format": "two-pair-counterexample-local-moduli-v1",
        "basis_order": [
            f"xi1^{i}*xi2^{4-i}*z1^{j}*z2^{4-j}"
            for i in range(5)
            for j in range(5)
        ],
        "effective_group": "PGL2",
        "sl2_stabilizer_groebner_basis": [
            "bb", "cc", "aa-dd", "(dd-1)*(dd+1)"
        ],
        "effective_stabilizer": "trivial",
        "orbit_dimension": 3,
        "I2_F": str(invariant_i2),
        "all_moment_tangent_dimension": 13,
        "linear_functionals_orders_0_through_11": matrix_strings(first_rows),
        "tangent_basis_columns": matrix_strings(tangent_basis),
        "tangent_coordinate_change_orbit_radial_first": matrix_strings(
            quotient_coordinate_basis
        ),
        "tail_identity_coordinates_orders_16_through_25": tail_coordinates,
        "quadratic_obstruction_orders": list(range(12, 19)),
        "quadratic_coordinate_monomial_order": [
            f"u{i}*u{j}" for i in range(13) for j in range(i, 13)
        ],
        "primitive_quadratic_obstruction_coefficients": matrix_strings(
            obstruction_vectors
        ),
        "quotient_quadratic_coordinate_monomial_order": [
            f"v{i}*v{j}" for i in range(9) for j in range(i, 9)
        ],
        "primitive_quotient_quadratic_obstruction_coefficients": (
            matrix_strings(quotient_obstruction_vectors)
        ),
        "quotient_quadratic_scheme": {
            "dimension": 5,
            "degree": 3,
            "apolar_adjoint_eigenspace_dimensions": {
                "even": reduced_plane_eigenspaces[0],
                "odd": reduced_plane_eigenspaces[1],
            },
            "radical_linear_forms": [
                str(polynomial) for polynomial in radical_linear_forms
            ],
        },
        "all_moment_fiber_quotient_tangent_eigenspaces": {
            "even": quotient_fiber_eigenspaces[0],
            "odd": quotient_fiber_eigenspaces[1],
        },
        "reduced_quadratic_cone_parameter_order": [
            str(parameter) for parameter in residual_parameters
        ],
        "polynomial_tangent_correction_for_third_order": [
            str(sp.factor(polynomial))
            for polynomial in polynomial_third_correction
        ],
        "third_order_correction_linear_system": {
            "equation_count": correction_matrix.rows,
            "unknown_coefficient_count": correction_matrix.cols,
            "rank": correction_rank,
            "replay_obstruction_orders": list(range(12, 26)),
        },
        "later_obstruction_coordinates_orders_19_through_25": (
            later_obstruction_coordinates
        ),
        "quadratic_obstruction_rank": 7,
        "family": (
            "F_ab=(a*R+b*Z)/2*"
            "(2*W*(a*R+b*Z)^2-2*a*b*R^3-b^2*R^2*Z)"
        ),
        "family_mixed_defect": (
            "(4*m+2)!*a^(2*m+1)*b^(m-1)*m!/(2*m+1)!!"
        ),
        "family_orbit_modulus_split": (
            "F_ab=a^2*b*(diag(t,t^-1).F), t^-2=b/a"
        ),
        "local_quotient_dimension_lower_bound": 1,
        "written_proof": (
            "extended-geometry/TWO_PAIR_COUNTEREXAMPLE_LOCAL_MODULI.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS local SIC2C4: effective stabilizer is trivial and orbit dimension is 3")
    print("PASS local SIC2C4: all-order moment tangent has dimension 13")
    print("PASS local SIC2C4: seven independent quadratic lifting obstructions")
    print(
        "PASS local SIC2C4: reduced quadratic cone is a five-plane "
        "with apolar split 1+4"
    )
    print("PASS local SIC2C4: every reduced direction lifts through third order")
    print("PASS local SIC2C4: F_ab is an all-order defect-preserving family")
    print(f"PASS local SIC2C4: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

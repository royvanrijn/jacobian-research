#!/usr/bin/env python3
"""Search nonlinear parent-preserving charts for direct unit Schur blocks.

HC4MCP6 found 54 reverse-order cubic--quadratic mixed Hamiltonian support
patterns with an exact constant-Hessian parent line.  On that line, the
cubic coefficient ``b`` is arbitrary nonzero and the quadratic coefficient
``a`` is one of ``+/-1/2`` or ``+/-1/4``.  This checker specializes

    b in {-2, -1, 1, 2}

and searches oblique constant directions in

    ({-1, 0, 1}^6 - {0}) / {+/-1}.

For a direction c, a direct scalar unit pivot means

    D_c^2 Phi in Q^x.

For a two-plane <c,d>, a polynomial 2-by-2 Schur block means that Phi is
jointly quadratic along the plane,

    D_c^3 Phi = D_c^2 D_d Phi = D_c D_d^2 Phi = D_d^3 Phi = 0,

and the restricted Hessian determinant

    (D_c^2 Phi)(D_d^2 Phi) - (D_c D_d Phi)^2

is a nonzero rational constant.  This is the exact condition that the two
critical equations are affine in the two pivot variables with a unimodular
polynomial coefficient matrix.  Its inverse is polynomial, so exact Schur
elimination would give a four-variable constant-Hessian potential and
would retain the transformed marked critical collision.

The default mode is the original finite-box search.  With
``--symbolic-classification``, the checker instead treats ``b`` and all
direction coordinates as indeterminates.  Scalar directions are handled by
an exact saturated ideal.  Two-planes are covered by the 15 standard affine
Pluecker charts on Gr(2,6); exact joint-cubic equations and determinant
equalities at fixed rational probe points give unit-ideal exclusion
certificates over Q.  A unit sampled ideal is sufficient because a constant
determinant must agree at every probe point.  Any sampled survivor falls back
to the full coefficient ideal.

Neither mode classifies longer symplectic words, different Hamiltonian
supports, coefficient-dependent nonlinear directions, or general polynomial
symplectomorphisms.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations, product
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

import search_hc4_mixed_canonical_pivots as base
import verify_hc4_symbolic_quadratic_cubic_words as symbolic_words


PRIME = 1_000_003
B_BOX = (-2, -1, 1, 2)
DIRECTION_VALUES = (-1, 0, 1)
SYMBOLIC_DIRECTION_VARIABLES = sp.symbols("c0:6")
GRASSMANN_VARIABLES = sp.symbols("r0:8")
SATURATION_VARIABLE = sp.Symbol("z_saturation")


@dataclass(frozen=True)
class PolynomialDerivatives:
    hessian: tuple[tuple[sp.Poly, ...], ...]
    third: tuple[tuple[tuple[sp.Poly, ...], ...], ...]


def exceptional_a_value(
    pattern: symbolic_words.Pattern,
) -> sp.Rational:
    """Return the exact parent-preserving coefficient from HC4MCP6."""

    if pattern.incidence == "reciprocal":
        value = sp.Rational(-pattern.epsilon1, 4)
    elif pattern.incidence == "h1_source_hits_h2_dual":
        value = sp.Rational(-pattern.epsilon1, 2)
    else:
        assert pattern.incidence == "h1_dual_hits_h2_source"
        value = sp.Rational(pattern.epsilon2, 2)
    return value


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def primitive_integer_expression(
    expression: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> sp.Expr:
    """Normalize a rational polynomial without changing its zero set."""

    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    _, cleared = polynomial.clear_denoms(convert=True)
    primitive = cleared.primitive()[1]
    if primitive.LC() < 0:
        primitive = -primitive
    return sp.expand(primitive.as_expr())


def normalized_equations(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Remove zero, duplicate, and rationally associated generators."""

    unique: dict[str, sp.Expr] = {}
    for equation in equations:
        if equation == 0:
            continue
        normalized = primitive_integer_expression(equation, variables)
        encoded = sp.sstr(normalized)
        unique.setdefault(encoded, normalized)
    return tuple(unique.values())


def run_singular_program(program: str, *, timeout: int) -> dict[str, str]:
    completed = subprocess.run(
        ["Singular", "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0 or any(
        line.lstrip().startswith("?")
        for line in completed.stdout.splitlines()
    ):
        raise RuntimeError(
            "Singular failed:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return {
        key: value
        for line in completed.stdout.splitlines()
        if "=" in line
        for key, value in [line.split("=", 1)]
    }


def singular_localized_unit_ideal(
    equations: Sequence[sp.Expr],
    *,
    parameter_variables: Sequence[sp.Symbol],
    localized_at: sp.Expr,
    timeout: int,
) -> dict[str, object]:
    """Check emptiness after localizing at ``localized_at`` over Q."""

    equations = normalized_equations(equations, parameter_variables)
    ring_variables = (SATURATION_VARIABLE, *parameter_variables)
    saturation_equation = primitive_integer_expression(
        SATURATION_VARIABLE * localized_at - 1,
        ring_variables,
    )
    lines = [
        "ring R=0,("
        + ",".join(str(variable) for variable in ring_variables)
        + "),dp;",
        "ideal I=0;",
    ]
    lines.extend(
        f"I=I,{singular_expression(equation)};"
        for equation in equations
    )
    lines.extend(
        [
            f"I=I,{singular_expression(saturation_equation)};",
            "ideal G=std(I);",
            'print("BASIS_SIZE="+string(size(G)));',
            'if (reduce(1,G)==0) { print("UNIT=1"); }'
            ' else { print("UNIT=0"); }',
        ]
    )
    program = "\n".join(lines)
    markers = run_singular_program(program, timeout=timeout)
    return {
        "unit_ideal": markers.get("UNIT") == "1",
        "basis_size": int(markers["BASIS_SIZE"]),
        "generators": len(equations),
        "singular_input_sha256": hashlib.sha256(
            program.encode()
        ).hexdigest(),
    }


def scalar_saturated_certificate(
    hessian: sp.Matrix,
    *,
    timeout: int,
) -> dict[str, object]:
    """Prove that a constant symbolic scalar pivot must be zero."""

    direction = sp.Matrix(SYMBOLIC_DIRECTION_VARIABLES)
    second = sp.expand((direction.T * hessian * direction)[0])
    polynomial = sp.Poly(second, *base.variables)
    constant_monomial = (0,) * len(base.variables)
    constant_coefficient = polynomial.coeff_monomial(constant_monomial)
    nonconstant_coefficients = [
        coefficient
        for monomial, coefficient in polynomial.terms()
        if monomial != constant_monomial
    ]
    parameter_variables = (
        symbolic_words.b,
        *SYMBOLIC_DIRECTION_VARIABLES,
    )
    equations = normalized_equations(
        nonconstant_coefficients, parameter_variables
    )
    normalized_constant = primitive_integer_expression(
        constant_coefficient, parameter_variables
    )
    lines = [
        'LIB "elim.lib";',
        "ring R=0,("
        + ",".join(str(variable) for variable in parameter_variables)
        + "),dp;",
        "ideal I=0;",
    ]
    lines.extend(
        f"I=I,{singular_expression(equation)};"
        for equation in equations
    )
    lines.extend(
        [
            f"ideal B={symbolic_words.b};",
            "ideal G=std(sat(I,B));",
            f"poly Q={singular_expression(normalized_constant)};",
            'print("BASIS_SIZE="+string(size(G)));',
            'print("DIMENSION="+string(dim(G)));',
            'if (reduce(Q,G)==0) { print("CONSTANT_ZERO=1"); }'
            ' else { print("CONSTANT_ZERO=0"); }',
        ]
    )
    program = "\n".join(lines)
    markers = run_singular_program(program, timeout=timeout)
    return {
        "nonconstant_coefficient_generators": len(equations),
        "saturated_basis_size": int(markers["BASIS_SIZE"]),
        "saturated_scheme_dimension": int(markers["DIMENSION"]),
        "constant_coefficient_normal_form_zero": (
            markers.get("CONSTANT_ZERO") == "1"
        ),
        "singular_input_sha256": hashlib.sha256(
            program.encode()
        ).hexdigest(),
    }


def symbolic_schur_probe_points() -> tuple[tuple[int, ...], ...]:
    """Fixed exact probes for the symbolic determinant certificates."""

    axes = []
    for index in range(6):
        for sign in (1, -1):
            point = [0] * 6
            point[index] = sign
            axes.append(tuple(point))
    pairs = []
    for left, right in combinations(range(6), 2):
        point = [0] * 6
        point[left] = 1
        point[right] = 1
        pairs.append(tuple(point))
    sign_corners = tuple(product((-1, 1), repeat=6))
    magnitude_corners = tuple(product((1, 2), repeat=6))
    return tuple((*axes, *pairs, *sign_corners, *magnitude_corners))


def grassmann_chart_basis(
    pivots: tuple[int, int],
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the normalized row basis on one standard Gr(2,6) chart."""

    left: list[sp.Expr] = []
    right: list[sp.Expr] = []
    parameter_index = 0
    for coordinate in range(6):
        if coordinate == pivots[0]:
            left.append(sp.Integer(1))
            right.append(sp.Integer(0))
        elif coordinate == pivots[1]:
            left.append(sp.Integer(0))
            right.append(sp.Integer(1))
        else:
            left.append(GRASSMANN_VARIABLES[parameter_index])
            right.append(GRASSMANN_VARIABLES[parameter_index + 1])
            parameter_index += 2
    assert parameter_index == len(GRASSMANN_VARIABLES)
    return sp.Matrix(left), sp.Matrix(right)


def symbolic_quadratic_form(
    hessian: Sequence[Sequence[sp.Poly]],
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
) -> sp.Poly:
    result = sp.Poly(0, *base.variables, domain=sp.EX)
    for row in range(6):
        for column in range(6):
            coefficient = left[row] * right[column]
            if coefficient != 0:
                result += hessian[row][column].mul_ground(coefficient)
    return result


def polynomial_directional_derivative(
    polynomial: sp.Poly,
    direction: Sequence[sp.Expr],
) -> sp.Poly:
    result = sp.Poly(0, *base.variables, domain=sp.EX)
    for index, variable in enumerate(base.variables):
        if direction[index] != 0:
            result += polynomial.diff(variable).mul_ground(
                direction[index]
            )
    return result


def evaluated_hessian(
    hessian: Sequence[Sequence[sp.Poly]],
    point: Sequence[int],
) -> sp.Matrix:
    rows = []
    for row in hessian:
        evaluated_row = []
        for entry in row:
            evaluated = entry
            for variable, value in zip(
                base.variables, point, strict=True
            ):
                evaluated = evaluated.eval(variable, value)
            evaluated_row.append(evaluated.as_expr())
        rows.append(evaluated_row)
    return sp.Matrix(rows)


def grassmann_chart_certificate(
    spatial_hessian: Sequence[Sequence[sp.Poly]],
    pivots: tuple[int, int],
    *,
    origin_hessian: sp.Matrix,
    probe_hessians: Sequence[tuple[tuple[int, ...], sp.Matrix]],
    timeout: int,
) -> dict[str, object]:
    """Exclude unit Schur blocks on one affine Pluecker chart."""

    left, right = grassmann_chart_basis(pivots)
    left_left_polynomial = symbolic_quadratic_form(
        spatial_hessian, left, left
    )
    left_right_polynomial = symbolic_quadratic_form(
        spatial_hessian, left, right
    )
    right_right_polynomial = symbolic_quadratic_form(
        spatial_hessian, right, right
    )
    joint_cubic_polynomials = (
        polynomial_directional_derivative(
            left_left_polynomial, left
        ),
        polynomial_directional_derivative(
            left_left_polynomial, right
        ),
        polynomial_directional_derivative(
            right_right_polynomial, left
        ),
        polynomial_directional_derivative(
            right_right_polynomial, right
        ),
    )
    joint_cubic_equations = [
        coefficient
        for polynomial in joint_cubic_polynomials
        for coefficient in polynomial.coeffs()
        if coefficient != 0
    ]
    determinant_at_origin = sp.expand(
        (left.T * origin_hessian * left)[0]
        * (right.T * origin_hessian * right)[0]
        - (left.T * origin_hessian * right)[0] ** 2
    )
    parameter_variables = (
        symbolic_words.b,
        *GRASSMANN_VARIABLES,
    )
    localized_at = symbolic_words.b * determinant_at_origin

    initial = singular_localized_unit_ideal(
        joint_cubic_equations,
        parameter_variables=parameter_variables,
        localized_at=localized_at,
        timeout=timeout,
    )
    if initial["unit_ideal"]:
        return {
            "pivots": list(pivots),
            "certificate_kind": (
                "joint_cubic_plus_origin_determinant"
            ),
            "joint_cubic_coefficient_generators": initial[
                "generators"
            ],
            "determinant_probe_generators": 0,
            **initial,
        }

    probe_equations = []
    for _, point_hessian in probe_hessians:
        point_left_left = sp.expand(
            (left.T * point_hessian * left)[0]
        )
        point_left_right = sp.expand(
            (left.T * point_hessian * right)[0]
        )
        point_right_right = sp.expand(
            (right.T * point_hessian * right)[0]
        )
        determinant_at_point = sp.expand(
            point_left_left * point_right_right - point_left_right**2
        )
        difference = sp.expand(
            determinant_at_point - determinant_at_origin
        )
        if difference != 0:
            probe_equations.append(difference)
    sampled = singular_localized_unit_ideal(
        [*joint_cubic_equations, *probe_equations],
        parameter_variables=parameter_variables,
        localized_at=localized_at,
        timeout=timeout,
    )
    if sampled["unit_ideal"]:
        return {
            "pivots": list(pivots),
            "certificate_kind": "sampled_determinant_unit_ideal",
            "joint_cubic_coefficient_generators": initial[
                "generators"
            ],
            "determinant_probe_generators": len(probe_equations),
            **sampled,
        }

    left_left = left_left_polynomial.as_expr()
    left_right = left_right_polynomial.as_expr()
    right_right = right_right_polynomial.as_expr()
    determinant = sp.expand(left_left * right_right - left_right**2)
    determinant_polynomial = sp.Poly(determinant, *base.variables)
    constant_monomial = (0,) * len(base.variables)
    determinant_coefficients = [
        coefficient
        for monomial, coefficient in determinant_polynomial.terms()
        if monomial != constant_monomial
    ]
    full = singular_localized_unit_ideal(
        [*joint_cubic_equations, *determinant_coefficients],
        parameter_variables=parameter_variables,
        localized_at=localized_at,
        timeout=timeout,
    )
    return {
        "pivots": list(pivots),
        "certificate_kind": "full_determinant_coefficient_ideal",
        "joint_cubic_coefficient_generators": initial["generators"],
        "determinant_probe_generators": len(probe_equations),
        "full_determinant_coefficient_generators": len(
            determinant_coefficients
        ),
        **full,
    }


def analyze_symbolic_pattern(
    pattern_index: int,
    *,
    singular_timeout: int,
) -> dict[str, object]:
    pattern = symbolic_words.patterns()[pattern_index]
    a_value = exceptional_a_value(pattern)
    potential = sp.expand(
        symbolic_words.transformed_potential(
            pattern, "cubic-quadratic"
        ).subs(symbolic_words.a, a_value)
    )
    hessian = sp.hessian(potential, base.variables)
    hessian_polynomials = tuple(
        tuple(
            sp.Poly(
                hessian[row, column],
                *base.variables,
                domain=sp.QQ.poly_ring(symbolic_words.b),
            )
            for column in range(6)
        )
        for row in range(6)
    )
    spatial_hessian = tuple(
        tuple(
            sp.Poly(
                hessian[row, column],
                *base.variables,
                domain=sp.EX,
            )
            for column in range(6)
        )
        for row in range(6)
    )
    origin_hessian = evaluated_hessian(
        hessian_polynomials, (0,) * 6
    )
    probe_hessians = tuple(
        (point, evaluated_hessian(hessian_polynomials, point))
        for point in symbolic_schur_probe_points()
    )
    scalar = scalar_saturated_certificate(
        hessian, timeout=singular_timeout
    )
    charts = [
        grassmann_chart_certificate(
            spatial_hessian,
            pivots,
            origin_hessian=origin_hessian,
            probe_hessians=probe_hessians,
            timeout=singular_timeout,
        )
        for pivots in combinations(range(6), 2)
    ]
    return {
        "pattern_index": pattern_index,
        "pattern_id": pattern.pattern_id,
        "poisson_incidence": pattern.incidence,
        "linear_pairing": pattern.kappa,
        "a": sp.sstr(a_value),
        "scalar": scalar,
        "grassmann_charts": charts,
    }


def analyze_symbolic_pattern_worker(
    arguments: tuple[int, int],
) -> dict[str, object]:
    pattern_index, singular_timeout = arguments
    return analyze_symbolic_pattern(
        pattern_index, singular_timeout=singular_timeout
    )


def primitive_directions() -> tuple[tuple[int, ...], ...]:
    """Return the declared projective direction box with a fixed sign."""

    result = []
    for direction in product(DIRECTION_VALUES, repeat=6):
        if not any(direction):
            continue
        first_nonzero = next(value for value in direction if value)
        if first_nonzero == 1:
            result.append(direction)
    assert len(result) == (3**6 - 1) // 2
    return tuple(result)


def deterministic_points() -> tuple[tuple[int, ...], ...]:
    """Small exact points used only for modular rejection witnesses."""

    axes = tuple(
        tuple(1 if index == active else 0 for index in range(6))
        for active in range(6)
    )
    return (
        (0, 0, 0, 0, 0, 0),
        *axes,
        (1, 1, 1, 1, 1, 1),
        (-1, 2, -2, 1, 2, -1),
        (2, -1, 1, -2, 1, 3),
        (1, -2, 3, -1, 2, -3),
    )


def rational_mod(value: sp.Expr) -> int:
    rational = sp.Rational(value)
    numerator = int(rational.p) % PRIME
    denominator = int(rational.q) % PRIME
    if denominator == 0:
        raise ZeroDivisionError("bad modular prime for a rational coefficient")
    return numerator * pow(denominator, -1, PRIME) % PRIME


class ModularPolynomial:
    """Sparse evaluator for a rational polynomial modulo ``PRIME``."""

    def __init__(self, polynomial: sp.Poly):
        self.terms = tuple(
            (monomial, rational_mod(coefficient))
            for monomial, coefficient in polynomial.terms()
        )

    def evaluate(self, point: Sequence[int]) -> int:
        value = 0
        for monomial, coefficient in self.terms:
            term = coefficient
            for coordinate, exponent in zip(point, monomial, strict=True):
                if exponent:
                    term = (
                        term
                        * pow(int(coordinate) % PRIME, exponent, PRIME)
                        % PRIME
                    )
            value = (value + term) % PRIME
        return value


def polynomial_derivatives(potential: sp.Poly) -> PolynomialDerivatives:
    gradient = tuple(
        potential.diff(variable) for variable in base.variables
    )
    hessian = tuple(
        tuple(entry.diff(variable) for variable in base.variables)
        for entry in gradient
    )
    third = tuple(
        tuple(
            tuple(entry.diff(variable) for variable in base.variables)
            for entry in row
        )
        for row in hessian
    )
    return PolynomialDerivatives(hessian=hessian, third=third)


def evaluate_derivatives(
    derivatives: PolynomialDerivatives,
    points: Sequence[Sequence[int]],
) -> tuple[
    tuple[
        tuple[tuple[int, ...], ...],
        tuple[tuple[tuple[int, ...], ...], ...],
    ],
    ...,
]:
    modular_hessian = tuple(
        tuple(ModularPolynomial(entry) for entry in row)
        for row in derivatives.hessian
    )
    modular_third = tuple(
        tuple(
            tuple(ModularPolynomial(entry) for entry in row)
            for row in matrix
        )
        for matrix in derivatives.third
    )
    result = []
    for point in points:
        hessian = tuple(
            tuple(entry.evaluate(point) for entry in row)
            for row in modular_hessian
        )
        third = tuple(
            tuple(
                tuple(entry.evaluate(point) for entry in row)
                for row in matrix
            )
            for matrix in modular_third
        )
        result.append((hessian, third))
    return tuple(result)


def matrix_vector(
    matrix: Sequence[Sequence[int]], direction: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        sum(
            int(matrix[row][column]) * int(direction[column])
            for column in range(6)
        )
        % PRIME
        for row in range(6)
    )


def third_square_vector(
    tensor: Sequence[Sequence[Sequence[int]]],
    direction: Sequence[int],
) -> tuple[int, ...]:
    """Return k -> T(c,c,e_k) modulo the good prime."""

    return tuple(
        sum(
            int(direction[left])
            * int(direction[right])
            * int(tensor[left][right][last])
            for left in range(6)
            for right in range(6)
        )
        % PRIME
        for last in range(6)
    )


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        int(a_value) * int(b_value)
        for a_value, b_value in zip(left, right, strict=True)
    ) % PRIME


def exact_origin_bilinear(
    derivatives: PolynomialDerivatives,
    left: Sequence[int],
    right: Sequence[int],
) -> sp.Rational:
    zero_monomial = (0,) * 6
    return sp.Rational(
        sum(
            int(left[row])
            * int(right[column])
            * derivatives.hessian[row][column].coeff_monomial(
                zero_monomial
            )
            for row in range(6)
            for column in range(6)
        )
    )


def quadratic_form_poly(
    derivatives: PolynomialDerivatives,
    left: Sequence[int],
    right: Sequence[int],
) -> sp.Poly:
    result = sp.Poly(0, *base.variables, domain=sp.QQ)
    for row in range(6):
        if left[row] == 0:
            continue
        for column in range(6):
            coefficient = left[row] * right[column]
            if coefficient:
                result += coefficient * derivatives.hessian[row][column]
    return result


def directional_derivative(
    polynomial: sp.Poly, direction: Sequence[int]
) -> sp.Poly:
    result = sp.Poly(0, *base.variables, domain=sp.QQ)
    for coefficient, variable in zip(
        direction, base.variables, strict=True
    ):
        if coefficient:
            result += coefficient * polynomial.diff(variable)
    return result


def plane_key(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...] | None:
    """Canonical primitive Pluecker coordinates for a rational two-plane."""

    minors = [
        int(left[i]) * int(right[j]) - int(left[j]) * int(right[i])
        for i in range(6)
        for j in range(i + 1, 6)
    ]
    common = 0
    for value in minors:
        common = gcd(common, abs(value))
    if common == 0:
        return None
    primitive = [value // common for value in minors]
    if next(value for value in primitive if value) < 0:
        primitive = [-value for value in primitive]
    return tuple(primitive)


def constant_nonzero(polynomial: sp.Poly) -> bool:
    return polynomial.total_degree() <= 0 and not polynomial.is_zero


def witness_record(polynomial: sp.Poly) -> dict[str, object]:
    encoded = sp.sstr(polynomial.as_expr())
    leading_monomial, leading_coefficient = polynomial.terms()[0]
    return {
        "total_degree": polynomial.total_degree(),
        "terms": len(polynomial.terms()),
        "leading_monomial": list(leading_monomial),
        "leading_coefficient": sp.sstr(leading_coefficient),
        "sha256": hashlib.sha256(encoded.encode()).hexdigest(),
    }


def scalar_pivots(
    directions: Sequence[tuple[int, ...]],
    derivatives: PolynomialDerivatives,
    evaluated: Sequence[
        tuple[
            Sequence[Sequence[int]],
            Sequence[Sequence[Sequence[int]]],
        ]
    ],
) -> tuple[list[dict[str, object]], int]:
    modular_candidates = []
    for direction in directions:
        values = []
        for hessian, _ in evaluated:
            hessian_times_direction = matrix_vector(hessian, direction)
            values.append(dot(direction, hessian_times_direction))
        origin_value = exact_origin_bilinear(
            derivatives, direction, direction
        )
        if len(set(values)) == 1 and origin_value != 0:
            modular_candidates.append(direction)

    exact = []
    for direction in modular_candidates:
        second = quadratic_form_poly(
            derivatives, direction, direction
        )
        if constant_nonzero(second):
            exact.append(
                {
                    "direction": list(direction),
                    "pivot_coefficient": sp.sstr(second.as_expr()),
                }
            )
    return exact, len(modular_candidates)


def cubic_null_directions(
    directions: Sequence[tuple[int, ...]],
    derivatives: PolynomialDerivatives,
    evaluated: Sequence[
        tuple[
            Sequence[Sequence[int]],
            Sequence[Sequence[Sequence[int]]],
        ]
    ],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    modular_candidates = []
    for direction in directions:
        if all(
            dot(direction, third_square_vector(third, direction)) == 0
            for _, third in evaluated
        ):
            modular_candidates.append(direction)

    exact = []
    for direction in modular_candidates:
        second = quadratic_form_poly(
            derivatives, direction, direction
        )
        if directional_derivative(second, direction).is_zero:
            exact.append(direction)
    return tuple(exact), len(modular_candidates)


def schur_blocks(
    cubic_directions: Sequence[tuple[int, ...]],
    derivatives: PolynomialDerivatives,
    evaluated: Sequence[
        tuple[
            Sequence[Sequence[int]],
            Sequence[Sequence[Sequence[int]]],
        ]
    ],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, int],
]:
    third_square = {
        direction: tuple(
            third_square_vector(third, direction)
            for _, third in evaluated
        )
        for direction in cubic_directions
    }
    hessian_vectors = {
        direction: tuple(
            matrix_vector(hessian, direction)
            for hessian, _ in evaluated
        )
        for direction in cubic_directions
    }

    distinct_planes: dict[
        tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]
    ] = {}
    for left, right in combinations(cubic_directions, 2):
        key = plane_key(left, right)
        if key is not None:
            distinct_planes.setdefault(key, (left, right))

    jointly_quadratic_modular = []
    determinant_modular = []
    for key, (left, right) in distinct_planes.items():
        if not all(
            dot(right, third_square[left][point_index]) == 0
            and dot(left, third_square[right][point_index]) == 0
            for point_index in range(len(evaluated))
        ):
            continue
        jointly_quadratic_modular.append((key, left, right))

        determinant_values = []
        for point_index in range(len(evaluated)):
            left_left = dot(left, hessian_vectors[left][point_index])
            right_right = dot(right, hessian_vectors[right][point_index])
            left_right = dot(left, hessian_vectors[right][point_index])
            determinant_values.append(
                (left_left * right_right - left_right * left_right)
                % PRIME
            )
        if len(set(determinant_values)) == 1:
            left_left_origin = exact_origin_bilinear(
                derivatives, left, left
            )
            left_right_origin = exact_origin_bilinear(
                derivatives, left, right
            )
            right_right_origin = exact_origin_bilinear(
                derivatives, right, right
            )
            origin_determinant = (
                left_left_origin * right_right_origin
                - left_right_origin**2
            )
            if origin_determinant != 0:
                determinant_modular.append((key, left, right))

    exact = []
    exact_jointly_quadratic = 0
    rejections = []
    for key, left, right in determinant_modular:
        left_left = quadratic_form_poly(derivatives, left, left)
        left_right = quadratic_form_poly(derivatives, left, right)
        right_right = quadratic_form_poly(derivatives, right, right)
        mixed_left = directional_derivative(left_left, right)
        mixed_right = directional_derivative(right_right, left)
        if not mixed_left.is_zero or not mixed_right.is_zero:
            witness = mixed_left if not mixed_left.is_zero else mixed_right
            rejections.append(
                {
                    "basis": [list(left), list(right)],
                    "pluecker_coordinates": list(key),
                    "reason": "not_jointly_quadratic",
                    "witness": witness_record(witness),
                }
            )
            continue
        exact_jointly_quadratic += 1
        determinant = left_left * right_right - left_right * left_right
        if constant_nonzero(determinant):
            exact.append(
                {
                    "basis": [list(left), list(right)],
                    "pluecker_coordinates": list(key),
                    "block_determinant": sp.sstr(determinant.as_expr()),
                }
            )
        else:
            rejections.append(
                {
                    "basis": [list(left), list(right)],
                    "pluecker_coordinates": list(key),
                    "reason": (
                        "zero_block_determinant"
                        if determinant.is_zero
                        else "nonconstant_block_determinant"
                    ),
                    "witness": (
                        None
                        if determinant.is_zero
                        else witness_record(determinant)
                    ),
                }
            )

    counts = {
        "distinct_candidate_planes": len(distinct_planes),
        "jointly_quadratic_modular_survivors": len(
            jointly_quadratic_modular
        ),
        "unit_determinant_modular_survivors": len(
            determinant_modular
        ),
        "exact_jointly_quadratic_near_misses": exact_jointly_quadratic,
    }
    return exact, rejections, counts


def analyze_chart(
    pattern: symbolic_words.Pattern,
    b_value: int,
    directions: Sequence[tuple[int, ...]],
    points: Sequence[Sequence[int]],
) -> dict[str, object]:
    a_value = exceptional_a_value(pattern)
    potential_expression = sp.expand(
        symbolic_words.transformed_potential(
            pattern, "cubic-quadratic"
        ).subs({symbolic_words.a: a_value, symbolic_words.b: b_value})
    )
    potential = sp.Poly(
        potential_expression, *base.variables, domain=sp.QQ
    )
    derivatives = polynomial_derivatives(potential)
    evaluated = evaluate_derivatives(derivatives, points)

    parent_determinants = [
        base.determinant_mod(hessian) for hessian, _ in evaluated
    ]
    expected_parent = base.BASE_HESSIAN_DETERMINANT % PRIME
    if any(value != expected_parent for value in parent_determinants):
        raise AssertionError("HC4MCP6 parent-Hessian regression failed")

    scalar_exact, scalar_modular_count = scalar_pivots(
        directions, derivatives, evaluated
    )
    cubic_exact, cubic_modular_count = cubic_null_directions(
        directions, derivatives, evaluated
    )
    block_exact, block_rejections, block_counts = schur_blocks(
        cubic_exact, derivatives, evaluated
    )

    return {
        "chart_id": f"{pattern.pattern_id}-b{b_value:+d}",
        "pattern_id": pattern.pattern_id,
        "poisson_incidence": pattern.incidence,
        "linear_pairing": pattern.kappa,
        "a": sp.sstr(a_value),
        "b": b_value,
        "potential_terms": len(potential.terms()),
        "parent_hessian_value_mod_p": expected_parent,
        "scalar_modular_survivors": scalar_modular_count,
        "scalar_unit_pivots": scalar_exact,
        "cubic_null_modular_survivors": cubic_modular_count,
        "cubic_null_directions": len(cubic_exact),
        **block_counts,
        "unit_schur_blocks": block_exact,
        "unit_determinant_near_miss_rejections": block_rejections,
    }


def sum_field(rows: Iterable[dict[str, object]], field: str) -> int:
    return sum(int(row[field]) for row in rows)


def run_search() -> dict[str, object]:
    directions = primitive_directions()
    points = deterministic_points()
    patterns = symbolic_words.patterns()
    rows = []
    total = len(patterns) * len(B_BOX)
    for pattern_index, pattern in enumerate(patterns):
        for b_index, b_value in enumerate(B_BOX):
            index = pattern_index * len(B_BOX) + b_index + 1
            print(
                f"progress={index}/{total} "
                f"pattern={pattern.pattern_id} b={b_value:+d}",
                flush=True,
            )
            rows.append(
                analyze_chart(
                    pattern, b_value, directions, points
                )
            )

    scalar_exact = sum(
        len(row["scalar_unit_pivots"]) for row in rows
    )
    block_exact = sum(
        len(row["unit_schur_blocks"]) for row in rows
    )
    cubic_null_total = sum_field(rows, "cubic_null_directions")
    plane_total = sum_field(rows, "distinct_candidate_planes")
    determinant_near_misses = sum_field(
        rows, "unit_determinant_modular_survivors"
    )
    exact_joint_near_misses = sum_field(
        rows, "exact_jointly_quadratic_near_misses"
    )
    assert len(patterns) == 54
    assert len(rows) == 216
    assert scalar_exact == 0
    assert block_exact == 0
    assert cubic_null_total == 4968
    assert plane_total == 32360
    assert determinant_near_misses == 72
    assert exact_joint_near_misses == 72
    return {
        "status": "exact_finite_box_search",
        "scope": {
            "base": "collision-centred foundational cubic Keller doubling",
            "dependency": (
                "HC4MCP6 exact reverse cubic--quadratic "
                "parent-preserving family"
            ),
            "composition_order": "cubic then quadratic",
            "patterns": "all 54 noncommuting shared-dual supports",
            "b_box": list(B_BOX),
            "a_values": ["-1/2", "1/2", "-1/4", "1/4"],
            "direction_box": "({-1,0,1}^6-{0})/{+/-1}",
            "directions": len(directions),
            "prime": PRIME,
            "modular_points": [list(point) for point in points],
            "exact_survivor_check": "symbolic polynomial identities over Q",
            "limitations": [
                "directions outside the declared coefficient box",
                "cubic coefficients b outside the declared box",
                "different Hamiltonian supports or degrees",
                "words of length at least three",
                "coefficient-dependent or nonlinear pivot directions",
                "general polynomial symplectomorphisms",
            ],
        },
        "summary": {
            "patterns": len(patterns),
            "b_values": len(B_BOX),
            "charts": len(rows),
            "directions_per_chart": len(directions),
            "scalar_direction_trials": len(rows) * len(directions),
            "scalar_modular_survivors": sum_field(
                rows, "scalar_modular_survivors"
            ),
            "scalar_unit_pivots": scalar_exact,
            "cubic_null_modular_survivors": sum_field(
                rows, "cubic_null_modular_survivors"
            ),
            "cubic_null_directions": sum_field(
                rows, "cubic_null_directions"
            ),
            "distinct_candidate_planes": sum_field(
                rows, "distinct_candidate_planes"
            ),
            "jointly_quadratic_modular_survivors": sum_field(
                rows, "jointly_quadratic_modular_survivors"
            ),
            "unit_determinant_modular_survivors": sum_field(
                rows, "unit_determinant_modular_survivors"
            ),
            "exact_jointly_quadratic_near_misses": sum_field(
                rows, "exact_jointly_quadratic_near_misses"
            ),
            "unit_schur_blocks": block_exact,
            "parent_hessian_regressions": len(rows),
        },
        "charts": rows,
    }


def run_symbolic_classification(
    *,
    jobs: int,
    singular_timeout: int,
    checkpoint_directory: Path | None,
) -> dict[str, object]:
    if shutil.which("Singular") is None:
        raise RuntimeError(
            "Singular is required for the symbolic classification"
        )
    if jobs < 1:
        raise ValueError("--jobs must be positive")

    patterns = symbolic_words.patterns()
    rows: list[dict[str, object] | None] = [None] * len(patterns)
    checker_sha256 = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    if checkpoint_directory is not None:
        checkpoint_directory.mkdir(parents=True, exist_ok=True)
        for pattern_index in range(len(patterns)):
            checkpoint_path = (
                checkpoint_directory / f"pattern-{pattern_index:02d}.json"
            )
            if not checkpoint_path.exists():
                continue
            checkpoint = json.loads(checkpoint_path.read_text())
            if checkpoint.get("checker_sha256") != checker_sha256:
                continue
            row = checkpoint.get("row")
            if not isinstance(row, dict):
                continue
            if int(row.get("pattern_index", -1)) != pattern_index:
                continue
            rows[pattern_index] = row
    arguments = [
        (pattern_index, singular_timeout)
        for pattern_index in range(len(patterns))
        if rows[pattern_index] is None
    ]
    loaded = len(patterns) - len(arguments)

    def store_checkpoint(row: dict[str, object]) -> None:
        if checkpoint_directory is None:
            return
        pattern_index = int(row["pattern_index"])
        checkpoint_path = (
            checkpoint_directory / f"pattern-{pattern_index:02d}.json"
        )
        checkpoint_path.write_text(
            json.dumps(
                {"checker_sha256": checker_sha256, "row": row},
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    if loaded:
        print(
            f"symbolic_checkpoints_loaded={loaded}/{len(patterns)}",
            flush=True,
        )
    if jobs == 1:
        for completed, arguments_item in enumerate(
            arguments, start=loaded + 1
        ):
            row = analyze_symbolic_pattern_worker(arguments_item)
            rows[int(row["pattern_index"])] = row
            store_checkpoint(row)
            print(
                f"symbolic_progress={completed}/{len(patterns)} "
                f"pattern={row['pattern_id']}",
                flush=True,
            )
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    analyze_symbolic_pattern_worker, arguments_item
                ): arguments_item[0]
                for arguments_item in arguments
            }
            for completed, future in enumerate(
                as_completed(futures), start=loaded + 1
            ):
                row = future.result()
                rows[int(row["pattern_index"])] = row
                store_checkpoint(row)
                print(
                    f"symbolic_progress={completed}/{len(patterns)} "
                    f"pattern={row['pattern_id']}",
                    flush=True,
                )

    completed_rows = [row for row in rows if row is not None]
    if len(completed_rows) != len(patterns):
        raise AssertionError("symbolic pattern census is incomplete")
    scalar_failures = [
        row["pattern_id"]
        for row in completed_rows
        if not row["scalar"][
            "constant_coefficient_normal_form_zero"
        ]
    ]
    chart_failures = [
        (row["pattern_id"], chart["pivots"])
        for row in completed_rows
        for chart in row["grassmann_charts"]
        if not chart["unit_ideal"]
    ]
    if scalar_failures:
        raise AssertionError(
            "scalar classification survivors: "
            + ", ".join(str(value) for value in scalar_failures)
        )
    if chart_failures:
        raise AssertionError(
            "two-plane classification survivors: "
            + ", ".join(str(value) for value in chart_failures)
        )

    certificate_counts: dict[str, int] = {}
    for row in completed_rows:
        for chart in row["grassmann_charts"]:
            kind = str(chart["certificate_kind"])
            certificate_counts[kind] = certificate_counts.get(kind, 0) + 1
    scalar_basis_sizes = sorted(
        {
            int(row["scalar"]["saturated_basis_size"])
            for row in completed_rows
        }
    )
    scalar_dimensions = sorted(
        {
            int(row["scalar"]["saturated_scheme_dimension"])
            for row in completed_rows
        }
    )
    return {
        "status": "exact_symbolic_family_classification",
        "scope": {
            "base": "collision-centred foundational cubic Keller doubling",
            "dependency": (
                "HC4MCP6 exact reverse cubic--quadratic "
                "parent-preserving family"
            ),
            "composition_order": "cubic then quadratic",
            "patterns": "all 54 noncommuting shared-dual supports",
            "a_values": ["-1/2", "1/2", "-1/4", "1/4"],
            "b": "arbitrary nonzero complex coefficient",
            "scalar_directions": "all constant c in C^6",
            "two_planes": (
                "all constant planes via the 15 standard affine "
                "Pluecker charts on Gr(2,6)"
            ),
            "determinant_probe_points": [
                list(point) for point in symbolic_schur_probe_points()
            ],
            "coefficient_field": "Q; unit ideals exclude the complex locus",
            "limitations": [
                "different Hamiltonian supports or degrees",
                "words of length at least three",
                "coefficient-dependent or nonlinear pivot directions",
                "general polynomial symplectomorphisms",
            ],
        },
        "summary": {
            "patterns": len(completed_rows),
            "scalar_symbolic_schemes": len(completed_rows),
            "scalar_nonzero_constant_pivots": 0,
            "scalar_saturated_basis_sizes": scalar_basis_sizes,
            "scalar_saturated_scheme_dimensions": scalar_dimensions,
            "grassmann_charts_per_pattern": 15,
            "grassmann_charts": sum(
                len(row["grassmann_charts"])
                for row in completed_rows
            ),
            "constant_two_plane_unit_schur_blocks": 0,
            "certificate_kinds": certificate_counts,
        },
        "patterns": completed_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete JSON census to this path",
    )
    parser.add_argument(
        "--symbolic-classification",
        action="store_true",
        help=(
            "classify arbitrary constant scalar directions and "
            "two-planes over Q[b] instead of running the finite box"
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="parallel pattern workers for --symbolic-classification",
    )
    parser.add_argument(
        "--singular-timeout",
        type=int,
        default=300,
        help="seconds allowed for each exact Singular standard basis",
    )
    parser.add_argument(
        "--checkpoint-directory",
        type=Path,
        help=(
            "optional resumable per-pattern checkpoint directory for "
            "--symbolic-classification"
        ),
    )
    args = parser.parse_args()

    if args.symbolic_classification:
        result = run_symbolic_classification(
            jobs=args.jobs,
            singular_timeout=args.singular_timeout,
            checkpoint_directory=args.checkpoint_directory,
        )
    else:
        result = run_search()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(f"artifact={args.output}")
        print(
            "artifact_sha256="
            f"{hashlib.sha256(encoded.encode()).hexdigest()}"
        )
    print("HC4_NONLINEAR_UNIT_SCHUR_BLOCK_SUMMARY")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    if args.symbolic_classification:
        print(
            "SCOPE: all constant scalar directions and two-planes in "
            "the 54 HC4MCP6 families for arbitrary b != 0; no "
            "classification of nonlinear polynomial symplectomorphisms"
        )
    else:
        print(
            "SCOPE: exact finite coefficient/direction box only; no "
            "classification of nonlinear polynomial symplectomorphisms"
        )


if __name__ == "__main__":
    main()

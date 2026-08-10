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

Evaluations modulo the good prime 1000003 are rejection certificates only.
Every modular survivor is checked symbolically over Q.  The result is a
finite-box computation, not a classification of arbitrary directions,
coefficients, longer symplectic words, or polynomial symplectomorphisms.
"""

from __future__ import annotations

import argparse
import hashlib
import json
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete JSON census to this path",
    )
    args = parser.parse_args()

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
    print(
        "SCOPE: exact finite coefficient/direction box only; no "
        "classification of nonlinear polynomial symplectomorphisms"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact classification and dependency certificate for Nagao linear sections.

The classification is deliberately scoped to polynomial sections with
``x=m*T+n`` and ``degree(y)<=3``.  It is not a classification of every
rational section of the elliptic surface.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable

import sympy as sp

from nagao_1994 import (
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
)
from nagao_linear_sections import (
    COMPANION_JACOBIAN_RELATIONS,
    LINEAR_COMPANION_SECTIONS,
)


Q = Fraction
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_nagao_linear_sections.py"
)


def _sympy_rational(value: Fraction) -> sp.Rational:
    value = Q(value)
    return sp.Rational(value.numerator, value.denominator)


def _interpolated_quartic(parameter: sp.Symbol, x: sp.Symbol) -> sp.Expr:
    coefficients: list[sp.Expr] = []
    for index in range(5):
        values = []
        for sample in range(1, 9):
            coefficient = primitive_quartic_coefficients(
                RANK13_CONSTRUCTION, Q(sample)
            )[index]
            values.append((sample, _sympy_rational(coefficient)))
        polynomial = sp.expand(sp.interpolate(values, parameter))
        check = primitive_quartic_coefficients(RANK13_CONSTRUCTION, Q(9))[index]
        if polynomial.subs(parameter, 9) != _sympy_rational(check):
            raise AssertionError("quartic interpolation failed its check value")
        coefficients.append(polynomial)
    return sp.expand(sum(value * x**index for index, value in enumerate(coefficients)))


def _monic(polynomial: sp.Poly) -> sp.Poly:
    return sp.Poly(polynomial.as_expr() / polynomial.LC(), *polynomial.gens)


def classify_abscissae() -> dict[str, Any]:
    """Eliminate the ordinate and classify the linear-abscissa ansatz."""

    parameter, x, slope, intercept = sp.symbols("T X m n")
    y0, y1, y2 = sp.symbols("y0 y1 y2")
    quartic = _interpolated_quartic(parameter, x)

    # Away from slope +/-1, comparison of the leading four coefficients
    # determines y3,y2,y1,y0.  The opposite choice of y3 negates the whole
    # ordinate and yields the same abscissa classification.
    y3 = 3 * (1 - slope**2)
    y2_formula = -6 * slope * (intercept - 75)
    y1_formula = -(
        105975 * slope**4
        + 9 * slope**2 * intercept**2
        - 1350 * slope**2 * intercept
        + 96825 * slope**2
        - 9 * intercept**2
        + 1350 * intercept
        - 79600
    ) / (3 * (slope - 1) * (slope + 1))
    y0_formula = -5 * slope * (
        42390 * slope**4 * intercept
        - 3210921 * slope**4
        - 84780 * slope**2 * intercept
        + 6337386 * slope**2
        - 6890 * intercept
        + 569535
    ) / (3 * (slope - 1) ** 2 * (slope + 1) ** 2)
    ordinate = (
        y3 * parameter**3
        + y2_formula * parameter**2
        + y1_formula * parameter
        + y0_formula
    )
    residual = sp.Poly(sp.together(quartic.subs(x, slope * parameter + intercept) - ordinate**2), parameter)
    equations = []
    for degree in (2, 1, 0):
        numerator = sp.factor(
            sp.together(residual.coeff_monomial(parameter**degree))
        ).as_numer_denom()[0]
        equations.append(numerator)

    resultants = (
        sp.Poly(sp.resultant(equations[0], equations[1], intercept), slope),
        sp.Poly(sp.resultant(equations[1], equations[2], intercept), slope),
        sp.Poly(sp.resultant(equations[0], equations[2], intercept), slope),
    )
    common = _monic(sp.gcd(sp.gcd(resultants[0], resultants[1]), resultants[2]))
    expected = _monic(
        sp.Poly(
            (3 * slope - 5)
            * (3 * slope + 5)
            * (15 * slope - 7)
            * (15 * slope + 7)
            * (15 * slope - 1)
            * (15 * slope + 1),
            slope,
        )
    )
    if common != expected:
        raise AssertionError("unexpected resultant gcd in the nonsingular branch")

    recovered: list[tuple[sp.Rational, sp.Rational]] = []
    for section in LINEAR_COMPANION_SECTIONS:
        candidate_slope = _sympy_rational(section.slope)
        specialized = [sp.Poly(eq.subs(slope, candidate_slope), intercept) for eq in equations]
        common_intercept = _monic(
            sp.gcd(sp.gcd(specialized[0], specialized[1]), specialized[2])
        )
        expected_intercept = sp.Poly(
            intercept - _sympy_rational(section.intercept), intercept
        )
        if common_intercept != expected_intercept:
            raise AssertionError(f"failed to recover {section.label}")
        recovered.append((candidate_slope, _sympy_rational(section.intercept)))

    # The leading cubic coefficient vanishes at slopes +/-1.  A separate
    # Groebner elimination gives exactly the six Mestre root intercepts.
    root_polynomial = sp.Poly(
        sp.prod(intercept - _sympy_rational(root) for root in RANK13_CONSTRUCTION.roots),
        intercept,
    )
    exceptional: dict[str, str] = {}
    for candidate_slope in (sp.Integer(-1), sp.Integer(1)):
        ordinate2 = y2 * parameter**2 + y1 * parameter + y0
        difference = sp.Poly(
            sp.expand(
                quartic.subs(x, candidate_slope * parameter + intercept)
                - ordinate2**2
            ),
            parameter,
        )
        equations2 = [difference.coeff_monomial(parameter**degree) for degree in range(5)]
        basis = sp.groebner(equations2, y2, y1, y0, intercept, order="lex")
        eliminants = [
            sp.Poly(item.as_expr(), intercept)
            for item in basis.polys
            if not any(item.as_expr().has(variable) for variable in (y2, y1, y0))
        ]
        if len(eliminants) != 1 or _monic(eliminants[0]) != _monic(root_polynomial):
            raise AssertionError("the slope +/-1 elimination missed the Mestre roots")
        exceptional[str(candidate_slope)] = str(sp.factor(eliminants[0].as_expr()))

    return {
        "nonsingular_slope_resultant_gcd": str(sp.factor(common.as_expr())),
        "recovered_companion_slope_intercept_pairs": [
            {"slope": str(item[0]), "intercept": str(item[1])}
            for item in recovered
        ],
        "slope_plus_or_minus_one_intercept_eliminants": exceptional,
        "classification_scope": (
            "all Q-rational polynomial sections with x=m*T+n and degree(y)<=3, "
            "up to negating y; this is not a classification of arbitrary rational sections"
        ),
    }


def _visible_symbolic_points(parameter: sp.Symbol) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    points: list[tuple[sp.Expr, sp.Expr]] = []
    for index in range(12):
        values = []
        for sample in range(1, 9):
            _, ordinate = primitive_visible_points(RANK13_CONSTRUCTION, Q(sample))[index]
            values.append((sample, _sympy_rational(ordinate)))
        ordinate_polynomial = sp.expand(sp.interpolate(values, parameter))
        _, check = primitive_visible_points(RANK13_CONSTRUCTION, Q(9))[index]
        if ordinate_polynomial.subs(parameter, 9) != _sympy_rational(check):
            raise AssertionError("visible ordinate interpolation failed")
        root = _sympy_rational(RANK13_CONSTRUCTION.roots[index // 2])
        sign = -1 if index % 2 == 0 else 1
        points.append((root + sign * parameter, ordinate_polynomial))
    return tuple(points)


def _covariant_map(
    quartic: sp.Expr, parameter: sp.Symbol, x_symbol: sp.Symbol, point: tuple[sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    polynomial = sp.Poly(quartic, x_symbol)
    e, d, c, b, a = (polynomial.coeff_monomial(x_symbol**degree) for degree in range(5))
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    g2 = c**2 / 12 - b * d / 8 - a * e
    g3 = c * d / 12 - b * e / 2
    g4 = d**2 / 16 - c * e / 6
    point_x, point_y = point
    g = g0 * point_x**4 + g1 * point_x**3 + g2 * point_x**2 + g3 * point_x + g4
    ux = 4 * a * point_x**3 + 3 * b * point_x**2 + 2 * c * point_x + d
    uy = b * point_x**3 + 2 * c * point_x**2 + 3 * d * point_x + 4 * e
    gx = 4 * g0 * point_x**3 + 3 * g1 * point_x**2 + 2 * g2 * point_x + g3
    gy = g1 * point_x**3 + 2 * g2 * point_x**2 + 3 * g3 * point_x + 4 * g4
    h = (ux * gy - uy * gx) / 8
    return (
        sp.cancel(36 * g / point_y**2),
        sp.cancel(108 * h / point_y**3),
        sp.expand(-27 * invariant_i),
        sp.expand(-27 * invariant_j),
    )


def _elliptic_add(
    left: tuple[sp.Expr, sp.Expr] | None,
    right: tuple[sp.Expr, sp.Expr] | None,
    coefficient_a: sp.Expr,
) -> tuple[sp.Expr, sp.Expr] | None:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if sp.cancel(x1 - x2) == 0:
        if sp.cancel(y1 + y2) == 0:
            return None
        slope = sp.cancel((3 * x1**2 + coefficient_a) / (2 * y1))
    else:
        slope = sp.cancel((y2 - y1) / (x2 - x1))
    x3 = sp.cancel(slope**2 - x1 - x2)
    y3 = sp.cancel(-y1 + slope * (x1 - x3))
    return x3, y3


def verify_symbolic_relations() -> list[dict[str, Any]]:
    """Prove the five companion dependencies in Q(T) by exact simplification."""

    parameter, x_symbol = sp.symbols("T X")
    quartic = _interpolated_quartic(parameter, x_symbol)
    visible_quartic = _visible_symbolic_points(parameter)
    visible_jacobian = []
    coefficient_a: sp.Expr | None = None
    coefficient_b: sp.Expr | None = None
    for point in visible_quartic:
        mapped_x, mapped_y, coefficient_a, coefficient_b = _covariant_map(
            quartic, parameter, x_symbol, point
        )
        visible_jacobian.append((mapped_x, mapped_y))
    assert coefficient_a is not None and coefficient_b is not None

    companion_jacobian: dict[str, tuple[sp.Expr, sp.Expr]] = {}
    for section in LINEAR_COMPANION_SECTIONS:
        point_x = _sympy_rational(section.slope) * parameter + _sympy_rational(section.intercept)
        point_y = sum(
            _sympy_rational(coefficient) * parameter**degree
            for degree, coefficient in enumerate(section.ordinate_coefficients)
        )
        if sp.expand(quartic.subs(x_symbol, point_x) - point_y**2) != 0:
            raise AssertionError(f"symbolic quartic identity failed for {section.label}")
        mapped_x, mapped_y, _, _ = _covariant_map(
            quartic, parameter, x_symbol, (point_x, point_y)
        )
        if sp.cancel(mapped_y**2 - mapped_x**3 - coefficient_a * mapped_x - coefficient_b) != 0:
            raise AssertionError(f"symbolic Jacobian identity failed for {section.label}")
        companion_jacobian[section.label] = (mapped_x, mapped_y)

    basis = tuple(visible_jacobian[index] for index in range(11)) + (
        companion_jacobian["plus-1/15"],
    )
    records: list[dict[str, Any]] = []
    for label, coefficients in COMPANION_JACOBIAN_RELATIONS.items():
        total: tuple[sp.Expr, sp.Expr] | None = None
        for point, coefficient in zip(basis, coefficients):
            if coefficient == 0:
                continue
            if coefficient not in (-1, 1):
                raise AssertionError("the pinned relations should have unit coefficients")
            summand = point if coefficient == 1 else (point[0], -point[1])
            total = _elliptic_add(total, summand, coefficient_a)
        if total is None:
            raise AssertionError(f"relation for {label} collapsed to infinity")
        target = companion_jacobian[label]
        if sp.cancel(total[0] - target[0]) != 0 or sp.cancel(total[1] - target[1]) != 0:
            raise AssertionError(f"symbolic Mordell--Weil relation failed for {label}")
        records.append(
            {
                "section": label,
                "basis_coefficients": list(coefficients),
                "exact_identity_in_Q_of_T": True,
            }
        )
    return records


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/generated-results/elliptic-curves/elliptic_nagao_linear_sections.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    classification = classify_abscissae()
    relations = verify_symbolic_relations()
    script_path = Path(__file__).resolve()
    formulas_path = script_path.with_name("nagao_linear_sections.py")
    artifact = {
        "schema_version": 1,
        "status": "exact symbolic classification and dependency certificate",
        "classification": classification,
        "companion_sections": [
            {
                "label": section.label,
                "slope": str(section.slope),
                "intercept": str(section.intercept),
                "ordinate_coefficients_ascending": [
                    str(value) for value in section.ordinate_coefficients
                ],
            }
            for section in LINEAR_COMPANION_SECTIONS
        ],
        "exact_jacobian_relations": relations,
        "frontier_consequence": (
            "five points previously labelled nonvisible by bounded searches are generic "
            "dependent companions; point-yield objectives must exclude them, while prior "
            "numerical height-rank values remain unchanged"
        ),
        "target_hit": False,
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "formulas_sha256": hashlib.sha256(formulas_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"classified {12 + len(LINEAR_COMPANION_SECTIONS)} linear sections")
    print(f"proved {len(relations)} exact companion relations in Q(T)")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

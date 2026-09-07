#!/usr/bin/env python3
"""Exact generic relation audit for the six D-square affine companions.

The split-infinity base change of the six-root tuple
``(0,25,95,143,168,205)`` has one previously certified nonvisible section.
This checker determines whether the remaining five displayed companions add
new generic directions.  It works in ``Q(u)`` and uses the compact Mestre
square recursion, primitive quartic normalization, the covariant map, and
exact short-Weierstrass addition.  No expanded two-section residual occurs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from search_mestre_dsquare_four import FAMILIES, base_parameter
from search_mestre_root_tuple_scale import (
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)


U, X = sp.symbols("u X")
BASE_CONSTANT = 39_146
ROOTS = (0, 25, 95, 143, 168, 205)
QUARTIC_SQUARE_SCALE = 6_600
QUARTIC_CONTENT = QUARTIC_SQUARE_SCALE**2
# Original-root coordinates: x=a+b*T, T=(39146-u^2)/(2u).
COMPANION_LINES = (
    (sp.Rational(-2375, 23), sp.Rational(37, 23)),
    (sp.Rational(-2375, 23), -sp.Rational(37, 23)),
    (sp.Rational(14575, 115), sp.Rational(13, 23)),
    (sp.Rational(14575, 115), -sp.Rational(13, 23)),
    (sp.Rational(17220, 115), sp.Rational(7, 23)),
    (sp.Rational(17220, 115), -sp.Rational(7, 23)),
)


def setup() -> tuple[object, object, list[object], object]:
    """Return Q(u), the base parameter, primitive quartic, and square root."""

    field = sp.QQ.frac_field(U)
    parameter = field.from_sympy((BASE_CONSTANT - U**2) / (2 * U))
    q = sp.Poly(sp.prod(X - root for root in ROOTS), X, domain=field)
    lower = sp.Poly(q.as_expr().subs(X, X - parameter.as_expr()), X, domain=field)
    upper = sp.Poly(q.as_expr().subs(X, X + parameter.as_expr()), X, domain=field)
    product = lower * upper
    approximant = sp.Poly(X**6, X, domain=field)
    for degree in range(5, -1, -1):
        target = product.coeff_monomial(X ** (6 + degree))
        square = (approximant * approximant).coeff_monomial(X ** (6 + degree))
        approximant += sp.Poly(
            (target - square) * X**degree / 2, X, domain=field
        )
    remainder = approximant * approximant - product
    if any(remainder.coeff_monomial(X**degree) for degree in range(5, 13)):
        raise AssertionError("the Mestre recursion did not leave a quartic remainder")
    quartic = [
        field.from_sympy(remainder.coeff_monomial(X**degree))
        / (parameter**2 * field.convert(QUARTIC_CONTENT))
        for degree in range(5)
    ]
    return field, parameter, quartic, approximant


def covariant_point(
    field: object, quartic: list[object], x_value: object, y_value: object
) -> tuple[object, object]:
    """Map a primitive quartic point to the primitive short Jacobian."""

    e, d, c, b, a = quartic
    g0 = b**2 / field.convert(16) - a * c / field.convert(6)
    g1 = b * c / field.convert(12) - a * d / field.convert(2)
    g2 = c**2 / field.convert(12) - b * d / field.convert(8) - a * e
    g3 = c * d / field.convert(12) - b * e / field.convert(2)
    g4 = d**2 / field.convert(16) - c * e / field.convert(6)
    g_value = g0 * x_value**4 + g1 * x_value**3 + g2 * x_value**2 + g3 * x_value + g4
    g_x = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
    g_y = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
    u_x = 4 * a * x_value**3 + 3 * b * x_value**2 + 2 * c * x_value + d
    u_y = b * x_value**3 + 2 * c * x_value**2 + 3 * d * x_value + 4 * e
    h_value = (u_x * g_y - u_y * g_x) / field.convert(8)
    return 36 * g_value / y_value**2, 108 * h_value / y_value**3


def short_add(
    coefficient_a: object,
    left: tuple[object, object] | None,
    right: tuple[object, object] | None,
) -> tuple[object, object] | None:
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    slope = (
        (3 * x_left**2 + coefficient_a) / (2 * y_left)
        if x_left == x_right
        else (y_right - y_left) / (x_right - x_left)
    )
    x_sum = slope**2 - x_left - x_right
    return x_sum, -y_left + slope * (x_left - x_sum)


def square_root_in_field(field: object, value: object) -> object:
    """Extract a rational-function square root after an exact factor audit."""

    numerator, denominator = map(
        sp.factor, sp.fraction(sp.cancel(value.as_expr()))
    )
    numerator_constant, numerator_factors = sp.factor_list(numerator)
    denominator_constant, denominator_factors = sp.factor_list(denominator)
    if any(exponent % 2 for _, exponent in numerator_factors + denominator_factors):
        raise AssertionError("a declared companion ordinate is not rational over Q(u)")
    root = sp.sqrt(numerator_constant / denominator_constant)
    root *= sp.prod(factor ** (exponent // 2) for factor, exponent in numerator_factors)
    root /= sp.prod(factor ** (exponent // 2) for factor, exponent in denominator_factors)
    return field.from_sympy(root)


def negate(point: tuple[object, object]) -> tuple[object, object]:
    return point[0], -point[1]


def evaluate_at_u(point: tuple[object, object], value: int) -> tuple[object, object]:
    return tuple(sp.cancel(coordinate.as_expr()).subs(U, value) for coordinate in point)


def replay() -> dict[str, object]:
    field, parameter, quartic, approximant = setup()
    family = FAMILIES[2]
    parameter_u = 197
    parameter_at_197 = base_parameter(family, parameter_u)
    quartic_at_197 = family.construction.primitive_quartic_coefficients(
        parameter_at_197
    )
    visible_at_197 = primitive_visible_points(family.construction, parameter_at_197)
    visible = {}
    for index, root in enumerate(ROOTS):
        for sign in (-1, 1):
            x_value = field.convert(root) + field.convert(sign) * parameter
            y_value = field.convert(approximant.eval(x_value)) / (
                parameter * field.convert(QUARTIC_SQUARE_SCALE)
            )
            if y_value**2 != sum(
                coefficient * x_value**degree
                for degree, coefficient in enumerate(quartic)
            ):
                raise AssertionError("a displayed visible ordinate failed its exact square identity")
            point = covariant_point(field, quartic, x_value, y_value)
            expected = quartic_point_to_jacobian(
                family.construction,
                parameter_at_197,
                visible_at_197[2 * index + (sign + 1) // 2],
            )
            specialization = evaluate_at_u(point, parameter_u)
            if specialization == expected:
                visible[index, sign] = point
            elif (specialization[0], -specialization[1]) == expected:
                visible[index, sign] = negate(point)
            else:
                raise AssertionError("a generic visible point did not match its u=197 orientation")

    companions = []
    for intercept, slope in COMPANION_LINES:
        x_value = field.from_sympy(intercept + slope * parameter.as_expr())
        value_on_quartic = sum(
            coefficient * x_value**degree
            for degree, coefficient in enumerate(quartic)
        )
        y_value = square_root_in_field(field, value_on_quartic)
        if y_value**2 != value_on_quartic:
            raise AssertionError("a declared companion ordinate failed its exact square identity")
        point = covariant_point(field, quartic, x_value, y_value)
        x_at_197 = sp.Rational(intercept) + sp.Rational(slope) * parameter_at_197
        y_at_197_squared = quartic_value(quartic_at_197, x_at_197)
        y_at_197 = sp.sqrt(sp.Rational(y_at_197_squared.numerator)) / sp.sqrt(
            sp.Rational(y_at_197_squared.denominator)
        )
        expected = quartic_point_to_jacobian(
            family.construction,
            parameter_at_197,
            (x_at_197, y_at_197),
        )
        specialization = evaluate_at_u(point, parameter_u)
        if specialization == expected:
            companions.append(point)
        elif (specialization[0], -specialization[1]) == expected:
            companions.append(negate(point))
        else:
            raise AssertionError("a generic companion did not match its u=197 orientation")

    # These support-five relations were found independently at u=197.  The
    # signs here are normalized by exact specialization before the generic
    # equality checks below.
    relations = (
        ((4, -1, 1), (4, 1, -1), (5, -1, -1), (5, 1, 1), ("P1", 0, -1)),
        ((0, 1, 1), (3, -1, 1), (4, 1, 1), (5, -1, 1), ("P1", 0, 1)),
        ((0, -1, -1), (3, 1, -1), (4, 1, -1), (5, -1, -1), ("P1", 0, -1)),
        ((1, 1, -1), (2, -1, -1), (4, -1, -1), (5, 1, -1), ("P1", 0, 1)),
        ((1, -1, 1), (2, 1, 1), (4, -1, 1), (5, 1, 1), ("P1", 0, -1)),
    )
    coefficient_a = -27 * (
        12 * quartic[4] * quartic[0]
        - 3 * quartic[3] * quartic[1]
        + quartic[2] ** 2
    )
    relation_text = (
        "P2=V(168,-)-V(168,+)-V(205,-)+V(205,+)-P1",
        "P3=V(0,+)+V(143,-)+V(168,+)+V(205,-)+P1",
        "P4=-V(0,-)-V(143,+)-V(168,+)-V(205,-)-P1",
        "P5=-V(25,+)-V(95,-)-V(168,-)-V(205,+)+P1",
        "P6=V(25,-)+V(95,+)+V(168,-)+V(205,+)-P1",
    )
    for target_index, relation in enumerate(relations, start=1):
        total = None
        for source, displayed_sign, group_sign in relation:
            if source == "P1":
                point = companions[0] if group_sign > 0 else negate(companions[0])
            else:
                point = visible[source, displayed_sign]
                if group_sign < 0:
                    point = negate(point)
            total = short_add(coefficient_a, total, point)
        if total != companions[target_index]:
            raise AssertionError(f"generic companion relation P{target_index + 1} failed")
        if evaluate_at_u(total, parameter_u) != evaluate_at_u(
            companions[target_index], parameter_u
        ):
            raise AssertionError(f"u=197 replay of P{target_index + 1} failed")
    return {
        "status": "exact generic D-square companion relations verified",
        "base_field": "Q(u), with T=(39146-u^2)/(2u)",
        "roots": list(ROOTS),
        "method": (
            "recursive Mestre square root, primitive quartic normalization, "
            "exact rational-factor ordinate extraction, covariant map, and short-Weierstrass addition"
        ),
        "expanded_two_section_residual_materialized": False,
        "first_companion": "P1=(-2375+37*T)/23",
        "remaining_companion_relations": list(relation_text),
        "conclusion": (
            "all six displayed affine companions lie in the generic subgroup generated by "
            "the twelve visible sections, split infinity, and P1; the other five do not "
            "supply a second generic direction on this D-square branch"
        ),
        "not_established": [
            "the rank or saturation of the displayed thirteen-generator subgroup",
            "a global component identity for the rank-six two-section moduli germ",
            "a Shioda Gram matrix or a generic rank upper bound",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

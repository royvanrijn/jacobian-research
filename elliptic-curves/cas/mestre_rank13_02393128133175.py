#!/usr/bin/env python3
"""A low-height six-root Mestre family with thirteen independent sections.

The six centers are ``(0,23,93,128,133,175)``.  Besides Mestre's twelve
paired-root points, its primitive quartic has the polynomial section

``x=(31*T+4333)/21``.

The leading quartic coefficient is ``9*(T^2+14406)``.  The conic base change

``T=(14406-u^2)/(2*u)``

therefore splits the two points at infinity.  This module supplies exact
formulas only.  The separate verifier proves that thirteen of the resulting
Jacobian sections are independent by a good specialization at ``u=1``.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)


Q = Fraction

ROOTS = (0, 23, 93, 128, 133, 175)
CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
BASE_CHANGE_CONSTANT = Q(14_406)


def primitive_quartic_coefficients(
    parameter_t: Fraction,
) -> tuple[Fraction, ...]:
    """Return the primitive quartic in ascending powers of ``x``."""

    parameter_t = Q(parameter_t)
    t = parameter_t
    coefficients = (
        9 * t**6 - 253_406 * t**4 + 1_434_086_185 * t**2
        + 7_050_150_764_944,
        12 * (276 * t**4 - 2_861_579 * t**2 - 31_373_984_992),
        -3 * (6 * t**4 - 158_072 * t**2 - 2_245_309_213),
        -36 * (92 * t**2 + 1_377_831),
        9 * (t**2 + BASE_CHANGE_CONSTANT),
    )
    if coefficients != CONSTRUCTION.primitive_quartic_coefficients(parameter_t):
        raise AssertionError("the explicit primitive quartic changed")
    return coefficients


def linear_extra_point(parameter_t: Fraction) -> tuple[Fraction, Fraction]:
    """Return one nonvisible polynomial section on the quartic."""

    parameter_t = Q(parameter_t)
    t = parameter_t
    point = (
        (31 * t + 4_333) / 21,
        (
            520 * t**3
            + 148_862 * t**2
            + 17_234_623 * t
            + 577_212_405
        )
        / 147,
    )
    if point[1] ** 2 != quartic_value(
        primitive_quartic_coefficients(parameter_t), point[0]
    ):
        raise AssertionError("the nonvisible linear section failed exactly")
    return point


def base_parameter(parameter_u: Fraction) -> Fraction:
    """Return ``T=(14406-u^2)/(2u)``."""

    parameter_u = Q(parameter_u)
    if parameter_u == 0:
        raise ValueError("the base change has a pole at u=0")
    return (BASE_CHANGE_CONSTANT - parameter_u**2) / (2 * parameter_u)


def leading_square(parameter_u: Fraction) -> Fraction:
    """Return a square root of the base-changed leading coefficient."""

    parameter_u = Q(parameter_u)
    if parameter_u == 0:
        raise ValueError("the base change has a pole at u=0")
    return 3 * (BASE_CHANGE_CONSTANT + parameter_u**2) / (2 * parameter_u)


def point_on_short_curve(
    coefficients: Sequence[Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    a1, a2, a3, a4, a6 = (Q(value) for value in coefficients)
    x_value, y_value = (Q(value) for value in point)
    return (
        y_value**2 + a1 * x_value * y_value + a3 * y_value
        == x_value**3 + a2 * x_value**2 + a4 * x_value + a6
    )


def split_infinity_jacobian_point(
    parameter_u: Fraction, *, sign: int = 1
) -> tuple[Fraction, Fraction]:
    """Return the covariant-map limit at a split quartic infinity."""

    if sign not in (-1, 1):
        raise ValueError("the infinity sign must be -1 or 1")
    parameter_u = Q(parameter_u)
    parameter_t = base_parameter(parameter_u)
    _, d, c, b, a = primitive_quartic_coefficients(parameter_t)
    square_root = leading_square(parameter_u)
    if square_root**2 != a:
        raise AssertionError("the leading quartic coefficient did not split")
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    point = (
        36 * g0 / a,
        sign * 54 * (a * g1 - b * g0) / square_root**3,
    )
    coefficients = short_jacobian_coefficients(CONSTRUCTION, parameter_t)
    if not point_on_short_curve(coefficients, point):
        raise AssertionError("the split-infinity limit missed the Jacobian")
    return point


def known_jacobian_points(
    parameter_u: Fraction,
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return twelve visible images, the extra image, and split infinity."""

    parameter_u = Q(parameter_u)
    parameter_t = base_parameter(parameter_u)
    quartic_points = primitive_visible_points(CONSTRUCTION, parameter_t) + (
        linear_extra_point(parameter_t),
    )
    affine = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, parameter_t, point)
        for point in quartic_points
    )
    answer = affine + (split_infinity_jacobian_point(parameter_u),)
    coefficients = short_jacobian_coefficients(CONSTRUCTION, parameter_t)
    if any(not point_on_short_curve(coefficients, point) for point in answer):
        raise AssertionError("a displayed Jacobian point failed exactly")
    return answer


def base_changed_short_jacobian_coefficients(
    parameter_u: Fraction,
) -> tuple[Fraction, ...]:
    return short_jacobian_coefficients(CONSTRUCTION, base_parameter(parameter_u))

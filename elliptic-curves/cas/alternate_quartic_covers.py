#!/usr/bin/env python3
"""Exact degree-two models of a short elliptic curve from rational points.

Let ``E: y^2=x^3+A*x+B`` and let ``Q=(xq,yq)`` be rational.  The function

``t = (y+yq)/(x-xq)``

has pole divisor ``O+Q``.  Eliminating ``x`` gives the binary-quartic chart

``v^2 = t^4 - 6*xq*t^2 - 8*yq*t - 3*xq^2 - 4*A``.

Different classes of ``Q`` modulo ``2E(Q)`` give the degree-two coordinate
systems that the rank search needs.  This module constructs and checks those
charts exactly; it makes no Selmer-completeness or rank claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm
from typing import Iterable, Sequence


Q = Fraction
AffinePoint = tuple[Fraction, Fraction]
EllipticPoint = AffinePoint | None
MobiusMatrix = tuple[int, int, int, int]


def point_on_short_curve(
    coefficients: Sequence[Fraction], point: AffinePoint
) -> bool:
    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        raise ValueError("expected a short Weierstrass coefficient vector")
    coefficient_a, coefficient_b = (Q(value) for value in coefficients[3:])
    x_value, y_value = (Q(value) for value in point)
    return y_value**2 == x_value**3 + coefficient_a * x_value + coefficient_b


def short_add(
    coefficients: Sequence[Fraction],
    left: EllipticPoint,
    right: EllipticPoint,
) -> EllipticPoint:
    """Exact group addition on a short Weierstrass curve."""

    if left is None:
        return right
    if right is None:
        return left
    if not point_on_short_curve(coefficients, left) or not point_on_short_curve(
        coefficients, right
    ):
        raise ValueError("a summand is not on the short curve")
    coefficient_a = Q(coefficients[3])
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right:
        if y_left == -y_right:
            return None
        if y_left == 0:
            return None
        slope = (3 * x_left**2 + coefficient_a) / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_sum = slope**2 - x_left - x_right
    y_sum = -y_left + slope * (x_left - x_sum)
    answer = (x_sum, y_sum)
    if not point_on_short_curve(coefficients, answer):
        raise AssertionError("exact group addition left the curve")
    return answer


def short_subset_sum(
    coefficients: Sequence[Fraction],
    points: Sequence[AffinePoint],
    indices: Iterable[int],
) -> EllipticPoint:
    """Add the zero-based selected points exactly."""

    answer: EllipticPoint = None
    for index in indices:
        if index < 0 or index >= len(points):
            raise IndexError("point index is outside the supplied basis")
        answer = short_add(coefficients, answer, points[index])
    return answer


def three_point_mobius_matrix(
    image_zero: Fraction,
    image_one: Fraction,
    image_infinity: Fraction,
) -> MobiusMatrix:
    """Return primitive integers for a map sending ``0,1,infinity`` as named.

    The returned matrix ``(a,b,c,d)`` represents
    ``t=(a*s+b)/(c*s+d)``.  The three images must be distinct.
    """

    t_zero, t_one, t_infinity = (
        Q(image_zero),
        Q(image_one),
        Q(image_infinity),
    )
    if len({t_zero, t_one, t_infinity}) != 3:
        raise ValueError("three distinct rational images are required")
    c_value = (t_one - t_zero) / (t_infinity - t_one)
    rational_entries = (t_infinity * c_value, t_zero, c_value, Q(1))
    common_denominator = 1
    for entry in rational_entries:
        common_denominator = lcm(common_denominator, entry.denominator)
    entries = tuple(
        int(entry * common_denominator) for entry in rational_entries
    )
    common_factor = 0
    for entry in entries:
        common_factor = gcd(common_factor, abs(entry))
    if common_factor == 0:
        raise AssertionError("the Mobius matrix vanished")
    entries = tuple(entry // common_factor for entry in entries)
    if next(entry for entry in entries if entry) < 0:
        entries = tuple(-entry for entry in entries)
    a_value, b_value, c_integer, d_value = entries
    if a_value * d_value - b_value * c_integer == 0:
        raise AssertionError("the Mobius matrix is singular")
    answer = (a_value, b_value, c_integer, d_value)
    if mobius_image(answer, Q(0)) != t_zero:
        raise AssertionError("the chart missed its zero image")
    if mobius_image(answer, Q(1)) != t_one:
        raise AssertionError("the chart missed its one image")
    if Q(a_value, c_integer) != t_infinity:
        raise AssertionError("the chart missed its infinity image")
    return answer


def mobius_image(matrix: Sequence[int], parameter: Fraction) -> Fraction | None:
    """Evaluate ``(a*s+b)/(c*s+d)``; return ``None`` at its pole."""

    if len(matrix) != 4:
        raise ValueError("a 2-by-2 matrix has four entries")
    a_value, b_value, c_value, d_value = (int(value) for value in matrix)
    if a_value * d_value - b_value * c_value == 0:
        raise ValueError("the Mobius matrix must be invertible")
    parameter = Q(parameter)
    denominator = c_value * parameter + d_value
    if denominator == 0:
        return None
    return Q(a_value * parameter + b_value, denominator)


def mobius_preimage(matrix: Sequence[int], value: Fraction) -> Fraction | None:
    """Invert an integral Mobius chart; return ``None`` for infinity."""

    if len(matrix) != 4:
        raise ValueError("a 2-by-2 matrix has four entries")
    a_value, b_value, c_value, d_value = (int(entry) for entry in matrix)
    if a_value * d_value - b_value * c_value == 0:
        raise ValueError("the Mobius matrix must be invertible")
    value = Q(value)
    denominator = value * c_value - a_value
    if denominator == 0:
        return None
    return Q(b_value - value * d_value, denominator)


@dataclass(frozen=True)
class AlternateQuarticCover:
    """The degree-two chart associated with ``O+base_point``."""

    coefficient_a: Fraction
    coefficient_b: Fraction
    base_point: AffinePoint

    def __post_init__(self) -> None:
        coefficients = self.short_coefficients
        if 4 * self.coefficient_a**3 + 27 * self.coefficient_b**2 == 0:
            raise ValueError("the short curve is singular")
        if not point_on_short_curve(coefficients, self.base_point):
            raise ValueError("the chart base point is not on the short curve")

    @property
    def short_coefficients(self) -> tuple[Fraction, ...]:
        return (Q(0), Q(0), Q(0), Q(self.coefficient_a), Q(self.coefficient_b))

    @property
    def coefficients(self) -> tuple[Fraction, ...]:
        """Quartic coefficients in ascending powers of ``t``."""

        x_base, y_base = self.base_point
        return (
            -3 * x_base**2 - 4 * self.coefficient_a,
            -8 * y_base,
            -6 * x_base,
            Q(0),
            Q(1),
        )

    def value(self, parameter: Fraction) -> Fraction:
        parameter = Q(parameter)
        answer = Q(0)
        for coefficient in reversed(self.coefficients):
            answer = answer * parameter + coefficient
        return answer

    def cover_point_to_curve(self, point: AffinePoint) -> AffinePoint:
        """Map an affine quartic point to ``E`` and verify it exactly."""

        parameter, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.value(parameter):
            raise ValueError("the supplied point is not on the quartic cover")
        x_base, y_base = self.base_point
        x_value = (parameter**2 - x_base + ordinate) / 2
        y_value = parameter * (x_value - x_base) - y_base
        answer = (x_value, y_value)
        if not point_on_short_curve(self.short_coefficients, answer):
            raise AssertionError("the exact quartic map left the elliptic curve")
        return answer

    def curve_point_to_cover(self, point: AffinePoint) -> AffinePoint:
        """Map an affine point other than the chart base point to the quartic."""

        if not point_on_short_curve(self.short_coefficients, point):
            raise ValueError("the supplied point is not on the elliptic curve")
        x_value, y_value = (Q(value) for value in point)
        x_base, y_base = self.base_point
        if x_value == x_base and y_value == y_base:
            raise ValueError("the chart base point maps to infinity")
        if x_value == x_base:
            if y_base == 0 or y_value != -y_base:
                raise ValueError("the exceptional affine point is not -Q")
            parameter = -(3 * x_base**2 + self.coefficient_a) / (2 * y_base)
        else:
            parameter = (y_value + y_base) / (x_value - x_base)
        ordinate = 2 * x_value + x_base - parameter**2
        answer = (parameter, ordinate)
        if ordinate**2 != self.value(parameter):
            raise AssertionError("the inverse chart map left the quartic")
        if self.cover_point_to_curve(answer) != (x_value, y_value):
            raise AssertionError("the quartic round trip changed the point")
        return answer


def alternate_cover(
    coefficients: Sequence[Fraction], base_point: AffinePoint
) -> AlternateQuarticCover:
    """Construct the exact chart from a short coefficient vector and point."""

    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        raise ValueError("expected a short Weierstrass coefficient vector")
    return AlternateQuarticCover(
        Q(coefficients[3]),
        Q(coefficients[4]),
        (Q(base_point[0]), Q(base_point[1])),
    )

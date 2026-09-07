"""Exact rational arithmetic for general Weierstrass models."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Iterable, Sequence


Point = tuple[Fraction, Fraction] | None


def _q(value: int | str | Fraction) -> Fraction:
    return Fraction(value)


@dataclass(frozen=True)
class EllipticCurve:
    """The curve ``y^2+a1*x*y+a3*y=x^3+a2*x^2+a4*x+a6``."""

    coefficients: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]

    def __init__(self, coefficients: Sequence[int | str | Fraction]):
        if len(coefficients) != 5:
            raise ValueError("a general Weierstrass model has five coefficients")
        object.__setattr__(self, "coefficients", tuple(map(_q, coefficients)))

    def is_on_curve(self, point: Point) -> bool:
        if point is None:
            return True
        x_value, y_value = map(_q, point)
        a1, a2, a3, a4, a6 = self.coefficients
        return (
            y_value * y_value + a1 * x_value * y_value + a3 * y_value
            == x_value**3 + a2 * x_value**2 + a4 * x_value + a6
        )

    def negate(self, point: Point) -> Point:
        if point is None:
            return None
        x_value, y_value = point
        a1, _a2, a3, _a4, _a6 = self.coefficients
        return x_value, -y_value - a1 * x_value - a3

    def add(self, left: Point, right: Point) -> Point:
        if left is None:
            return right
        if right is None:
            return left
        x1, y1 = map(_q, left)
        x2, y2 = map(_q, right)
        a1, a2, a3, a4, a6 = self.coefficients
        if x1 == x2:
            if y2 == -y1 - a1 * x1 - a3:
                return None
            denominator = 2 * y1 + a1 * x1 + a3
            if denominator == 0:
                return None
            slope = (3 * x1**2 + 2 * a2 * x1 + a4 - a1 * y1) / denominator
            intercept = (-x1**3 + a4 * x1 + 2 * a6 - a3 * y1) / denominator
        else:
            slope = (y2 - y1) / (x2 - x1)
            intercept = (y1 * x2 - y2 * x1) / (x2 - x1)
        x3 = slope**2 + a1 * slope - a2 - x1 - x2
        y3 = -(slope + a1) * x3 - intercept - a3
        answer = Fraction(x3), Fraction(y3)
        if not self.is_on_curve(answer):
            raise ArithmeticError("elliptic addition produced an off-curve point")
        return answer

    def multiply(self, point: Point, multiplier: int) -> Point:
        multiplier = int(multiplier)
        if multiplier < 0:
            return self.multiply(self.negate(point), -multiplier)
        answer: Point = None
        addend = point
        while multiplier:
            if multiplier & 1:
                answer = self.add(answer, addend)
            multiplier >>= 1
            if multiplier:
                addend = self.add(addend, addend)
        return answer

    def linear_combination(
        self, points: Sequence[Point], coordinates: Sequence[int]
    ) -> Point:
        if len(points) != len(coordinates):
            raise ValueError("point and coordinate counts differ")
        answer: Point = None
        for point, coefficient in zip(points, coordinates):
            if coefficient:
                answer = self.add(answer, self.multiply(point, int(coefficient)))
        return answer


def point_complexity(point: Point) -> dict[str, int | bool]:
    """Return exact, model-independent-of-format rational-size metadata."""

    if point is None:
        return {
            "infinity": True,
            "integral": False,
            "x_numerator_bits": 0,
            "x_denominator_bits": 0,
            "y_numerator_bits": 0,
            "y_denominator_bits": 0,
            "total_bits": 0,
        }
    x_value, y_value = map(Fraction, point)
    fields = {
        "infinity": False,
        "integral": x_value.denominator == y_value.denominator == 1,
        "x_numerator_bits": abs(x_value.numerator).bit_length(),
        "x_denominator_bits": x_value.denominator.bit_length(),
        "y_numerator_bits": abs(y_value.numerator).bit_length(),
        "y_denominator_bits": y_value.denominator.bit_length(),
    }
    fields["total_bits"] = sum(
        int(fields[key])
        for key in (
            "x_numerator_bits",
            "x_denominator_bits",
            "y_numerator_bits",
            "y_denominator_bits",
        )
    )
    return fields


def denominator_square_root(point: Point) -> int | None:
    """Return ``d`` for ``den(x)=d^2``, checking the elliptic denominator law."""

    if point is None:
        return None
    denominator = Fraction(point[0]).denominator
    root = int(denominator**0.5)
    while (root + 1) ** 2 <= denominator:
        root += 1
    while root * root > denominator:
        root -= 1
    if root * root != denominator:
        raise ArithmeticError("x-coordinate denominator is not a square")
    return root

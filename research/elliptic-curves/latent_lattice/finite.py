"""Exact finite-reduction quotient codes for general Weierstrass models."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .elliptic import EllipticCurve, Point


FinitePoint = tuple[int, int] | None


def _reduce(value: Fraction | int, prime: int) -> int:
    value = Fraction(value)
    if value.denominator % prime == 0:
        raise ValueError(f"denominator is not a unit modulo {prime}")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def discriminant(coefficients: Sequence[Fraction | int]) -> Fraction:
    """Return the exact discriminant of a general Weierstrass equation."""

    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4**3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


def finite_negate(point: FinitePoint, coefficients: Sequence[int], prime: int) -> FinitePoint:
    if point is None:
        return None
    x_value, y_value = point
    a1, _a2, a3, _a4, _a6 = coefficients
    return x_value, (-y_value - a1 * x_value - a3) % prime


def finite_add(
    left: FinitePoint,
    right: FinitePoint,
    coefficients: Sequence[int],
    prime: int,
) -> FinitePoint:
    """Add points on a nonsingular general Weierstrass curve over F_p."""

    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    a1, a2, a3, a4, a6 = coefficients
    if x1 == x2:
        if (y2 + y1 + a1 * x1 + a3) % prime == 0:
            return None
        denominator = (2 * y1 + a1 * x1 + a3) % prime
        if denominator == 0:
            return None
        slope = (
            (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1)
            * pow(denominator, -1, prime)
        ) % prime
        intercept = (
            (-x1**3 + a4 * x1 + 2 * a6 - a3 * y1)
            * pow(denominator, -1, prime)
        ) % prime
    else:
        inverse = pow((x2 - x1) % prime, -1, prime)
        slope = (y2 - y1) * inverse % prime
        intercept = (y1 * x2 - y2 * x1) * inverse % prime
    x3 = (slope * slope + a1 * slope - a2 - x1 - x2) % prime
    y3 = (-(slope + a1) * x3 - intercept - a3) % prime
    return x3, y3


def finite_multiply(
    point: FinitePoint,
    multiplier: int,
    coefficients: Sequence[int],
    prime: int,
) -> FinitePoint:
    multiplier = int(multiplier)
    if multiplier < 0:
        return finite_multiply(
            finite_negate(point, coefficients, prime), -multiplier, coefficients, prime
        )
    answer: FinitePoint = None
    addend = point
    while multiplier:
        if multiplier & 1:
            answer = finite_add(answer, addend, coefficients, prime)
        multiplier >>= 1
        if multiplier:
            addend = finite_add(addend, addend, coefficients, prime)
    return answer


def finite_curve_points(coefficients: Sequence[int], prime: int) -> tuple[FinitePoint, ...]:
    a1, a2, a3, a4, a6 = coefficients
    points: list[FinitePoint] = [None]
    for x_value in range(prime):
        right = (x_value**3 + a2 * x_value**2 + a4 * x_value + a6) % prime
        for y_value in range(prime):
            left = (y_value**2 + a1 * x_value * y_value + a3 * y_value) % prime
            if left == right:
                points.append((x_value, y_value))
    return tuple(points)


@dataclass(frozen=True)
class FiniteQuotientBlock:
    reduction_prime: int
    relation_prime: int
    group_order: int
    multiple_subgroup_order: int
    quotient_dimension: int
    rows: tuple[tuple[int, ...], ...]

    def vector_class(self, coordinates: Sequence[int]) -> tuple[int, ...]:
        if self.rows and len(coordinates) != len(self.rows[0]):
            raise ValueError("coordinate width differs from quotient-code width")
        return tuple(
            sum(value * coefficient for value, coefficient in zip(row, coordinates))
            % self.relation_prime
            for row in self.rows
        )

    def to_record(self) -> dict[str, object]:
        return {
            "reduction_prime": self.reduction_prime,
            "relation_prime": self.relation_prime,
            "group_order": self.group_order,
            "multiple_subgroup_order": self.multiple_subgroup_order,
            "quotient_dimension": self.quotient_dimension,
            "rows": [list(row) for row in self.rows],
        }


def finite_quotient_block(
    curve: EllipticCurve,
    points: Sequence[Point],
    reduction_prime: int,
    relation_prime: int,
) -> FiniteQuotientBlock:
    """Compute the exact map of supplied points into E(F_p)/ell E(F_p)."""

    prime = int(reduction_prime)
    ell = int(relation_prime)
    if prime <= 2 or ell <= 1 or prime == ell:
        raise ValueError("use distinct primes, with reduction prime greater than two")
    coefficients = tuple(_reduce(value, prime) for value in curve.coefficients)
    if _reduce(discriminant(curve.coefficients), prime) == 0:
        raise ValueError("the chosen prime has bad reduction")
    if any(point is None for point in points):
        raise ValueError("finite quotient input must contain finite points")
    reduced = tuple(
        (_reduce(point[0], prime), _reduce(point[1], prime))  # type: ignore[index]
        for point in points
    )
    group = finite_curve_points(coefficients, prime)
    multiples = {
        finite_multiply(point, ell, coefficients, prime) for point in group
    }
    representatives: list[FinitePoint] = [None]
    representative_coordinates: list[tuple[int, ...]] = [()]
    quotient_basis: list[FinitePoint] = []
    for point in group:
        if any(
            finite_add(
                point,
                finite_negate(representative, coefficients, prime),
                coefficients,
                prime,
            )
            in multiples
            for representative in representatives
        ):
            continue
        quotient_basis.append(point)
        old_representatives = tuple(representatives)
        old_coordinates = tuple(representative_coordinates)
        representatives = []
        representative_coordinates = []
        for scalar in range(ell):
            multiple = finite_multiply(point, scalar, coefficients, prime)
            for representative, coordinate in zip(
                old_representatives, old_coordinates
            ):
                representatives.append(
                    finite_add(representative, multiple, coefficients, prime)
                )
                representative_coordinates.append(coordinate + (scalar,))
    if len(representatives) * len(multiples) != len(group):
        raise ArithmeticError("quotient representatives do not cover E(F_p)")
    rows = [[0] * len(points) for _ in quotient_basis]
    for point_index, point in enumerate(reduced):
        coordinate_index = next(
            (
                index
                for index, representative in enumerate(representatives)
                if finite_add(
                    point,
                    finite_negate(representative, coefficients, prime),
                    coefficients,
                    prime,
                )
                in multiples
            ),
            None,
        )
        if coordinate_index is None:
            raise ArithmeticError("a rational point missed every quotient coset")
        for row, value in zip(rows, representative_coordinates[coordinate_index]):
            row[point_index] = value
    return FiniteQuotientBlock(
        reduction_prime=prime,
        relation_prime=ell,
        group_order=len(group),
        multiple_subgroup_order=len(multiples),
        quotient_dimension=len(quotient_basis),
        rows=tuple(tuple(row) for row in rows),
    )

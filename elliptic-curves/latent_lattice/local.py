"""Exact component codes at declared multiplicative reduction places."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .elliptic import EllipticCurve, Point


def valuation(value: Fraction | int, prime: int) -> int:
    value = Fraction(value)
    if value == 0:
        return 10**9
    numerator = abs(value.numerator)
    denominator = value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def singular_reduction(curve: EllipticCurve, point: Point, prime: int) -> bool:
    if point is None:
        return False
    x_value, y_value = point
    if valuation(x_value, prime) < 0 or valuation(y_value, prime) < 0:
        return False
    a1, a2, a3, a4, _a6 = curve.coefficients
    derivative_x = a1 * y_value - 3 * x_value**2 - 2 * a2 * x_value - a4
    derivative_y = 2 * y_value + a1 * x_value + a3
    return valuation(derivative_x, prime) > 0 and valuation(derivative_y, prime) > 0


def multiplicative_depth(
    curve: EllipticCurve, point: Point, prime: int, fibre_order: int
) -> int:
    if not singular_reduction(curve, point, prime):
        return 0
    assert point is not None
    x_value, y_value = point
    a1, _a2, a3, _a4, _a6 = curve.coefficients
    return min(
        valuation(2 * y_value + a1 * x_value + a3, prime),
        fibre_order // 2,
    )


@dataclass(frozen=True)
class ComponentBlock:
    prime: int
    kodaira_symbol: str
    modulus: int
    classes: tuple[int, ...]
    orientation: str

    def vector_class(self, coordinates: Sequence[int]) -> int:
        if len(coordinates) != len(self.classes):
            raise ValueError("coordinate width differs from component-code width")
        return sum(a * b for a, b in zip(coordinates, self.classes)) % self.modulus

    def to_record(self) -> dict[str, object]:
        return {
            "prime": self.prime,
            "kodaira_symbol": self.kodaira_symbol,
            "modulus": self.modulus,
            "classes": list(self.classes),
            "orientation": self.orientation,
        }


def multiplicative_component_block(
    curve: EllipticCurve,
    points: Sequence[Point],
    *,
    prime: int,
    fibre_order: int,
    split: bool,
) -> ComponentBlock:
    """Resolve a homomorphic component code for a declared I_n fibre.

    The caller is responsible for supplying a locally minimal model and
    independently certified Kodaira type.  For split reduction the first
    ambiguous supplied point fixes the unavoidable global sign.  Every pair
    sum is replayed before the code is returned.
    """

    if fibre_order < 1:
        raise ValueError("fibre_order must be positive")
    if not split:
        modulus = 2 if fibre_order % 2 == 0 else 1
        values = tuple(
            int(singular_reduction(curve, point, prime)) % modulus
            if modulus > 1
            else 0
            for point in points
        )
        if modulus > 1:
            for left, left_value in zip(points, values):
                for right, right_value in zip(points, values):
                    observed = int(singular_reduction(curve, curve.add(left, right), prime))
                    if (left_value + right_value) % modulus != observed:
                        raise ArithmeticError("nonsplit component pair-sum replay failed")
        return ComponentBlock(
            prime, f"I{fibre_order}", modulus, values, "canonical nonsplit code"
        )

    modulus = fibre_order
    depths = [
        multiplicative_depth(curve, point, prime, fibre_order) for point in points
    ]
    options = [
        {depth % modulus, (-depth) % modulus} for depth in depths
    ]
    values: list[int | None] = [
        next(iter(option)) if len(option) == 1 else None for option in options
    ]
    anchor = next((i for i, option in enumerate(options) if len(option) == 2), None)
    if anchor is not None:
        values[anchor] = depths[anchor]
    changed = True
    while changed:
        changed = False
        for index, option in enumerate(options):
            if values[index] is not None:
                continue
            valid = []
            for candidate in option:
                if all(
                    other_value is None
                    or min(
                        (candidate + other_value) % modulus,
                        (-candidate - other_value) % modulus,
                    )
                    == multiplicative_depth(
                        curve, curve.add(points[index], points[other]), prime, fibre_order
                    )
                    for other, other_value in enumerate(values)
                ):
                    valid.append(candidate)
            if len(valid) == 1:
                values[index] = valid[0]
                changed = True
    if any(value is None for value in values):
        raise ArithmeticError("component orientations were not resolved")
    oriented = tuple(int(value) for value in values)
    for left, left_value in enumerate(oriented):
        for right, right_value in enumerate(oriented):
            residue = (left_value + right_value) % modulus
            expected = min(residue, (-residue) % modulus)
            observed = multiplicative_depth(
                curve, curve.add(points[left], points[right]), prime, fibre_order
            )
            if expected != observed:
                raise ArithmeticError("split component pair-sum replay failed")
    return ComponentBlock(
        prime,
        f"I{fibre_order}",
        modulus,
        oriented,
        "first ambiguous supplied point positive; all pair sums replayed",
    )

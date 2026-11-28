"""Directed rational enclosure for the four-cell cone in the FR4 benchmark."""

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    @classmethod
    def exact(cls, value):
        return cls(Fraction(value), Fraction(value))

    def __add__(self, other):
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        products = [x * y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    def __truediv__(self, other):
        assert other.lo * other.hi > 0
        return self * Interval(1 / other.hi, 1 / other.lo)

    def square(self):
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return Interval(Fraction(lower), max(self.lo**2, self.hi**2))

    def sqrt(self):
        assert self.lo >= 0
        denominator = 1 << 64

        def round_down(value):
            integer = (value.numerator * denominator**2) // value.denominator
            return Fraction(isqrt(integer), denominator)

        return Interval(round_down(self.lo), round_down(self.hi) + Fraction(1, denominator))


I = Interval.exact


def small_atan(argument):
    """Bounds using 28 terms of the alternating series at each endpoint."""

    def series(value):
        assert 0 <= value <= Fraction(3, 4)
        term = value
        partial = Fraction(0)
        for degree in range(28):
            partial += term / (2 * degree + 1)
            term *= -(value**2)
        return partial, partial + term / 57

    return Interval(series(argument.lo)[0], series(argument.hi)[1])


PI = I(16) * small_atan(I(Fraction(1, 5))) - I(4) * small_atan(I(Fraction(1, 239)))


def atan(argument):
    if argument.lo >= 1:
        return PI / I(2) - small_atan(I(1) / argument)
    return small_atan(argument)


def third_weight(weights, first_moments):
    def mixed_coefficient(r, s):
        variance = I(1 + r * r)
        root = variance.sqrt()
        normalizer = (I(2) * PI).sqrt()
        sector = I(Fraction(1, 4)) - atan(I(s) / root) / (I(2) * PI)
        correction = I(s) * normalizer / (I(2) * PI * variance * (variance + I(s * s)))
        radial_second = normalizer * sector / (variance * root) - correction
        return -I(r) * radial_second / PI

    total = I(0)
    for index, weight in enumerate(weights):
        other_weights = [value for other, value in enumerate(weights) if other != index]
        r, s = (weight / other for other in other_weights)
        left = mixed_coefficient(r, s)
        right = mixed_coefficient(s, r)
        pure = -(first_moments[index] + left + right)
        total += I(Fraction(2, 3)) * pure.square() + I(2) * (left.square() + right.square())
    return total


def main():
    inverse_root_two = I(Fraction(1, 2)).sqrt()
    inverse_root_five = I(Fraction(1, 5)).sqrt()
    density_zero = I(1) / (I(2) * PI).sqrt()
    tangent = I(Fraction(1, 2)) - atan(inverse_root_five) / PI
    z = I(2) * density_zero * inverse_root_five * tangent
    axis_difference = atan(inverse_root_two) / (I(2) * PI * PI.sqrt())
    x = z + axis_difference
    cone_first = I(4) * (I(2) * x.square() + z.square())
    cone_defect = I(8) * axis_difference.square()

    tied_mass = I(Fraction(1, 4)) + atan(inverse_root_two / I(2)) / (I(2) * PI)
    tetra_centroid = tied_mass / PI.sqrt()
    tetra_first = I(12) * tetra_centroid.square()
    shared_second = I(1) / PI.square()
    tetra_third = third_weight((Fraction(1),) * 3, (tetra_centroid,) * 3)
    cone_third = third_weight((Fraction(1), Fraction(1), Fraction(1, 2)), (x, x, z))
    tetra_fourth = I(Fraction(10, 27)) * shared_second

    first_gap = tetra_first - cone_first
    third_gain = cone_third - tetra_third
    variance = I(Fraction(3, 4))
    remainder_two = variance - tetra_first - shared_second - tetra_third
    remainder_three = remainder_two - tetra_fourth
    margin_two = remainder_two / I(4) - I(Fraction(3, 2)) * first_gap
    margin_three = remainder_three / I(8) - I(Fraction(7, 4)) * first_gap + third_gain / I(4)

    bounds = {
        "degree_one_defect_lower": (cone_defect, Fraction(1, 50)),
        "degree_two_bracket_at_half_lower": (margin_two, Fraction(11, 250)),
        "degree_three_bracket_at_half_lower": (margin_three, Fraction(7, 500)),
        "cone_minus_tetra_degree_three_weight_lower": (third_gain, Fraction(13, 2500)),
        "tetra_remainder_after_four_lower": (remainder_three, Fraction(1, 5)),
    }
    for interval, bound in bounds.values():
        assert interval.lo > bound

    print(json.dumps({
        "arithmetic": "exact Fraction arithmetic with 64-bit outward square roots and alternating arctangent bounds",
        "correlation_domain": "1/2 <= rho < 1",
        "reflection_tolerance_domain": "0 < epsilon <= 1/50",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "strict_bounds": {name: str(bound) for name, (_, bound) in bounds.items()},
        "status": "verified",
    }, indent=2))


if __name__ == "__main__":
    main()

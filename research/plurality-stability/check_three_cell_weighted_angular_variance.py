"""Exact rational checks for the three lower bounds in WV5."""

import json
from fractions import Fraction
from math import comb, factorial


def scaled_cosine_integral(frequency):
    """Return pi squared times the signed Fourier integral in SV4."""
    if frequency == 0:
        return Fraction(0)
    cosine = Fraction(1) if frequency % 3 == 0 else Fraction(-1, 2)
    return (9 * cosine - 1 - 8 * (-1)**frequency) / frequency**2


def scaled_moment(order):
    return sum(
        comb(order, j) * scaled_cosine_integral(abs(order - 2 * j))
        for j in range(order + 1)
    ) / 2**order


def minimum_bernstein_coefficient(coefficients, left, right):
    degree = len(coefficients) - 1
    shifted = [
        (right - left)**j * sum(
            coefficients[k] * comb(k, j) * left**(k - j)
            for k in range(j, degree + 1)
        )
        for j in range(degree + 1)
    ]
    return min(
        sum(shifted[j] * Fraction(comb(i, j), comb(degree, j)) for j in range(i + 1))
        for i in range(degree + 1)
    )


def verify_centered_partner():
    pi_lo, pi_hi = Fraction("3.14159"), Fraction("3.14160")
    root2_lo, root2_hi = Fraction("1.41421"), Fraction("1.41422")
    root3_lo, root3_hi = Fraction("1.73205"), Fraction("1.73206")
    cosine_pi8_hi = Fraction("0.92389")
    assert root2_lo**2 < 2 < root2_hi**2
    assert root3_lo**2 < 3 < root3_hi**2
    assert 4 * cosine_pi8_hi**2 > 2 + root2_hi

    h = Fraction(5, 12)
    negative = 2 * (-(1 + root2_lo) / (2 * pi_hi**2) + h * root3_hi / (2 * pi_lo)) / h**2
    positive_half = 2 * ((cosine_pi8_hi - Fraction(1, 2)) / pi_lo**2
                         - h * root3_lo / (4 * pi_hi)) / (h/2)**2
    positive = -Fraction(12, 25)
    ratio = max((2 * negative + positive_half) / 3, (negative + positive) / 2, positive)
    assert (-Fraction(3, 2) / pi_lo**2 - ratio) / 4 > Fraction(3, 100)


def verify():
    polynomial = [scaled_moment(order) / factorial(order) for order in range(1, 17)]
    assert polynomial[0] == Fraction(5, 2)

    # Bounds after division by a and multiplication by pi squared.
    # In each cubic coefficient of WV5, replace pi squared by 10.
    bounds = (
        (Fraction(48, 25), Fraction(1, 500)),
        (Fraction(3, 2), Fraction(21, 200)),
        (Fraction(7, 5), Fraction(14, 125)),
    )
    floors = []
    for linear, cubic in bounds:
        coefficients = polynomial.copy()
        coefficients[0] -= linear
        coefficients[2] -= cubic
        minimum = min(
            minimum_bernstein_coefficient(coefficients, Fraction(i, 16), Fraction(i + 1, 16))
            for i in range(128)
        )
        assert minimum > 0
        floors.append((10**6 * minimum).__floor__())

    exponential_at_eight = sum(Fraction(28, 5)**j / factorial(j) for j in range(13))
    third_bound_at_eight = Fraction(7, 5) * (Fraction(8, 9) + Fraction(8**3, 125))
    assert (exponential_at_eight - 1) / 32 > third_bound_at_eight
    for linear, cubic in (
        (Fraction(48, 25), Fraction(1, 5000)),
        (Fraction(3, 2), Fraction(21, 2000)),
    ):
        assert Fraction(14, 1250) > cubic
        assert (Fraction(14, 1250) - cubic) * 8**2 > (linear - Fraction(7, 5)) / 9
    verify_centered_partner()

    return {
        "status": "pass",
        "polynomial_degree": 15,
        "intervals": 128,
        "interval_width": "1/16",
        "minimum_bernstein_floors_in_millionths": floors,
        "tail_boundary": 8,
        "one_sided_positive_taylor_coefficient": "3/100",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))

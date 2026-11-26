"""Exact coefficient certificate for the final third-band bilinear tie."""

from collections import defaultdict
from fractions import Fraction
from math import comb


def dyadic_valuation(value):
    magnitude = abs(value)
    assert magnitude > 0
    return (magnitude & -magnitude).bit_length() - 1


def binomial(n, r):
    return comb(n, r) if 0 <= r <= n else 0


def add_shifted_power(coefficients, scale, shift, exponent):
    """Add scale * Z**shift * (1-Z)**exponent."""
    for degree in range(exponent + 1):
        coefficients[shift + degree] += scale * (-1) ** degree * comb(exponent, degree)


def large_order_numerator(j):
    """Return the coefficients of the numerator in equation (5)."""
    numerator = defaultdict(int)
    for r in range(j + 1):
        scalar = comb(j + r, r)
        add_shifted_power(numerator, scalar * binomial(j - 1, j - r), j - r + 2, r)
        add_shifted_power(numerator, scalar * binomial(j + 1, j - r), j - r + 2, r + 1)
    for r in range(j + 2):
        scalar = 2 * (j + 1) * comb(j + 1, r) * comb(j + 1 + r, r)
        add_shifted_power(numerator, 21 * scalar, j + 1 - r, r)
        add_shifted_power(numerator, -scalar, j + 2 - r, r)
    return numerator


def moment_coefficient(n, j, order):
    """Return [T**n] F_j(T) in the Z[q] basis for q**order = -1."""
    coefficients = [0] * order
    for r in range(n - j + 1):
        quotient, remainder = divmod(r, order)
        value = comb(j + r, j) * comb(n - r, j)
        coefficients[remainder] += (-1) ** quotient * value
    return coefficients


def combined_coefficient(k, j, order, calibration=20):
    terminal = moment_coefficient(k, j, order)
    middle = moment_coefficient(k - 1, j + 1, order)
    q_middle = [-middle[-1], *middle[:-1]]
    return [
        terminal[i] - (j + 1) * (calibration * middle[i] + q_middle[i])
        for i in range(order)
    ]


def cyclotomic_valuation(coefficients, order):
    """Change to epsilon=1-q, an Eisenstein uniformizer of value 1/order."""
    valuations = []
    for degree in range(order):
        coefficient = sum(
            (-1) ** degree * comb(i, degree) * coefficients[i]
            for i in range(degree, order)
        )
        if coefficient:
            valuations.append(Fraction(dyadic_valuation(coefficient)) + Fraction(degree, order))
    return min(valuations) if valuations else None


def scalar_weight(j, order):
    return Fraction(j, 8) - j.bit_count() + Fraction(j, order)


def main():
    coefficient_controls = 0
    for j in range(1, 32):
        for degree, coefficient in large_order_numerator(j).items():
            if coefficient == 0:
                continue
            intercept = Fraction(dyadic_valuation(coefficient)) + Fraction(j, 8) - j.bit_count()
            slope = degree - j - 3
            assert intercept > 0, (j, degree, intercept)
            assert intercept + Fraction(slope, 64) > 0, (j, degree, intercept, slope)
            coefficient_controls += 1

    cyclotomic_controls = 0
    for order, last_cross_count in ((16, 15), (32, 31)):
        for k in range(order, 512, 2 * order):
            for j in range(1, last_cross_count + 1):
                valuation = cyclotomic_valuation(combined_coefficient(k, j, order), order)
                assert valuation is None or valuation + scalar_weight(j, order) > 0, (
                    order, k, j, valuation
                )
                cyclotomic_controls += 1

    wrong = cyclotomic_valuation(combined_coefficient(16, 6, 16, calibration=19), 16)
    assert wrong is not None and wrong + scalar_weight(6, 16) < 0
    print(f"{coefficient_controls} polynomial coefficient bounds, "
          f"{cyclotomic_controls} cyclotomic controls, one rejected calibration")


if __name__ == "__main__":
    main()

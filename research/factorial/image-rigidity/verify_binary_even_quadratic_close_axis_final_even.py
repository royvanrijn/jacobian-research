#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial


def valuation(value):
    if not value:
        return float("inf")
    ratio = Fraction(value)
    numerator, denominator = abs(ratio.numerator), ratio.denominator
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


def eta(order, k, d):
    return Fraction(k, d) + (order - k).bit_count() - order.bit_count()


def check_radial_blocks():
    shift_checks = 0
    coefficient_checks = 0
    for e in range(2, 8):
        d = 2**e
        for high in range(1, 8):
            n = 2 * d * high
            j = n - 1
            minimum_odd = 2 - Fraction(1, d) + valuation(high)
            for k in range(1, n + 1):
                value = eta(n, k, d)
                assert value == valuation(n) + Fraction(1, d) - 1 + eta(j, k - 1, d)
                assert value >= 1 + valuation(high)
                assert (value == 1 + valuation(high)) == (k in (d, 2 * d))
                if k % 2:
                    assert value >= minimum_odd
                    if k in (d - 1, 2 * d - 1):
                        assert value == minimum_odd
                    else:
                        assert value >= minimum_odd + Fraction(1, 2)
                depth = valuation(k)
                if 1 <= depth < e:
                    assert value >= 2 - Fraction(2**depth, d) + valuation(high)
                if (k - 1) % (2 * d) < d:
                    assert value == eta(n, k + d, d)
                shift_checks += 1
            radial = Fraction(
                2 * comb(n, 2 * d - 1) * factorial(4 * d - 1),
                4**d * comb(n, d - 1) * factorial(2 * d - 1),
            )
            axis = Fraction(
                2 * comb(j, 2 * d - 1) * factorial(4 * d - 2),
                4**d * comb(j, d - 1) * factorial(2 * d - 2),
            )
            assert valuation(radial) == valuation(axis) == 0
            assert valuation(radial - axis) >= 2
            coefficient_checks += 1
    return shift_checks, coefficient_checks


def check_pair_counts():
    checks = 0
    for e in range(2, 7):
        d = 2**e
        for k in range(1, d):
            for m in range(1, d, 2):
                for curvature in range(m // 2 + 1):
                    linear = m - 2 * curvature
                    def coefficient(total):
                        if total < curvature:
                            return 0
                        return comb(total, curvature) * comb(total - curvature, linear)
                    a, b = coefficient(k), coefficient(k + d)
                    depth = valuation(k)
                    assert valuation(a) >= depth and valuation(b) >= depth
                    assert valuation(a - b) >= depth + 1
                    if m < d // 2:
                        assert (a - b) % 4 == 0
                    checks += 1
    return checks


def beta(m, h):
    return m * h - valuation(m + 1)


def check_beta():
    checks = 0
    for e in range(2, 10):
        d = 2**e
        r = Fraction(1, d)
        if d <= 16:
            h = 2 * r
        elif d == 32:
            h = Fraction(1, 16)
        else:
            h = (Fraction(3, 2) + r) / 31
        for m in range(1, 8 * d):
            value = beta(m, h)
            if d <= 16:
                assert value >= -Fraction(17, 8)
                if m >= d:
                    assert value > -(2 - r)
            elif d == 32:
                assert value >= -Fraction(49, 16)
                if 16 <= m < 32 and m != 31:
                    assert value >= -2
                if m >= 32 and m != 63:
                    assert value > -(2 - r)
            else:
                assert value >= -(Fraction(7, 2) - r)
                if m >= d // 2:
                    assert value > -(3 - r)
                if m >= d:
                    assert value > -(2 - r)
            checks += 1
    return checks


def main():
    shifts, coefficients = check_radial_blocks()
    pairs = check_pair_counts()
    beta_checks = check_beta()
    print(
        f"{shifts:,} factorial shifts, {coefficients} coefficient-ratio checks, "
        f"{pairs:,} paired counts, {beta_checks:,} dyadic inequalities"
    )


if __name__ == "__main__":
    main()

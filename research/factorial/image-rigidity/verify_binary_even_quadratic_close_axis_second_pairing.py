#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial


def valuation(value):
    if not value:
        return float("inf")
    ratio = Fraction(value)
    numerator, denominator = abs(ratio.numerator), ratio.denominator
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


def check_axis_pairs():
    crossing_checks = 0
    gaps = 0
    radial_checks = 0
    for e in range(3, 8):
        d = 2**e
        for tail in (1, 2, 3):
            n = 2 * d * tail
            j = n - 1
            axes = [comb(j, k) * factorial(2 * k) for k in range(j + 1)]
            u_star = Fraction(2 * axes[2 * d - 1], 4**d * axes[d - 1])
            a0, b0 = d - 1, d // 2 - 1
            for k in range(j + 1):
                block, remainder = divmod(k, 2 * d)
                if remainder >= d:
                    continue
                unit = Fraction(2 * axes[k + d], 4**d * axes[k])
                assert valuation(unit) == valuation(u_star) == 0
                if block == 0:
                    expected = u_star if k >= d // 2 else -u_star
                    assert valuation(unit - expected) >= 2
                    crossing_checks += 1
                axis_gap = (
                    Fraction(k - a0, d)
                    + (j - k).bit_count()
                    - (j - a0).bit_count()
                )
                error = axis_gap + valuation(unit - u_star)
                if k == b0:
                    assert error == Fraction(3, 2)
                else:
                    assert error >= 2
                gaps += 1
            for k in (a0, b0):
                b_ratio = Fraction(2 * axes[k + d], 4**d * axes[k])
                c_ratio = Fraction(
                    2 * comb(n, k + d) * factorial(2 * (k + d) + 1),
                    4**d * comb(n, k) * factorial(2 * k + 1),
                )
                assert valuation(b_ratio - c_ratio) >= 2
                radial_checks += 1
            radial_ratio = Fraction(
                comb(n, b0) * factorial(2 * b0 + 1),
                comb(n, a0) * factorial(2 * a0 + 1),
            )
            coefficient_ratio = Fraction(axes[b0], axes[a0])
            assert valuation(radial_ratio / coefficient_ratio - 1) >= 2
            radial_checks += 1
    return crossing_checks, gaps, radial_checks


def check_small_counts():
    count = 0
    for e in range(3, 7):
        d = 2**e
        for k in range(2, d, 4):
            for m in range(1, d // 2, 2):
                for curvature in range(m // 2 + 1):
                    linear = m - 2 * curvature
                    def coefficient(total):
                        if total < curvature:
                            return 0
                        return comb(total, curvature) * comb(total - curvature, linear)
                    assert (coefficient(k + d) - coefficient(k)) % 8 == 0
                    count += 1
        for m in range(1, d // 2, 2):
            for curvature in range(m // 2 + 1):
                linear = m - 2 * curvature
                def coefficient(total):
                    return comb(total, curvature) * comb(total - curvature, linear)
                assert (coefficient(d - 1) - coefficient(d // 2 - 1)) % 2 == 0
                count += 1
    return count


def beta(m, h):
    return m * h - valuation(m + 1)


def check_strengths():
    thresholds = 0
    frontiers = 0
    for e in range(3, 10):
        d = 2**e
        r = Fraction(1, d)
        for tau in range(e + 2):
            for c in (Fraction(7, 2), 4):
                strengths = [
                    (s - c + r - tau) / (2**s - 1)
                    for s in range(1, tau + 13)
                ]
                maximum = (5 - c + r) / (2**(tau + 5) - 1)
                assert max(strengths) == maximum
                assert strengths.index(maximum) == tau + 4
                thresholds += 1
            if e - tau >= 5:
                special = (e - tau - 3 + r) / (d - 1)
                assert special > (e - tau - 1 + r) / (2 * d - 1)
                assert (2 * r < special < 4 * r) == (e - tau in (5, 6))
                if e - tau in (5, 6):
                    assert special < 4 * r
                    assert (1 + r) / (2**(tau + 5) - 1) < 4 * r
                thresholds += 1
        for h in (Fraction(9, 4) * r, 3 * r, Fraction(15, 4) * r):
            target = beta(d - 1, h)
            for m in range(1, 8 * d):
                if m < d // 2:
                    gain = Fraction(3, 2) if h < 3 * r else 2
                    assert beta(m, h) + gain > target + 1
                elif m < d and m != d - 1:
                    assert beta(m, h) > target
                elif m >= d:
                    assert beta(m, h) > target + 1
                frontiers += 1
    return thresholds, frontiers


def main():
    crossing, gaps, radial = check_axis_pairs()
    counts = check_small_counts()
    thresholds, frontiers = check_strengths()
    print(
        f"{crossing:,} coefficient crossings, {gaps:,} axis gaps, "
        f"{radial} radial ratios, {counts:,} paired counts, "
        f"{thresholds} thresholds and {frontiers:,} frontier checks"
    )


if __name__ == "__main__":
    main()

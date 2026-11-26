"""Exact replay of the two escaping-root translations at the negative specialization."""

from fractions import Fraction
from math import comb, factorial, gcd, prod


def moment_coefficients(rows, order):
    return [
        comb(order, k) * prod(factorial((order - k) * a + k * b) for a, b in rows)
        for k in range(order + 1)
    ]


def rising(j, parameter, degree, repetitions, remove_zero=False):
    factors = [j * parameter + index for index in range(1, degree + 1)]
    if remove_zero:
        assert factors.count(0) == 1
        factors.remove(0)
        factors.append(j)
    else:
        assert 0 not in factors
    return prod(factors) ** repetitions


def extracted_translations(n, upper_order, h, r, s, lower, upper):
    first = Fraction(-1, upper_order)
    second = Fraction(-h, n)
    width = n - h + 1

    def factor(j, k, remove_first=False, remove_second=False):
        return (
            rising(j, first, k, r, remove_first)
            * rising(j, second, k, s, remove_second)
        )

    left = lower[h - 1] * factor(n, h - 1)
    left_neighbor = lower[h - 2] * factor(n, h - 2) if h >= 2 else Fraction(0)
    right = lower[n] * factor(n, n, remove_second=True)
    right_neighbor = lower[n - 1] * factor(n, n - 1, remove_second=True)
    top = upper[upper_order] * factor(upper_order, upper_order, remove_first=True)
    top_neighbor = upper[upper_order - 1] * factor(upper_order, upper_order - 1, remove_first=True)
    return (left_neighbor / left - right_neighbor / right) / width, -top_neighbor / (upper_order * top)


def closed_translations(n, upper_order, h, r, s, lower, upper):
    width = n - h + 1
    low_ratio = Fraction(0)
    if h >= 2:
        low_ratio = (-1) ** s * Fraction(lower[h - 2], lower[h - 1]) / (
            h - 1 - Fraction(n, upper_order)
        ) ** r
    high_ratio = Fraction(lower[n - 1], lower[n]) / (
        Fraction(n * (upper_order - 1), upper_order) ** r * (width - 1) ** s
    )
    high_upper = -Fraction(upper[upper_order - 1], upper_order * upper[upper_order]) / (
        (upper_order - 1)**r * Fraction(upper_order * (width - 1), n) ** s
    )
    return (low_ratio - high_ratio) / width, high_upper


if __name__ == "__main__":
    count = 0
    for rows in ((), ((2, 3), (0, 1)), ((3, 2), (4, 6)), ((8, 1), (4, 2))):
        degree_pairs = [(n, n + 1) for n in range(2, 14)]
        degree_pairs.extend([(3, 8), (4, 9), (5, 12), (7, 10), (8, 13), (4, 3), (10, 3), (12, 5), (13, 6)])
        for n, upper_order in degree_pairs:
            lower = moment_coefficients(rows, n)
            upper = moment_coefficients(rows, upper_order)
            for h in range(1, n):
                for r, s in ((1, 1), (2, 3), (4, 2), (6, 5)):
                    extracted = extracted_translations(n, upper_order, h, r, s, lower, upper)
                    assert extracted == closed_translations(n, upper_order, h, r, s, lower, upper)
                    if gcd(n - h + 1, upper_order) > 1:
                        if s % 2 and n < upper_order and all(a <= b for a, b in rows):
                            assert extracted[0] < extracted[1] < 0
                        if not r % 2 and not s % 2 and upper_order % 2:
                            assert extracted[0] > 0 > extracted[1]
                    count += 1
    print(f"Verified both escape translations in {count:,} exact coefficient cases.")

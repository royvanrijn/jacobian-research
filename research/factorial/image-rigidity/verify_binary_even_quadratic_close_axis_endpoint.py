#!/usr/bin/env python3
from fractions import Fraction
from itertools import product
from math import comb, factorial


def valuation(integer):
    if not integer:
        return float("inf")
    return (integer & -integer).bit_length() - 1


def multiply(left, right):
    coefficients = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            coefficients[i + j] += first * second
    return coefficients


def radial(a, b, c, order):
    coefficients = [1]
    for _ in range(order):
        coefficients = multiply(coefficients, [a, b - 2 * a, a - b + c])
    return factorial(2 * order + 1) * sum(
        Fraction(value, i + 1) for i, value in enumerate(coefficients)
    )


def native(a, b, c, order):
    result = 0
    for i in range(order + 1):
        for ell in range(order - i + 1):
            j = order - i - ell
            result += (
                comb(order, i) * comb(order - i, ell)
                * factorial(2 * i + ell) * factorial(2 * j + ell)
                * a**i * b**ell * c**j
            )
    return result


def check_radial():
    count = 0
    for a, b, c in product((-2, 0, 1), repeat=3):
        for order in range(9):
            assert radial(a, b, c, order) == native(a, b, c, order)
            count += 1
    return count


def eta(order, degree, d):
    return Fraction(degree, d) + (order - degree).bit_count() - order.bit_count()


def check_blocks():
    count = 0
    for exponent in range(2, 7):
        d = 2**exponent
        for high in range(5):
            j = 2 * d * (high + 1) - 1
            selected = eta(j, 2 * d - 1, d)
            for degree in range(j + 1):
                q, w = divmod(degree, 2 * d)
                cost = eta(j, degree, d) - selected
                assert cost >= 2 * q - q.bit_count()
                assert (cost == 0) == (degree in (d - 1, 2 * d - 1))
                if w < d:
                    assert eta(j, degree + d, d) - eta(j, degree, d) == 0
                if q:
                    assert (2 * degree + 1).bit_length() - 1 - (exponent + 1) <= (
                        2 * q - q.bit_count()
                    )
                count += 1
    return count


def check_frontier():
    count = 0
    for exponent in range(2, 6):
        d = 2**exponent
        for high in range(3):
            j = 2 * d * (high + 1) - 1
            scale = 4 * d
            for strength in (2, 4, 5, 7, 8):
                h = Fraction(strength, scale)
                paired = h >= Fraction(1, d)
                gain = min(Fraction(1), (d - 1) * h) if paired else 0
                valuations = [int(scale * eta(j, k, d)) for k in range(j + 1)]
                powers = [
                    m * strength - scale * valuation(m + 1)
                    + (int(scale * gain) if paired and m < d else 0)
                    for m in range(2 * j + 1)
                ]
                target = valuations[2 * d - 1] + powers[2 * d - 1]
                least = set()
                for degree, ordinate in enumerate(valuations):
                    for power in range(2 * degree + 1):
                        value = ordinate + powers[power]
                        assert value >= target
                        if value == target:
                            least.add((degree, power))
                        count += 1
                expected = {(2 * d - 1, 2 * d - 1)}
                if strength == 8:
                    expected.update(
                        (k, m)
                        for k in (d - 1, 2 * d - 1)
                        for m in (d // 2 - 1, d - 1)
                    )
                assert least == expected
    return count


def check_count_parity():
    count = 0
    for exponent in range(2, 6):
        d = 2**exponent
        for k in range(2 * d):
            for m in range(d):
                for ell in range(m // 2 + 1):
                    linear = m - 2 * ell
                    def coefficient(n):
                        return comb(n, ell) * comb(n - ell, linear) if n >= ell else 0
                    assert (coefficient(k) - coefficient(k + d)) % 2 == 0
                    count += 1
        for m, expected in ((d // 2 - 1, 1), (d - 1, 3)):
            numerator = comb(2 * d - 1, m)
            denominator = comb(d - 1, m)
            assert (numerator - expected * denominator) % 4 == 0
            count += 1
    return count


def main():
    radial = check_radial()
    blocks = check_blocks()
    frontier = check_frontier()
    counts = check_count_parity()
    print(
        f"{radial} radial identities, {blocks:,} block checks, "
        f"{frontier:,} frontier terms, {counts:,} count and ratio checks"
    )


if __name__ == "__main__":
    main()

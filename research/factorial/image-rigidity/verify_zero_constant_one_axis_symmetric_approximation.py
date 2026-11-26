#!/usr/bin/env python3
from fractions import Fraction
from itertools import product

from verify_zero_constant_one_axis_unit_pole_removal import axis, exact_moment


def valuation(value):
    if not value:
        return float("inf")
    value = Fraction(value)
    n, d = abs(value.numerator), value.denominator
    return (n & -n).bit_length() - (d & -d).bit_length()


def check_base():
    count = 0
    for gamma, e in product((-2, 0, 1, 2), repeat=2):
        moments = [exact_moment(j, (1, 1, gamma, e, e)) for j in range(10)]
        boundary = [axis(j, gamma, e) for j in range(10)]
        assert moments[1] == 2 + gamma + 2 * e
        count += 1
        for j in range(2, 10):
            assert moments[j] == (
                2 * j * moments[j - 1] + e * e * moments[j - 2]
                + boundary[j] + e * boundary[j - 1]
            )
            count += 1
        if not e:
            continue
        normalized = [
            Fraction(moment, e**j if j % 2 == 0 else 2 * e ** (j - 1))
            for j, moment in enumerate(moments)
        ]
        for j in range(2, 10):
            force = boundary[j] + e * boundary[j - 1]
            if j % 2:
                expected = (
                    j * normalized[j - 1] + normalized[j - 2]
                    + Fraction(force, 2 * e ** (j - 1))
                )
            else:
                expected = (
                    Fraction(4 * j, e * e) * normalized[j - 1] + normalized[j - 2]
                    + Fraction(force, e**j)
                )
            assert normalized[j] == expected
            count += 1
    return count


def multiply_polynomial(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return result


def check_filtered_numerator():
    count = 0
    for sigma in range(1, 5):
        for w in range(sigma + 1, sigma + 7):
            r = 2**sigma
            s = r + 2**w
            a, a0 = Fraction(1, 1 + r * r), Fraction(1, 1 + s * s)
            u, u0 = 1 - a, 1 - a0
            d, e = r * a, s * a0
            assert valuation(d - e) == w
            for sign in (-1, 1):
                beta = 1 + sign * 2**w
                numerator = [1, beta + d + e, d * e]
                base = [1, 1 + 2 * e, e * e]
                first = multiply_polynomial(
                    multiply_polynomial([1, a0], numerator), [1, u0]
                )
                second = multiply_polynomial(
                    multiply_polynomial([1, a], base), [1, u]
                )
                assert first[0] == second[0]
                for x, y, bound in zip(
                    first[1:], second[1:], (w, w, w + sigma, w + 3 * sigma)
                ):
                    assert valuation(x - y) >= bound
                    count += 1
    return count


def field_product(a, b, modulus):
    result = 0
    top_bit = 1 << (modulus.bit_length() - 1)
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & top_bit:
            a ^= modulus
        b >>= 1
    return result


def field_power(a, exponent, modulus):
    size = 1 << (modulus.bit_length() - 1)
    exponent %= size - 1
    value = 1
    for _ in range(exponent):
        value = field_product(value, a, modulus)
    return value


def check_geometric_sums():
    count = 0
    for modulus in (0b111, 0b1011, 0b10011, 0b100101):
        size = 1 << (modulus.bit_length() - 1)
        for xi in range(2, size):
            g = xi ^ field_power(xi, -1, modulus)
            assert g
            even, odd_c_sum, even_p_sum = 1, g, 1
            for m in range(61):
                if m:
                    even ^= field_power(xi, 2 * m, modulus)
                    even ^= field_power(xi, -2 * m, modulus)
                    odd_c_sum ^= field_power(xi, 2 * m + 1, modulus)
                    odd_c_sum ^= field_power(xi, -2 * m - 1, modulus)
                    even_p_sum ^= even
                assert (even == 0) == (field_power(xi, 2 * m + 1, modulus) == 1)
                count += 1
                if field_power(xi, 2 * m + 3, modulus) == 1:
                    assert odd_c_sum == 1
                    assert field_product(g, even_p_sum, modulus) == 1
                    count += 2
    return count


def main():
    base = check_base()
    filtered = check_filtered_numerator()
    geometric = check_geometric_sums()
    print(
        f"{base} base and normalized recurrences, "
        f"{filtered} filtered coefficient bounds, {geometric:,} residue sums"
    )


if __name__ == "__main__":
    main()

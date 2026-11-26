#!/usr/bin/env python3
from fractions import Fraction
from itertools import combinations


def check_axis_minima():
    count = 0
    for e in range(2, 7):
        d = 2**e
        for order in range(4 * d + 5):
            levels = [
                Fraction(i, d) + (order - i).bit_count() - order.bit_count()
                for i in range(order + 1)
            ]
            high, low = divmod(order, d)
            least = {i for i, value in enumerate(levels) if value == min(levels)}
            assert least == ({low, low + d} if high % 2 else {low})
            count += 1
    return count


def check_even_wraps():
    axis_count = 0
    cross_count = 0
    for e in range(2, 7):
        d = 2**e
        for high in range(1, 10):
            order = 2 * d * high
            levels = [
                Fraction(i, d) + (order - i).bit_count() - order.bit_count()
                for i in range(order + 1)
            ]
            low_terms = {i for i in range(1, order + 1) if levels[i] <= 1}
            assert low_terms == ({d, 2 * d} if high % 2 else set())
            assert all(levels[i] == 1 for i in low_terms)
            axis_count += 1
            j = order - 1
            cross = [
                Fraction(k, d) + (j - k).bit_count()
                - j.bit_count() - k.bit_count()
                for k in range(j + 1)
            ]
            low_terms = {k for k, value in enumerate(cross) if value == min(cross)}
            assert low_terms == ({2 * d - 1} if high % 2 else {2 * d - 1, 4 * d - 1})
            cross_count += 1
    return axis_count, cross_count


def multiply(a, b, modulus):
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


def power(value, exponent, modulus):
    size = 1 << (modulus.bit_length() - 1)
    exponent %= size - 1
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value, modulus)
        value = multiply(value, value, modulus)
        exponent >>= 1
    return result


def divide(a, b, modulus):
    assert b
    return multiply(a, power(b, -1, modulus), modulus)


def check_residues():
    odd_wrap = 0
    even_wrap = 0
    denominator = 0
    for modulus in (0b111, 0b1011, 0b10011, 0b100101, 0b1000011):
        size = 1 << (modulus.bit_length() - 1)
        assert all(multiply(x, power(x, -1, modulus), modulus) == 1
                   for x in range(1, size))
        for e in range(2, 7):
            d = 2**e
            axis = {x: power(x, d, modulus) for x in range(1, size)}
            odd = {
                x: multiply(power(x, 1 - d, modulus), 1 ^ axis[x], modulus)
                for x in range(2, size)
            }
            even = {
                x: divide(power(x, 1 - d, modulus), 1 ^ axis[x], modulus)
                for x in range(2, size)
            }
            for a, c in combinations(range(2, size), 2):
                total = a ^ c
                if odd[a] == odd[c]:
                    mixed = divide(
                        multiply(a, 1 ^ axis[a], modulus)
                        ^ multiply(c, 1 ^ axis[c], modulus),
                        axis[a] ^ axis[c],
                        modulus,
                    )
                    assert mixed == odd[a]
                    outer = divide(power(a, d - 1, modulus)
                                   ^ power(c, d - 1, modulus), total, modulus)
                    inner = divide(power(a, d - 2, modulus)
                                   ^ power(c, d - 2, modulus), total, modulus)
                    assert outer ^ multiply(mixed, inner, modulus) == power(
                        multiply(a, c, modulus), -1, modulus
                    )
                    odd_wrap += 1
                if even[a] == even[c]:
                    assert total != 1
                    even_wrap += 1
            for total in range(1, size):
                residue = 1 ^ axis[total] ^ power(total, 2 * d, modulus)
                assert (residue == 0) == (total != 1 and power(total, 3, modulus) == 1)
                denominator += 1
    return odd_wrap, even_wrap, denominator


def main():
    axis = check_axis_minima()
    even_axis, cross = check_even_wraps()
    odd_wrap, even_wrap, denominator = check_residues()
    print(
        f"{axis} axis minima, {even_axis} even-wrap axis checks, "
        f"{cross} pure-cross minima, {odd_wrap} odd-wrap residue pairs, "
        f"{even_wrap} even-wrap residue pairs, {denominator} denominator checks"
    )


if __name__ == "__main__":
    main()

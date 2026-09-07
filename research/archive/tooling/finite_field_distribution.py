#!/usr/bin/env python3
"""Enumerate the finite-field fibers and check the proved histograms.

Prime fields are checked directly. Binary extensions use a small polynomial-
basis implementation, keeping this certificate dependency-free.
"""

from collections import Counter
from itertools import product


def image_point(x, y, z, add, neg, mul, one):
    """Evaluate the integral formula using field operations."""
    def n(a, integer):
        result = 0
        for _ in range(abs(integer)):
            result = add(result, a)
        return neg(result) if integer < 0 else result

    def power(a, exponent):
        result = one
        while exponent:
            if exponent & 1:
                result = mul(result, a)
            a = mul(a, a)
            exponent //= 2
        return result

    u = add(one, mul(x, y))
    xy = mul(x, y)
    four_plus_3xy = add(n(one, 4), n(xy, 3))
    a = add(mul(power(u, 3), z), mul(mul(power(y, 2), u), four_plus_3xy))
    b = add(
        add(y, n(mul(mul(x, power(u, 2)), z), 3)),
        n(mul(mul(x, power(y, 2)), four_plus_3xy), 3),
    )
    c = add(add(n(x, 2), n(mul(power(x, 2), y), -3)), n(mul(power(x, 3), z), -1))
    return a, b, c


def histogram(elements, add, neg, mul, one):
    fibers = Counter(
        image_point(x, y, z, add, neg, mul, one)
        for x, y, z in product(elements, repeat=3)
    )
    result = Counter(fibers.values())
    result[0] = len(elements) ** 3 - len(fibers)
    return dict(sorted(result.items()))


def expected_odd(q, characteristic):
    if characteristic == 3:
        n3 = q * q * (q - 1) // 6
        return {0: 2 * n3, 1: q * q * (q + 1) // 2, 3: n3}
    n3 = (q - 1) * (q * q + 2) // 6
    return {0: 2 * n3, 1: (q**3 + q**2 - 2 * q + 2) // 2, 3: n3}


def expected_binary(q):
    counts = Counter({
        0: 2 * (q - 1) ** 2,
        1: q**3 - 2 * q**2 + 2 * q - 1,
        q - 1: 2 * (q - 1),
        2 * q - 1: 1,
    })
    # Counter construction overwrites rather than adds coincident keys.
    if q == 2:
        counts[1] = q**3 - 2 * q**2 + 2 * q - 1 + 2 * (q - 1)
    return dict(sorted(counts.items()))


def check_prime(p):
    elements = range(p)
    add = lambda a, b: (a + b) % p
    neg = lambda a: (-a) % p
    mul = lambda a, b: (a * b) % p
    actual = histogram(elements, add, neg, mul, 1)
    expected = expected_binary(p) if p == 2 else expected_odd(p, p)
    assert actual == expected, (p, actual, expected)
    print(f"PASS F_{p}: {actual}")


# Irreducible binary polynomials, encoded with the leading bit included.
BINARY_MODULI = {
    2: 0b111,       # X^2 + X + 1
    3: 0b1011,      # X^3 + X + 1
    4: 0b10011,     # X^4 + X + 1
    5: 0b100101,    # X^5 + X^2 + 1
    6: 0b1000011,   # X^6 + X + 1
}


def check_binary_extension(degree):
    modulus = BINARY_MODULI[degree]
    q = 1 << degree

    def mul(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a >> degree:
                a ^= modulus
        return result

    actual = histogram(range(q), lambda a, b: a ^ b, lambda a: a, mul, 1)
    expected = expected_binary(q)
    assert actual == expected, (q, actual, expected)
    print(f"PASS F_{q}: {actual}")


TERNARY_MODULI = {
    2: (1, 0, 1),       # X^2 + 1
    3: (1, 2, 0, 1),    # X^3 + 2X + 1
}


def check_ternary_extension(degree):
    """Check GF(3^degree) in a polynomial basis."""
    p = 3
    q = p**degree
    modulus = TERNARY_MODULI[degree]

    def digits(value):
        result = []
        for _ in range(degree):
            result.append(value % p)
            value //= p
        return result

    def encode(coefficients):
        return sum(coefficient * p**i for i, coefficient in enumerate(coefficients))

    def add(a, b):
        return encode([(x + y) % p for x, y in zip(digits(a), digits(b))])

    def neg(a):
        return encode([(-x) % p for x in digits(a)])

    def mul(a, b):
        left, right = digits(a), digits(b)
        coefficients = [0] * (2 * degree - 1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                coefficients[i + j] = (coefficients[i + j] + x * y) % p
        for power in range(2 * degree - 2, degree - 1, -1):
            leading = coefficients[power]
            for j in range(degree):
                coefficients[power - degree + j] = (
                    coefficients[power - degree + j] - leading * modulus[j]
                ) % p
        return encode(coefficients[:degree])

    actual = histogram(range(q), add, neg, mul, 1)
    expected = expected_odd(q, 3)
    assert actual == expected, (q, actual, expected)
    print(f"PASS F_{q}: {actual}")


if __name__ == "__main__":
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        check_prime(prime)
    for extension_degree in range(2, 7):
        check_binary_extension(extension_degree)
    for extension_degree in (2, 3):
        check_ternary_extension(extension_degree)

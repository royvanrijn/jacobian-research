#!/usr/bin/env python3
"""Finite exact checks of the common-marker integer family."""

from fractions import Fraction
from math import factorial


# A polynomial is a dict from exponent tuples to integer or rational coefficients.
def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponents, coefficient in polynomial.items():
            result[exponents] = result.get(exponents, 0) + coefficient
    return {exponents: c for exponents, c in result.items() if c}


def scale(factor, polynomial):
    return {exponents: factor * c for exponents, c in polynomial.items() if factor * c}


def multiply(left, right):
    result = {}
    for m, a in left.items():
        for n, b in right.items():
            exponents = tuple(x + y for x, y in zip(m, n))
            result[exponents] = result.get(exponents, 0) + a * b
    return {exponents: c for exponents, c in result.items() if c}


def power(polynomial, exponent):
    result = {(0,) * len(next(iter(polynomial))): 1}
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def evaluate(polynomial, point):
    total = 0
    for exponents, coefficient in polynomial.items():
        for exponent, coordinate in zip(exponents, point):
            coefficient *= coordinate**exponent
        total += coefficient
    return total


def original_contraction(polynomial):
    # Variable order: xi_1, xi_2, z_1, z_2.
    assert all(a + b == c + d for a, b, c, d in polynomial)
    return sum(
        value * factorial(c) * factorial(d)
        for (a, b, c, d), value in polynomial.items()
        if a == c and b == d
    )


def simple_root_specialization(polynomial):
    # Set xi_1=z_1=1, z_2=0, xi_2=u.
    result = {}
    for (_, exponent, _, z_2), coefficient in polynomial.items():
        if z_2 == 0:
            result[(exponent,)] = result.get((exponent,), 0) + coefficient
    return {exponents: c for exponents, c in result.items() if c}


def marked_interval_value(a, m):
    return Fraction(
        2 ** ((a + 2) * m) * factorial((a + 1) * m) * factorial(m),
        factorial((a + 2) * m + 1),
    )


def verify_original():
    radial = {(1, 0, 1, 0): 1, (0, 1, 0, 1): 1}
    x = {(1, 0, 0, 1): 1}
    y = {(0, 1, 1, 0): 2}
    t = {(1, 0, 1, 0): 1, (0, 1, 0, 1): -1}
    linear = add(radial, x)
    quadratic = add(scale(2, multiply(y, radial)), scale(-1, power(t, 2)))
    cubic = add(
        scale(2, multiply(y, power(radial, 2))),
        scale(-1, multiply(add(scale(2, radial), x), power(t, 2))),
    )
    f4 = multiply(linear, cubic)
    w = multiply(t, linear)
    g = multiply(linear, quadratic)
    shifted_w = add(power(radial, 2), w)
    f7 = multiply(f4, g)
    assert multiply(radial, f7) == add(power(f4, 2), multiply(f4, power(w, 2)))

    sample_count = 0
    for a, b, max_m in ((0, 1, 5), (0, 0, 4), (1, 0, 4), (2, 0, 3), (1, 1, 3), (0, 2, 3)):
        degree = 4 + 2 * a + 3 * b
        phase = multiply(multiply(f4, power(shifted_w, a)), power(g, b))
        assert {(i + j, k + l) for i, j, k, l in phase} == {(degree, degree)}
        assert evaluate(phase, (1, 1, 1, -1)) == -4 * (-2)**a * 4**b
        expected = scale(2**a, multiply({(0,): -2, (1,): 4}, power({(0,): -1, (1,): 4}, b)))
        actual = simple_root_specialization(phase)
        assert actual == expected
        assert evaluate(actual, (Fraction(1, 2),)) == 0
        derivative = {(e - 1,): e * c for (e,), c in actual.items() if e}
        assert evaluate(derivative, (Fraction(1, 2),)) != 0

        phase_power = {(0, 0, 0, 0): 1}
        for m in range(1, max_m + 1):
            phase_power = multiply(phase_power, phase)
            assert original_contraction(phase_power) == 0
            assert original_contraction(multiply(x, phase_power)) == (
                factorial(degree * m + 2) * marked_interval_value(a, m)
            ) != 0
            sample_count += 1

    corrupted = add(f7, power(radial, 7))
    assert original_contraction(corrupted) == factorial(8) != 0
    return sample_count


def verify_laurent():
    # Literal multiplication on the sphere in the two Laurent exponents (angle, height).
    x = {(1, 0): 1}
    t = {(0, 1): 1}
    y = {(-1, 0): Fraction(1, 2), (-1, 2): Fraction(-1, 2)}
    unit_plus_x = add({(0, 0): 1}, x)
    e = multiply(
        unit_plus_x,
        add(scale(2, y), scale(-1, multiply(add({(0, 0): 2}, x), power(t, 2)))),
    )
    g = multiply(unit_plus_x, add(scale(2, y), scale(-1, power(t, 2))))
    phase = multiply(e, g)

    def interval_constant_term(polynomial):
        return sum(
            Fraction(c, height + 1)
            for (angle, height), c in polynomial.items()
            if angle == 0 and height % 2 == 0
        )

    phase_power = {(0, 0): 1}
    for m in range(1, 9):
        phase_power = multiply(phase_power, phase)
        assert interval_constant_term(phase_power) == 0
        assert interval_constant_term(multiply(x, phase_power)) == marked_interval_value(0, m)


def main():
    if not __debug__:
        raise SystemExit("Run the replay without -O so its assertions are enabled.")
    count = verify_original()
    verify_laurent()
    print(f"Verified {count} pure and {count} marked original contractions.")
    print("Verified eight Laurent power/marker pairs, radial evaluations, and simple roots.")
    print("Negative control passed: adding R^7 to the septic yields first contraction 8!.")


if __name__ == "__main__":
    main()

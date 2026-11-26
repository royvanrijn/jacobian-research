#!/usr/bin/env python3
from fractions import Fraction
from itertools import product
from math import comb, factorial


def moment_terms(order):
    for x_square in range(order + 1):
        for cross in range(order - x_square + 1):
            for y_square in range(order - x_square - cross + 1):
                remaining = order - x_square - cross - y_square
                for x_linear in range(remaining + 1):
                    y_linear = remaining - x_linear
                    x_degree = 2 * x_square + cross + x_linear
                    y_degree = 2 * y_square + cross + y_linear
                    denominator = (
                        2 ** (x_square + y_square)
                        * factorial(x_square) * factorial(cross)
                        * factorial(y_square) * factorial(x_linear) * factorial(y_linear)
                    )
                    multiplier, remainder = divmod(
                        factorial(x_degree) * factorial(y_degree), denominator
                    )
                    assert remainder == 0
                    yield multiplier, (x_square, cross, y_square, x_linear, y_linear)


def exact_moment(j, values):
    answer = 0
    for coefficient, counts in moment_terms(j):
        term = coefficient
        for count, value in zip(counts, values):
            term *= value**count
        answer += term
    return answer


def axis(j, alpha, linear):
    return sum(
        comb(j, square) * factorial(j + square)
        // (2**square * factorial(j)) * alpha**square * linear ** (j - square)
        for square in range(j + 1)
    )


def check_interior():
    count = 0
    for beta, gamma, d, e in product((-1, 0, 2), repeat=4):
        values = (1, beta, gamma, d, e)
        moments = [exact_moment(j, values) for j in range(8)]
        x_axis = [axis(j, 1, d) for j in range(8)]
        y_axis = [axis(j, gamma, e) for j in range(8)]
        delta = gamma - beta * beta
        boundary = 1 + gamma - 2 * beta
        linear = gamma * d * d - 2 * beta * d * e + e * e
        for j in range(2, 8):
            assert x_axis[j] == (2 * j - 1) * x_axis[j - 1] + d * d * x_axis[j - 2]
            assert y_axis[j] == gamma * (2 * j - 1) * y_axis[j - 1] + e * e * y_axis[j - 2]
            assert boundary * moments[j] == (
                2 * delta * j * moments[j - 1] + linear * moments[j - 2]
                + (gamma - beta) * y_axis[j]
                + (gamma * d - beta * e) * y_axis[j - 1]
                + (1 - beta) * x_axis[j] + (e - beta * d) * x_axis[j - 1]
            )
            count += 3
    return count


def check_slopes():
    count = 0
    for d in (2, 4, -4, 8, -8, 12):
        a = [axis(j, 1, d) for j in range(19)]
        ratios = [Fraction(a[j], a[j - 1]) for j in range(1, len(a))]
        slopes = [second - first for first, second in zip(ratios, ratios[1:])]
        assert slopes[0] == 2 - Fraction(d, 1 + d)
        count += 1
        for j in range(len(slopes) - 1):
            assert slopes[j + 1] == 2 - d * d * slopes[j] / ratios[j] / ratios[j + 1]
            count += 1
    return count


def multiply(a, b, precision):
    answer = 0
    while b:
        if b & 1:
            answer ^= a
        a <<= 1
        b >>= 1
    return answer & ((1 << precision) - 1)


def power(a, n, precision):
    answer = 1
    for _ in range(n):
        answer = multiply(answer, a, precision)
    return answer


def matching(variables, maximum, precision):
    alpha, beta, gamma, d, e = variables
    mult = lambda first, second: multiply(first, second, precision)
    d2, e2 = mult(d, d), mult(e, e)
    denominator = (
        1, alpha ^ gamma, d2 ^ mult(alpha, gamma) ^ e2,
        mult(alpha, e2) ^ mult(d2, gamma), mult(d2, e2),
    )
    numerator = (1, beta ^ d ^ e, mult(d, e))
    coefficients = []
    for j in range(maximum + 1):
        coeff = numerator[j] if j < len(numerator) else 0
        for i in range(1, min(j, 4) + 1):
            coeff ^= mult(denominator[i], coefficients[j - i])
        coefficients.append(coeff)
    return coefficients


def check_ramified_matching():
    count = 0
    maximum = 7
    terms = [
        [counts for coefficient, counts in moment_terms(j) if coefficient % 2]
        for j in range(maximum + 1)
    ]
    for precision in range(2, 7):
        for beta, gamma, d, e in product((0, 2, 4, 6), repeat=4):
            if max(beta, gamma, d, e) >= 2**precision:
                continue
            variables = (1, beta | 1, gamma, d, e)
            expected = matching(variables, maximum, precision)
            for j in range(maximum + 1):
                actual = 0
                for counts in terms[j]:
                    term = 1
                    for c, p in zip(variables, counts):
                        term = multiply(term, power(c, p, precision), precision)
                    actual ^= term
                assert actual == expected[j]
                count += 1
    return count


def main():
    interior = check_interior()
    slopes = check_slopes()
    matching_count = check_ramified_matching()
    print(
        f"{interior:,} interior and axis checks, {slopes} slope checks, "
        f"{matching_count:,} ramified matching coefficients"
    )


if __name__ == "__main__":
    main()

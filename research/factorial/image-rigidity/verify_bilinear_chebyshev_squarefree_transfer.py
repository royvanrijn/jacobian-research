#!/usr/bin/env python3
from fractions import Fraction
from math import comb, factorial


def next_chebyshev(current, previous):
    size = max(len(current) + 1, len(previous))
    answer = [0] * size
    for i, value in enumerate(current):
        answer[i + 1] += 2 * value
    for i, value in enumerate(previous):
        answer[i] -= value
    return answer


def evaluate(polynomial, x):
    value = 0
    for coefficient in reversed(polynomial):
        value = coefficient + x * value
    return value


def derivative(polynomial):
    return [i * c for i, c in enumerate(polynomial) if i]


def transform(polynomial, x):
    coefficients = []
    while polynomial:
        coefficients.append(evaluate(polynomial, x))
        polynomial = derivative(polynomial)
    return coefficients


def check_identity():
    checks = 0
    u = [[1], [0, 2]]
    t = [[1], [0, 1]]
    for _ in range(30):
        u.append(next_chebyshev(u[-1], u[-2]))
        t.append(next_chebyshev(t[-1], t[-2]))
    for n in range(1, 31):
        for x in (0, 1, -1, Fraction(1, 3), Fraction(-2, 5)):
            next_poly = transform(u[n], x)
            previous_poly = transform(u[n - 1], x)
            rhs = [0] * (n + 2)
            for k, coefficient in enumerate(next_poly):
                rhs[k] += x * coefficient
                rhs[k + 1] += (k - n) * coefficient
            rhs[0] -= evaluate(t[n + 1], x)
            assert rhs == previous_poly + [0, 0]
            checks += 1
    return checks


def check_opposite():
    checks = 0
    for m in range(1, 41):
        n = 2 * m
        for r in range(m + 1):
            j = n - 2 * r
            cross_expansion = (-1)**r * factorial(j) * comb(j + r, r)
            transformed = (-1)**r * factorial(2 * m - r) // factorial(r)
            robin = (-1)**r * comb(m, m - r) * factorial(2 * m - r) * factorial(m - r) // factorial(m)
            assert cross_expansion == transformed == robin
            checks += 1
    return checks


def main():
    print(f"{check_identity()} consecutive Chebyshev identities, {check_opposite()} Robin coefficients")


if __name__ == "__main__":
    main()

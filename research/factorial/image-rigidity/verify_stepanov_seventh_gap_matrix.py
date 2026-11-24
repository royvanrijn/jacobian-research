#!/usr/bin/env python3
"""Independent direct-product/matrix verification of the seven gap-seven lifts."""

import json
from math import comb, prod


MODULUS = 7 ** 4
IDENTITY = ((1, 0), (0, 1))
THETA = ((0, -3), (1, 1))


def add(left, right):
    return tuple(tuple((left[i][j] + right[i][j]) % MODULUS for j in range(2))
                 for i in range(2))


def multiply(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(2)) % MODULUS
                       for j in range(2)) for i in range(2))


def scalar(number):
    return ((number % MODULUS, 0), (0, number % MODULUS))


def element(a, b):
    return ((a % MODULUS, -3 * b % MODULUS),
            (b % MODULUS, (a + b) % MODULUS))


def coefficients(order):
    gamma = 3 * order + 1
    values = []
    for k in range(28):
        # Independent of the generator's block-by-block division.
        numerator = prod(range(gamma, gamma + 7 * k))
        denominator = 7 ** k
        assert numerator % denominator == 0
        values.append(comb(order, k) * (numerator // denominator) % MODULUS)
    return values


def evaluate(poly, argument):
    result = scalar(0)
    for coefficient in reversed(poly):
        result = add(multiply(result, argument), scalar(coefficient))
    return result


def valuation(integer):
    assert integer
    exponent = 0
    while integer % 7 == 0:
        exponent += 1
        integer //= 7
    return exponent


def main():
    assert add(add(multiply(THETA, THETA), element(0, -1)), scalar(3)) == scalar(0)
    # Every nonzero quadratic residue class is tested explicitly for the
    # claimed irreducibility of x^2-x+3 over F_7.
    assert all((x * x - x + 3) % 7 for x in range(7))
    lifts = [(1517, 2053), (1860, 1710), (2203, 1367), (145, 1024),
             (488, 681), (831, 338), (1174, 2396)]
    residues = [(4, 6), (0, 4), (3, 2), (6, 0), (2, 5), (5, 3), (1, 1)]
    records = []
    for branch, (lift, residue) in enumerate(zip(lifts, residues)):
        terminal = 446 + 2401 * branch
        lower = coefficients(terminal - 1)
        upper = coefficients(terminal)
        assert lower == coefficients(terminal - 1 + 7 ** 5)
        assert upper == coefficients(terminal + 7 ** 5)
        argument = element(*lift)
        assert evaluate(lower, argument) == scalar(0)
        actual = evaluate(upper, argument)
        expected = element(343 * residue[0], 343 * residue[1])
        assert actual == expected and actual != scalar(0)
        # Check the other residue-field conjugate without regenerating a root.
        conjugate = element(lift[0] + lift[1], -lift[1])
        assert evaluate(lower, conjugate) == scalar(0)
        conjugate_upper = evaluate(upper, conjugate)
        assert conjugate_upper == element(343 * (residue[0] + residue[1]),
                                           -343 * residue[1])
        assert conjugate_upper != scalar(0)
        derivative = [k * lower[k] % MODULUS for k in range(1, len(lower))]
        derivative_matrix = evaluate(derivative, argument)
        determinant = ((derivative_matrix[0][0] * derivative_matrix[1][1]
                        - derivative_matrix[0][1] * derivative_matrix[1][0]) % 7)
        assert determinant != 0
        records.append({"N": terminal, "root_z": lift,
                        "upper_over_343_mod7": residue,
                        "derivative_norm_mod7": determinant,
                        "lower_coefficients_mod2401": lower,
                        "upper_coefficients_mod2401": upper})
    binomial_period_guards = [valuation(comb(7 ** 5, j)) for j in range(1, 28)]
    assert min(binomial_period_guards) >= 4
    print(json.dumps({"status": "PASS", "method": "direct rising products and 2x2 matrix ring; no generator imports",
                      "scope": "complete seven-class certificate once the written period and tail lemmas are supplied",
                      "modulus": MODULUS, "period": 7 ** 5,
                      "binomial_shift_valuations_j1_to27": binomial_period_guards,
                      "records": records}, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    main()

#!/usr/bin/env python3
from fractions import Fraction
from itertools import product

from verify_zero_constant_one_axis_symmetric_approximation import field_power, field_product
from verify_zero_constant_one_axis_unit_pole_removal import moment_terms


DEGREE = 8


def monomial(coefficient, exponent=0):
    quotient, remainder = divmod(exponent, DEGREE)
    values = [Fraction(0)] * DEGREE
    values[remainder] = coefficient * Fraction(2) ** quotient
    return tuple(values)


def add(*values):
    return tuple(sum(at_index) for at_index in zip(*values))


def multiply(first, second):
    result = monomial(0)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            if a and b:
                result = add(result, monomial(a * b, i + j))
    return result


def residue(value):
    for coefficient in value:
        assert coefficient.denominator % 2
    return int(value[0].numerator % 2)


def matrix_product(first, second):
    return [
        [add(*(multiply(first[i][l], second[l][j]) for l in range(2))) for j in range(2)]
        for i in range(2)
    ]


def matrix_vector(matrix, vector):
    return [add(*(multiply(row[i], vector[i]) for i in range(2))) for row in matrix]


def transition(m, delta):
    b = monomial(m, -delta)
    q = monomial(2 * m + 1)
    return [[monomial(1), b], [q, add(monomial(1), multiply(q, b))]]


def forcing(m, c, d):
    return [c, add(multiply(monomial(2 * m + 1), c), d)]


def check_matrices():
    count = 0
    for delta in range(4):
        z, inverse = monomial(1, delta), monomial(1, -delta)
        for m in range(2, 44, 2):
            even = transition(m, delta)
            odd = transition(m - 1, delta)
            combined = matrix_product(even, odd)
            combined[0][1] = multiply(z, combined[0][1])
            combined[1][0] = multiply(inverse, combined[1][0])
            assert [[residue(x) for x in row] for row in combined] == [[1, 1], [0, 1]]
            count += 4
            for prev, curr in product(range(2), repeat=2):
                first = forcing(m - 1, monomial(prev), monomial(1, 3))
                second = forcing(m, monomial(curr), monomial(1, 4))
                out = [add(a, b) for a, b in zip(second, matrix_vector(even, first))]
                out[0] = multiply(z, out[0])
                assert [residue(x) for x in out] == [prev ^ curr if delta == 0 else 0, curr]
                count += 2
            for left, right, curr in product(range(2), repeat=3):
                state = [monomial(left, -delta), monomial(right)]
                out = matrix_vector(transition(m + 1, delta), state)
                raw_force = forcing(m + 1, monomial(curr), monomial(1, 3))
                result = [residue(multiply(z, add(a, b))) for a, b in zip(out, raw_force)]
                inter = left ^ right ^ (curr if delta == 0 else 0)
                assert result == [inter, inter ^ (right if delta == 0 else 0)]
                count += 2
    return count


def check_residue_sums():
    count = 0
    for modulus in (0b111, 0b1011, 0b10011, 0b100101):
        size = 1 << (modulus.bit_length() - 1)
        for xi in range(2, size):
            g = xi ^ field_power(xi, -1, modulus)
            prefix = 0
            for r in range(121):
                v = 0
                for s in range(-r, r + 1):
                    v ^= field_power(xi, 4 * s, modulus)
                assert (not v) == (field_power(xi, 2 * r + 1, modulus) == 1)
                count += 1
                if not v:
                    assert field_product(prefix, field_product(g, g, modulus), modulus) == 1
                    forcing_sum = 0
                    for s in range(-2 * r, 2 * r + 1):
                        forcing_sum ^= field_power(xi, 2 * s, modulus)
                    assert forcing_sum == 1
                    count += 2
                prefix ^= v
    return count


def check_jets():
    count = 0
    for j in range(1, 12):
        beta, linear = 0, 0
        for coeff, counts in moment_terms(j):
            _, cross, y_square, x_linear, y_linear = counts
            if y_square == y_linear == 0:
                if x_linear == 0:
                    beta += coeff * cross
                elif x_linear == 1:
                    linear += coeff
        assert beta % 2 == linear % 2 == 1
        count += 2
    return count


def main():
    matrices = check_matrices()
    sums = check_residue_sums()
    jets = check_jets()
    print(f"{matrices:,} transfer checks, {sums:,} residue sums, {jets} first-order checks")


if __name__ == "__main__":
    main()

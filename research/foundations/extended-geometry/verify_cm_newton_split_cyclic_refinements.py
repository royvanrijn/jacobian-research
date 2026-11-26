#!/usr/bin/env python3
"""Check the fractional-part table and split-prime Hasse–Witt ranks."""

from fractions import Fraction
from math import comb


def fractional_part(number):
    return number - number.numerator // number.denominator


def refinement_bits(x):
    return tuple(-1 + fractional_part(Fraction(t, 3)-x) + fractional_part(3*x)
                 + fractional_part(-Fraction(t, 3)-2*x) for t in range(3))


def check_table():
    endpoints = tuple(map(Fraction, ("0", "1/6", "1/3", "1/2", "2/3", "5/6", "1")))
    expected = ((1, 0, 0), (1, 0, 1), (0, 1, 0),
                (1, 1, 0), (0, 0, 1), (0, 1, 1))
    breakpoints = {Fraction(0), Fraction(1)}
    for t in range(3):
        for constant, slope in ((Fraction(t, 3), -1), (Fraction(0), 3),
                                (-Fraction(t, 3), -2)):
            for integer in range(-3, 4):
                x = (integer-constant) / slope
                if 0 < x < 1:
                    breakpoints.add(x)
    assert tuple(sorted(breakpoints)) == endpoints
    # The three slopes in the fractional-part sum add to zero on each interval.
    assert tuple(refinement_bits((left+right)/2) for left, right in
                 zip(endpoints, endpoints[1:])) == expected


def original_matrix(n, prime, i):
    size = 1 + (6*i)//n - 2*((3*i)//n)
    z_exponent = prime * ((n-3)*i % n) // n
    polynomial_exponent = prime * (3*i % n) // n

    def coefficient(height):
        index, remainder = divmod(height-z_exponent, 3)
        if remainder or not 0 <= index <= polynomial_exponent:
            return 0
        return ((-1)**(polynomial_exponent-index) * comb(polynomial_exponent, index)) % prime

    return tuple(tuple(coefficient(prime*source-target) for source in range(1, size+1))
                 for target in range(1, size+1))


def stable_rank(matrix, prime):
    if len(matrix) == 1:
        return int(matrix[0][0] != 0)
    (a, b), (c, d) = matrix
    if (a*d-b*c) % prime:
        return 2
    return int(any(entry % prime for entry in (a*a+b*c, a*b+b*d, c*a+d*c, c*b+d*d)))


def check_original_cover(n, prime):
    assert prime % n == 1 and prime % 3 == 2
    matrices = {i: original_matrix(n, prime, i) for i in range(1, n, 2)}
    ranks = {}
    for i, matrix in matrices.items():
        fixed, first, second = refinement_bits(Fraction(i, n))
        assert len(matrix) == fixed+first+second
        ranks[i] = stable_rank(matrix, prime)
        assert ranks[i] == fixed+2*first*second
    assert sum(ranks.values()) == n//4 + 2*((n+6)//12)
    order = (prime-1)//n
    assert matrices[1] == ((comb(3*order, order) % prime,),)
    if n == 16:
        assert matrices[1] == ((3,),)
        assert matrices[5] == ((11, 0), (0, 0))
        assert matrices[11] == ((0,),)
    print(f"n={n}, p={prime}: selected scalar {matrices[1][0][0]}, primitive p-rank {sum(ranks.values())}")


if __name__ == "__main__":
    check_table()
    check_original_cover(16, 17)
    check_original_cover(64, 257)

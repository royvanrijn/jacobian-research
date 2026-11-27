#!/usr/bin/env python3
"""Exact finite controls for MCT21, not a replay of MCT16 or a BSD proof.

Run from the repository root with Python 3.12 (standard library only):
    python3 research/arithmetic/verify_mordell_prime_density.py

Exhaustively enumerate dimensions 2 and 4; compare with two general
counting formulas through r=12. Bound the infinite product with rational
arithmetic. The unrestricted recurrence is proved in the canonical note.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations, product
from math import prod


def rank_mod_three(matrix):
    rows = [[x % 3 for x in row] for row in matrix]
    rank = 0
    for column in range(len(rows)):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][column], -1, 3)
        rows[rank] = [(x * inv) % 3 for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank:
                scalar = rows[i][column]
                rows[i] = [(x - scalar * y) % 3 for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def enumerate_ranks(dimension):
    positions = tuple(combinations(range(dimension), 2))
    counts = {rank: 0 for rank in range(0, dimension + 1, 2)}
    for entries in product(range(3), repeat=len(positions)):
        matrix = [[0] * dimension for _ in range(dimension)]
        for (i, j), value in zip(positions, entries):
            matrix[i][j] = value
            matrix[j][i] = -value
        rank = rank_mod_three(matrix)
        assert rank in counts, (dimension, entries, rank)
        counts[rank] += 1
    assert sum(counts.values()) == 3 ** len(positions)
    return counts


def nonsingular_count(r):
    return 3 ** (r * (r - 1)) * prod(3 ** (2 * j - 1) - 1 for j in range(1, r + 1))


def orbit_count(r):
    """GL_(2r)(F3) / Sp_(2r)(F3), counted independently of the row recurrence."""
    n = 2 * r
    gl = prod(3**n - 3**i for i in range(n))
    sp = 3 ** (r * r) * prod(3 ** (2 * j) - 1 for j in range(1, r + 1))
    quotient, remainder = divmod(gl, sp)
    assert remainder == 0
    return quotient


def density(r):
    return Fraction(1, 12) * prod(
        1 - Fraction(1, 3 ** (2 * j - 1)) for j in range(1, r + 1)
    )


def limit_interval(truncation):
    upper = density(truncation)
    tail_sum = Fraction(1, 3 ** (2 * truncation + 1)) / (1 - Fraction(1, 9))
    return upper * (1 - tail_sum), upper


def decimal_percent(value):
    with localcontext() as context:
        context.prec = 35
        return Decimal(100 * value.numerator) / Decimal(value.denominator)


def main():
    assert rank_mod_three([[0, 0], [0, 0]]) == 0
    assert rank_mod_three([[0, 1], [-1, 0]]) == 2
    for r in (1, 2):
        counts = enumerate_ranks(2 * r)
        assert counts[2 * r] == nonsingular_count(r) == orbit_count(r)
        assert Fraction(counts[2 * r], 12 * sum(counts.values())) == density(r)
        print(f"dimension {2*r}: rank counts {counts}")
    previous = 1
    for r in range(1, 13):
        recurrence = (3 ** (2 * r - 1) - 1) * 3 ** (2 * r - 2) * previous
        assert recurrence == nonsingular_count(r) == orbit_count(r)
        assert Fraction(recurrence, 12 * 3 ** (r * (2 * r - 1))) == density(r)
        assert density(r) < density(r - 1)
        previous = recurrence
    assert [density(r) for r in (1, 2, 3)] == [
        Fraction(1, 18), Fraction(13, 243), Fraction(3146, 59049)
    ]
    for r in (1, 2, 3):
        print(f"delta_{r} = {density(r)} = {decimal_percent(density(r))}%")
    lower, upper = limit_interval(20)
    assert lower > Fraction(23, 432)
    assert upper - lower < Fraction(1, 10**20)
    for n in range(1, 20):
        old_lower, old_upper = limit_interval(n)
        new_lower, new_upper = limit_interval(n + 1)
        assert old_lower < new_lower < new_upper < old_upper
    print(f"limit percentage lies between {decimal_percent(lower)} and {decimal_percent(upper)}")
    print("PASS: finite counts, recurrence/orbit controls, and rational product bounds")
    print("Scope: no replay of governing fields, prime equidistribution, converse, or BSD proofs")


if __name__ == "__main__":
    main()

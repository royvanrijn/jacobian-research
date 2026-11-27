"""Enumerate the multiplicity-two outer cusp budgets of the (12,42) transfer."""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb


def profiles(max_delta=99):
    """Return the numerical solutions for (m, cusp denominators)."""
    admissible = {}
    for number_of_cusps in range(2, 5):
        for denominators in combinations_with_replacement(
            range(2, max_delta + 2), number_of_cusps
        ):
            if sum(denominator - 1 for denominator in denominators) > max_delta:
                continue
            rate = sum(
                (1 - Fraction(1, denominator) for denominator in denominators),
                Fraction(),
            )
            if rate <= Fraction(9, 5):
                continue
            candidate = rate / (5 * rate - 9)
            if candidate.denominator != 1 or candidate < 1:
                continue
            m = candidate.numerator
            ramification_index = 5 * m - 1
            if any(ramification_index % denominator for denominator in denominators):
                continue
            admissible.setdefault(m, []).append(denominators)
    return admissible


def main():
    admissible = profiles()
    expected = [
        1, 2, 3, 4, 5, 8, 10, 11, 17, 20, 23, 29, 31, 38,
        47, 65, 80, 101, 119, 173, 227, 335,
    ]
    assert sorted(admissible) == expected
    assert sum(map(len, admissible.values())) == 29
    assert admissible[1] == [(4, 4, 4), (2, 2, 2, 4)]
    assert admissible[2] == [(3, 3, 3)]
    assert admissible[335] == [(6, 31)]

    for m, options in admissible.items():
        e = 5 * m - 1
        for denominators in options:
            assert sum(e - e // denominator for denominator in denominators) == 9 * m
            assert sum(denominator - 1 for denominator in denominators) <= 99

    for degree in range(2, 26):
        boundary_constant = sum(
            Fraction((-1) ** j * comb(degree - 1, j), 2 * j + 1)
            for j in range(degree)
        )
        assert boundary_constant > 0

    print("Verified 29 multiplicity-two cusp budgets at 22 values of m.")


if __name__ == "__main__":
    main()

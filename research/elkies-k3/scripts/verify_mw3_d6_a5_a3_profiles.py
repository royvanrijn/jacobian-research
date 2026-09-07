#!/usr/bin/env python3
"""Exact Shioda-height replay for the D6+A5+A3 MW3 profile."""

from fractions import Fraction


def a_correction(rank: int, i: int, j: int) -> Fraction:
    """Inverse-Cartan pairing for labels i,j on A_rank."""
    if i == 0 or j == 0:
        return Fraction(0)
    return Fraction(min(i, j) * (rank + 1 - max(i, j)), rank + 1)


# (D6 discriminant class, I6 label, I4 label, P.O)
profiles = [
    ("0", 2, 1, 0),
    ("v", 5, 0, 0),
    ("v", 2, 1, 1),
]


def d6_correction(c1: str, c2: str) -> Fraction:
    # The chosen D6 vector class has inverse-Cartan norm one.
    return Fraction(1) if c1 == c2 == "v" else Fraction(0)


def local_correction(p, q) -> Fraction:
    return (
        d6_correction(p[0], q[0])
        + a_correction(5, p[1], q[1])
        + a_correction(3, p[2], q[2])
    )


gram = [[Fraction(0) for _ in profiles] for _ in profiles]
for i, p in enumerate(profiles):
    gram[i][i] = Fraction(4 + 2 * p[3]) - local_correction(p, p)
    for j in range(i):
        q = profiles[j]
        section_intersection = 1
        value = (
            Fraction(2 + p[3] + q[3] - section_intersection)
            - local_correction(p, q)
        )
        gram[i][j] = gram[j][i] = value

expected_scaled = [
    [23, 8, -1],
    [8, 26, 8],
    [-1, 8, 35],
]
assert [[12 * x for x in row] for row in gram] == expected_scaled

print("D6A5A3PROFILE|status=PASS|scaled_gram=" + str(expected_scaled))

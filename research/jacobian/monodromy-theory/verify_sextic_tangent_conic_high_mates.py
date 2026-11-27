#!/usr/bin/env python3
"""Replay the finite monomial-support comparisons in the tangent-conic reduction."""

from math import gcd


def monomial_cases():
    rows = {}
    for m in range(4, 13):
        for n in range(3, m):
            divisor = gcd(m, n)
            cases = []
            for b in range(-24, 7):
                a, c = b*(m//divisor), b*(n//divisor)
                p_degrees = tuple(
                    i for i in range(4)
                    if a + i >= -(m//2) and a + 4*i <= 6 - m
                )
                q_degrees = tuple(
                    i for i in range(4)
                    if c + i >= -(n//2) and c + 4*i <= 6 - n
                )
                if p_degrees and q_degrees:
                    cases.append((a, c, p_degrees, q_degrees))
            if cases:
                rows[m, n] = tuple(cases)
    return rows


def constant_case():
    return ((0, 0, (0,), (0,)),)


assert monomial_cases() == {
    (4, 3): constant_case(),
    (5, 3): constant_case(),
    (5, 4): constant_case(),
    (6, 3): ((-4, -2, (1,), (1,)), (-2, -1, (0,), (0, 1))) + constant_case(),
    (6, 4): ((-3, -2, (0,), (0, 1)),) + constant_case(),
    (6, 5): constant_case(),
    (8, 4): ((-4, -2, (0,), (0, 1)), (-2, -1, (0,), (0,))),
    (8, 6): ((-4, -3, (0,), (0,)),),
    (9, 3): ((-3, -1, (0,), (0, 1)),),
    (9, 6): ((-3, -2, (0,), (0,)),),
    (10, 4): ((-5, -2, (0,), (0, 1)),),
    (10, 5): ((-4, -2, (0,), (0,)),),
    (10, 6): ((-5, -3, (0,), (0,)),),
    (10, 8): ((-5, -4, (0,), (0,)),),
    (12, 4): ((-6, -2, (0,), (0, 1)),),
    (12, 6): ((-6, -3, (0,), (0,)),),
    (12, 8): ((-6, -4, (0,), (0,)),),
    (12, 10): ((-6, -5, (0,), (0,)),),
}
print("Tangent-conic high mates: exact monomial-support table verified.")

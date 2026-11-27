#!/usr/bin/env python3
"""Replay the remaining degree pairs and exceptional quadratic bracket identity."""

from math import gcd

from sympy import Poly, expand, symbols


def degree_pairs(removed):
    """Maximal degree of the common power root after removing fibre factors."""
    result = {}
    for m in range(3, 6):
        for n in range(2, m):
            divisor = gcd(m, n)
            max_degree = min(
                (6 - m)//(m//divisor),
                (6 - n - removed)//(n//divisor),
            )
            if max_degree >= 1:
                result[m, n] = max_degree
    return result


assert degree_pairs(0) == {(3, 2): 1, (4, 2): 1}
assert degree_pairs(3) == {(4, 2): 1}

y, r, p, q, f, K, Df, DK = symbols("y r p q f K Df DK")
P_y = 3*p*y**2 + f
Q_y = 2*q*y
DP = Df*y
DQ = DK
bracket = Poly(expand(DP*Q_y - P_y*DQ), y)
assert bracket.nth(0) == -f*DK
assert bracket.nth(2) == expand(2*q*Df - 3*p*DK)

print("Sextic cubic component: degree pairs and centered bracket verified.")

"""Exact finite-order checks of the universal radial cancellation formulas.

The proof in the paper is symbolic for arbitrary m. This script independently
checks the combinatorial sums over a broad integer parameter grid.
"""

from math import comb, factorial
from fractions import Fraction


def general_mixed_direct(m, r, s, u, v, a, d, b, c, mu):
    a, d, b, c = map(Fraction, (a, d, b, c))
    total = Fraction(0)
    for p in range(m + 1):
        for j in range(m + 1):
            if j + r != j + s or m - j + u != p + v:
                continue
            total += (
                comb(m, p)
                * comb(m, j)
                * a ** (m - p)
                * d**p
                * b ** (m - j)
                * c**j
                * mu(j + r, m - j + u)
            )
    return total


def general_mixed_formula(m, r, s, u, v, a, d, b, c, mu):
    a, d, b, c = map(Fraction, (a, d, b, c))
    if r != s:
        return Fraction(0)
    delta = u - v
    x = a * c
    total = Fraction(0)
    for j in range(m + 1):
        shifted = j - delta
        if 0 <= shifted <= m:
            total += (
                (-1) ** j
                * comb(m, j)
                * comb(m, shifted)
                * mu(j + r, m - j + u)
            )
    return a ** (-delta) * d**delta * (-1) ** m * x**m * total


def rising(alpha, n):
    out = Fraction(1)
    for j in range(n):
        out *= alpha + j
    return out


def pure_sum(m, a, d, b, c, rho=Fraction(1)):
    total = Fraction(0)
    for j in range(m + 1):
        total += (
            comb(m, j) ** 2
            * a**j
            * d ** (m - j)
            * b ** (m - j)
            * c**j
            * rho
            * factorial(j)
            * factorial(m - j)
        )
    return total


def mixed_sum(m, q, a, d, b, c, rho=Fraction(1)):
    total = Fraction(0)
    if q > m:
        return total
    for j in range(m - q + 1):
        total += (
            comb(m, j)
            * comb(m, j + q)
            * a ** (j + q)
            * d ** (m - j - q)
            * b ** (m - j)
            * c**j
            * rho
            * factorial(j)
            * factorial(m - j)
        )
    return total


checks = 0

# The general phase-invariant formula is checked against raw expansion for
# a deliberately nonradial moment table.
mu = lambda r, u: Fraction((r + 1) ** 2 + 3 * (u + 1))
for a, d, b in [(1, 2, 3), (2, -3, 5), (-2, 5, 3)]:
    c = Fraction(-d * b, a)
    for m in range(1, 9):
        for r in range(3):
            for s in range(3):
                for u in range(3):
                    for v in range(3):
                        assert general_mixed_direct(
                            m, r, s, u, v, a, d, b, c, mu
                        ) == general_mixed_formula(
                            m, r, s, u, v, a, d, b, c, mu
                        )
                        checks += 1

# Within the two-shape Gamma--Dirichlet family, (1,1) passes every tested
# circuit equation, while the first two equations exclude other sample
# shapes exactly as in the proof.
for m in range(1, 21):
    theta = sum(
        (-1) ** j
        * comb(m, j) ** 2
        * rising(Fraction(1), j)
        * rising(Fraction(1), m - j)
        for j in range(m + 1)
    )
    assert theta == 0
    checks += 1
for alpha, beta in [
    (Fraction(2), Fraction(2)),
    (Fraction(1), Fraction(2)),
    (Fraction(3, 2), Fraction(3, 2)),
]:
    theta1 = beta - alpha
    theta2 = beta * (beta + 1) - 4 * alpha * beta + alpha * (alpha + 1)
    assert theta1 != 0 or theta2 != 0
    checks += 1

for a, d, b in [(1, 1, 1), (2, -3, 5), (-2, 5, 3), (3, 4, -2)]:
    c = Fraction(-d * b, a)
    for m in range(1, 31):
        assert pure_sum(m, a, d, b, c) == 0
        for q in range(1, m + 2):
            expected = (
                0
                if q > m
                else factorial(m)
                * a**q
                * d ** (m - q)
                * b**m
                * comb(m - 1, q - 1)
            )
            assert mixed_sum(m, q, a, d, b, c) == expected
            checks += 1

for a, d, b, c in [(1, 2, 3, 4), (2, -1, 5, 7), (-3, 2, -4, 1)]:
    for m in range(1, 31):
        expected = factorial(m) * (d * b + a * c) ** m
        assert pure_sum(m, a, d, b, c) == expected
        checks += 1

print("exact_checks", checks)
print("max_m", 30)
print("status", "PASS")

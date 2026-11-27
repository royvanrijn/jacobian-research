"""Exact low-gap crossing checks with no differential-order cutoff."""

from fractions import Fraction
from math import gcd


def crossings_at_gap(gap):
    rows = set()
    for difference in range(1, gap):
        for f1 in range(2, gap - difference + 2):
            u, remainder = divmod(f1 * gap, difference)
            if remainder:
                continue
            v = u + gap
            f2 = f1 + difference
            for r in range(1, gap + 1):
                for t in range(r):
                    if (f1 - 1) * (r - t) + difference * r != gap:
                        continue
                    rows.add(((u, v), (r, t), (f1, f2)))
    return rows


def scan_at_gap(gap):
    rows = set()
    for u in range(3, gap * gap + 1):
        v = u + gap
        divisor = gcd(u, v)
        for k in range(1, divisor):
            f1, f2 = k * u // divisor, k * v // divisor
            if f1 < 2:
                continue
            for r in range(1, u):
                for t in range(r):
                    if (f2 - 1) * (u - r) == (f1 - 1) * (v - t):
                        rows.add(((u, v), (r, t), (f1, f2)))
    return rows


def check_gaps():
    survivors = set()
    for gap in range(1, 6):
        rows = crossings_at_gap(gap)
        assert rows == scan_at_gap(gap)
        for C0, C1, (f1, f2) in rows:
            _, v = C0
            r, t = C1
            divisor = gcd(f1 - 1, f2 - 1)
            rho, sigma = (f2 - 1) // divisor, (1 - f1) // divisor
            delta = Fraction(v - t, rho) - t
            x, y = r - delta * Fraction(sigma, rho), t + delta
            if divisor == 1 and y >= 2 and x == y - Fraction(1, rho):
                continue
            survivors.add((gap, C0, C1, (f1, f2), (rho, sigma)))
    assert survivors == {
        (4, (6, 10), (1, 0), (3, 5), (2, -1)),
        (4, (8, 12), (2, 0), (2, 3), (2, -1)),
        (5, (10, 15), (3, 1), (2, 3), (2, -1)),
    }


if __name__ == "__main__":
    check_gaps()
    print("Three crossings before local successor exclusions at diagonal gaps through five.")

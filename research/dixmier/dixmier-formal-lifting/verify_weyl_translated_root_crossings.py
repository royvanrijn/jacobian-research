"""Exact controls for translated-root crossings with first coordinate at most five."""

from fractions import Fraction
from math import gcd


def solve_weight_equations():
    rows = set()
    for u in range(3, 6):
        for f1 in range(2, u):
            for r in range(1, u):
                for t in range(r):
                    denominator = u - f1 * r
                    if denominator <= 0:
                        continue
                    v = Fraction(u * (u - r - (f1 - 1) * t), denominator)
                    f2 = f1 * v / u
                    if v > u and v.denominator == f2.denominator == 1:
                        rows.add(((u, int(v)), (r, t), (f1, int(f2))))
    return rows


def scan_weight_equations():
    rows = set()
    for u in range(3, 6):
        for v in range(u + 1, u * (u - 1) + 1):
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


def check_crossings():
    rows = solve_weight_equations()
    assert rows == scan_weight_equations()
    assert rows == {
        ((3, 6), (1, 0), (2, 4)),
        ((4, 6), (1, 0), (2, 3)),
        ((4, 12), (1, 0), (3, 9)),
        ((5, 10), (2, 1), (2, 4)),
        ((5, 10), (1, 0), (3, 6)),
        ((5, 15), (2, 0), (2, 6)),
        ((5, 20), (1, 0), (4, 16)),
    }

    remaining = set()
    rejected = {}
    for C0, C1, (f1, f2) in rows:
        _, v = C0
        r, t = C1
        divisor = gcd(f2 - 1, f1 - 1)
        rho, sigma = (f2 - 1) // divisor, (1 - f1) // divisor
        delta = Fraction(v - t, rho) - t
        C2 = (r - delta * Fraction(sigma, rho), t + delta)
        if divisor == 1 and C2[1] >= 2 and C2[0] == C2[1] - Fraction(1, rho):
            rejected[C0, C1] = C2
        else:
            remaining.add((C0, C1, (rho, sigma)))
    assert rejected == {
        ((3, 6), (1, 0)): (Fraction(5, 3), 2),
        ((4, 6), (1, 0)): (Fraction(5, 2), 3),
        ((5, 10), (2, 1)): (Fraction(8, 3), 3),
        ((5, 10), (1, 0)): (Fraction(9, 5), 2),
    }
    assert remaining == {
        ((4, 12), (1, 0), (4, -1)),
        ((5, 15), (2, 0), (5, -1)),
        ((5, 20), (1, 0), (5, -1)),
    }


def check_order_ten_scalar_successors():
    candidates = set()
    for multiplicity in (2, 3, 4):
        beta = 4 * multiplicity - 5
        for total in range(5, 2 * beta + 1, 2):
            scaled_slope = 5 - Fraction(beta * total, total * multiplicity - 1)
            denominator = scaled_slope.denominator
            mate_quotient = total - 2
            if denominator <= 2 * multiplicity and (
                2 * multiplicity % denominator == 0
                or mate_quotient * multiplicity % denominator == 0
            ):
                candidates.add((multiplicity, mate_quotient, scaled_slope / 5))
    assert candidates == {(2, 3, Fraction(2, 3)), (3, 3, Fraction(1, 2))}


def check_gcd_five_quotients():
    # The written denominator argument gives the unbounded families.
    # Here an independent divisibility scan checks their parameterization.
    limit = 200
    scanned = set()
    for total in range(5, limit + 1):
        for a in range(2, total // 2 + 1):
            b = total - a
            if gcd(a, b) != 1:
                continue
            for multiplicity in (2, 3, 4):
                beta = 4 * multiplicity - 5
                denominator = (5 - Fraction(beta * total, multiplicity * total - 1)).denominator
                if a % denominator == 0 or b % denominator == 0:
                    scanned.add((multiplicity, a, b))

    table = set()
    for j in range(limit):
        rows = [
            (2, 2 * j + 3, j + 2),
            (3, 3 * j + 2, 4 * j + 3),
            (4, 4 * j + 5, 7 * j + 9),
        ]
        if j >= 2 and j % 2 == 0:
            rows.append((3, 6 * j + 4, j + 1))
        if j % 2 == 1:
            rows.append((4, 8 * j + 10, 3 * j + 4))
        for multiplicity, a, b in rows:
            if a + b <= limit:
                table.add((multiplicity, min(a, b), max(a, b)))
    assert scanned == table
    assert {b for _, a, b in table if a == 2} == {3}
    assert {b if a == 3 else a for _, a, b in table if 3 in (a, b)} == {2, 5, 16}


if __name__ == "__main__":
    check_crossings()
    check_order_ten_scalar_successors()
    check_gcd_five_quotients()
    print("Seven integer crossings, four forbidden successors, and gcd-five quotient controls passed.")

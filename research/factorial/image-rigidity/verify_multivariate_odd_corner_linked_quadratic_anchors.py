#!/usr/bin/env python3
"""Exact checks for the linked quadratic anchors at odd-terminal corners."""

from collections import Counter
from fractions import Fraction
from math import comb, gcd, lcm, prod
import json


def coefficients(rows, order, height):
    height = Fraction(height)
    num = height.numerator
    den = height.denominator
    prefixes = []
    for t, d, sign in rows:
        prefix = [1]
        for r in range(1, order*d+1):
            prefix.append(prefix[-1] * (order*num + den*(r-order*t)))
        prefixes.append((d, sign, prefix))
    output = []
    for k in range(order + 1):
        lengths = [(d * (order-k if sign > 0 else k), prefix) for d, sign, prefix in prefixes]
        integer = comb(order, k) * prod(prefix[length] for length, prefix in lengths)
        output.append(Fraction(integer, den**sum(length for length, _ in lengths)))
    return output


def prime_divisors(value):
    out = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            out.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        out.append(value)
    return out


def integer_valuation(value, prime):
    assert value
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def valuation(value, prime):
    if value == 0:
        return float("inf")
    value = Fraction(value)
    return integer_valuation(value.numerator, prime) - integer_valuation(value.denominator, prime)


def leading_ratio(rows, degree, offset):
    value = Fraction(comb(degree, 2))
    for t, d, sign in rows:
        if sign < 0:
            value *= prod(range(degree * t + offset - 2 * d, degree * t + offset))
        else:
            start = degree * (t - d) + offset
            value /= prod(range(start, start + 2 * d))
    return value


def increment(rows, degree, offset, prime):
    result = 0
    for t, d, sign in rows:
        start = degree * (t - d if sign > 0 else t) + offset
        result += valuation(start, prime) - valuation(start + sign * 2 * d, prime)
    return result


def remainder(quadratic, other):
    work = list(map(Fraction, other))
    a, b, c = map(Fraction, quadratic)
    for degree in range(len(work) - 1, 1, -1):
        factor = work[degree] / c
        work[degree - 2] -= factor * a
        work[degree - 1] -= factor * b
    return work[:2]


def sylvester(marked, other):
    m = len(marked) - 1
    n = len(other) - 1
    columns = [
        [0] * shift + values + [0] * (m + n - shift - len(values))
        for values, shifts in ((marked, n), (other, m))
        for shift in range(shifts)
    ]
    return [list(row) for row in zip(*columns)]


def determinant(matrix):
    a = [row.copy() for row in matrix]
    size = len(a)
    old = 1
    sign = 1
    for i in range(size - 1):
        pivot = next((k for k in range(i, size) if a[k][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            sign = -sign
        for j in range(i + 1, size):
            for k in range(i + 1, size):
                a[j][k] = (a[j][k] * a[i][i] - a[j][i] * a[i][k]) // old
        old = a[i][i]
        for j in range(i + 1, size):
            a[j][i] = 0
    return sign * a[-1][-1]


def minor(matrix, row, column):
    return [
        entries[:column] + entries[column + 1:]
        for index, entries in enumerate(matrix)
        if index != row
    ]


def integer_scale(matrix):
    denominator = lcm(*(Fraction(value).denominator for row in matrix for value in row))
    return denominator, [[int(value * denominator) for value in row] for row in matrix]


def minor_content(matrix):
    scale, integer = integer_scale(matrix)
    dimension = len(integer)
    content = 0
    for row in range(dimension):
        for column in range(dimension):
            content = gcd(content, determinant(minor(integer, row, column)))
            if content == 1:
                return 1
    return content // gcd(content, scale**(dimension-1))


def profiles(m):
    q = m // 2
    corner = ((1, 1, 1), (1, q, -1))
    return (
        corner,
        corner * 3,
        corner + ((5, 4, 1), (7, 3*q, -1), (2, 2*q, -1)),
        corner + ((2, 2, 1), (3, 3, 1), (6, 4*q, -1), (29, 27, 1)),
        corner + ((1009, 17, 1), (997, 2*q, -1), (1, 1, 1)),
    )


def verify_anchors(counts):
    for m in (6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 30, 32, 40):
        n = m + 1
        q = m // 2
        p = next((prime for prime in prime_divisors(m) if prime > 2), 2)
        pair = (p-2, p-1) if p > 2 else (1, 3)
        for rows in profiles(m):
            delta = sum(d * sign for _, d, sign in rows)
            chosen_ratios = [leading_ratio(rows, m, h) for h in pair]
            values = [valuation(value, p) for value in chosen_ratios]
            assert values[0] > values[1]
            assert values[1] - values[0] == sum(increment(rows, m, h, p) for h in range(*pair))
            assert any(value != 2 * valuation(m, p) * delta for value in values)
            counts["strict_linked_anchors"] += 1
            for s in (m, n):
                other = m+n-s
                for h in range(1, 3*q-s+1):
                    anchor = -Fraction(h, s)
                    P = coefficients(rows, s, anchor)
                    Q = coefficients(rows, other, anchor)
                    assert all(P[:3]) and not any(P[3:]) and all(Q)
                    assert P[2]/P[0] == leading_ratio(rows, s, h)
                    counts["quadratic_specializations"] += 1
                    if s == m and h in pair and valuation(P[2]/P[0], p) != 2 * valuation(m, p) * delta:
                        counts["nonzero_selected_prefilters"] += 1
                        if m <= 16:
                            assert any(remainder(P[:3], Q))
                            counts["nonzero_selected_remainders"] += 1
                    if s == m and h in pair:
                        target = valuation(m, p) * delta
                        lower = valuation(Q[0], p)
                        assert all(
                            valuation(c, p) == lower + k * target + valuation(comb(n, k), p)
                            for k, c in enumerate(Q)
                        )
                        counts["unmarked_affine_valuations"] += 1


def verify_minors(counts):
    for m in (6, 8):
        n = m + 1
        rows = profiles(m)[0]
        X = max(t for t, _, _ in rows) + 17
        for s in (m, n):
            other = m+n-s
            for h in (1, 3*(m//2)-s):
                anchor = -Fraction(h, s)
                modulus = (s*X+h)//gcd(s, h)
                P = coefficients(rows, s, anchor)
                Q = coefficients(rows, other, anchor)
                matrix = sylvester(P, Q)
                B0 = coefficients(rows, s, X)
                B1 = coefficients(rows, other, X)
                ordinary = sylvester(B0, B1)
                scale, integer = integer_scale(matrix)
                size = len(matrix)
                for row, column in ((0, 0), (0, size-1), (size-1, 0), (size-1, size-1), (2, 3)):
                    d0 = determinant(minor(ordinary, row, column))
                    d1 = determinant(minor(integer, row, column))
                    assert (d0 * scale**(size-1) - d1) % modulus == 0
                    counts["minor_congruences"] += 1
                D = minor_content(matrix)
                R = remainder(P[:3], Q)
                assert (D == 0) == (not any(R))
                counts["minor_rank_checks"] += 1
                for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29):
                    if gcd(prime, scale) == 1 and valuation(P[2], prime) == valuation(Q[-1], prime) == 0:
                        assert valuation(D, prime) == min(valuation(c, prime) for c in R)
                        counts["local_remainder_ideals"] += 1

    for prime, depth in ((3, 2), (5, 3), (7, 4)):
        P = [2, 1, 1]
        Q = [2 + prime**depth, 5 + prime**(depth+1), 5, 3, 1]
        R = remainder(P, Q)
        matrix = sylvester(P + [0]*4, Q)
        D = minor_content(matrix)
        assert valuation(D, prime) == min(valuation(c, prime) for c in R) == depth
        counts["deep_remainder_ideals"] += 1


def run():
    counts = Counter()
    verify_anchors(counts)
    verify_minors(counts)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()

#!/usr/bin/env python3
"""Exact replay of simultaneous odd-anchor divisibility and private witnesses."""

from collections import Counter
from math import comb, gcd, isqrt, prod
import json

from verify_multivariate_odd_corner_linked_quadratic_anchors import (
    determinant,
    leading_ratio,
    minor,
    sylvester,
    valuation,
)


def factorization(value):
    factors = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        factors[value] = factors.get(value, 0) + 1
    return factors


def is_prime(value):
    return value > 1 and all(value % divisor for divisor in range(2, isqrt(value) + 1))


def hits(row, order, modulus):
    a, d, _ = row
    return order * (a + d) // modulus != order * a // modulus


def neutral(row, witnesses):
    return all(
        not hits(row, order, prime)
        for terminal, prime, _ in witnesses
        for order in (terminal-1, terminal)
    )


def spaced(row, rows, spacing):
    return all(abs(row[0] - earlier[0]) > spacing for earlier in rows)


def construction():
    B = 30064
    m = 6
    n = m+1
    q = m//2
    cutoff = 9
    rows = [(B-1, 1, 1)]*(2*q) + [(B-1, q, -1)]*2
    witnesses = []
    for terminal in range(2, cutoff):
        sign = 1 if terminal % 2 else -1
        for a in range((B+7)//8, B//6+1):
            row = (a, 1, sign)
            prime = terminal*a+1
            if not spaced(row, rows, 2*cutoff) or not neutral(row, witnesses):
                continue
            if any(hits(old, order, prime) for old in rows for order in (terminal-1, terminal)):
                continue
            if is_prime(prime):
                witnesses.append((terminal, prime, len(rows)))
                rows.append(row)
                break
        else:
            raise AssertionError("required witness unavailable")
    balance = sum(d*sign for _, d, sign in rows)
    if balance:
        assert abs(balance) == 1
        for a in range((B+7)//8, B//6+1):
            row = (a, 1, -balance)
            if spaced(row, rows, 2*cutoff) and neutral(row, witnesses):
                rows.append(row)
                break
        else:
            raise AssertionError("balancing witness unavailable")
    fan = [(s, h) for s in (m, n) for h in range(1, 3*q-s+1)]
    targets = []
    for s, h in fan:
        t = m+n-s
        V = (s*B+h)//gcd(s,h)
        for p, E in factorization(V).items():
            power = p**(E+valuation(t, p))
            targets.append((power, s, h, t, p, E))
    targets.sort(reverse=True)
    for power, _, _, t, _, _ in targets:
        for sign in (1, -1):
            for a in range((B+3)//4, B//3+1):
                row = (a, t, sign)
                if (t*(a+1))//power == (t*a)//power:
                    continue
                if spaced(row, rows, 2*cutoff) and neutral(row, witnesses):
                    rows.append(row)
                    break
            else:
                raise AssertionError("required residue row unavailable")
    assert sum(d*sign for _, d, sign in rows) == 0
    assert all(a+d <= B for a, d, sign in rows if sign > 0)
    assert all(2*d <= m*(B-a) for a, d, sign in rows if sign < 0)
    return B, m, rows, witnesses, targets


def primitive_representation(rows, order):
    cache = []
    common = None
    scalar_content = 0
    for k in range(order+1):
        line_powers = Counter()
        constant = comb(order, k)
        for a, d, sign in rows:
            draw = order-k if sign > 0 else k
            for h in range(1, d*draw+1):
                y = order*a+h
                line_powers[y] += 1
                constant *= gcd(order, y)
        common = line_powers if common is None else common & line_powers
        scalar_content = gcd(scalar_content, constant)
        cache.append((constant, line_powers))
    return [(constant//scalar_content, line_powers-common) for constant, line_powers in cache]


def evaluate_modulo(representation, order, shift, modulus):
    result = []
    for constant, lines in representation:
        value = constant % modulus
        for y, exponent in lines.items():
            common = gcd(order, y)
            primitive = ((order//common)*shift+y//common) % modulus
            value = value * pow(primitive, exponent, modulus) % modulus
        result.append(value)
    return result


def verify_crt_depth():
    m = 6
    pairs = ((6, 1), (6, 2), (6, 3), (7, 1), (7, 2))
    fixed_primes = (101, 103, 107, 109, 113)
    base, modulus = 0, 1
    for (s, h), prime in zip(pairs, fixed_primes):
        common = gcd(s, h)
        q = prime**2
        residue = (prime-h//common)*pow(s//common, -1, q) % q
        base += modulus*((residue-base)*pow(modulus, -1, q) % q)
        modulus *= q
    assert all(
        valuation((s*base+h)//gcd(s, h), prime) == 1
        for (s, h), prime in zip(pairs, fixed_primes)
    )
    return len(fixed_primes)


def run():
    B, m, rows, witnesses, targets = construction()
    counts = Counter({"crt_depths": verify_crt_depth(), "ordinary_row_count": len(rows)})
    for terminal, prime, marked in witnesses:
        hit_rows = [
            (order, index)
            for order in (terminal-1, terminal)
            for index, row in enumerate(rows)
            if hits(row, order, prime)
        ]
        assert hit_rows == [(terminal, marked)]
        at = rows[marked]
        assert at[1] == 1 and prime == terminal*at[0]+1
        counts["private_prime_certificates"] += 1

    representations = {j: primitive_representation(rows, j) for j in (m, m+1)}
    relative = [(B-a, d, sign) for a, d, sign in rows]
    assert valuation(leading_ratio(relative, m, 1), 3) > valuation(leading_ratio(relative, m, 2), 3)
    for _, s, h, t, p, E in targets:
        small_modulus = p**E
        assert not any(evaluate_modulo(representations[t], t, 0, small_modulus))
        counts["primitive_coefficient_layers"] += t+1
        modulus = p**((s-1)*E)
        common = gcd(s, h)
        shift = (-B - (h//common)*pow(s//common, -1, modulus)) % modulus
        marked = evaluate_modulo(representations[s], s, shift, modulus)
        other = evaluate_modulo(representations[t], t, shift, modulus)
        assert not any(value % small_modulus for value in other)
        assert not any(marked[3:])
        matrix = sylvester(marked, other)
        dimension = len(matrix)
        for row, column in ((0, 0), (0, dimension-1), (dimension-1, 0), (dimension-1, dimension-1), (3, 4)):
            assert determinant(minor(matrix, row, column)) % modulus == 0
            counts["primitive_minor_prime_power_checks"] += 1
        counts["joint_fan_prime_powers"] += 1

    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()

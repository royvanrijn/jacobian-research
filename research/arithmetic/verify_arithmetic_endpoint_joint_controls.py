#!/usr/bin/env python3
"""Finite exact controls for (365.J2), (675.W1)--(675.W3), and the 2376 return.

Run from the repository root with Python 3, standard library only:
    python3 research/arithmetic/verify_arithmetic_endpoint_joint_controls.py

The Pell range is 2 <= D <= 5000, squarefree D only. The CRT cover is
exhausted for 1 <= N <= 12 and 1 <= t <= 500. These are controls of the
normalization and finite cover, not proofs of the asymptotic claims or
of the unrestricted endpoints. No floating point arithmetic is used.
"""

import json
from collections import Counter
from math import gcd, isqrt, prod


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def primes_through(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"[:limit + 1]
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            for n in range(p * p, limit + 1, p):
                flags[n] = 0
    return [p for p in range(2, limit + 1) if flags[p]]


def squarefree_sieve(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0] = 0
    for p in primes_through(isqrt(limit)):
        for n in range(p * p, limit + 1, p * p):
            flags[n] = 0
    return flags


def pell(D):
    """Continued-fraction convergents, stopping at the first norm-one unit."""
    root = isqrt(D)
    require(root * root != D, ("nonsquare", D))
    m, denominator, a = 0, 1, root
    p0, p1, q0, q1 = 0, 1, 1, 0
    while True:
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
        if p1 * p1 - D * q1 * q1 == 1:
            return p1, q1
        m = denominator * a - m
        denominator = (D - m * m) // denominator
        a = (root + m) // denominator


def check_pell_normalization():
    flags = squarefree_sieve(5000)
    branches = Counter()
    for D in range(2, 5001):
        if not flags[D]:
            continue
        T, U = pell(D)
        H = gcd(D, U)
        R, Q = D // H, U // H
        h = H // 2 if H % 2 == 0 else H
        A, B = (T - 1, T + 1) if T % 2 == 0 else ((T - 1) // 2, (T + 1) // 2)
        d, e = gcd(h, A), h // gcd(h, A)
        f = gcd(A, D // h)
        g = (D // h) // f
        require(A % (d**3 * f) == B % (e**3 * g) == 0, ("divisibility", D))
        u2, v2 = A // (d**3 * f), B // (e**3 * g)
        u, v = isqrt(u2), isqrt(v2)
        require(u * u == u2 and v * v == v2, ("square extraction", D))
        require(d * e * f * g == D and gcd(d, e) == 1, ("kernel", D))
        require(gcd(f * u2, g * v2) == 1, ("content", D))
        require(g * v2 * e**3 - f * u2 * d**3 == B - A, ("Thue", D))
        if T % 2 == 0:
            key, expected = "T_even", (R, Q, 2)
        elif D % 2:
            require(Q % 2 == 0, ("even Q", D))
            key, expected = "T_odd_D_odd", (R, Q // 2, 1)
        else:
            key, expected = "T_odd_D_even", (2 * R, Q, 1)
        require((f * g, u * v, B - A) == expected, ("prime two table", D))
        branches[key] += 1
    require(len(branches) == 3, "all parity branches must occur")
    return dict(branches)


def joint_class(positive, assignment):
    """None rejects an overlap conflict or a hit on a positive site."""
    used = {}
    for b, q in assignment.items():
        residue = (-b) % (q * q)
        if q in used and used[q] != residue:
            return None
        if any((a - b) % (q * q) == 0 for a in positive):
            return None
        used[q] = residue
    r, M = 0, 1
    for q, residue in sorted(used.items()):
        m = q * q
        r += M * (((residue - r) * pow(M, -1, m)) % m)
        M *= m
    require(0 <= r < M, ("canonical CRT representative", r, M))
    return r, M


def class_count(X, r, M):
    return (X - r) // M - (-r) // M


def check_finite_cover():
    X = 500
    flags = squarefree_sieve(X + 12)
    assignment_count = 0
    return_count = 0
    for N in range(1, 13):
        positive = [a for a in range(1, N + 1) if flags[a]]
        negative = [a for a in range(1, N + 1) if not flags[a]]
        primes = primes_through(isqrt(X + N))
        covered, passing = set(), set()
        count_sum = 0

        def visit(index, assignment):
            nonlocal assignment_count, count_sum
            if index < len(negative):
                b = negative[index]
                for q in primes:
                    assignment[b] = q
                    if joint_class(positive, assignment) is not None:
                        visit(index + 1, assignment)
                assignment.pop(b)
                return
            assignment_count += 1
            r, M = joint_class(positive, assignment)
            candidates = set(range(r if r else M, X + 1, M))
            require(len(candidates) == class_count(X, r, M), ("floor formula", N, r, M))
            count_sum += len(candidates)
            covered.update(candidates)
            passing.update(t for t in candidates if all(flags[t + a] for a in positive))

        visit(0, {})
        direct = {t for t in range(1, X + 1)
                  if all(flags[a] == flags[t + a] for a in range(1, N + 1))}
        require(direct <= covered and passing == direct, ("exact finite cover", N))
        require(len(direct) <= count_sum, ("union bound", N))
        return_count += len(direct)
    require(joint_class([], {4: 2, 8: 2}) == (0, 4), "legal overlap")
    require(joint_class([], {4: 2, 9: 2}) is None, "reject incompatible overlap")
    require(joint_class([1], {5: 2}) is None, "reject positive-site hit")
    require(class_count(3, 0, 4) == 0 and class_count(4, 0, 4) == 1, "exclude t=0")
    return {"prefixes": 12, "height": X, "safe_assignments": assignment_count,
            "returns_summed_over_prefixes": return_count}


def check_exchange():
    t, N = 2376, 49
    flags = squarefree_sieve(t + N + 1)
    require(all(flags[a] == flags[t + a] for a in range(1, N + 1)), "return through 49")
    require(not flags[50] and flags[t + 50], "failure at 50")
    positive = [a for a in range(1, N + 1) if flags[a]]
    negative = [a for a in range(1, N + 1) if not flags[a]]
    primes = primes_through(isqrt(t + N))
    assignment = {b: next(q for q in primes if (t + b) % (q * q) == 0)
                  for b in negative}
    require(assignment[25] == 7 and assignment[49] == 5, "witness exchange")
    require(joint_class(positive, assignment) == (2376, 44100), "small joint representative")
    require(prod(q * q for q in set(assignment.values())) == 44100, "squarefree prime modulus")
    return {"prefix": N, "shift": t, "modulus": 44100, "first_failure": 50}


if __name__ == "__main__":
    result = {"status": "PASS", "pell_squarefree_radicals_through": 5000,
              "pell_parity_branches": check_pell_normalization(),
              "finite_crt_cover": check_finite_cover(), "witness_exchange": check_exchange(),
              "asymptotic_or_endpoint_verification": False}
    print(json.dumps(result, indent=2, sort_keys=True))

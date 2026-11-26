"""Finite diagnostic for rStation Paper 125 Lemma 2 on one abc fibre.

Run: .venv/bin/python research/arithmetic/verify_abc_determinant_ratio_fibre.py
This counts only the stated finite parameter box; it proves no asymptotic bound.
"""

from collections import Counter, defaultdict
from math import gcd
from sympy import isprime

N = 10**10
p, q = 10000019, 10000169
assert isprime(p) and isprime(q) and p != q
A, Q, d = N, 3 * N, p * q


def egcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        k = old_r // r
        old_r, r = r, old_r - k * r
        old_s, s = s, old_s - k * s
        old_t, t = t, old_t - k * t
    return old_r, old_s, old_t


rows = []
for u in range((3 * N + p - 1) // p, (4 * N - 1) // p + 1):
    x = p * u
    if x % 72 != 1:
        continue
    for v in range((N + q - 1) // q, (2 * N - 1) // q + 1):
        y = q * v
        if y % 72 != 12 or gcd(u, v) != 1 or not x < Q + y:
            continue
        assert gcd(x, y) == 1
        rows.append((y, x, u, v))

assert rows
profiles = defaultdict(Counter)
for a, b, _, _ in rows:
    for c, e, _, _ in rows:
        D = abs(a * e - b * c)
        if not D:
            continue
        assert D % d == 0 and D < 6 * A * Q
        g, z, w = egcd(c, e)
        assert g == 1
        h = (z * a + w * b) % d
        assert gcd(h, d) == 1
        profiles[D][h] += 1

pair_count = sum(sum(cs.values()) for cs in profiles.values())
ratio_collision_count = sum(sum(v * (v - 1) for v in cs.values()) for cs in profiles.values())
det_collision_count = sum(n * (n - 1) for n in (sum(cs.values()) for cs in profiles.values()))
max_fibre = max(max(cs.values()) for cs in profiles.values())
max_lemma_rhs = max(2 * (D + 8 * A * A) // d for D in profiles)
assert all(max(cs.values()) <= 2 * (D + 8 * A * A) / d for D, cs in profiles.items())

print('N,p,q,A,Q,d', N, p, q, A, Q, d)
print('rows, ordered nonzero-determinant pairs', len(rows), pair_count)
print('distinct D, distinct (D,h) cells', len(profiles), sum(map(len, profiles.values())))
print('ordered same-D pair collisions', det_collision_count)
print('ordered same-(D,h) pair collisions', ratio_collision_count)
print('largest actual cell, largest Lemma 2 RHS', max_fibre, max_lemma_rhs)
print('weight per abc triple: 1; sum_h cell masses = N_D for every D')
print('ratio is not fixed by the p,q record: h mod p = v/v_prime, h mod q = u/u_prime')

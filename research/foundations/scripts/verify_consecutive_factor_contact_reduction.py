#!/usr/bin/env python3
"""Exact checks for the consecutive-factor contact reduction."""

from itertools import product
import sympy as sp

L = sp.Symbol("L")


def projective(n):
    return sum(L**i for i in range(n + 1))


def coprime_class(a, b):
    if a == 0 or b == 0:
        return projective(a + b)
    lower = sum(
        projective(d) * coprime_class(a - d, b - d)
        for d in range(1, min(a, b) + 1)
    )
    return sp.expand(projective(a) * projective(b) - lower)


for a in range(1, 7):
    for b in range(1, 7):
        assert sp.expand(coprime_class(a, b) - L ** (a + b - 1) * (L + 1)) == 0

P2 = projective(2)
complete_hyperplane = sp.expand(P2**2 + L**3)
positive_gcd = sp.expand(P2**2 - L**4)
section = sp.expand(complete_hyperplane - positive_gcd)
assert section == L**4 + L**3
assert sp.expand(coprime_class(2, 3) - section) == L**5 - L**3

# Closed p=2 formula (17). In the rank-one osculating case,
# K=Q=T=P^1, O_3=A^1, and O_0=pt.
rank_one_class = sp.expand(
    L**5
    - L**3
    + L * (projective(1) - projective(1))
    + L**2 * L
    + L**3 * (1 - projective(1))
)
assert rank_one_class == L**5 - L**4


def projective_points(n, q):
    for first in range(n + 1):
        for tail in product(range(q), repeat=n - first):
            yield (0,) * first + (1,) + tail


def trim(f, q):
    f = list(f)
    while f and f[0] % q == 0:
        f.pop(0)
    return f


def remainder(f, g, q):
    f, g = trim(f, q), trim(g, q)
    while f and len(f) >= len(g):
        c = f[0] * pow(g[0], -1, q) % q
        for i, value in enumerate(g):
            f[i] = (f[i] - c * value) % q
        f = trim(f, q)
    return f


def coprime(f, g, q):
    if f[0] == g[0] == 0:
        return False
    f, g = trim(f, q), trim(g, q)
    while g:
        f, g = g, remainder(f, g, q)
    return len(f) == 1


def multiply(f, g, q):
    answer = [0] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        for j, y in enumerate(g):
            answer[i + j] = (answer[i + j] + x * y) % q
    return answer


q = 7
hyperplanes = {
    (5,): (0, 0, 0, 0, 0, 1),
    (4, 1): (0, 0, 0, 0, 4, 1),
    (3, 2): (0, 0, 0, 5, 1, 1),
}
counts = dict.fromkeys(hyperplanes, 0)
total = 0
for quadratic in projective_points(2, q):
    for cubic in projective_points(3, q):
        if not coprime(quadratic, cubic, q):
            continue
        total += 1
        quintic = multiply(quadratic, cubic, q)
        for partition, functional in hyperplanes.items():
            if sum(x * y for x, y in zip(functional, quintic)) % q:
                counts[partition] += 1

assert total == q**4 * (q + 1)
assert counts == {(5,): 14406, (4, 1): 16464, (3, 2): 16513}
assert all(value != q**5 for value in counts.values())

# Two squarefree quintics with five rational roots over F_11.  Their contact
# partitions agree, but their configurations are not PGL_2-equivalent and
# the complement counts differ.
q = 11
simple_hyperplanes = (
    (0, 7, 6, 9, 9, 1),  # roots 0,1,2,3,4
    (0, 6, 6, 3, 0, 1),  # roots 0,1,2,3,5
)
simple_counts = [0, 0]
cubics = tuple(projective_points(3, q))
for quadratic in projective_points(2, q):
    for cubic in cubics:
        if not coprime(quadratic, cubic, q):
            continue
        quintic = multiply(quadratic, cubic, q)
        for index, functional in enumerate(simple_hyperplanes):
            if sum(x * y for x, y in zip(functional, quintic)) % q:
                simple_counts[index] += 1
assert simple_counts == [158620, 159742]

print("PASS consecutive contact reduction: universal coprime class")
print("PASS consecutive contact reduction: p=2 tangent class L^5-L^3")
print("PASS consecutive contact reduction: rank-one class L^5-L^4")
print("PASS representative F_7 contact counts")
print("PASS F_11 cross-ratio sensitivity within the squarefree contact stratum")

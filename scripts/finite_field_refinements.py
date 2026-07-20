#!/usr/bin/env python3
"""Check refined fiber, factorization, discriminant, and collision counts."""

from collections import Counter, defaultdict
from itertools import product


def evaluate(point, p):
    x, y, z = point
    u = (1 + x * y) % p
    return (
        (u**3 * z + y**2 * u * (4 + 3 * x * y)) % p,
        (y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) % p,
        (2 * x - 3 * x**2 * y - x**3 * z) % p,
    )


def factorization_type(a, b, c, p):
    """Factorization partition of c*T^3-2*T^2+b*T-2*a for c != 0."""
    values = [
        (c * t**3 - 2 * t**2 + b * t - 2 * a) % p
        for t in range(p)
    ]
    roots = [t for t, value in enumerate(values) if value == 0]
    derivative = lambda t: (3 * c * t**2 - 4 * t + b) % p
    if len(roots) == 3:
        return "1+1+1"
    if len(roots) == 2:
        return "2+1"
    if len(roots) == 1:
        r = roots[0]
        if derivative(r) != 0:
            return "1+2"
        return "3"
    return "3-irreducible"


def discriminant_class(a, b, c, p):
    # Disc_T(c*T^3 - 2*T^2 + b*T - 2*a) = -4Q.
    qpoly = (27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b) % p
    disc = (-4 * qpoly) % p
    if disc == 0:
        return "zero"
    return "square" if pow(disc, (p - 1) // 2, p) == 1 else "nonsquare"


def expected(p):
    split = (p - 1) * (p - 2) // 6
    return {
        "c=0 fibers": {1: p * (p + 1) // 2, 3: p * (p - 1) // 2},
        "c!=0 fibers": {
            0: (p - 1) * (p * p + 2) // 3,
            1: (p - 1) ** 2 * (p + 2) // 2,
            3: (p - 1) ** 2 * (p - 2) // 6,
        },
        "factorization per c": {
            "1+1+1": split,
            "1+2": p * (p - 1) // 2,
            "2+1": p - 1,
            "3": 1,
            "3-irreducible": (p * p - 1) // 3,
        },
        "discriminant per c": {
            "nonsquare": p * (p - 1) // 2,
            "square": p * (p - 1) // 2,
            "zero": p,
        },
        "off-diagonal double": (p - 1) * (p * p + 2),
        "distinct triple": (p - 1) * (p * p + 2),
    }


def check(p):
    fibers = defaultdict(list)
    for source in product(range(p), repeat=3):
        fibers[evaluate(source, p)].append(source)

    strata = {"c=0 fibers": Counter(), "c!=0 fibers": Counter()}
    factorization = Counter()
    discriminants_by_c = defaultdict(Counter)
    off_diagonal_double = 0
    distinct_triple = 0

    for target in product(range(p), repeat=3):
        a, b, c = target
        size = len(fibers[target])
        strata["c=0 fibers" if c == 0 else "c!=0 fibers"][size] += 1
        discriminants_by_c[c][discriminant_class(a, b, c, p)] += 1
        if c != 0:
            factorization[factorization_type(a, b, c, p)] += 1
        off_diagonal_double += size * (size - 1)
        distinct_triple += size * (size - 1) * (size - 2)

    wanted = expected(p)
    assert dict(sorted(strata["c=0 fibers"].items())) == wanted["c=0 fibers"]
    assert dict(sorted(strata["c!=0 fibers"].items())) == wanted["c!=0 fibers"]

    # Scaling the nonzero root sum makes the factorization law identical for
    # every nonzero c; compare after dividing the aggregate by p-1.
    factorization_per_c = {key: value // (p - 1) for key, value in factorization.items()}
    assert factorization_per_c == wanted["factorization per c"]
    for counts in discriminants_by_c.values():
        assert dict(sorted(counts.items())) == wanted["discriminant per c"]
    assert off_diagonal_double == wanted["off-diagonal double"]
    assert distinct_triple == wanted["distinct triple"]
    print(f"PASS F_{p}: strata, factorization, discriminants, and fiber products")


if __name__ == "__main__":
    for prime in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        check(prime)

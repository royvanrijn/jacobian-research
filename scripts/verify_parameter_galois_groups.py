#!/usr/bin/env python3
"""Exact Galois-group classification for C24 parameter degree at most six."""
from __future__ import annotations

import math

import sympy as sp

from master_cancellation import parameter_polynomial


q = sp.symbols("q")

# name, expected order
EXPECTED = {
    (1, 1): ("1", 1),
    (1, 2): ("S2", 2),
    (1, 3): ("S3", 6),
    (1, 4): ("S4", 24),
    (1, 5): ("S5", 120),
    (1, 6): ("S6", 720),
    (2, 1): ("S2", 2),
    (2, 2): ("D4", 8),
    (2, 3): ("PGL2(5)", 120),
    (3, 1): ("S3", 6),
    (3, 2): ("S6", 720),
    (4, 1): ("S4", 24),
    (5, 1): ("S5", 120),
    (6, 1): ("PGL2(5)", 120),
}


def cycle_type(permutation: sp.combinatorics.Permutation, degree: int) -> tuple[int, ...]:
    lengths = [len(cycle) for cycle in permutation.cyclic_form]
    lengths.extend([1] * (degree - sum(lengths)))
    return tuple(sorted(lengths, reverse=True))


for pair, (name, expected_order) in EXPECTED.items():
    m, r = pair
    degree = m * r
    polynomial = sp.Poly(parameter_polynomial(m, r, q), q, domain=sp.QQ)
    assert polynomial.degree() == degree
    assert polynomial.is_irreducible

    if degree == 1:
        order = 1
        in_alternating = True
        cycle_types = {(1,)}
    else:
        group, in_alternating = sp.polys.numberfields.galois_group(
            polynomial.as_expr(), q
        )
        order = int(group.order())
        assert group.is_transitive()
        cycle_types = {cycle_type(element, degree) for element in group.elements}

    assert order == expected_order
    if name == "1":
        assert degree == 1 and order == 1
    elif name.startswith("S"):
        assert order == math.factorial(degree)
    elif name == "D4":
        assert degree == 4 and (4,) in cycle_types and (2, 2) in cycle_types
    else:
        assert degree == 6 and (5, 1) in cycle_types and (2, 2, 2) in cycle_types

    if degree >= 2:
        assert not in_alternating
    print(f"PASS (m,r)={pair}: degree={degree}, Gal={name}, order={order}")

assert set(EXPECTED) == {
    (m, r)
    for m in range(1, 7)
    for r in range(1, 7)
    if m * r <= 6
}
print("PASS: complete C24 parameter-Galois classification through degree six")

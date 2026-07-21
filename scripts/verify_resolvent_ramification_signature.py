#!/usr/bin/env python3
"""Exact small-case regressions for the resolvent--ramification signature.

The all-(m,r) statements are proved in the accompanying note.  This file only
checks the four displayed cases, including the infinity polynomial and the
distinguished cancellation branch.  It must not be cited as an all-degree
proof.
"""
from __future__ import annotations

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup

from master_cancellation import fiber_antiderivative, parameter_polynomial


T, P, Q, R, q, u, v = sp.symbols("T P Q R q u v")


def infinity_polynomial(m: int, r: int) -> sp.Expr:
    """Return u^(-r-1) integral_0^u v^r(1-v)^(mr) dv."""
    integral = sp.integrate(v**r * (1 - v) ** (m * r), (v, 0, u))
    return sp.cancel(integral / u ** (r + 1))


def star_cycles(m: int, r: int) -> list[Permutation]:
    """A minimal transitive factorization by m+1 cycles of length r+1."""
    cycles = []
    for block in range(m + 1):
        support = [0] + list(
            range(1 + block * r, 1 + (block + 1) * r)
        )
        cycles.append(Permutation(*support, size=r * (m + 1) + 1))
    return cycles


expected_arithmetic_boundary_degrees = {
    (1, 1): [],
    (1, 2): [1],
    (2, 1): [1],
    (2, 2): [1, 2],
}

for m, r in expected_arithmetic_boundary_degrees:
    degree = r * (m + 1) + 1
    anti = fiber_antiderivative(m, r, T, P, Q)
    derivative = sp.factor(sp.diff(anti, T))
    critical = 1 - T * (Q - P * T) ** m
    assert sp.expand(derivative - critical**r) == 0
    assert sp.Poly(anti, T).degree() == degree
    assert sp.Poly(critical, T).degree() == m + 1

    # The projective T=infinity chart over P=0.  Its nonzero roots are
    # simple; one is the finite A=0 cancellation branch.
    K = infinity_polynomial(m, r)
    assert sp.Poly(K, u).degree() == m * r
    assert sp.gcd(sp.Poly(K, u), sp.Poly(sp.diff(K, u), u)).degree() == 0
    modulus = parameter_polynomial(m, r, q)
    distinguished = -q / (1 - q)
    numerator = sp.together(K.subs(u, distinguished)).as_numer_denom()[0]
    assert sp.rem(sp.Poly(numerator, q), sp.Poly(modulus, q)).is_zero

    # The other infinity roots are the other cancellation parameters after
    # the same Mobius change.  Their residue degrees over QQ(q) are therefore
    # the nontrivial orbit sizes of a point stabilizer in Gal(M).  This avoids
    # heuristic numerical factorization in an algebraic number field.
    if sp.Poly(modulus, q).degree() == 1:
        degrees = []
    else:
        parameter_group, _ = sp.polys.numberfields.galois_group(modulus, q)
        stabilizer_orbits = parameter_group.stabilizer(0).orbits()
        degrees = sorted(len(orbit) for orbit in stabilizer_orbits
                         if 0 not in orbit)
    assert degrees == expected_arithmetic_boundary_degrees[(m, r)]

    # Finite branch cycles have length r+1.  The star model is a regression
    # for the group-theoretic hypertree lemma used in the prose proof.
    group = PermutationGroup(star_cycles(m, r))
    expected_order = sp.factorial(degree)
    expected_name = f"S_{degree}"
    if r % 2 == 0:
        expected_order //= 2
        expected_name = f"A_{degree}"
    assert group.order() == expected_order

    print(
        f"PASS (m,r)=({m},{r}): degree={degree}, monodromy={expected_name}, "
        f"critical partition={r}^{m + 1}, P=0 loss={m*r - 1}, "
        f"arithmetic boundary degrees={degrees}"
    )

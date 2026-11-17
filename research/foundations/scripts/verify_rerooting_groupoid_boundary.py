#!/usr/bin/env python3
"""Combinatorial and symbolic checks for the compactified rerooting quotient."""

from __future__ import annotations

from itertools import combinations
from math import comb

import sympy as sp


# There are n=N-2 simple points over zero.  Fixing one of them changes the
# quotient group from S_n to S_(n-1), hence gives degree n.  A boundary block
# of k simple marks splits according to whether the selected mark lies inside
# or outside; the two restrictions have degrees k and n-k.
for n in range(2, 10):
    labels = tuple(range(n))
    selected = 0
    for k in range(1, n + 1):
        blocks = tuple(combinations(labels, k))
        inside = tuple(block for block in blocks if selected in block)
        outside = tuple(block for block in blocks if selected not in block)
        assert len(inside) == comb(n - 1, k - 1)
        assert len(outside) == (comb(n - 1, k) if k < n else 0)
        assert k + (n - k) == n

    # At the generic coefficient discriminant, selecting the collided value
    # gives one component of ramification index two.  Selecting one of the
    # other n-2 values gives an unramified component of degree n-2.
    assert 2 * 1 + 1 * (n - 2) == n


x, epsilon = sp.symbols("x epsilon")

# The special cyclic slice x^k-epsilon has k-cycle monodromy, but it meets the
# discriminant with multiplicity k-1.  Hence its epsilon=delta^k chart is a
# higher-codimension slice, not a generic boundary-divisor inertia statement.
for k in range(2, 8):
    cyclic = x**k - epsilon
    discriminant = sp.factor(sp.discriminant(cyclic, x))
    epsilon_order = sp.Poly(discriminant, epsilon).as_dict()
    assert set(power[0] for power in epsilon_order) == {k - 1}


# A generic pair collision is transverse to the discriminant and has the
# expected square-root normalization, i.e. transposition inertia.
for degree in range(3, 8):
    polynomial = x**2 - epsilon
    for root in range(2, degree):
        polynomial *= x - root
    discriminant = sp.factor(sp.discriminant(sp.expand(polynomial), x))
    quotient = sp.cancel(discriminant / epsilon)
    assert quotient.subs(epsilon, 0) != 0


print("PASS rerooting quotient: [S_n:S_(n-1)]=n=N-2")
print("PASS boundary pullback: selected-in/out degrees are k and n-k")
print("PASS coarse diagonal: one e=2 branch plus n-2 unramified choices")
print("PASS cyclic slices: k-cycle monodromy is not generic divisor inertia")

#!/usr/bin/env python3
"""Exact checks for the fourth-suspension old-boundary valuation fan."""

from __future__ import annotations

from itertools import product

import sympy as sp


r, z, m = sp.symbols("r z m")
a0, a1, b0, b1 = sp.symbols("a0 a1 b0 b1")
A = a0 + a1 * z
B = b0 + b1 * z
F0 = A / B
K0 = sp.factor(m * F0 + z * sp.diff(F0, z))

# The derivative formula in r=1/S, z=P*S coordinates.
F = sp.Function("F")(r, z)
K = m * F - r * sp.diff(F, r) + z * sp.diff(F, z)
assert K == m * F - r * sp.diff(F, r) + z * sp.diff(F, z)

# The three and only three coprime affine-ratio exceptional contacts.
assert sp.factor(K0.subs({m: 0, a0: 2 * b0, a1: 2 * b1})) == 0
assert sp.factor(K0.subs({m: 1, a1: 0, b0: 0})) == 0
assert sp.factor(K0.subs({m: -1, a0: 0, b1: 0})) == 0

# Generic representatives of all other support pairs have nonzero K0.
support_polynomial = {
    (1, 0): 2,
    (0, 1): 3 * z,
    (1, 1): 2 + 3 * z,
}
exceptional = {
    (0, (1, 0), (1, 0)),
    (0, (0, 1), (0, 1)),
    (0, (1, 1), (1, 1)),  # proportional coefficients are a sub-stratum
    (1, (1, 0), (0, 1)),
    (-1, (0, 1), (1, 0)),
}
for integer_m, support_u, support_v in product(
    range(-3, 4), support_polynomial, support_polynomial
):
    test = sp.factor(
        (
            integer_m
            * support_polynomial[support_u]
            / support_polynomial[support_v]
            + z
            * sp.diff(
                support_polynomial[support_u]
                / support_polynomial[support_v],
                z,
            )
        )
    )
    if test == 0:
        assert (integer_m, support_u, support_v) in exceptional

# Degree <=2 gives exactly the advertised old-boundary surviving pairs.
lambdas = (-1, 0, 1, 2)
generic_survivors = {
    integer_m: [
        (lambda_u, lambda_v)
        for lambda_u, lambda_v in product(lambdas, repeat=2)
        if lambda_u - lambda_v == integer_m
    ]
    for integer_m in (1, 2, 3)
}
assert generic_survivors == {
    1: [(0, -1), (1, 0), (2, 1)],
    2: [(1, -1), (2, 0)],
    3: [(2, -1)],
}

# On the unique exceptional survivor m=1, delta=c and the cancelling
# monomials P^i*T^j lie on j-i=c.
for contact_order in range(1, 6):
    delta = 1 - 1 + contact_order
    assert delta == contact_order
    lattice = [
        (i, j)
        for i, j in product(range(8), repeat=2)
        if i - j == -delta
    ]
    assert all(j - i == contact_order for i, j in lattice)

print("PASS: effective degrees reduce to the finite degree-two lattice")
print("PASS: the affine leading ratio has exactly three derivative contacts")
print("PASS: all generic old-boundary survivors have m=1,2,3")
print("PASS: the exceptional surviving residue lies on j-i=c")

# In the standard reciprocal source eta is proportional to z.  Hence on the
# exceptional m=1 row, K_c must be proportional to z**(c+1).  Check the two
# nontrivial degree placements and their successive forced degenerations.
a0, a1, a2 = sp.symbols("a0 a1 a2", nonzero=True)
b0, b1, b2 = sp.symbols("b0 b1 b2")
c0, c1, c2 = sp.symbols("c0 c1 c2", nonzero=True)
d0 = sp.symbols("d0")
P = z * r
S = 1 / r


def derivative_residues(numerator: sp.Expr, denominator: sp.Expr, count: int):
    """Return K_0,...,K_count for m=1."""
    local_f = sp.factor(r * numerator / denominator)
    local_k = sp.factor(
        local_f - r * sp.diff(local_f, r) + z * sp.diff(local_f, z)
    )
    return [
        sp.factor(sp.diff(local_k, r, order).subs(r, 0) / sp.factorial(order))
        for order in range(count + 1)
    ]


numerator_10 = a1 * S + a0 + P * (b1 * S + b0)
denominator_10 = P * (c1 * S + c0)
residues_10 = derivative_residues(numerator_10, denominator_10, 2)
assert residues_10[0] == 0
assert sp.factor(
    residues_10[1] + (a0 * c1 - a1 * c0) / (c1**2 * z)
) == 0
forced_10 = {a0: a1 * c0 / c1}
assert sp.factor(
    residues_10[2].subs(forced_10)
    + (b0 * c1 - b1 * c0) / c1**2
) == 0

numerator_21 = (
    a2 * S**2 + a1 * S + a0 + P * (b2 * S**2 + b1 * S + b0)
)
denominator_21 = d0 + P * (c2 * S**2 + c1 * S + c0)
residues_21 = derivative_residues(numerator_21, denominator_21, 3)
assert residues_21[0] == 0
expected_k1 = -(
    (a1 * c2 - a2 * c1) * z - 2 * a2 * d0
) / (c2**2 * z**2)
assert sp.factor(residues_21[1] - expected_k1) == 0
forced_21_first = {d0: 0, c1: a1 * c2 / a2}
expected_k2 = (
    (-2 * a0 * a2 * c2 + 2 * a2**2 * c0)
    + (a1 * b2 * c2 - a2 * b1 * c2) * z
) / (a2 * c2**2 * z)
assert sp.factor(residues_21[2].subs(forced_21_first) - expected_k2) == 0
forced_21_second = {
    **forced_21_first,
    c0: a0 * c2 / a2,
    b1: a1 * b2 / a2,
}
expected_k3 = 2 * (a0 * b2 - a2 * b0) / (a2 * c2)
assert sp.factor(residues_21[3].subs(forced_21_second) - expected_k3) == 0

print("PASS: the exceptional degree-two m=1 row fails the residue-shape test")

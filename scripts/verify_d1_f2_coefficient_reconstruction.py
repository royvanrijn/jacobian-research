#!/usr/bin/env python3
"""Symbolic function-field checks for the coefficient-space D1/F2 proof."""

from __future__ import annotations

import sympy as sp


W, rho, a = sp.symbols("W rho a", nonzero=True)


def twice_primitive(polynomial: sp.Expr, degree: int) -> sp.Expr:
    expanded = sp.Poly(polynomial, W)
    return sp.expand(
        sum(
            expanded.coeff_monomial(W**i) * W ** (i + 2) / ((i + 1) * (i + 2))
            for i in range(degree + 1)
        )
    )


for N in (4, 5, 6, 8):
    n = N - 2
    coefficients = sp.symbols(f"k{N}_0:{n}")

    # Work over Q(k_0,...,k_(n-1),rho) and solve the root-incidence equation
    # J_K(rho)=0 for the leading coefficient k_n.
    partial = sum(
        coefficients[i] * W**i for i in range(n)
    )
    partial_J = twice_primitive(partial, n - 1)
    leading = sp.cancel(
        -N * (N - 1) * partial_J.subs(W, rho) / rho**N
    )
    K = sp.expand(partial + leading * W**n)
    J = twice_primitive(K, n)

    assert sp.cancel(J.subs(W, rho)) == 0
    assert sp.expand(sp.diff(J, W, 2) - K) == 0
    assert J.subs(W, 0) == 0
    assert sp.diff(J, W).subs(W, 0) == 0

    seed = sp.cancel(-J.subs(W, rho * W) / (rho * sp.diff(J, W).subs(W, rho)))
    assert sp.cancel(seed.subs(W, 0)) == 0
    assert sp.cancel(sp.diff(seed, W).subs(W, 0)) == 0
    assert sp.cancel(seed.subs(W, 1)) == 0
    assert sp.cancel(sp.diff(seed, W).subs(W, 1) + 1) == 0
    assert sp.cancel(
        sp.diff(seed, W, 2)
        + rho * K.subs(W, rho * W) / sp.diff(J, W).subs(W, rho)
    ) == 0

    # Equivariance: K(W)->K(aW), rho->rho/a leaves the reconstructed seed
    # unchanged.  The identity for J is coefficient-independent.
    scaled_K = sp.expand(K.subs(W, a * W))
    scaled_J = twice_primitive(scaled_K, n)
    assert sp.expand(scaled_J - J.subs(W, a * W) / a**2) == 0
    scaled_rho = rho / a
    scaled_seed = sp.cancel(
        -scaled_J.subs(W, scaled_rho * W)
        / (scaled_rho * sp.diff(scaled_J, W).subs(W, scaled_rho))
    )
    assert sp.cancel(scaled_seed - seed) == 0


print("PASS D1/F2: symbolic root-incidence inverse in degrees 4, 5, 6, and 8")
print("PASS D1/F2: reconstruction is invariant under coefficient and coordinate scaling")
print("PASS D1/F2: the marked inverse is scheme-level over the symbolic function field")

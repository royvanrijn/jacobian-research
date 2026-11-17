#!/usr/bin/env python3
"""Exact bounded regression for the constant-term Gaussian--Lagrange lemma.

The all-order proof is in
extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md.  This script deliberately
uses a nonlinear two-variable Phi with two nonzero constant terms.  It checks
the fixed point and compares the circular-Wick derivative sum with

    A(g(u)) / det(I-u JPhi(g(u)))

through an exact finite order.  This is a regression, not a replacement for
the written formal-residue proof.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


u, z1, z2 = sp.symbols("u z1 z2")
CUTOFF = 6


def truncate(expr: sp.Expr) -> sp.Expr:
    """Return the exact u-series through CUTOFF."""
    return sp.expand(sp.series(expr, u, 0, CUTOFF + 1).removeO())


def main() -> None:
    # Both entries have nonzero constant term, so the inverse zero section
    # moves: g_1=u+O(u^2), g_2=-2u+O(u^2).
    phi = sp.Matrix(
        [
            1 + 2 * z1 - z2 + z1 * z2 + z2**2,
            -2 + z1 + 3 * z2 + z1**2 - z1 * z2,
        ]
    )
    A = 1 + 3 * z1 - 2 * z2 + z1**2 + z1 * z2 + 2 * z2**2 + z1**3
    assert phi.subs({z1: 0, z2: 0}) == sp.Matrix([1, -2])

    g = sp.Matrix([sp.Integer(0), sp.Integer(0)])
    for _ in range(CUTOFF + 1):
        substituted = phi.subs({z1: g[0], z2: g[1]}, simultaneous=True)
        g = sp.Matrix([truncate(u * substituted[i]) for i in range(2)])

    fixed_error = g - u * phi.subs({z1: g[0], z2: g[1]}, simultaneous=True)
    assert all(truncate(entry) == 0 for entry in fixed_error)
    assert sp.expand(g[0]).coeff(u, 1) == 1
    assert sp.expand(g[1]).coeff(u, 1) == -2

    # Circular Wick says E(W_1^a W_2^b C(Z)) =
    # partial_1^a partial_2^b C(0).  After expanding exp(u W.Phi),
    # this is exactly the following coefficientwise finite sum.
    wick_series = sp.Integer(0)
    for n in range(CUTOFF + 1):
        coefficient = sp.Integer(0)
        for a in range(n + 1):
            b = n - a
            differentiated = sp.diff(A * phi[0] ** a * phi[1] ** b, z1, a, z2, b)
            coefficient += differentiated.subs({z1: 0, z2: 0}) / (
                factorial(a) * factorial(b)
            )
        wick_series += sp.expand(coefficient) * u**n

    jacobian = phi.jacobian([z1, z2])  # rows: Phi components; columns: z variables
    jacobian_at_g = jacobian.subs({z1: g[0], z2: g[1]}, simultaneous=True)
    denominator = sp.det(sp.eye(2) - u * jacobian_at_g)
    numerator = A.subs({z1: g[0], z2: g[1]}, simultaneous=True)
    lagrange_series = truncate(numerator / denominator)

    assert truncate(wick_series - lagrange_series) == 0
    assert truncate(denominator).subs(u, 0) == 1
    print(
        "PASS formal Gaussian--Lagrange regression: nonzero constant-term case, "
        f"exact through u^{CUTOFF}"
    )


if __name__ == "__main__":
    main()

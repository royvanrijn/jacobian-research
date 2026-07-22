#!/usr/bin/env python3
"""Exact regression for recovery of h from weighted-bridge mixed moments.

For g=u*h(g), the bridge gives

    [u^m]g = E(Q P^m)/m!.

This script constructs those coefficients by univariate Lagrange inversion,
reverts g without using the identity being tested, and verifies h=z/eta.
The written argument in WEIGHTED_GAUSSIAN_BRIDGE.md is the all-order proof;
these finite symbolic and rational calculations are regressions.
"""

from __future__ import annotations

import sympy as sp


u, z = sp.symbols("u z")


def truncate(expr: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, variable, 0, order).removeO())


def mixed_generating_series(h: sp.Expr, degree: int) -> sp.Expr:
    """Return g through u^(degree+1) from the exact mixed-moment formula."""
    result = sp.Integer(0)
    for m in range(1, degree + 2):
        coefficient = sp.expand(h**m).coeff(z, m - 1) / sp.Integer(m)
        result += coefficient * u**m
    return sp.expand(result)


def compositional_inverse(g: sp.Expr, order: int) -> sp.Expr:
    """Revert g=u+O(u^2) through z^(order-1), coefficient by coefficient."""
    assert sp.expand(g).coeff(u, 1) == 1
    eta = z
    for n in range(2, order):
        trial_coefficient = sp.symbols(f"c{n}")
        trial = eta + trial_coefficient * z**n
        composed = truncate(g.subs(u, trial), z, n + 1)
        equation = sp.expand(composed - z).coeff(z, n)
        solution = sp.solve(sp.Eq(equation, 0), trial_coefficient)
        assert len(solution) == 1
        eta = sp.expand(trial.subs(trial_coefficient, solution[0]))
    return eta


def recover(h: sp.Expr, degree: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    g = mixed_generating_series(h, degree)
    eta = compositional_inverse(g, degree + 2)
    recovered = truncate(z / eta, z, degree + 1)
    return g, eta, recovered


def main() -> None:
    # A universal degree-four audit: M_1,...,M_5 recover every coefficient.
    a1, a2, a3, a4 = sp.symbols("a1 a2 a3 a4")
    symbolic_h = 1 + a1 * z + a2 * z**2 + a3 * z**3 + a4 * z**4
    _, symbolic_eta, symbolic_recovered = recover(symbolic_h, 4)
    assert sp.expand(symbolic_recovered - symbolic_h) == 0
    assert truncate(symbolic_eta * symbolic_h - z, z, 6) == 0
    print("PASS Gaussian fingerprint: symbolic degree four recovered from M_1,...,M_5")

    # A concrete degree-five weighted seed with rational scale.
    seed = z**4 * (1 - z)
    normalized_seed = sp.expand(-seed / sp.diff(seed, z).subs(z, 1))
    scale = sp.Rational(2, 3)
    assert sp.diff(normalized_seed, z).subs(z, 1) == -1
    concrete_h = 1 + scale * normalized_seed
    concrete_g, concrete_eta, concrete_recovered = recover(concrete_h, 5)
    assert sp.expand(concrete_recovered - concrete_h) == 0
    assert truncate(concrete_g.subs(u, concrete_eta) - z, z, 7) == 0
    recovered_scale = -sp.diff(concrete_recovered, z).subs(z, 1)
    recovered_seed = sp.cancel((concrete_recovered - 1) / recovered_scale)
    assert recovered_scale == scale
    assert sp.expand(recovered_seed - normalized_seed) == 0
    print("PASS Gaussian fingerprint: weighted degree five recovered from M_1,...,M_6")


if __name__ == "__main__":
    main()

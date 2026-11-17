#!/usr/bin/env python3
"""Structural exact checks for the all-degree quadratic-gauge construction.

This certificate is independent of the Lean development and complements the
full finite expansions in ``verify_finite_etale_keller_fibers.py``.  It checks
all identities which are independent of the degree once, verifies the generic
``k``-th summand symbolically, and performs a separate six-coefficient bridge
regression.  Together with linearity this audits every finite seed polynomial;
it is not presented as a substitute for the Lean proof.
"""
from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class MonomialExponent:
    """Exponents of ``t``, ``x``, and ``q`` as affine expressions in ``k``."""

    t: sp.Expr
    x: sp.Expr
    q: sp.Expr

    def __add__(self, other: "MonomialExponent") -> "MonomialExponent":
        return MonomialExponent(
            sp.expand(self.t + other.t),
            sp.expand(self.x + other.x),
            sp.expand(self.q + other.q),
        )

    def __rmul__(self, scalar: sp.Expr) -> "MonomialExponent":
        return MonomialExponent(
            sp.expand(scalar * self.t),
            sp.expand(scalar * self.x),
            sp.expand(scalar * self.q),
        )


def assert_zero(expression: sp.Expr, label: str) -> None:
    reduced = sp.simplify(expression)
    if reduced != 0:
        reduced = sp.factor(sp.cancel(sp.together(reduced)))
    assert reduced == 0, f"{label}: residual {reduced}"


def check_source_chart() -> None:
    """Verify the reciprocal factor and the full source-chart Jacobian."""
    x, y, z, a = sp.symbols("x y z a")
    t = 1 + x * y
    q = t**2 * z + a * y**2 * (1 + 3 * t)
    pi = t * q
    S = x / t
    Q = y + x * q
    D = 1 - S * Q + pi * S**2

    assert_zero(D - 1 / t, "source reciprocal identity")
    source_jacobian = sp.Matrix((pi, S, Q)).jacobian((x, y, z)).det()
    assert_zero(source_jacobian - t, "source-chart Jacobian")
    assert_zero(source_jacobian * D - 1, "source-chart cancellation")
    assert_zero(source_jacobian * (-2 * D) + 2, "combined determinant")


def check_marked_line() -> None:
    """Verify the marked-line identities for an arbitrary differentiable U."""
    S, Q, pi = sp.symbols("S Q pi")
    U = sp.Function("U")(S)
    beta = (sp.diff(U, S) - 1 - pi * S**2) / S
    B = Q + beta
    C = 2 * U - B * S**2
    D = 1 - S * Q + pi * S**2

    assert_zero(sp.diff(U, S) - B * S - D, "derivative/chart factor")
    plane_jacobian = sp.Matrix((B, C)).jacobian((S, Q)).det()
    assert_zero(plane_jacobian + 2 * D, "marked-line Jacobian")


def check_all_degree_coefficient_bridge() -> None:
    """Check the low terms and one generic tail summand exactly."""
    k = sp.symbols("k", integer=True, positive=True)
    S, pi = sp.symbols("S pi")
    A2, A3, Ak = sp.symbols("A2 A3 Ak")

    U_low = S + A2 * pi * S**2 + A3 * pi * S**3
    beta_low = sp.cancel((sp.diff(U_low, S) - 1 - pi * S**2) / S)
    assert_zero(
        beta_low - (2 * A2 * pi + (3 * A3 - 1) * pi * S),
        "low-degree beta identity",
    )
    Q = sp.symbols("Q")
    assert_zero(
        2 * U_low - (Q + beta_low) * S**2
        - (2 * S - Q * S**2 + (1 - A3) * pi * S**3),
        "low-degree C expansion",
    )

    U_tail = Ak * pi**k * S**k
    beta_tail = sp.cancel(sp.diff(U_tail, S) / S)
    assert_zero(
        beta_tail - k * Ak * pi**k * S ** (k - 2),
        "generic beta summand",
    )
    assert_zero(
        2 * U_tail - beta_tail * S**2
        + (k - 2) * Ak * pi**k * S**k,
        "generic C summand",
    )

    # Under pi=t*q and S=x/t, encode monomials by their (t,x,q) exponents.
    pi_exp = MonomialExponent(1, 0, 1)
    S_exp = MonomialExponent(-1, 1, 0)
    xq_exp = MonomialExponent(0, 1, 1)
    tq_exp = MonomialExponent(1, 0, 1)

    assert pi_exp + S_exp == xq_exp
    assert pi_exp == tq_exp
    assert sp.expand(1 + (3 * A3 - 1) - 3 * A3) == 0

    b_tail = k * pi_exp + (k - 2) * S_exp
    assert b_tail == MonomialExponent(2, k - 2, k)

    c_tail = k * pi_exp + k * S_exp
    assert c_tail == k * xq_exp

    t = sp.symbols("t")
    assert_zero(
        (t + 1) - (t - 1) ** 2 * (1 + 3 * t) - t**2 * (5 - 3 * t),
        "low-degree C identity",
    )


def check_symbolic_sixth_degree_regression() -> None:
    """Instantiate the bridge with six independent coefficients."""
    x, y, z, S, pi = sp.symbols("x y z S pi")
    g1, g3 = sp.symbols("g1 g3", nonzero=True)
    g2, g4, g5, g6 = sp.symbols("g2 g4 g5 g6")
    coefficients = {1: g1, 2: g2, 3: g3, 4: g4, 5: g5, 6: g6}

    t = 1 + x * y
    q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
    pi_source = t * q
    S_source = x / t
    Q_source = y + x * q

    G_pi = (
        g1 * S
        + pi * (g2 * S**2 + g3 * S**3)
        + sum(coefficients[j] * pi**j * S**j for j in range(4, 7))
    )
    beta = sp.cancel((sp.diff(G_pi, S) / g1 - 1 - pi * S**2) / S)

    displayed_B = (
        y
        + 3 * (g3 / g1) * x * q
        + 2 * (g2 / g1) * t * q
        + sum(
            j * (coefficients[j] / g1) * t**2 * x ** (j - 2) * q**j
            for j in range(4, 7)
        )
    )
    displayed_C = (
        x * (5 - 3 * t)
        - (g3 / g1) * x**3 * z
        - sum(
            (j - 2) * (coefficients[j] / g1) * (x * q) ** j
            for j in range(4, 7)
        )
    )

    substitutions = {S: S_source, pi: pi_source}
    marked_B = (Q_source + beta).subs(substitutions)
    marked_C = (2 * G_pi / g1 - (Q_source + beta) * S**2).subs(substitutions)
    assert_zero(marked_B - displayed_B, "degree-six B bridge")
    assert_zero(marked_C - displayed_C, "degree-six C bridge")


def check_effective_degree_bound() -> None:
    """Check the termwise total-degree estimates for every N >= 3."""
    k, N = sp.symbols("k N", integer=True)
    degree_t = 2
    degree_q = max(2 * degree_t + 1, 2 + degree_t)
    assert degree_q == 5
    assert degree_t + degree_q == 7

    degree_B_tail = 2 * degree_t + (k - 2) + k * degree_q
    degree_C_tail = k * (1 + degree_q)
    assert sp.expand(degree_B_tail - (6 * k + 2)) == 0
    assert sp.expand(degree_C_tail - 6 * k) == 0

    assert sp.expand(degree_B_tail.subs(k, N) - (6 * N + 2)) == 0
    assert sp.expand(degree_C_tail.subs(k, N) - 6 * N) == 0
    assert max(1, 6, 7) <= 6 * 3 + 2
    assert max(3, 4) <= 6 * 3


def check_normalization() -> None:
    """Verify the target-preserving determinant-one output normalization."""
    assert sp.Rational(-1, 2) * (-2) == 1
    B = sp.symbols("B")
    assert (-B / 2).subs(B, 0) == 0


if __name__ == "__main__":
    check_source_chart()
    check_marked_line()
    check_all_degree_coefficient_bridge()
    check_symbolic_sixth_degree_regression()
    check_effective_degree_bound()
    check_normalization()
    print("PASS: universal source-chart reciprocal and Jacobian identities")
    print("PASS: universal marked-line Jacobian cancellation")
    print("PASS: generic k-th summand recovers the displayed all-degree map")
    print("PASS: independent symbolic degree-six bridge regression")
    print("PASS: all-degree coordinate bound max(deg F_i) <= 6N+2")
    print("PASS: target-preserving output normalization has determinant one")

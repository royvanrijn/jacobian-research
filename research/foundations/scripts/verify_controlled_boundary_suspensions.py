#!/usr/bin/env python3
"""Exact checks for the boundary-cancelled incidence lemma and its examples.

The symbolic identities certify the three family ledgers and the stated
propositions.  The finite
exponent box is only a regression against transcription errors; Proposition
4.1 itself is proved by the three coefficient equations in the note.
"""
from __future__ import annotations

import sympy as sp


def check_weighted_ledger() -> None:
    x, gamma, b0, c = sp.symbols("x gamma b0 c", nonzero=True)
    C = x * gamma
    jac_alpha = b0 * x**3 * gamma**2
    jac_core = -c**2 * gamma
    jac_beta = -c * C**3
    jac_F = sp.cancel(jac_core * jac_alpha / jac_beta)
    assert jac_F == b0 * c


def check_cancellation_ledger() -> None:
    A, C = sp.symbols("A C", nonzero=True)
    r = sp.symbols("r", positive=True, integer=True)
    D_on_source = A**-1
    jac_alpha = -A**r
    jac_core = C * D_on_source**r
    assert sp.cancel(jac_alpha * jac_core) == -C


def check_quadratic_gauge_ledger() -> None:
    P, S, Q = sp.symbols("P S Q")
    D = 1 - S * Q + P * S**2
    # The foundational normalized seed G_P=S+P*S^3 already contains the
    # universal marked-line ledger; coefficient decorations do not change it.
    Y = 2 * (S + P * S**3)
    B = Q + 2 * P * S
    C = Y - B * S**2
    jac_core = sp.factor(
        sp.det(sp.Matrix([P, B, C]).jacobian((P, S, Q)))
    )
    assert sp.factor(jac_core + 2 * D) == 0
    assert sp.cancel(jac_core / D) == -2


def check_simple_section_normal_form() -> None:
    w, q, u = sp.symbols("w q u", nonzero=True)
    for degree in range(1, 7):
        coefficients = sp.symbols(f"h0:{degree + 1}")
        h = sum(coefficients[j] * w**j for j in range(degree + 1))
        H = sp.integrate(h, w)
        T = u * (w * q - H)
        jacobian = sp.det(sp.Matrix([q, T]).jacobian([w, q]))
        assert sp.expand(jacobian - u * (h - q)) == 0


def check_critical_normalizations() -> None:
    w, q, P, Y = sp.symbols("w q P Y", nonzero=True)
    h = w**4 + 2 * w + 3
    weighted_D = q - h
    assert sp.expand(weighted_D.subs(q, h)) == 0

    for m in range(1, 7):
        s = Y**-m
        Q = Y + P * Y**-m
        cancellation_D = 1 - s * (Q - P * s) ** m
        assert sp.cancel(cancellation_D) == 0


def check_two_boundary_jacobian_formula() -> None:
    x, y, z = sp.symbols("x y z")
    f1 = y**2 + y + 1
    f2 = y**3-y+2
    A1 = 1 + x * f1
    A2 = 1 + x * f2

    samples = (
        (1, 2, 2, 1, -1, 0),
        (2, 1, 1, 3, 0, -1),
        (1, 1, 2, 2, -1, -1),
        (3, 2, 1, 1, 1, -2),
    )
    for a1, a2, b1, b2, d1, d2 in samples:
        g = x * y + y**2
        B = A1**b1 * A2**b2 * z + g
        P = A1**a1 * A2**a2 * B
        s = x * A1**d1 * A2**d2
        Q = y + s * P
        direct = sp.factor(sp.det(sp.Matrix([s, P, Q]).jacobian([x, y, z])))
        extra = (1 + d1 + d2) * A1 * A2 - d1 * A2 - d2 * A1
        formula = -A1 ** (a1 + b1 + d1 - 1) * A2 ** (
            a2 + b2 + d2 - 1
        ) * extra
        assert sp.cancel(direct - formula) == 0


def check_third_divisor_classification() -> None:
    # N has coefficients (1+d1+d2, -d2, -d1) on
    # (A1*A2, A1, A2).  It is a monomial iff exactly one is nonzero.
    exceptional: set[tuple[int, int]] = set()
    for d1 in range(-12, 13):
        for d2 in range(-12, 13):
            coefficients = (1 + d1 + d2, -d2, -d1)
            if sum(coefficient != 0 for coefficient in coefficients) == 1:
                exceptional.add((d1, d2))
    assert exceptional == {(0, 0), (-1, 0), (0, -1)}

    # The coefficient equations also prove exhaustion over all integers.
    d1, d2 = sp.symbols("d1 d2")
    cases = [
        sp.solve((sp.Eq(d1, 0), sp.Eq(d2, 0)), (d1, d2), dict=True),
        sp.solve((sp.Eq(1 + d1 + d2, 0), sp.Eq(d2, 0)), (d1, d2), dict=True),
        sp.solve((sp.Eq(1 + d1 + d2, 0), sp.Eq(d1, 0)), (d1, d2), dict=True),
    ]
    assert cases == [[{d1: 0, d2: 0}], [{d1: -1, d2: 0}], [{d1: 0, d2: -1}]]


def main() -> None:
    check_weighted_ledger()
    check_cancellation_ledger()
    check_quadratic_gauge_ledger()
    check_simple_section_normal_form()
    check_critical_normalizations()
    check_two_boundary_jacobian_formula()
    check_third_divisor_classification()
    print(
        "PASS boundary-cancelled incidence ledgers for weighted, "
        "cancellation, and quadratic-gauge families"
    )
    print("PASS controlled-boundary plane forms and two-boundary obstruction")


if __name__ == "__main__":
    main()

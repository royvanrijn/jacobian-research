#!/usr/bin/env python3
"""Verify the exact implicit quintic and conductor pullback on F2 k=1."""

from __future__ import annotations

import sympy as sp


t, u, P, Q = sp.symbols("t u P Q")
a, b, c, d = sp.symbols("a b c d")


def expected_implicit_quintic() -> sp.Expr:
    return sp.expand(
        P**5
        + (a * b + b**3 + 3 * c) * P**4
        + 3 * b * P**3 * Q
        + (-a * b * c + 4 * a * d + 3 * b**2 * d + 3 * c**2) * P**3
        + (-5 * a**2 - 4 * a * b**2 + 3 * b * c + 3 * d) * P**2 * Q
        + (
            a**3 * c
            + a**2 * b**2 * c
            - a**2 * b * d
            - 2 * a * b * c**2
            + 5 * a * c * d
            + 3 * b * d**2
            + c**3
        )
        * P**2
        - 5 * a * P * Q**2
        + (a**3 * b - 3 * a**2 * c - 5 * a * b * d + 3 * c * d) * P * Q
        + (
            a**4 * d
            + a**3 * b**2 * d
            - 2 * a**2 * b * c * d
            + 2 * a**2 * d**2
            + a * c**2 * d
            + d**3
        )
        * P
        - Q**3
        + (2 * a**2 * b - 2 * a * c) * Q**2
        + (
            -a**5
            - a**4 * b**2
            + 2 * a**3 * b * c
            - 2 * a**3 * d
            - a**2 * c**2
            - a * d**2
        )
        * Q
    )


def implicitization_audit() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    implicit = sp.expand(-sp.resultant(P - p, Q - q, t))
    assert sp.expand(implicit - expected_implicit_quintic()) == 0
    assert sp.expand(implicit.subs({P: p, Q: q})) == 0

    polynomial = sp.Poly(implicit, P, Q)
    assert polynomial.total_degree() == 5
    assert polynomial.degree(P) == 5
    assert polynomial.degree(Q) == 3
    assert {
        monomial
        for monomial, coefficient in polynomial.terms()
        if coefficient != 0
    } == {
        (5, 0),
        (4, 0),
        (3, 1),
        (3, 0),
        (2, 1),
        (2, 0),
        (1, 2),
        (1, 1),
        (1, 0),
        (0, 3),
        (0, 2),
        (0, 1),
    }
    top = sum(
        coefficient * P**monomial[0] * Q**monomial[1]
        for monomial, coefficient in polynomial.terms()
        if sum(monomial) == 5
    )
    assert top == P**5
    return p, q, implicit


def conductor_gradient_audit(p: sp.Expr, q: sp.Expr, implicit: sp.Expr) -> sp.Expr:
    collision_quartic = (
        u**4
        + b * u**3
        + a * u**2
        + (2 * a * b - c) * u
        - (a**2 + d)
    )
    pair_polynomial = t**2 - u * t + (u**2 + a)
    conductor = sp.expand(sp.resultant(collision_quartic, pair_polynomial, u))
    assert sp.Poly(conductor, t).degree() == 8

    # Compare the expanded identities directly.  Asking a multivariate
    # factorizer to rediscover the displayed factors is substantially slower
    # and nondeterministic across SymPy runs.
    gradient_p = sp.expand(sp.diff(implicit, P).subs({P: p, Q: q}))
    gradient_q = sp.expand(sp.diff(implicit, Q).subs({P: p, Q: q}))
    assert sp.expand(gradient_p - sp.diff(q, t) * conductor) == 0
    assert sp.expand(gradient_q + sp.diff(p, t) * conductor) == 0
    assert sp.expand(gradient_p * sp.diff(p, t) + gradient_q * sp.diff(q, t)) == 0
    return conductor


def generic_witness_audit(conductor: sp.Expr) -> None:
    witness = {a: 1, b: 0, c: 0, d: 0}
    witness_conductor = sp.Poly(conductor.subs(witness), t)
    assert witness_conductor.as_expr() == (
        t**8 + 3 * t**6 + 4 * t**4 + 2 * t**2 + 1
    )
    assert sp.discriminant(witness_conductor.as_expr(), t) == 4_000_000
    p_witness = t**3 + t
    q_witness = t**5
    assert sp.gcd(sp.diff(p_witness, t), sp.diff(q_witness, t)) == 1


def main() -> None:
    p, q, implicit = implicitization_audit()
    conductor = conductor_gradient_audit(p, q, implicit)
    generic_witness_audit(conductor)
    print(
        "PASS: the F2 k=1 normal form has an exact 12-support implicit "
        "quintic; its pulled-back gradient is the degree-eight nodal "
        "conductor times (q',-p')"
    )


if __name__ == "__main__":
    main()

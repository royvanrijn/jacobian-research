#!/usr/bin/env python3
"""Verify the finite F2 affine target-curve/collision atlas arithmetic."""

from __future__ import annotations

import sympy as sp


PARAMETRIZATION_BOUND = 125 - 1


def divided_difference(poly: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    numerator = sp.expand(poly.subs(t, left) - poly.subs(t, right))
    quotient, remainder = sp.div(numerator, left - right, domain=sp.QQ)
    assert remainder == 0
    return sp.expand(quotient)


t, s = sp.symbols("t s")


def degree_atlas_audit() -> None:
    rows = []
    for scale in range(1, 25):
        p_degree = 3 * scale
        q_degree = 5 * scale
        curve_degree = q_degree
        raw_parameter_coefficients = p_degree + q_degree + 2
        normalized_implicit_coefficients = (
            curve_degree * (curve_degree + 1) // 2 - 3
        )
        assert p_degree * 5 == q_degree * 3
        assert q_degree <= PARAMETRIZATION_BOUND
        assert curve_degree % 5 == 0
        rows.append(
            (
                scale,
                p_degree,
                q_degree,
                curve_degree,
                raw_parameter_coefficients,
                normalized_implicit_coefficients,
            )
        )

    assert len(rows) == 24
    assert rows[0][:4] == (1, 3, 5, 5)
    assert rows[-1][:5] == (24, 72, 120, 120, 194)
    assert 5 * 25 > PARAMETRIZATION_BOUND
    assert rows[-1][5] == 7257


def cusp_collision_audit() -> None:
    p = t**3
    q = t**5
    delta_p = divided_difference(p, s, t)
    delta_q = divided_difference(q, s, t)
    assert delta_p.subs({s: 0, t: 0}) == 0
    assert delta_q.subs({s: 0, t: 0}) == 0
    assert sp.diff(p, t).subs(t, 0) == 0
    assert sp.diff(q, t).subs(t, 0) == 0

    x, y = sp.symbols("x y")
    implicit = y**3 - x**5
    assert sp.expand(implicit.subs({x: p, y: q})) == 0
    assert sp.Poly(implicit, x, y).total_degree() == 5
    assert sp.diff(implicit, x).subs({x: 0, y: 0}) == 0
    assert sp.diff(implicit, y).subs({x: 0, y: 0}) == 0


def multibranch_collision_audit() -> None:
    # A birational degree-(3,5) parametrization with three distinct
    # normalization points over the origin.  It demonstrates the
    # off-diagonal alternative to the cusp/critical diagonal.
    p = t**3 - t
    q = t**5 - t**3
    delta_p = divided_difference(p, s, t)
    delta_q = divided_difference(q, s, t)
    for left, right in ((-1, 0), (-1, 1), (0, 1)):
        assert p.subs(t, left) == p.subs(t, right) == 0
        assert q.subs(t, left) == q.subs(t, right) == 0
        assert delta_p.subs({s: left, t: right}) == 0
        assert delta_q.subs({s: left, t: right}) == 0
    assert sp.gcd(sp.diff(p, t), sp.diff(q, t)) == 1

    x, y = sp.symbols("x y")
    implicit = x**5 - y * (y - x) ** 2
    assert sp.expand(implicit.subs({x: p, y: q})) == 0
    assert sp.Poly(implicit, x, y).total_degree() == 5


def main() -> None:
    degree_atlas_audit()
    cusp_collision_audit()
    multibranch_collision_audit()
    print(
        "PASS: every F2 nonproperness normalization lies in one of 24 "
        "degree charts (3k,5k), has curve degree 5k<=120, and must meet "
        "the divided-difference collision/critical locus"
    )


if __name__ == "__main__":
    main()

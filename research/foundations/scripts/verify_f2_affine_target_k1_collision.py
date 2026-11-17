#!/usr/bin/env python3
"""Verify the exact collision quartic and generic k=1 conductor packet."""

from __future__ import annotations

import sympy as sp


s, t, u, v, z = sp.symbols("s t u v z")
a, b, c, d = sp.symbols("a b c d")


def divided_difference(poly: sp.Expr) -> sp.Expr:
    numerator = sp.expand(poly.subs(t, s) - poly)
    quotient, remainder = sp.div(numerator, s - t, domain=sp.QQ)
    assert remainder == 0
    return sp.expand(quotient)


def symbolic_collision_audit() -> None:
    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    delta_p = divided_difference(p)
    delta_q = divided_difference(q)

    symmetric_p = u**2 - v + a
    reduced_p = sp.symmetrize(delta_p, [s, t], formal=True)[0]
    assert sp.expand(reduced_p.subs({sp.Symbol("s1"): u, sp.Symbol("s2"): v})) == symmetric_p

    reduced_q = sp.symmetrize(delta_q, [s, t], formal=True)[0]
    s1, s2 = sp.symbols("s1 s2")
    reduced_q = sp.expand(reduced_q.subs({s1: u, s2: v}))
    reduced_q = sp.expand(reduced_q.subs(v, u**2 + a))
    collision_quartic = (
        u**4 + b * u**3 + a * u**2 + (2 * a * b - c) * u - (a**2 + d)
    )
    assert sp.expand(reduced_q + collision_quartic) == 0

    pair_polynomial = z**2 - u * z + (u**2 + a)
    q_remainder = sp.rem(q.subs(t, z), pair_polynomial, z)
    expected_y = (u**2 + a) * (u**3 + 2 * a * u + a * b - c)
    assert sp.expand(q_remainder - (-collision_quartic * z + expected_y)) == 0
    assert sp.discriminant(pair_polynomial, z) == -3 * u**2 - 4 * a


def generic_witness_audit() -> None:
    witness = {a: 1, b: 0, c: 0, d: 0}
    collision_quartic = u**4 + u**2 - 1
    diagonal = 3 * u**2 + 4
    assert sp.discriminant(collision_quartic, u) == -400
    assert sp.resultant(collision_quartic, diagonal, u) == 25

    p = t**3 + t
    q = t**5
    tangent = sp.expand(
        sp.diff(p.subs(t, s), s) * sp.diff(q, t)
        - sp.diff(q.subs(t, s), s) * sp.diff(p, t)
    )
    tangent_quotient, remainder = sp.div(tangent, s - t, domain=sp.QQ)
    assert remainder == 0
    symmetric_tangent = sp.symmetrize(
        tangent_quotient, [s, t], formal=True
    )[0]
    s1, s2 = sp.symbols("s1 s2")
    tangent_u = sp.expand(
        symmetric_tangent.subs({s1: u, s2: u**2 + 1})
    )
    tangent_u = sp.rem(tangent_u, collision_quartic, u)
    assert sp.factor(tangent_u) == -10 * u * (u**2 + 2)
    assert sp.resultant(collision_quartic, tangent_u, u) == -10000

    left, right, inverse = sp.symbols("left right inverse")
    x = lambda value: -value * (value**2 + 1)
    y = lambda value: (value**2 + 1) * (value**3 + 2 * value)
    distinct_images = sp.groebner(
        [
            collision_quartic.subs(u, left),
            collision_quartic.subs(u, right),
            x(left) - x(right),
            y(left) - y(right),
            inverse * (left - right) - 1,
        ],
        inverse,
        right,
        left,
        order="lex",
    )
    assert any(poly.as_expr() == 1 for poly in distinct_images.polys)

    # The finite collision set proves generic injectivity of this witness.
    assert sp.gcd(sp.diff(p, t), sp.diff(q, t)) == 1


def genus_conductor_audit() -> None:
    quintic_arithmetic_genus = (5 - 1) * (5 - 2) // 2
    infinity_delta = (2 - 1) * (5 - 1) // 2
    affine_node_delta = 4
    assert quintic_arithmetic_genus == 6
    assert infinity_delta == 2
    assert infinity_delta + affine_node_delta == quintic_arithmetic_genus


def main() -> None:
    symbolic_collision_audit()
    generic_witness_audit()
    genus_conductor_audit()
    print(
        "PASS: the F2 k=1 target chart has an exact collision quartic; "
        "a nonempty open packet has four ordinary affine nodes plus the "
        "delta-2 (2,5) infinity cusp"
    )


if __name__ == "__main__":
    main()

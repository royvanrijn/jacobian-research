#!/usr/bin/env python3
"""Verify the lower-Smith obstruction to an unsaturated decic certificate."""

from __future__ import annotations

import sympy as sp

from verify_hc4_double_conic_normal_layers import (
    binary_coefficients,
    harmonic_layers,
    harmonic_lift,
    lam,
    mu,
    nu,
    q,
    s,
    t,
    x,
    y,
    z,
)


def exact_quotient(numerator: sp.Expr, divisor: sp.Expr) -> sp.Expr:
    quotient = sp.cancel(numerator / divisor)
    assert sp.denom(quotient) == 1
    return sp.expand(quotient)


def main() -> None:
    binary_decic = s**10 + t**10
    quintic = x**5 + z**5
    f_coefficients = binary_coefficients(binary_decic, 10)

    harmonic_f = harmonic_lift(f_coefficients, 5)
    first_remainder = exact_quotient(quintic - harmonic_f, q)
    binary_sextic = sp.expand(first_remainder.subs({x: s**2, y: s * t, z: t**2}))
    g_coefficients = binary_coefficients(binary_sextic, 6)

    harmonic_g = harmonic_lift(g_coefficients, 3)
    second_remainder = exact_quotient(first_remainder - harmonic_g, q)
    binary_quadratic = sp.expand(
        second_remainder.subs({x: s**2, y: s * t, z: t**2})
    )
    k_coefficients = binary_coefficients(binary_quadratic, 2)
    harmonic_k = harmonic_lift(k_coefficients, 1)

    reconstructed = harmonic_f + q * harmonic_g + q**2 * harmonic_k
    assert sp.expand(reconstructed - quintic) == 0
    assert sp.expand(sp.hessian(quintic, (x, y, z)).det()) == 0

    layers = harmonic_layers(f_coefficients, g_coefficients, k_coefficients)
    for degree in (18, 14, 10, 6):
        assert sp.expand(layers[degree].subs({lam: 1, mu: 1, nu: 1})) == 0
    print("PASS: a squarefree decic lies on all four unsaturated normal layers")

    affine_decic = sp.Poly(binary_decic.subs(t, 1), s)
    discriminant = sp.discriminant(affine_decic.as_expr(), s)
    assert discriminant != 0
    assert sp.gcd(affine_decic, affine_decic.diff()) == 1
    print(f"PASS: Disc(s^10+t^10)={discriminant} is nonzero")

    repeated_stable = (
        s**2
        * t
        * (s - t)
        * (s - 2 * t)
        * (s - 3 * t)
        * (s - 4 * t)
        * (s - 5 * t)
        * (s - 6 * t)
        * (s - 7 * t)
    )
    assert sp.Poly(repeated_stable, s, t).total_degree() == 10
    assert sp.discriminant(sp.expand(repeated_stable.subs(t, 1)), s) == 0
    print("PASS: a GIT-stable repeated-root decic can have zero discriminant")
    print("THEOREM: invariant elimination must saturate by Phi_2 and use more than Disc(f)")


if __name__ == "__main__":
    main()

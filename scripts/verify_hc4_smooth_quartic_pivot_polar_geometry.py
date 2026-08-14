#!/usr/bin/env python3
"""Verify the polar-resultant identity for the HC4NHM16 pivot."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    s, t = sp.symbols("s t")
    p, q, r, tau, epsilon = sp.symbols("p q r tau epsilon")

    cubic = s**3 + t**3
    quadratic = p * s**2 - q * s * t / 3 + r * t**2
    resultant = sp.factor(sp.resultant(cubic.subs(t, 1), quadratic.subs(t, 1), s))
    hesse_cubic = p**3 + q**3 / 27 + r**3 - p * q * r
    assert sp.expand(resultant - hesse_cubic) == 0

    first = tau**5 + 6 * tau**2
    middle = tau**4 + tau
    last = 6 * tau**3 + 1
    direction = first * s**2 + 3 * middle * s * t + last * t**2

    varied_resultant = hesse_cubic.subs(
        {
            p: p + epsilon * first,
            q: q - 9 * epsilon * middle,
            r: r + epsilon * last,
        },
        simultaneous=True,
    )
    polar = sp.Poly(sp.expand(varied_resultant), epsilon).coeff_monomial(epsilon)
    pivot = (
        (3 * p**2 - q * r) * tau**5
        + (9 * p * r - q**2) * tau**4
        + (18 * r**2 - 6 * p * q) * tau**3
        + (18 * p**2 - 6 * q * r) * tau**2
        + (9 * p * r - q**2) * tau
        + (3 * r**2 - p * q)
    )
    assert sp.expand(polar - pivot) == 0
    print("PASS: Delta is the first polar of Res(s^3+t^3,H) in direction K_tau")

    direction_resultant = sp.factor(
        sp.resultant(cubic.subs(t, 1), direction.subs(t, 1), s)
    )
    expected = (
        (tau + 1)
        * (tau**2 - tau + 1)
        * (tau**4 - 4 * tau**3 + 10 * tau**2 - 4 * tau + 1)
        * (
            tau**8
            + 4 * tau**7
            + 6 * tau**6
            + 32 * tau**5
            + 83 * tau**4
            + 32 * tau**3
            + 6 * tau**2
            + 4 * tau
            + 1
        )
    )
    assert sp.expand(direction_resultant - expected) == 0
    assert sp.gcd(expected, sp.diff(expected, tau)) == 1
    print("PASS: the 15 finite degeneration parameters are squarefree and explicit")

    h1, h2, h3, k1, k2, k3 = sp.symbols("h1 h2 h3 k1 k2 k3")
    root_polar = k1 * h2 * h3 + k2 * h1 * h3 + k3 * h1 * h2
    conic_hessian = sp.hessian(root_polar, (h1, h2, h3))
    assert sp.expand(conic_hessian.det() - 2 * k1 * k2 * k3) == 0
    assert sp.factor(root_polar.subs(k1, 0)) == h1 * (h2 * k3 + h3 * k2)
    print("PASS: the generic polar is a smooth conic and each degeneration is two lines")
    print("THEOREM: the visible smooth-quartic pivot divisor is a polar-resultant divisor")


if __name__ == "__main__":
    main()

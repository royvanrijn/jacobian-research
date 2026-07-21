"""Exact regressions for all-parameter cancellation construction thick boundary intersections."""

from __future__ import annotations

import sympy as sp


Y, Q, P, R, C, u = sp.symbols("Y Q P R C u")


for m in range(1, 4):
    for r in range(1, 4):
        F = sp.expand(P - (Q - Y) * Y**m)
        integrand = sp.expand(
            (Y**m - u * (Q + (Y - Q) * u) ** m) ** r
        )
        I = sp.integrate(integrand, (u, 0, 1))
        G = sp.expand(C * I - R * Y ** (m * (r + 1)))

        beta = sp.factor(
            (-1) ** r
            * Q ** (m * r)
            * sp.factorial(r)
            * sp.factorial(m * r)
            / sp.factorial(m * r + r + 1)
        )
        assert sp.factor(I.subs(Y, 0) - beta) == 0
        assert sp.factor(I.subs(Y, Q) - Q ** (m * r) / (r + 1)) == 0

        resultant_trace = sp.factor(
            sp.resultant(F.subs(P, 0), G, Y)
        )
        expected = Q ** (m * r * (m + 1)) * ((r + 1) * R * Q**m - C)
        quotient = sp.factor(resultant_trace / expected)
        assert quotient != 0 and not quotient.has(Q, R)

        exponent = m * r * (m + 1)
        finite_factor = (r + 1) * R * Q**m - C
        assert sp.gcd(
            sp.Poly(Q**exponent, Q, R), sp.Poly(finite_factor, Q, R)
        ) == 1

        label = "coordinate trace exponent" if (m, r) == (1, 1) else "nilpotency index"
        print(f"PASS (m,r)=({m},{r}): {label} {exponent}")

print("PASS: all bounded cancellation construction scheme-boundary traces match the uniform formula")

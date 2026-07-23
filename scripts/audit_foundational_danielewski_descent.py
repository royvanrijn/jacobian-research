#!/usr/bin/env python3
"""Test natural Danielewski descents of the foundational Keller map.

For a source ordering (a,b,c), the open immersion associated with
P(s)=s Q(s) has image ring

    k[a,c,s=bc,w=b Q(s)] inside k[a,b,c].

Writing a polynomial in the monomial normal form obtained by replacing
b^p c^q with b^(p-q)s^q or c^(q-p)s^p, membership requires the coefficient
of every positive b^i to be divisible by Q(s)^i.  This script computes the
largest squarefree nonconstant Q allowed simultaneously by all three
coordinates of the foundational Keller map, for every permutation of its
source variables.
"""

from __future__ import annotations

import itertools

import sympy as sp


x, y, z = sp.symbols("x y z")
u = 1 + x * y
foundational = (
    u**3 * z + y**2 * u * (4 + 3 * x * y),
    y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)

a, b, c, s = sp.symbols("a b c s")


def positive_b_coefficients(polynomial: sp.Expr) -> dict[int, sp.Expr]:
    """Return h_i(a,s) in the b-positive part of monomial normal form."""

    result: dict[int, sp.Expr] = {}
    expanded = sp.Poly(sp.expand(polynomial), a, b, c)
    for (a_power, b_power, c_power), coefficient in expanded.terms():
        if b_power <= c_power:
            continue
        residual_b = b_power - c_power
        contribution = coefficient * a**a_power * s**c_power
        result[residual_b] = sp.expand(
            result.get(residual_b, sp.Integer(0)) + contribution
        )
    return result


def coefficient_gcd(polynomial: sp.Expr) -> sp.Expr:
    """GCD in k[s] of the coefficients with respect to a."""

    coefficients = sp.Poly(polynomial, a).all_coeffs()
    gcd = sp.Poly(coefficients[0], s)
    for coefficient in coefficients[1:]:
        gcd = sp.gcd(gcd, sp.Poly(coefficient, s))
    return sp.monic(gcd).as_expr()


candidate_results: list[sp.Expr] = []
for suspension, positive, divisor in itertools.permutations((x, y, z)):
    substitution = {suspension: a, positive: b, divisor: c}
    required_valuations: dict[sp.Expr, int] | None = None
    b_profiles: list[dict[int, sp.Expr]] = []
    for coordinate in foundational:
        renamed = sp.expand(coordinate.xreplace(substitution))
        profile = positive_b_coefficients(renamed)
        b_profiles.append(profile)
        for power, coefficient in profile.items():
            gcd = coefficient_gcd(coefficient)
            factorization = sp.factor_list(gcd, s)[1]
            valuations = {
                factor: exponent // power
                for factor, exponent in factorization
                if factor.subs(s, 0) != 0 and exponent // power > 0
            }
            if required_valuations is None:
                required_valuations = valuations
            else:
                required_valuations = {
                    factor: min(exponent, valuations.get(factor, 0))
                    for factor, exponent in required_valuations.items()
                    if valuations.get(factor, 0) > 0
                }

    candidate = sp.Integer(1)
    for factor, exponent in (required_valuations or {}).items():
        # Q must be squarefree for sQ(s) to define a smooth Danielewski
        # surface, so only the presence of a factor matters.
        if exponent > 0:
            candidate *= factor
    candidate = sp.factor(candidate)
    candidate_results.append(candidate)
    print(
        "ORDER",
        f"a={suspension}",
        f"b={positive}",
        f"c={divisor}",
        "candidate_Q=" + str(candidate),
        "b_powers=" + str(tuple(tuple(sorted(p)) for p in b_profiles)),
    )

assert len(candidate_results) == 6
assert all(candidate == 1 for candidate in candidate_results)
print("PASS: no coordinate permutation admits a nontrivial smooth Danielewski descent")

#!/usr/bin/env python3
"""Exact bounded audit of the fixed-r branch-at-infinity mechanism.

This does not attempt the residual all-r wedge.  It constructs the limiting
endpoint system for r=5,6,7 and checks whether the algebraic ingredients used
in the r=6 proof survive in the first subsequent column: a squarefree branch
polynomial away from c=0, no branch at z=infinity, and unique reconstruction
of z by the linear subresultant.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import QQ


c, z = sp.symbols("c z")


def beta_limit(k: int) -> sp.Rational:
    return sp.factorial(k) / sp.Integer(k) ** (k + 1)


def tail_limit(k: int) -> sp.Expr:
    return sp.expand(
        sum(
            sp.binomial(k, j)
            * (-c) ** (k - j)
            * sp.factorial(j)
            / sp.Integer(k) ** (j + 1)
            for j in range(k + 1)
        )
    )


def limiting_system(r: int) -> tuple[sp.Expr, sp.Expr]:
    endpoint = c**r * z ** (r - 1)
    for k in range(1, r):
        endpoint += (
            (-1) ** k
            * sp.binomial(r, k)
            * z ** (r - k - 1)
            * (beta_limit(k) - z**k * tail_limit(k))
            * (-c) ** (r - k - 1)
        )
    moment = beta_limit(r) - z**r * tail_limit(r)
    return sp.factor(endpoint), sp.factor(moment)


def audit(r: int, expected_valuation: int, expected_degree: int) -> None:
    endpoint, moment = limiting_system(r)
    assert (sp.degree(endpoint, z), sp.degree(moment, z)) == (r - 1, r)

    resultant = sp.Poly(sp.resultant(endpoint, moment, z), c, domain=QQ)
    valuation = next(j for j in range(resultant.degree() + 1) if resultant.nth(j))
    assert valuation == expected_valuation
    branch = sp.Poly(resultant.as_expr() / c**valuation, c, domain=QQ)
    branch = branch.clear_denoms(convert=True)[1].primitive()[1]
    assert branch.degree() == expected_degree
    assert branch.eval(0) != 0
    assert sp.gcd(branch, branch.diff()).degree() == 0
    assert sp.gcd(branch, sp.Poly(tail_limit(r), c, domain=QQ)).degree() == 0

    subresultants = sp.subresultants(endpoint, moment, z)
    assert [sp.degree(member, z) for member in subresultants] == list(
        range(r, -1, -1)
    )
    linear = sp.Poly(subresultants[-2], z, domain=QQ.poly_ring(c))
    coefficient = sp.Poly(linear.coeff_monomial(z), c, domain=QQ)
    assert linear.degree() == 1
    assert sp.gcd(branch, coefficient).degree() == 0
    print(
        f"PASS fixed-r limiting schema r={r}: "
        f"c^{valuation} times a squarefree degree-{branch.degree()} branch "
        "polynomial with unique finite z"
    )


audit(5, expected_valuation=5, expected_degree=20)
audit(6, expected_valuation=7, expected_degree=29)
audit(7, expected_valuation=7, expected_degree=42)
print("SCOPE: fixed-r limiting systems only; no uniform continuation in r")

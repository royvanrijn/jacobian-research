#!/usr/bin/env python3
"""Exact checks for the fixed-r Newton-ramification extraction."""
from __future__ import annotations

import math

import sympy as sp


x = sp.symbols("x")


def derivative_polynomial(m: int, r: int) -> sp.Poly:
    degree = m * r
    return sp.Poly(
        sum(sp.binomial(r + j, r) * x**j for j in range(degree + 1)),
        x,
        domain=sp.ZZ,
    )


def reciprocal_numerator(m: int, r: int) -> sp.Poly:
    degree = m * r
    polynomial = derivative_polynomial(m, r).as_expr()
    expression = sp.expand(
        (x - 1) ** (r + 1) * x**degree * polynomial.subs(x, 1 / x)
    )
    return sp.Poly(expression, x, domain=sp.ZZ)


def scaled_remainder_is_unit(
    remainder: sp.Poly, cyclotomic: sp.Poly, prime: int, scale: int = 1
) -> bool:
    expression = sum(
        (int(coefficient) // scale) * x ** exponent[0]
        for exponent, coefficient in remainder.as_dict().items()
    )
    reduced = sp.Poly(expression, x, modulus=prime)
    cyclotomic_reduced = sp.Poly(cyclotomic.as_expr(), x, modulus=prime)
    return not reduced.is_zero and sp.gcd(reduced, cyclotomic_reduced).degree() == 0


def assert_newton_edge(m: int, r: int, prime: int) -> None:
    degree = m * r
    near_degree = degree + 1
    root_order = near_degree // prime
    assert near_degree == prime * root_order
    assert root_order > 1 and math.gcd(prime, root_order) == 1
    assert prime > r and prime <= degree - 3

    polynomial = reciprocal_numerator(m, r)
    cyclotomic = sp.Poly(sp.cyclotomic_poly(root_order, x), x, domain=sp.ZZ)

    for j in range(prime + 1):
        taylor_coefficient = sp.Poly(
            sp.diff(polynomial.as_expr(), x, j) / math.factorial(j),
            x,
            domain=sp.ZZ,
        )
        remainder = taylor_coefficient.rem(cyclotomic)
        coefficients = [int(coefficient) for coefficient in remainder.all_coeffs()]

        if j < prime:
            assert all(coefficient % prime == 0 for coefficient in coefficients)
        if j == 0:
            assert any(coefficient % (prime * prime) != 0 for coefficient in coefficients)
            assert scaled_remainder_is_unit(remainder, cyclotomic, prime, prime)
        if j == prime:
            assert scaled_remainder_is_unit(remainder, cyclotomic, prime)


# Structural derivative and reciprocal-numerator identities.
for r in range(1, 7):
    for m in range(1, 9):
        degree = m * r
        geometric_degree = (m + 1) * r
        geometric = sum(x**j for j in range(geometric_degree + 1))
        derivative = sp.Poly(
            sp.diff(geometric, x, r) / math.factorial(r), x, domain=sp.ZZ
        )
        expected = derivative_polynomial(m, r)
        assert derivative == expected

        numerator = reciprocal_numerator(m, r)
        assert numerator.degree() == geometric_degree + 1
        assert numerator.LC() == 1
        assert sp.rem(numerator.as_expr(), (x - 1) ** (r + 1), x) == 0

        # The BFLT congruence for a=d+1, including prime-power modulus p^ell.
        near_degree = degree + 1
        comparison = sp.Poly(
            numerator.as_expr() - (x**near_degree - 1) * x**r,
            x,
            domain=sp.ZZ,
        )
        for prime, exponent in sp.factorint(near_degree).items():
            prime = int(prime)
            if prime <= r:
                continue
            modulus = prime ** int(exponent)
            assert all(int(coefficient) % modulus == 0 for coefficient in comparison.all_coeffs())

print("PASS: derivative, reciprocal-numerator, and BFLT congruence identities")


# One exact unramified cyclotomic-cluster certificate for each r in 1,...,8.
# Entries are (m,p); a=mr+1=p*a_0 and Phi_(a_0) supplies the residue roots.
EDGE_CERTIFICATES = {
    1: (5, 2),
    2: (7, 3),
    3: (7, 11),
    4: (8, 11),
    5: (4, 7),
    6: (9, 11),
    7: (3, 11),
    8: (4, 11),
}

for r, (m, prime) in EDGE_CERTIFICATES.items():
    assert_newton_edge(m, r, prime)
    degree = m * r
    root_order = (degree + 1) // prime
    print(
        f"PASS r={r}, m={m}: Phi_{root_order} cluster has edge "
        f"(0,1)->({prime},0)"
    )

print("PASS: fixed-r Newton-ramification certificates through r=8")

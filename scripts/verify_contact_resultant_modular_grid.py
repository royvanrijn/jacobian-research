#!/usr/bin/env python3
"""Single-prime certificates for a finite grid of contact resultants.

The endpoint identity

    K_{m,k}(1-y) = beta_{m,k} * sum_{j=0}^{mk} binomial(j+k,k) y^j

and the triangular formula

    L_{m,r}(1-y) = sum_{k=0}^r (-1)^k binomial(r,k)
                         y^(m(r-k)) K_{m,k}(1-y)

construct the two rational polynomials directly, without integrations or a
Sylvester resultant.  Reduction modulo PRIME is legitimate because every
displayed denominator and both leading coefficients are units.  A modular
monic gcd equal to one then certifies coprimality over QQ.

This is finite evidence only; it does not prove the uniform r >= 6 theorem.
"""

from __future__ import annotations

import math

import sympy as sp


PRIME = 1_000_003
GRID = {
    5: 50,
    6: 40,
    7: 30,
    8: 25,
    9: 20,
    10: 16,
    11: 12,
    12: 10,
}

y = sp.Symbol("y")


def inverse_mod(value: int) -> int:
    """Return the inverse of a certified unit modulo PRIME."""

    residue = value % PRIME
    assert residue != 0
    return pow(residue, -1, PRIME)


def endpoint_moment_coefficients(m: int, k: int) -> list[int]:
    """Ascending coefficients of K_{m,k}(1-y), reduced modulo PRIME."""

    denominator = math.prod(m * k + j for j in range(1, k + 2))
    beta = math.factorial(k) * inverse_mod(denominator) % PRIME
    return [beta * math.comb(j + k, k) % PRIME for j in range(m * k + 1)]


def endpoint_pair(m: int, r: int) -> tuple[sp.Poly, sp.Poly]:
    """Return K_{m,r}(1-y), L_{m,r}(1-y) over GF(PRIME)."""

    degree = m * r
    moments = [endpoint_moment_coefficients(m, k) for k in range(r + 1)]
    k_coefficients = moments[r]
    l_coefficients = [0] * (degree + 1)
    for k, moment in enumerate(moments):
        shift = m * (r - k)
        multiplier = (-1) ** k * math.comb(r, k)
        for index, coefficient in enumerate(moment):
            l_coefficients[index + shift] = (
                l_coefficients[index + shift] + multiplier * coefficient
            ) % PRIME

    K = sp.Poly.from_list(k_coefficients[::-1], gens=y, modulus=PRIME)
    L = sp.Poly.from_list(l_coefficients[::-1], gens=y, modulus=PRIME)
    assert K.degree() == degree
    assert L.degree() == degree
    assert K.LC() % PRIME != 0
    assert L.LC() % PRIME != 0
    return K, L


certificate_count = 0
for r, maximum_m in GRID.items():
    for m in range(1, maximum_m + 1):
        K, L = endpoint_pair(m, r)
        assert sp.gcd(K, L).monic() == sp.Poly(1, y, modulus=PRIME)
        certificate_count += 1


print(
    "PASS contact resultant finite grid: "
    f"{certificate_count} monic gcd certificates modulo {PRIME}"
)
for r, maximum_m in GRID.items():
    print(f"PASS contact resultant finite grid: r={r}, 1<=m<={maximum_m}")
print("SCOPE: finite evidence only; the uniform continuation begins at r=6")

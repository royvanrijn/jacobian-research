#!/usr/bin/env python3
"""Exact algebra behind the irreducibility-to-contact corollary."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.boundary import (  # noqa: E402
    cancellation_branch_polynomial,
    cancellation_critical_coefficient,
)


w, q = sp.symbols("w q")


# K and L have the same degree and constant coefficient, but different linear
# coefficients.  Therefore irreducibility of K implies gcd(K,L)=1: a common
# root would make K divide L, and equal degrees/constants would then force
# K=L, contradicting the linear term.
for m in range(1, 7):
    for r in range(1, 7):
        K = sp.Poly(cancellation_branch_polynomial(m, r, w), w)
        L = sp.Poly(cancellation_critical_coefficient(m, r, w), w)
        degree = m * r

        assert K.degree() == L.degree() == degree
        assert K.coeff_monomial(1) == L.coeff_monomial(1) == sp.Rational(1, r + 1)
        assert K.coeff_monomial(w) == -sp.Rational(m * r, r + 2)
        assert L.coeff_monomial(w) == -sp.Rational(
            m * r * (r + 3), (r + 1) * (r + 2)
        )
        assert K.coeff_monomial(w) != L.coeff_monomial(w)

        # The fractional-linear transform identifies K with the cancellation
        # parameter polynomial, so every proved parameter-irreducibility range
        # transfers to K.
        parameter = sum(
            (-1) ** j
            * sp.binomial(degree + r + 1, j)
            * q ** (degree - j)
            for j in range(degree + 1)
        )
        transformed = sp.Poly(
            sum(
                K.coeff_monomial(w**power)
                * (-q) ** power
                * (1 - q) ** (degree - power)
                for power in range(degree + 1)
            ),
            q,
        )
        quotient, remainder = sp.div(sp.Poly(parameter, q), transformed)
        assert remainder.is_zero
        assert quotient.degree() == 0


print("PASS contact resultant: irreducible K cannot share a root with L")
print("PASS contact resultant: K is the parameter polynomial up to a Mobius transform")
print("PASS contact resultant: arithmetic irreducibility ranges transfer unchanged")

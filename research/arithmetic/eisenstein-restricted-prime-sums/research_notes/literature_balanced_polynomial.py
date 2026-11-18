#!/usr/bin/env python3
"""Exact elementary certificate for a proposed balanced-summand sieve weight.

This certifies the polynomial inequalities only.  It is NOT a proof of the
balanced Eisenstein theorem.  The required archimedean/twisted spectral
estimates are described separately and have not been independently audited.
"""

from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sympy as s


z, t = s.symbols("z t")
P = s.Poly(18000 * (s.Rational(1, 9) - z) * (z + s.Rational(1, 36))
           * (z - s.Rational(3, 10))**2
           * (z - s.Rational(9, 16))**2
           * (z - s.Rational(4, 5))**2
           * (z - s.Rational(24, 25))**2, z)


def fraction(x):
    return F(int(s.numer(x)), int(s.denom(x)))


def bernstein_coefficients(poly, left, right):
    """Power coefficients on [left,right] -> degree-d Bernstein coefficients."""
    d = poly.degree()
    a = [fraction(poly.nth(j)) for j in range(d + 1)]
    width = right - left
    shifted = [sum(a[j] * comb(j, k) * left**(j-k) * width**k
                   for j in range(k, d + 1)) for k in range(d + 1)]
    return [sum(shifted[k] * F(comb(i, k), comb(d, k))
                for k in range(i + 1)) for i in range(d + 1)]


def positive_certificate(poly, left, right, pieces):
    minima = []
    for i in range(pieces):
        a = left + (right-left) * F(i, pieces)
        b = left + (right-left) * F(i+1, pieces)
        minimum = min(bernstein_coefficients(poly, a, b))
        assert minimum >= 0, (i, minimum)
        minima.append(str(minimum))
    return {"interval": [str(left), str(right)], "equal_pieces": pieces,
            "minimum_Bernstein_coefficient_per_piece": minima,
            "nonnegative_coefficients_verified": True}


def main():
    upper = positive_certificate(s.Poly(1-P.as_expr(), z), F(0), F(1,9), 64)
    lower = positive_certificate(s.Poly(P.as_expr()+s.Rational(1,5), z), F(0), F(1), 256)
    # The factor (1/9-z) and squared factors prove P<=0 for z>=1/9.
    even_polynomial = P.as_expr().subs(z, t*t)
    coefficients = []
    for j in range(0,21,2):
        c = s.integrate(s.expand(even_polynomial*s.legendre(j,t)), (t,-1,1)) * s.Rational(2*j+1,2)
        coefficients.append(fraction(c))
    mean = coefficients[0]
    nonconstant_l1 = sum(abs(c) for c in coefficients[1:])
    assert mean > F(7,32)
    assert nonconstant_l1 < F(7,2)
    # A rigorous elementary majorant for the sum over odd inert primes.
    odd_tail_bound = sum(F(2,(6*k-1)*(6*k-2)) for k in range(1,11)) + F(1,3*(60-2))
    assert odd_tail_bound < F(1,7)
    out = {
        "scope": "Polynomial and rational-series inequalities only; not a proof of balanced existence",
        "P_factorization": str(s.factor(P.as_expr())),
        "variable": "z=t^2, t=(a-b)/(a+b)",
        "minorant": "P(t^2) <= indicator(|t|<=1/3)",
        "global_lower_bound": "P(t^2)>=-1/5 for -1<=t<=1",
        "upper_certificate": upper,
        "lower_certificate": lower,
        "Legendre_even_coefficients_degrees_0_to_20": [str(c) for c in coefficients],
        "mean": str(mean), "mean_decimal": float(mean),
        "mean_lower_bound": "7/32",
        "nonconstant_Legendre_l1": str(nonconstant_l1),
        "nonconstant_Legendre_l1_decimal": float(nonconstant_l1),
        "nonconstant_Legendre_l1_upper_bound": "7/2",
        "odd_inert_main_sum_bound": str(odd_tail_bound),
        "odd_inert_main_sum_upper_bound": "1/7",
    }
    path = Path(__file__).with_name("literature_balanced_polynomial_receipt.json")
    path.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if not k.endswith("certificate")}, indent=2))


if __name__ == "__main__":
    main()

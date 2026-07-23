#!/usr/bin/env python3
"""Exact audit of constructive adelic Keller-fiber engineering."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import (
    constructive_weak_approximation_lift,
    rational_good_reduction_certificate,
)


w = sp.symbols("w")
warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


def factor_degree_partition(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    """Return the squarefree factor-degree partition over F_p."""
    _, factors = sp.factor_list(poly.as_expr(), w, modulus=prime)
    assert all(multiplicity == 1 for _, multiplicity in factors)
    degrees = tuple(
        sorted((sp.Poly(factor, w, modulus=prime).degree() for factor, _ in factors),
               reverse=True)
    )
    assert sum(degrees) == poly.degree()
    return degrees


def reduce_rational_poly(poly: sp.Poly, prime: int) -> sp.Poly:
    """Reduce a Q-polynomial whose coefficient denominators are p-units."""
    expression = 0
    for (exponent,), coefficient in poly.terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        assert denominator % prime
        reduced = numerator * pow(denominator, -1, prime) % prime
        expression += reduced * w**exponent
    return sp.Poly(expression, w, modulus=prime)


# Lift through a small real box near the known zero-root chamber.  The
# resulting rational point is certified below by an exact Sturm count.  The
# two local witnesses are (s,t)=(0,1) modulo 7 and 11, with different types.
residue_data = {7: (0, 1), 11: (0, 1)}
intervals = (
    (Fraction(-301, 100), Fraction(-299, 100)),
    (Fraction(199, 100), Fraction(201, 100)),
)
s_fraction, t_fraction = constructive_weak_approximation_lift(
    residue_data, intervals
)
assert (s_fraction, t_fraction) == (Fraction(-308, 103), Fraction(617, 309))

for prime, residue_vector in residue_data.items():
    for value, residue in zip((s_fraction, t_fraction), residue_vector):
        assert value.denominator % prime
        reduction = value.numerator * pow(value.denominator, -1, prime) % prime
        assert reduction == residue

s = sp.Rational(s_fraction.numerator, s_fraction.denominator)
t = sp.Rational(t_fraction.numerator, t_fraction.denominator)
H4 = sp.Poly(w**4 / 4 - 3 * w**3 / 2 + 5 * w**2 / 4, w, domain=sp.QQ)
E = sp.Poly(H4.as_expr() - s * w + t, w, domain=sp.QQ)
primitive_E = sp.Poly(1236 * E.as_expr(), w, domain=sp.ZZ)
assert primitive_E == sp.Poly(
    309 * w**4 - 1854 * w**3 + 1545 * w**2 + 3696 * w + 2468,
    w,
    domain=sp.ZZ,
)
assert sp.gcd_list(primitive_E.all_coeffs()) == 1
a0 = sp.Rational(-5, 3)
certificate = rational_good_reduction_certificate(
    H4.as_expr(), w, c=1, b0=1, a0=a0
)
assert all(certificate["bad_integer"] % prime for prime in residue_data)

assert E.gcd(E.diff()).degree() == 0
assert E.count_roots(-sp.oo, sp.oo) == 0
assert factor_degree_partition(E, 7) == (4,)
assert factor_degree_partition(E, 11) == (2, 1, 1)

assert reduce_rational_poly(E, 7).monic() == sp.Poly(
    w**4 + w**3 - 2 * w**2 - 3, w, modulus=7
)
assert reduce_rational_poly(E, 11).monic() == sp.Poly(
    (w + 2) * (w - 4) * (w**2 - 4 * w + 5), w, modulus=11
)

# Irreducibility modulo 7 proves irreducibility over Q after primitive
# denominator clearing.  SymPy's exact Q-factorization independently audits it.
assert sp.Poly(E, w, domain=sp.QQ).is_irreducible
assert len(sp.factor_list(E.as_expr())[1]) == 1

A, B, C = t, s, sp.Integer(1)
assert sp.Poly(H4.as_expr() - B * C * w + A * C**2, w) == E

print("PASS: the denominator-one CRT grid meets the prescribed real box")
print("PASS: the quartic fiber has signature (0,2) and is globally irreducible")
print("PASS: 7 is inert and 11 has unramified factor degrees (2,1,1)")
print(f"PASS: target (A,B,C)=({A},{B},{C}) reconstructs the audited pencil")

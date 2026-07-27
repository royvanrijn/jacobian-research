#!/usr/bin/env python3
"""Exact checks for split and connected full fibers in one Keller map.

For every N >= 3, the theorem uses the determinant-one quadratic gauge
attached to G_N(S) = prod_{j=0}^{N-1}(S-j).  Its target line (1,0,c) has
inverse polynomial G_N(S)-g_1*c/2, irreducible over QQ(c).  Hilbert
irreducibility supplies infinitely many connected specializations.

The general scheme-theoretic reconstruction is proved and checked separately
in verified/FINITE_ETALE_KELLER_FIBERS.md and
scripts/verify_finite_etale_keller_fibers.py.  This script checks:

* the displayed map has Jacobian one;
* its two inverse polynomials are the claimed split and field polynomials;
* the second polynomial is irreducible modulo 5 (hence over Q);
* both inverse polynomials are squarefree.
* the all-degree line formula and admissibility through degree twelve.
"""
from __future__ import annotations

import sympy as sp

x, y, z, S = sp.symbols("x y z S")

t = 1 + x * y
q = t**2 * z + sp.Rational(24, 35) * y**2 * (1 + 3 * t)

# The second coordinate is -1/2 times the determinant-minus-two coordinate
# in the normalized quadratic-gauge theorem.
F = (
    t * q,
    -sp.Rational(1, 2) * y
    - sp.Rational(35, 16) * x * q
    + sp.Rational(25, 12) * t * q
    + sp.Rational(5, 6) * t**2 * x**2 * q**4
    - sp.Rational(5, 48) * t**2 * x**3 * q**5,
    x * (5 - 3 * t)
    - sp.Rational(35, 24) * x**3 * z
    + sp.Rational(5, 6) * (x * q) ** 4
    - sp.Rational(1, 8) * (x * q) ** 5,
)

assert sp.factor(sp.Matrix(F).jacobian((x, y, z)).det()) == 1


def inverse_polynomial(pi: sp.Rational, u: sp.Rational, c: sp.Rational) -> sp.Expr:
    """Inverse polynomial for the displayed determinant-one coordinates."""
    return sp.expand(
        pi**5 * S**5
        - 10 * pi**4 * S**4
        + 35 * pi * S**3
        + (-50 * pi + 24 * u) * S**2
        + 24 * S
        - 12 * c
    )


# Derive the inverse relation directly from the displayed map on the rational
# chart S=x/t.  Targets with pi=1 lie entirely in this chart because pi=tq.
assert sp.factor(sp.cancel(inverse_polynomial(*F).subs(S, x / t))) == 0

y_split = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
y_field = (sp.Integer(1), sp.Integer(0), sp.Integer(1))

P_split = inverse_polynomial(*y_split)
P_field = inverse_polynomial(*y_field)

assert sp.expand(P_split - sp.prod(S - root for root in range(5))) == 0
assert P_field == S**5 - 10 * S**4 + 35 * S**3 - 50 * S**2 + 24 * S - 12

for polynomial in (P_split, P_field):
    poly = sp.Poly(polynomial, S, domain=sp.QQ)
    assert sp.gcd(poly, poly.diff()).degree() == 0

# Modulo 5 this is the Artin--Schreier polynomial S^5-S-2.  Rabin's
# criterion gives a short machine-checkable irreducibility certificate.
p = 5
field_modulus = sp.Poly(P_field, S, modulus=p)
assert field_modulus == sp.Poly(S**5 - S - 2, S, modulus=p)
frobenius_1 = sp.Poly(S**p - S, S, modulus=p)
assert sp.gcd(field_modulus, frobenius_1).degree() == 0

power = sp.Poly(S, S, modulus=p)
for _ in range(5):
    power = (power**p).rem(field_modulus)
assert (power - sp.Poly(S, S, modulus=p)).rem(field_modulus).is_zero
assert field_modulus.is_irreducible

# All-degree structural regression.  The mathematical proof of generic
# irreducibility is the degree-one-in-c argument in the theorem note; this
# bounded loop is deliberately only an executable regression.
c_parameter = sp.symbols("c")
irreducible_witness_primes = {
    3: 2,
    4: 5,
    5: 5,
    6: 29,
    7: 7,
    8: 11,
    9: 47,
    10: 11,
    11: 11,
    12: 13,
}
for degree, witness_prime in irreducible_witness_primes.items():
    seed = sp.expand(sp.prod(S - root for root in range(degree)))
    seed_poly = sp.Poly(seed, S, domain=sp.QQ)
    g1 = seed_poly.coeff_monomial(S)
    g3 = seed_poly.coeff_monomial(S**3)
    assert g1 == (-1) ** (degree - 1) * sp.factorial(degree - 1)
    assert g3 != 0

    line_inverse = sp.expand(seed - g1 * c_parameter / 2)
    assert sp.diff(line_inverse, c_parameter) == -g1 / 2
    assert sp.Poly(line_inverse, S, domain=sp.QQ.frac_field(c_parameter)).is_irreducible

    witness = sp.Poly(line_inverse.subs(c_parameter, 1), S, domain=sp.QQ)
    assert sp.Poly(witness, S, modulus=witness_prime).is_irreducible

print("PASS: the displayed quadratic-gauge map has Jacobian one")
print("PASS: its displayed inverse polynomial vanishes identically at S=x/t")
print("PASS: y_split=(1,0,0) has inverse prod(S-j), 0<=j<=4")
print("PASS: y_field=(1,0,1) has the claimed squarefree quintic")
print("PASS: the field quintic is Artin--Schreier irreducible modulo 5")
print("PASS: the all-degree target-line formula regresses through degree twelve")

#!/usr/bin/env python3
"""Exact audit of local-to-global Keller-fiber synthesis.

The abstract theorem uses monogenicity, Hensel factor stability, Krasner's
lemma, and weak approximation.  This checker certifies all exact hypotheses
for the ramified quintic witness and reuses the quadratic-gauge verifier for
the final scheme-theoretic compilation.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.local_global import (
    automatic_local_stability_certificate,
    certify_polynomial_in_local_stability_ball,
    constructive_coefficient_crt,
    rational_valuation,
    synthesize_monic_polynomial,
)
from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber
from verify_finite_etale_keller_fibers import (
    S,
    check_scheme_reconstruction,
    x,
    y,
    z,
)


T = sp.symbols("T")


def valuation(integer: int, prime: int) -> int:
    """Return v_p(integer) for a nonzero integer."""
    assert integer
    value = abs(integer)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def coefficient_congruence(
    polynomial: sp.Poly, local_model: sp.Poly, prime: int, exponent: int
) -> None:
    """Check congruence of rational coefficients with p-unit denominators."""
    modulus = prime**exponent
    for degree in range(polynomial.degree() + 1):
        difference = polynomial.coeff_monomial(T**degree) - local_model.coeff_monomial(
            T**degree
        )
        numerator, denominator = map(int, sp.fraction(difference))
        assert denominator % prime
        if numerator:
            assert valuation(numerator, prime) >= exponent


def reduced_poly(polynomial: sp.Poly, prime: int) -> sp.Poly:
    """Reduce a rational polynomial with p-integral coefficients."""
    expression = 0
    for (degree,), coefficient in polynomial.terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        assert denominator % prime
        expression += numerator * pow(denominator, -1, prime) * T**degree
    return sp.Poly(expression, T, modulus=prime)


P = sp.Poly(
    T**5
    + sp.Rational(225, 1261) * T**4
    + sp.Rational(5765, 1261) * T**3
    + sp.Rational(190, 1261) * T**2
    + sp.Rational(4854, 1261) * T
    + sp.Rational(294, 1261),
    T,
    domain=sp.QQ,
)
primitive_P = sp.Poly(
    1261 * T**5 + 225 * T**4 + 5765 * T**3 + 190 * T**2 + 4854 * T + 294,
    T,
    domain=sp.ZZ,
)
assert sp.Poly(1261 * P.as_expr(), T, domain=sp.ZZ) == primitive_P
assert sp.gcd_list(primitive_P.all_coeffs()) == 1

ramified_cubic_2 = sp.Poly(T**3 - 2, T, domain=sp.ZZ)
unramified_quadratic_2 = sp.Poly(T**2 + T + 1, T, domain=sp.ZZ)
local_2 = ramified_cubic_2 * unramified_quadratic_2

ramified_quadratic_3 = sp.Poly(T**2 - 3, T, domain=sp.ZZ)
unramified_cubic_3 = sp.Poly(T**3 - T + 1, T, domain=sp.ZZ)
local_3 = ramified_quadratic_3 * unramified_cubic_3

local_5 = sp.Poly(T**5 - T - 1, T, domain=sp.QQ)
local_7 = sp.Poly(T * (T**2 + 1) * (T**2 + T + 3), T, domain=sp.QQ)
coefficient_box = (
    (sp.Rational(-1, 2), sp.Rational(1, 2)),
    (sp.Rational(7, 2), sp.Rational(9, 2)),
    (sp.Rational(-1, 2), sp.Rational(1, 2)),
    (sp.Rational(9, 2), sp.Rational(11, 2)),
    (sp.Rational(-1, 2), sp.Rational(1, 2)),
)
synthesis = synthesize_monic_polynomial(
    {
        2: (2, local_2),
        3: (2, local_3),
        5: (1, local_5),
        7: (1, local_7),
    },
    coefficient_box,
    T,
)
assert synthesis.polynomial == P
assert synthesis.certificate.base_denominator == 1
assert synthesis.certificate.crt_modulus == 4 * 9 * 5 * 7
assert synthesis.certificate.multiplier == 1
assert synthesis.certificate.common_denominator == 1261

certificate_2 = automatic_local_stability_certificate(local_2, T, 2)
certificate_3 = automatic_local_stability_certificate(local_3, T, 3)
assert (certificate_2.discriminant_valuation, certificate_2.coefficient_precision) == (
    2,
    5,
)
assert (certificate_3.discriminant_valuation, certificate_3.coefficient_precision) == (
    1,
    3,
)
automatic_synthesis = synthesize_monic_polynomial(
    {
        2: (None, local_2),
        3: (None, local_3),
        5: (None, local_5),
        7: (None, local_7),
    },
    coefficient_box,
    T,
)
automatic_P = sp.Poly(
    T**5
    - sp.Rational(9855, 30241) * T**4
    + sp.Rational(163265, 30241) * T**3
    + sp.Rational(190, 30241) * T**2
    + sp.Rational(113214, 30241) * T
    - sp.Rational(7266, 30241),
    T,
    domain=sp.QQ,
)
assert automatic_synthesis.polynomial == automatic_P
assert automatic_synthesis.certificate.crt_modulus == 2**5 * 3**3 * 5 * 7
assert automatic_synthesis.certificate.common_denominator == 30241
assert {
    certificate.prime: certificate.coefficient_precision
    for certificate in automatic_synthesis.automatic_local_certificates
} == {2: 5, 3: 3, 5: 1, 7: 1}
for certificate in automatic_synthesis.automatic_local_certificates:
    certify_polynomial_in_local_stability_ball(automatic_P, certificate)
assert automatic_P.count_roots(-sp.oo, sp.oo) == 1
assert automatic_P.is_irreducible
assert reduced_poly(automatic_P, 5).is_irreducible
assert reduced_poly(automatic_P, 7) == sp.Poly(local_7.as_expr(), T, modulus=7)

# Denominators divisible by the selected primes are absorbed into the
# prime-power CRT moduli rather than silently treated as p-adic units.
fractional_center_lift = constructive_coefficient_crt(
    {
        2: (2, (Fraction(1, 2),)),
        3: (2, (Fraction(1, 3),)),
    },
    ((Fraction(-1), Fraction(1)),),
)
assert fractional_center_lift.base_denominator == 6
assert fractional_center_lift.crt_modulus == 2**3 * 3**3
fractional_coefficient = fractional_center_lift.coefficients[0]
assert rational_valuation(fractional_coefficient - Fraction(1, 2), 2) == 2
assert rational_valuation(fractional_coefficient - Fraction(1, 3), 3) == 2

coefficient_congruence(P, local_2, 2, 2)
coefficient_congruence(P, local_3, 3, 2)

# Unit resultants make the two local factors Hensel-stable independently.
assert int(sp.resultant(ramified_cubic_2.as_expr(), unramified_quadratic_2.as_expr(), T)) % 2
assert int(sp.resultant(ramified_quadratic_3.as_expr(), unramified_cubic_3.as_expr(), T)) % 3
assert reduced_poly(unramified_quadratic_2, 2).is_irreducible
assert reduced_poly(unramified_cubic_3, 3).is_irreducible

# Tame binomial Krasner certificates.  A coefficient perturbation of
# valuation k gives a root displacement of valuation at least k-v(f'(alpha)).
local_certificates = (
    # p, ramification degree, perturbation valuation
    (2, 3, 2),
    (3, 2, 2),
)
for prime, ramification_degree, perturbation_valuation in local_certificates:
    derivative_valuation = Fraction(ramification_degree - 1, ramification_degree)
    conjugate_distance_valuation = Fraction(1, ramification_degree)
    root_displacement_valuation = (
        Fraction(perturbation_valuation) - derivative_valuation
    )
    assert Fraction(perturbation_valuation) > 2 * derivative_valuation
    assert root_displacement_valuation > conjugate_distance_valuation

# The selected unramified Frobenius fingerprints.
mod_5 = sp.Poly(local_5.as_expr(), T, modulus=5)
mod_7_factors = (
    sp.Poly(T, T, modulus=7),
    sp.Poly(T**2 + 1, T, modulus=7),
    sp.Poly(T**2 + T + 3, T, modulus=7),
)
assert reduced_poly(P, 5) == mod_5
assert mod_5.is_irreducible
assert reduced_poly(P, 7) == sp.Poly(
    sp.prod(factor.as_expr() for factor in mod_7_factors), T, modulus=7
)
assert all(factor.is_irreducible for factor in mod_7_factors)
assert len({factor.as_expr() for factor in mod_7_factors}) == 3

# The inert prime proves global irreducibility; Sturm gives the signature.
assert P.is_irreducible
assert P.gcd(P.diff()).degree() == 0
assert P.count_roots(-sp.oo, sp.oo) == 1

# Compile P with translation a=0.
assert P.diff().eval(0) == sp.Rational(4854, 1261)
assert P.diff((T, 3)).eval(0) == 6 * sp.Rational(5765, 1261)
compilation = compile_polynomial_to_keller_fiber(
    P,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
)
mapping_minus_two = compilation.determinant_minus_two_map
mapping_one = compilation.determinant_one_map
assert sp.factor(sp.Matrix(mapping_minus_two).jacobian((x, y, z)).det()) == -2
assert sp.factor(sp.Matrix(mapping_one).jacobian((x, y, z)).det()) == 1
assert compilation.coordinate_degrees == (7, 32, 30)
assert compilation.geometric_degree == 5
assert compilation.target == (sp.Integer(1), sp.Integer(0), sp.Rational(-98, 809))
assert compilation.inverse_polynomial == P.as_expr().subs(T, S)
check_scheme_reconstruction(P.as_expr().subs(T, S), sp.Integer(0))

# Automatic translation skips a critical center and retains the translated
# quotient algebra.
automatic_compilation = compile_polynomial_to_keller_fiber(
    T**3 + T**2 + 1,
    T,
    inverse_variable=S,
    source_variables=(x, y, z),
)
assert automatic_compilation.translation == 1
assert automatic_compilation.inverse_polynomial == sp.expand(
    (S + 1) ** 3 + (S + 1) ** 2 + 1
)
automatic_local_compilation = compile_polynomial_to_keller_fiber(
    automatic_P,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
)
assert automatic_local_compilation.target == (
    sp.Integer(1),
    sp.Integer(0),
    sp.Rational(2422, 18869),
)
assert automatic_local_compilation.coordinate_degrees == (7, 32, 30)

print("PASS: coefficient CRT realizes the certified Q_2 and Q_3 algebras")
print("PASS: the generic prime-power coefficient synthesizer reproduces P")
print("PASS: discriminant radii automatically certify a second global quintic")
print("PASS: the quintic has signature (1,2), is inert at 5, and has type (2,2,1) at 7")
print("PASS: the compiled map has determinant one, geometric degree five, and degrees (7,32,30)")
print("PASS: target (1,0,-98/809) has complete fiber Spec(Q[T]/(P))")

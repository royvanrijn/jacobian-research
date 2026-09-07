#!/usr/bin/env python3
"""Exact audit of stratified seed-and-target adelic engineering."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.chebotarev import (  # noqa: E402
    constructive_weak_approximation_lift,
    factorization_type_mod_prime,
    rational_good_reduction_certificate,
)


warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)
w, x = sp.symbols("w x")

# The degree-five type-(3,2) omitted seed.
H = sp.Poly(
    sp.Rational(1, 27) * w**2 * (1 - w) * (w**2 + 8 * w + 18),
    w,
    domain=sp.QQ,
)
assert H.eval(0) == 0
assert H.diff().eval(0) == 0
assert H.eval(1) == 0
assert H.diff().eval(1) == -1
kappa = H.diff().diff().eval(1)
assert kappa == -sp.Rational(128, 27)
assert kappa != -2

omitted = sp.Poly(H.as_expr() + w - 1, w, domain=sp.QQ)
assert omitted == sp.Poly(
    -sp.Rational(1, 27) * (w + 3) ** 3 * (w - 1) ** 2,
    w,
    domain=sp.QQ,
)
assert sp.factor_list(omitted.as_expr())[1] == [(w - 1, 2), (w + 3, 3)]

# The reduced cubic Hessian divisor has trivial affine stabilizer. After
# translating its barycenter to zero, both the linear and constant terms are
# nonzero, excluding respectively mu_3 and mu_2 symmetry.
hessian_core = sp.Poly((w + 3) * (5 * w**2 + 6 * w - 3), w, domain=sp.QQ)
assert hessian_core.discriminant() != 0
barycenter = -sp.Rational(7, 5)
centered_hessian = sp.Poly(
    sp.expand(hessian_core.as_expr().subs(w, x + barycenter)),
    x,
    domain=sp.QQ,
)
assert centered_hessian == sp.Poly(
    5 * x**3 - sp.Rational(72, 5) * x - sp.Rational(64, 25),
    x,
    domain=sp.QQ,
)
assert centered_hessian.coeff_monomial(x) != 0
assert centered_hessian.coeff_monomial(1) != 0

residue_data = {7: (0, 4), 11: (1, 1)}
target_data = {
    5: (
        (
            (Fraction(-157, 100), Fraction(-156, 100)),
            (Fraction(46, 100), Fraction(47, 100)),
        ),
        (Fraction(-161, 103), Fraction(48, 103)),
    ),
    3: (
        (
            (Fraction(-151, 100), Fraction(-149, 100)),
            (Fraction(-1, 10), Fraction(1, 10)),
        ),
        (Fraction(-1435, 963), Fraction(-82, 963)),
    ),
    1: (
        (
            (Fraction(-151, 100), Fraction(-149, 100)),
            (Fraction(-21, 10), Fraction(-19, 10)),
        ),
        (Fraction(-1435, 963), Fraction(-223, 107)),
    ),
}

audited_pencils = {}
for real_roots, (box, expected_target) in target_data.items():
    lifted = constructive_weak_approximation_lift(residue_data, box)
    assert lifted == expected_target
    slope_fraction, intercept_fraction = lifted
    slope = sp.Rational(slope_fraction.numerator, slope_fraction.denominator)
    intercept = sp.Rational(
        intercept_fraction.numerator, intercept_fraction.denominator
    )
    E = sp.Poly(H.as_expr() - slope * w + intercept, w, domain=sp.QQ)
    assert E.gcd(E.diff()).degree() == 0
    assert E.count_roots(-sp.oo, sp.oo) == real_roots
    for prime, (s_residue, t_residue) in residue_data.items():
        assert (
            int(slope.p) % prime * pow(int(slope.q), -1, prime) % prime
            == s_residue
        )
        assert (
            int(intercept.p) % prime * pow(int(intercept.q), -1, prime) % prime
            == t_residue
        )
    assert factorization_type_mod_prime(E.as_expr(), w, 0, 0, 7) == (5,)
    assert factorization_type_mod_prime(E.as_expr(), w, 0, 0, 11) == (2, 2, 1)
    assert E.is_irreducible
    audited_pencils[real_roots] = E

E = audited_pencils[5]

a0 = -sp.cancel((1 + kappa) / (2 + kappa))
certificate = rational_good_reduction_certificate(
    H.as_expr(), w, c=1, b0=1, a0=a0
)
assert all(certificate["bad_integer"] % prime for prime in residue_data)

primitive = sp.Poly(-2781 * E.as_expr(), w, domain=sp.ZZ)
assert primitive == sp.Poly(
    103 * w**5
    + 721 * w**4
    + 1030 * w**3
    - 1854 * w**2
    - 4347 * w
    - 1296,
    w,
    domain=sp.ZZ,
)
assert sp.gcd_list(primitive.all_coeffs()) == 1
factor_7 = sp.Poly(w**5 + 3 * w**3 + 3 * w**2 - 3, w, modulus=7)
assert factor_7.is_irreducible
factor_11 = sp.Poly(
    (w - 1) * (w**2 - 2 * w - 1) * (w**2 - w - 5),
    w,
    modulus=11,
)
assert factor_11.gcd(factor_11.diff()).degree() == 0
assert tuple(
    sorted(
        (
            factor.degree()
            for factor, multiplicity in factor_11.factor_list()[1]
            if multiplicity == 1
        ),
        reverse=True,
    )
) == (2, 2, 1)

print("PASS: the seed is normalized admissible with exact omitted type (3,2)")
print("PASS: the reduced Hessian divisor has trivial affine stabilizer")
print("PASS: the lifted targets have exactly five, three, and one real sheets")
print("PASS: Frobenius cycle types are (5) at 7 and (2,2,1) at 11")
print("PASS: the complete fibers are quintic fields of all three signatures")

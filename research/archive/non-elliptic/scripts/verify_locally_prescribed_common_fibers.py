#!/usr/bin/env python3
"""Exact audit of locally prescribed fibers in one fixed Keller-map pair."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.common_fibers import (  # noqa: E402
    fixed_common_fiber_polynomial,
    synthesize_fixed_common_fiber_parameter,
)
from jcsearch.local_global import rational_valuation  # noqa: E402
from jcsearch.weighted import WeightedSeedModel, w  # noqa: E402


T = sp.symbols("T")
degree = 6

# The local centers define genuinely ramified algebras.  At 2 the polynomial
# factors exactly into Q_2, Q_2(i), and the unramified cubic.  At 3 the
# Newton segment of length two and slope -1/2 gives one ramified quadratic;
# the two other quadratic factors are unramified.
local_2 = fixed_common_fiber_polynomial(degree, -1, T)
assert sp.factor(local_2.as_expr()) == (
    (T - 1) * (T**2 + 1) * (T**3 + T**2 + 1)
)
assert sp.Poly(T**3 + T**2 + 1, T, modulus=2).is_irreducible
assert sp.discriminant(local_2.as_expr(), T) == 4464
assert rational_valuation(Fraction(4464), 2) == 4

local_3 = fixed_common_fiber_polynomial(degree, 2, T)
assert sp.discriminant(local_3.as_expr(), T) == 2136417
assert rational_valuation(Fraction(2136417), 3) == 1
translated_3 = sp.Poly(sp.expand(local_3.as_expr().subs(T, 1 + T)), T)
ascending_valuations = tuple(
    rational_valuation(
        Fraction(translated_3.coeff_monomial(T**index)),
        3,
    )
    for index in range(7)
)
assert ascending_valuations == (1, 1, 0, 1, 1, 1, 0)
assert sp.Poly(T**2 + 1, T, modulus=3).is_irreducible
assert sp.Poly(T**2 + 2 * T + 2, T, modulus=3).is_irreducible
assert sp.Poly(local_3.as_expr(), T, modulus=3).as_expr() == sp.Poly(
    (T + 2) ** 2 * (T**2 + 1) * (T**2 + 2 * T + 2),
    T,
    modulus=3,
).as_expr()

# The auxiliary local center at 5 is an unramified degree-six field and
# therefore proves connectedness of every rational lift in its residue ball.
local_5 = fixed_common_fiber_polynomial(degree, 1, T)
assert sp.Poly(local_5.as_expr(), T, modulus=5).is_irreducible

synthesis = synthesize_fixed_common_fiber_parameter(
    degree,
    {2: -1, 3: 2, 5: 1},
    (Fraction(1, 2), Fraction(3, 2)),
    T,
)
assert synthesis.parameter == Fraction(95231, 69121)
assert synthesis.crt_certificate.base_denominator == 1
assert synthesis.crt_certificate.crt_modulus == 2**9 * 3**3 * 5
assert synthesis.crt_certificate.multiplier == 1
assert synthesis.crt_certificate.common_denominator == 69121
assert {
    certificate.prime: certificate.stability.coefficient_precision
    for certificate in synthesis.local_certificates
} == {2: 9, 3: 3, 5: 1}

parameter = sp.Rational(95231, 69121)
P = fixed_common_fiber_polynomial(degree, parameter, T)
assert synthesis.polynomial == P
for prime, center, precision in ((2, -1, 9), (3, 2, 3), (5, 1, 1)):
    value = rational_valuation(Fraction(parameter - center), prime)
    assert value is not None and value >= precision

assert P.gcd(P.diff()).degree() == 0
assert P.count_roots(-sp.oo, sp.oo) == 2
assert sp.Poly(P.as_expr(), T, domain=sp.QQ).is_irreducible

def reduce_rational_polynomial(polynomial: sp.Poly, prime: int) -> sp.Poly:
    expression = 0
    for (exponent,), coefficient in polynomial.terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        assert denominator % prime
        expression += numerator * pow(denominator, -1, prime) * T**exponent
    return sp.Poly(expression, T, modulus=prime)


assert reduce_rational_polynomial(P, 5) == sp.Poly(
    local_5.as_expr(), T, modulus=5
)

# Both fixed maps recover the same inverse polynomial at their transported
# targets.  Fixed target scalings normalize determinants -5 and -2 to one.
H = T**degree + T**3 - 2 * T**2
c = 1 - degree
weighted = WeightedSeedModel(sp.diff(H, T).subs(T, w), c=c)
weighted_target = (parameter / c, sp.Integer(-1), sp.Integer(1))
assert sp.expand(
    weighted.inverse_polynomial(*weighted_target).subs(w, T) - P.as_expr()
) == 0
weighted_target_one = (
    sp.cancel(-weighted_target[0] / 5),
    weighted_target[1],
    weighted_target[2],
)
assert weighted_target_one == (
    sp.Rational(95231, 1728025),
    -1,
    1,
)

quadratic_seed = H + T
quadratic_target = (sp.Integer(1), sp.Integer(0), -2 * parameter)
assert sp.expand(
    quadratic_seed
    - sp.Rational(1, 2) * quadratic_target[2]
    - P.as_expr()
) == 0
assert quadratic_target == (1, 0, sp.Rational(-190462, 69121))

print("PASS: automatic parameter radii are 2^9, 3^3, and 5")
print("PASS: u=95231/69121 preserves the prescribed ramified local algebras")
print("PASS: the common sextic field has signature (2,2) and is inert at 5")
print("PASS: two fixed stably inequivalent determinant-one maps share this fiber")

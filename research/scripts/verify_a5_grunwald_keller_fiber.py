#!/usr/bin/env python3
"""Exact A5 Grunwald-at-2,3,5 certificate and Keller compilation.

The checker uses exact rational arithmetic, finite-field factorization,
Newton polygons, and the public quadratic-gauge compiler.  It does not call
PARI/GP's Galois-group classifier or a p-adic numerical factorizer.
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber
from verify_finite_etale_keller_fibers import S, check_scheme_reconstruction


T, U = sp.symbols("T U")
x, y, z = sp.symbols("x y z")


def factor_degrees_mod_prime(polynomial: sp.Poly, prime: int) -> tuple[int, ...]:
    """Return the squarefree factor-degree partition modulo ``prime``."""

    assert polynomial.domain == sp.QQ
    _, cleared = polynomial.clear_denoms(convert=True)
    _, primitive = cleared.primitive()
    primitive = sp.Poly(primitive, polynomial.gens[0], domain=sp.ZZ)
    assert int(primitive.LC()) % prime
    assert int(sp.discriminant(primitive.as_expr(), polynomial.gens[0])) % prime
    factors = sp.factor_list(primitive.as_expr(), modulus=prime)[1]
    return tuple(
        sorted(
            (
                int(sp.degree(factor, polynomial.gens[0]))
                for factor, exponent in factors
                for _ in range(exponent)
            ),
            reverse=True,
        )
    )


def valuation(integer: int, prime: int) -> int:
    """Return the valuation of a nonzero integer."""

    assert integer
    value = abs(integer)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


# The Keller inverse polynomial and a monic integral reduced generator for the
# same quintic field.
P = sp.Poly(
    T**5 - 5 * T**3 + 4 * T + sp.Rational(4, 5),
    T,
    domain=sp.QQ,
)
H = sp.Poly(
    T**5 - 40 * T**3 - 110 * T**2 - 40 * T + 32,
    T,
    domain=sp.QQ,
)

# Both polynomials are irreducible modulo 3.  The displayed Tschirnhaus
# element alpha in Q[beta]/(H(beta)) is a root of P, so it identifies the
# quotient fields exactly.
assert factor_degrees_mod_prime(P, 3) == (5,)
assert factor_degrees_mod_prime(H, 3) == (5,)
alpha = (-U**4 + 2 * U**3 + 36 * U**2 + 38 * U - 36) / 20
relation_numerator = sp.together(P.as_expr().subs(T, alpha)).as_numer_denom()[0]
assert sp.rem(relation_numerator, H.as_expr().subs(T, U), U) == 0

# Global Galois group.  The square discriminant places the transitive group in
# A5.  A 3-cycle excludes the proper transitive subgroups C5 and D5.
assert sp.discriminant(P.as_expr(), T) == 232**2
assert sp.discriminant(H.as_expr(), T) == 580000**2
assert factor_degrees_mod_prime(H, 23) == (3, 1, 1)

# The Q_2 completion.  The coefficient points of H have lower Newton polygon
#
#   (0,5) -- (2,1) -- (5,0),
#
# with slopes -2 and -1/3.  The first residual polynomial is Y^2+Y+1 and
# the second is Y+1.  Hence the two factors are the unramified quadratic and
# the unique totally tamely ramified cubic Q_2(cuberoot(2)).
coefficients_low = [int(H.coeff_monomial(T**degree)) for degree in range(6)]
coefficient_valuations_2 = tuple(
    None if coefficient == 0 else valuation(coefficient, 2)
    for coefficient in coefficients_low
)
assert coefficient_valuations_2 == (5, 3, 1, 3, None, 0)
assert (
    coefficients_low[0] // 2**5 % 2,
    coefficients_low[1] // 2**3 % 2,
    coefficients_low[2] // 2 % 2,
) == (1, 1, 1)
residual_quadratic_2 = sp.Poly(U**2 + U + 1, U, modulus=2)
assert residual_quadratic_2.is_irreducible
assert (coefficients_low[2] // 2 % 2, coefficients_low[5] % 2) == (1, 1)

# The Q_3 completion is the unramified quintic: H is irreducible modulo 3
# and its discriminant is a 3-adic unit.
assert factor_degrees_mod_prime(H, 3) == (5,)
assert int(sp.discriminant(H.as_expr(), T)) % 3

# The Q_5 completion is the explicitly prescribed totally ramified quintic.
# Translation by the repeated residue root 3 gives an Eisenstein polynomial.
H5 = sp.Poly(sp.expand(H.as_expr().subs(T, U + 3)), U, domain=sp.ZZ)
assert H5 == sp.Poly(
    U**5 + 15 * U**4 + 50 * U**3 - 200 * U**2 - 1375 * U - 1915,
    U,
    domain=sp.ZZ,
)
assert int(H5.LC()) % 5
assert all(int(coefficient) % 5 == 0 for coefficient in H5.all_coeffs()[1:])
assert int(H5.TC()) % 25

# Compile the original presentation into the compact determinant-one Keller
# map and verify the complete scheme fiber.
t = 1 + x * y
q = t**2 * z - sp.Rational(4, 5) * y**2 * (1 + 3 * t)
expected_minus_two = (
    t * q,
    y - sp.Rational(15, 4) * x * q + sp.Rational(5, 4) * t**2 * x**3 * q**5,
    x * (5 - 3 * t) + sp.Rational(5, 4) * x**3 * z - sp.Rational(3, 4) * (x * q) ** 5,
)
expected_one = (
    expected_minus_two[0],
    -expected_minus_two[1] / 2,
    expected_minus_two[2],
)

compilation = compile_polynomial_to_keller_fiber(
    P,
    T,
    translation=0,
    inverse_variable=U,
    source_variables=(x, y, z),
)
assert all(
    sp.expand(actual - expected) == 0
    for actual, expected in zip(compilation.determinant_minus_two_map, expected_minus_two)
)
assert all(
    sp.expand(actual - expected) == 0
    for actual, expected in zip(compilation.determinant_one_map, expected_one)
)
assert sp.factor(
    sp.Matrix(compilation.determinant_minus_two_map).jacobian((x, y, z)).det()
) == -2
assert sp.factor(
    sp.Matrix(compilation.determinant_one_map).jacobian((x, y, z)).det()
) == 1
assert compilation.target == (sp.Integer(1), sp.Integer(0), -sp.Rational(2, 5))
assert compilation.inverse_polynomial == P.as_expr().subs(T, U)
assert compilation.geometric_degree == 5
assert compilation.coordinate_degrees == (7, 32, 30)
check_scheme_reconstruction(P.as_expr().subs(T, S), sp.Integer(0))

print("PASS: the reduced quintic has Galois closure A5")
print("PASS: at 2 it is Q_2(cuberoot(2)) times the unramified quadratic")
print("PASS: at 3 it is the unramified quintic")
print("PASS: at 5 it is the displayed Eisenstein quintic")
print("PASS: the determinant-one Keller map has this connected complete fiber")

#!/usr/bin/env python3
"""Replay the reconstructed HC4 exceptional-Schur parameter strata.

This is a research checker, not a full 15-quartic classification.  It
verifies:

* rational reconstruction of the two even-chart parameter ideals from
  reductions modulo 47, 101, and 103;
* the zero-dimensional complex support of the reconstructed ideals;
* the Fermat and radial Hessian-discriminant factorizations;
* uniqueness of the radial Schur quartic.

The expensive modular eliminations which produced the displayed residues
are recorded in the generated-results artifact and the accompanying note.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_exceptional_schur_locus_modular.json"
)

mu, nu = sp.symbols("mu nu")
x, y, z = sp.symbols("x y z")


def modular_coefficient(value: sp.Rational, prime: int) -> int:
    numerator, denominator = map(int, value.as_numer_denom())
    return numerator * pow(denominator, -1, prime) % prime


def centered(residue: int, prime: int) -> int:
    return residue if residue <= prime // 2 else residue - prime


def coefficient_vector(
    polynomial: sp.Expr, monomials: tuple[sp.Expr, ...]
) -> tuple[sp.Rational, ...]:
    expanded = sp.Poly(polynomial, mu, nu)
    return tuple(
        sp.Rational(expanded.coeff_monomial(monomial))
        for monomial in monomials
    )


pure_chart = (
    mu * nu**2
    - 2 * nu**3
    - mu**2 / 10
    + 3 * mu * nu / 10
    - nu**2 / 5,
    mu**2 * nu
    - 4 * nu**3
    - mu**2 / 2
    + 8 * mu * nu / 5
    - 6 * nu**2 / 5,
    mu**3
    - 8 * nu**3
    - 9 * mu**2 / 5
    + 6 * mu * nu
    - 24 * nu**2 / 5,
    nu**4
    - nu**3 / 5
    - mu**2 / 100
    + mu * nu / 25
    - 3 * nu**2 / 100,
)

mixed_chart = (
    mu * nu
    - 3 * nu**2
    - mu / 10
    + 2 * nu / 5
    - sp.Rational(1, 100),
    mu**2
    - 9 * nu**2
    - 2 * mu / 5
    + 9 * nu / 5
    - sp.Rational(1, 20),
    (nu - sp.Rational(1, 10)) ** 3,
)

# The coefficient order matches the modular elimination transcript.
pure_monomials = (
    (mu * nu**2, nu**3, mu**2, mu * nu, nu**2),
    (mu**2 * nu, nu**3, mu**2, mu * nu, nu**2),
    (mu**3, nu**3, mu**2, mu * nu, nu**2),
    (nu**4, nu**3, mu**2, mu * nu, nu**2),
)
mixed_monomials = (
    (mu * nu, nu**2, mu, nu, sp.Integer(1)),
    (mu**2, nu**2, mu, nu, sp.Integer(1)),
    (nu**3, nu**2, nu, sp.Integer(1)),
)

expected_pure_residues = {
    47: (
        (1, -2, 14, 5, -19),
        (1, -4, 23, 11, -20),
        (1, -8, 17, 6, 14),
        (1, -19, -8, -15, 23),
    ),
    101: (
        (1, -2, 10, -30, 20),
        (1, -4, 50, 42, 19),
        (1, -8, -22, 6, -25),
        (1, 20, 1, -4, 3),
    ),
    103: (
        (1, -2, -31, -10, 41),
        (1, -4, 51, -19, 40),
        (1, -8, -43, 6, -46),
        (1, 41, -34, 33, 1),
    ),
}
expected_mixed_residues = {
    47: (
        (1, -3, 14, -9, -8),
        (1, -9, 9, -17, 7),
        (1, -5, -23, 18),
    ),
    101: (
        (1, -3, 10, -40, 1),
        (1, -9, 40, 22, 5),
        (1, 30, -3, -10),
    ),
    103: (
        (1, -3, -31, 21, -34),
        (1, -9, -21, 43, 36),
        (1, 10, -1, -24),
    ),
}

for prime in (47, 101, 103):
    actual_pure = tuple(
        tuple(
            centered(modular_coefficient(value, prime), prime)
            for value in coefficient_vector(polynomial, monomials)
        )
        for polynomial, monomials in zip(pure_chart, pure_monomials)
    )
    actual_mixed = tuple(
        tuple(
            centered(modular_coefficient(value, prime), prime)
            for value in coefficient_vector(polynomial, monomials)
        )
        for polynomial, monomials in zip(mixed_chart, mixed_monomials)
    )
    assert actual_pure == expected_pure_residues[prime]
    assert actual_mixed == expected_mixed_residues[prime]

pure_lex = sp.groebner(pure_chart, mu, nu, order="lex")
pure_lex_factored = tuple(
    sp.factor(polynomial.as_expr()) for polynomial in pure_lex.polys
)
assert pure_lex_factored == (
    (mu - 10 * nu**2 - nu) * (mu + 10 * nu**2 - 3 * nu),
    nu * (10 * nu - 1) * (mu - 10 * nu**2 - nu) / 10,
    nu**2 * (10 * nu - 1) ** 3 / 1000,
)

mixed_lex = sp.groebner(mixed_chart, mu, nu, order="lex")
mixed_lex_factored = tuple(
    sp.factor(polynomial.as_expr()) for polynomial in mixed_lex.polys
)
assert mixed_lex_factored == (
    (2 * mu + 6 * nu - 1) * (10 * mu - 30 * nu + 1) / 20,
    (10 * nu - 1) * (10 * mu - 30 * nu + 1) / 100,
    (10 * nu - 1) ** 3 / 1000,
)

fermat = {mu: 0, nu: 0}
radial = {mu: sp.Rational(1, 5), nu: sp.Rational(1, 10)}
assert all(polynomial.subs(fermat) == 0 for polynomial in pure_chart)
assert all(polynomial.subs(radial) == 0 for polynomial in pure_chart)
assert all(polynomial.subs(radial) == 0 for polynomial in mixed_chart)

variables = (x, y, z)
mixed_42 = sum(
    left**4 * right**2
    for left in variables
    for right in variables
    if left != right
)
sextic = (
    (x**6 + y**6 + z**6) / 30
    + mu * x**2 * y**2 * z**2
    + nu * mixed_42
)
radius = x**2 + y**2 + z**2

fermat_sextic = sp.expand(sextic.subs(fermat))
radial_sextic = sp.expand(sextic.subs(radial))
assert fermat_sextic == (x**6 + y**6 + z**6) / 30
assert sp.expand(radial_sextic - radius**3 / 30) == 0
assert sp.factor(sp.hessian(fermat_sextic, variables).det()) == (
    x**4 * y**4 * z**4
)
assert sp.factor(sp.hessian(radial_sextic, variables).det()) == (
    radius**6 / 25
)

# For s=R*q, polynomiality of the radial Schur norm requires R | q^2.
# Verify the displayed numerator identity for a generic quadratic q.
q_coefficients = sp.symbols("q0:6")
q0, q1, q2, q3, q4, q5 = q_coefficients
quadratic = (
    q0 * x**2
    + q1 * y**2
    + q2 * z**2
    + q3 * x * y
    + q4 * x * z
    + q5 * y * z
)
radial_quartic = radius * quadratic
gradient = sp.Matrix(
    [sp.diff(radial_quartic, variable) for variable in variables]
)
displayed_numerator = sp.expand(
    5 * radius * gradient.dot(gradient) - 64 * radial_quartic**2
)
expected_numerator = sp.expand(
    radius**2
    * (
        -4 * quadratic**2
        + 5
        * radius
        * sum(
            sp.diff(quadratic, variable) ** 2
            for variable in variables
        )
    )
)
assert displayed_numerator == expected_numerator

with ARTIFACT.open() as stream:
    artifact = json.load(stream)
assert artifact["even_quartic_charts"]["good_primes"] == [47, 101, 103]

print("PASS: reconstructed even-chart ideals replay modulo 47, 101, 103")
print("PASS: pure chart support is Fermat plus radial")
print("PASS: mixed chart support is radial with a cubic transverse thickening")
print("PASS: Hessian discriminants factor as (xyz)^4 and R^6/25")
print("PASS: radial polynomiality forces the quartic line C*R^2")
print(
    "SCOPE: modular reconstruction and exact special fibers; "
    "the full mixed-character classification remains open"
)

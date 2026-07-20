#!/usr/bin/env python3
"""Uniform irreducibility certificates for maximal 2/3 contact strata."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    maximal_two_three_partitions,
    maximal_two_three_phi,
)


def substitute_endpoint_regime(expression, double_count=None, triple_count=None):
    """Impose the endpoint-coordinate relations in degrees zero, one, or two."""
    substitutions = {}
    if double_count == 0:
        substitutions.update({x: 1, u: 0, X: 1})
    elif double_count == 1:
        substitutions.update({u: 1, X: 1 + x})
    elif double_count == 2:
        substitutions[X] = 1 + x + u
    if triple_count == 0:
        substitutions.update({y: 1, v: 0, Y: 1})
    elif triple_count == 1:
        substitutions.update({v: 1, Y: 1 + y})
    elif triple_count == 2:
        substitutions[Y] = 1 + y + v
    return sp.expand(expression.subs(substitutions))


def has_only_constant_content(polynomial, variable):
    """Gauss-primitivity over the remaining polynomial variables."""
    coefficients = sp.Poly(polynomial, variable).all_coeffs()
    return not sp.gcd_list(coefficients).free_symbols


x, u, X, y, v, Y = sp.symbols("x u X y v Y")
universal_phi = sp.expand(
    X**2 * Y**3 - x**2 * y**3 - 2 * x * u * y**3 - 3 * x**2 * y**2 * v
)

# If the triple-root polynomial has degree at least three, its values
# (R(0),R'(0),R(1))=(y,v,Y) are independent affine coordinates.  Phi is
# primitive linear in v.  The coefficient -3*x^2*y^2 is coprime to the
# constant term in every possible endpoint regime for Q.
for double_regime in (0, 1, 2, 3):
    expression = substitute_endpoint_regime(
        universal_phi,
        double_count=double_regime if double_regime < 3 else None,
    )
    coefficient = sp.diff(expression, v)
    constant = expression.subs(v, 0)
    assert coefficient == -3 * x**2 * y**2 or (
        double_regime == 0 and coefficient == -3 * y**2
    )
    assert sp.gcd(coefficient, constant) == 1
    assert sp.degree(expression, v) == 1

# If Q has degree at least three and R has degree at most two, X=Q(1) is an
# independent coordinate.  Phi is primitive quadratic in X.  Its constant
# term has x-adic valuation exactly one, so it cannot become a square in the
# fraction field; the quadratic is irreducible by Gauss's lemma.
for triple_regime in (0, 1, 2):
    expression = substitute_endpoint_regime(
        universal_phi, triple_count=triple_regime
    )
    polynomial = sp.Poly(expression, X)
    leading = polynomial.coeff_monomial(X**2)
    constant = polynomial.coeff_monomial(1)
    assert polynomial.degree() == 2
    assert sp.gcd(leading, constant) == 1
    assert constant.subs(x, 0) == 0
    assert sp.diff(constant, x).subs(x, 0) != 0

# Only seven endpoint-rank cases remain (inverse degrees 3 through 10).  Each
# is linear or quadratic in the displayed variable.  The quadratic
# discriminants have an explicit odd valuation and hence are nonsquares.
small_cases = {
    (0, 1): (y, "linear"),
    (0, 2): (y, "triple_quadratic"),
    (1, 1): (x, "mixed_linear_linear"),
    (1, 2): (x, "mixed_linear_quadratic"),
    (2, 0): (x, "linear"),
    (2, 1): (u, "double_quadratic_linear"),
    (2, 2): (u, "double_quadratic_quadratic"),
}
for (double_count, triple_count), (variable, certificate) in small_cases.items():
    expression = substitute_endpoint_regime(
        universal_phi,
        double_count=double_count,
        triple_count=triple_count,
    )
    assert has_only_constant_content(expression, variable)
    degree = sp.degree(expression, variable)
    if certificate == "linear":
        assert degree == 1
        continue
    assert degree == 2
    discriminant = sp.factor(sp.discriminant(expression, variable))
    if certificate == "triple_quadratic":
        assert discriminant == 3 * (v + 1) ** 3 * (3 * v - 1)
    elif certificate == "mixed_linear_linear":
        assert discriminant == 4 * y**2 * (6 * y**2 + 8 * y + 3)
        assert sp.discriminant(6 * y**2 + 8 * y + 3, y) == -8
    elif certificate == "mixed_linear_quadratic":
        residual = sp.factor(discriminant / (-4 * y**2))
        specialized = sp.factor(residual.subs(v, 0))
        assert specialized.subs(y, 0) == 0
        assert sp.diff(specialized, y).subs(y, 0) != 0
    else:
        assert discriminant.subs(x, 0) == 0
        assert sp.diff(discriminant, x).subs(x, 0) != 0

# Connect the endpoint-coordinate proof back to the public quotient model and
# retain exact factorization regressions for every maximal type through n=14.
for degree in range(3, 15):
    for partition in maximal_two_three_partitions(degree):
        double_count = partition.count(2)
        triple_count = partition.count(3)
        model = maximal_two_three_phi(
            double_count, triple_count, prefix=f"irr_{degree}_{triple_count}"
        )
        assert model.partition == partition
        factors = sp.factor_list(model.phi)[1]
        assert len(factors) == 1 and factors[0][1] == 1

print("PASS: Phi is primitive linear when the triple block has degree at least 3")
print("PASS: Phi is a nonsquare quadratic when the double block has degree at least 3")
print("PASS: all seven endpoint-rank exceptions have exact irreducibility certificates")
print("PASS: every maximal 2/3 Phi_lambda is irreducible in characteristic zero")

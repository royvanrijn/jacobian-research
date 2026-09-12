#!/usr/bin/env python3
"""Verify the exceptional radial Schur pair on the HC4 two-parameter surface.

This checker certifies only the degree-fourteen Schur divisibility identity.
It does not assert that the pair extends through the lower collision faces.
"""

from __future__ import annotations

import sympy as sp


x, y, z, mu, nu = sp.symbols("x y z mu nu")
variables = (x, y, z)

radius = x**2 + y**2 + z**2
pair_sum = x**2 * y**2 + x**2 * z**2 + y**2 * z**2
triple_product = x**2 * y**2 * z**2
mixed_42 = sum(
    left**4 * right**2
    for left in variables
    for right in variables
    if left != right
)
sextic_family = (
    (x**6 + y**6 + z**6) / 30
    + mu * x**2 * y**2 * z**2
    + nu * mixed_42
)

exceptional_parameters = {mu: sp.Rational(1, 5), nu: sp.Rational(1, 10)}
exceptional_sextic = sp.expand(sextic_family.subs(exceptional_parameters))
quartic = radius**2
quadratic_quotient = 16 * radius

radial_difference = sp.Poly(
    sextic_family - radius**3 / 30, x, y, z
)
assert radial_difference.coeff_monomial(
    x**2 * y**2 * z**2
) == mu - sp.Rational(1, 5)
assert radial_difference.coeff_monomial(
    x**4 * y**2
) == nu - sp.Rational(1, 10)
assert sp.expand(exceptional_sextic - radius**3 / 30) == 0

hessian = sp.hessian(exceptional_sextic, variables)
hessian_determinant = sp.factor(hessian.det())
assert hessian_determinant == radius**6 / 25

quartic_gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in variables]
)
schur_numerator = sp.expand(
    (quartic_gradient.T * hessian.adjugate() * quartic_gradient)[0]
)
assert sp.expand(
    schur_numerator - hessian_determinant * quadratic_quotient
) == 0
assert quartic != 0
assert hessian_determinant != 0

# In the invariant coordinates
#
#   h6 = R^3/30 + A*R*P2 + B*P3,
#
# check that the fixed radial quartic R^2 has a polynomial Schur norm only
# at A=B=0.  Independent sign changes and permutations force its quadratic
# quotient to be a scalar multiple q*R.
A, B, q = sp.symbols("A B q")
invariant_sextic = (
    radius**3 / 30
    + A * radius * pair_sum
    + B * triple_product
)
invariant_hessian = sp.hessian(invariant_sextic, variables)
invariant_remainder = sp.expand(
    (
        quartic_gradient.T
        * invariant_hessian.adjugate()
        * quartic_gradient
    )[0]
    - invariant_hessian.det() * q * radius
)
invariant_equations = list(
    dict.fromkeys(
        sp.Poly(invariant_remainder, *variables).coeffs()
    )
)
assert len(invariant_equations) == 8
invariant_basis = sp.groebner(
    invariant_equations, q, A, B, order="lex"
)
assert [
    sp.expand(polynomial.as_expr())
    for polynomial in invariant_basis.polys
] == [q - 16, A, B]

print("PASS: (mu,nu)=(1/5,1/10) gives h6=(x^2+y^2+z^2)^3/30")
print("PASS: det(Hess(h6))=(x^2+y^2+z^2)^6/25 is nonzero")
print("PASS: s4=(x^2+y^2+z^2)^2 has Schur quotient 16*(x^2+y^2+z^2)")
print("PASS: this radial quartic occurs only at the radial parameter point")
print("SCOPE: exceptional Schur pair only; lower collision faces are not checked")

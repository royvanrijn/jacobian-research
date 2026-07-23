#!/usr/bin/env python3
"""Exact audit of the Danielewski Poisson-contraction frontier."""

from __future__ import annotations

from itertools import product

import sympy as sp


a, b, c, w, x = sp.symbols("a b c w x")
alpha, beta, gamma = sp.symbols("alpha beta gamma", nonzero=True)

P = x * (x**2 + 1)
P_prime = sp.diff(P, x)
relation = c * w - P


def bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    """Ambient formula for the quotient Poisson bracket."""

    return sp.expand(
        c
        * (
            sp.diff(first, c) * sp.diff(second, x)
            - sp.diff(first, x) * sp.diff(second, c)
        )
        + P_prime
        * (
            sp.diff(first, c) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, c)
        )
        + w
        * (
            sp.diff(first, x) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, x)
        )
    )


# The bracket descends to the hypersurface and has the stated generators.
assert bracket(c, x) == c
assert bracket(c, w) == P_prime
assert bracket(x, w) == w
for generator in (c, x, w):
    assert sp.expand(bracket(relation, generator)) == 0
print("PASS: the polynomial bracket descends to cw=P(x)")


# Pullback of dF wedge dG is bracket(F,G) times dc wedge dx/c.
# Test the identity after eliminating w=P/c for generic polynomials.
F_test = c**2 + c * x + w + x**2
G_test = c * w + c + x * w
F_chart = sp.cancel(F_test.subs(w, P / c))
G_chart = sp.cancel(G_test.subs(w, P / c))
chart_coefficient = sp.cancel(
    c
    * (
        sp.diff(F_chart, c) * sp.diff(G_chart, x)
        - sp.diff(F_chart, x) * sp.diff(G_chart, c)
    )
)
bracket_chart = sp.cancel(bracket(F_test, G_test).subs(w, P / c))
assert sp.cancel(chart_coefficient - bracket_chart) == 0

# The source residue pullback contributes the remaining constant sign.
source_x = b * c
jacobian_acx = sp.Matrix((a, c, source_x)).jacobian((a, b, c)).det()
assert jacobian_acx == -c
assert sp.cancel(jacobian_acx / c) == -1
print("PASS: a unit-bracket pair produces an ordinary determinant -1 map")


# Hamiltonian vector field of a general affine-linear first coordinate.
F_linear = alpha * c + beta * x + gamma * w
field = tuple(sp.expand(bracket(F_linear, z)) for z in (c, x, w))
expected_field = (
    -beta * c - gamma * P_prime,
    alpha * c - gamma * w,
    alpha * P_prime + beta * w,
)
assert all(
    sp.expand(component - expected) == 0
    for component, expected in zip(field, expected_field)
)

# In the fully mixed beta != 0 case, the advertised critical point works
# modulo H=alpha*gamma*(P')^2-beta^2*P.
mixed_substitution = {
    c: -gamma * P_prime / beta,
    w: -alpha * P_prime / beta,
}
critical_polynomial = alpha * gamma * P_prime**2 - beta**2 * P
for component in field:
    assert sp.cancel(component.subs(mixed_substitution)) == 0
assert sp.cancel(
    relation.subs(mixed_substitution)
    - critical_polynomial / beta**2
) == 0
assert sp.degree(critical_polynomial, x) == 4
print("PASS: the mixed affine-linear critical-point formula is exact")


# Edge cases at the simple root x=0.
root_derivative = P_prime.subs(x, 0)
assert root_derivative == 1

gamma_zero_point = {
    x: 0,
    c: 0,
    w: -alpha * root_derivative / beta,
    gamma: 0,
}
alpha_zero_point = {
    x: 0,
    c: -gamma * root_derivative / beta,
    w: 0,
    alpha: 0,
}
for substitution in (gamma_zero_point, alpha_zero_point):
    assert sp.cancel(relation.subs(substitution)) == 0
    assert all(sp.cancel(component.subs(substitution)) == 0 for component in field)
print("PASS: all beta-nonzero affine-linear edge cases have critical points")


# Pure c is locally nilpotent: delta_c(x)=c and
# delta_c^j(w)=c^(j-1) P^(j)(x), ending after deg(P)+1 steps.
def delta_c(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(bracket(c, polynomial))


assert delta_c(c) == 0
assert delta_c(x) == c
assert delta_c(w) == P_prime
iterate = w
for _ in range(sp.degree(P, x) + 1):
    iterate = delta_c(iterate)
assert iterate == 0

# delta_x is the weight derivation for wt(c,x,w)=(-1,0,1).
assert bracket(x, c) == -c
assert bracket(x, x) == 0
assert bracket(x, w) == w
for c_power, x_power, w_power in product(range(4), repeat=3):
    monomial = c**c_power * x**x_power * w**w_power
    weight = w_power - c_power
    assert sp.expand(bracket(x, monomial) - weight * monomial) == 0
print("PASS: the pure-coordinate LND and grading obstructions are exact")

# On a vanishing sphere joining 0 to i, take x(t)=i*t and
# c=r(t)*exp(i*theta).  The dr/r term wedges to zero with dx, leaving
# i*dtheta wedge dx and the nonzero period 2*pi*i*(i-0)=-2*pi.
t, theta = sp.symbols("t theta", real=True)
x_path = sp.I * t
period_density = sp.I * sp.diff(x_path, t)
period = sp.integrate(
    sp.integrate(period_density, (theta, 0, 2 * sp.pi)),
    (t, 0, 1),
)
assert sp.simplify(period + 2 * sp.pi) == 0
assert period != 0
print("PASS: the residue form has a nonzero vanishing-sphere period")


# For general A,F,G in B[a], expansion of dA wedge dF wedge dG relative
# to da wedge omega gives the coupled coefficient below.  Check it on a
# nontrivial polynomial triple by an ordinary chart Jacobian.
A_test = a**2 + c * x + w
F_coupled = a * c + x**2 + w
G_coupled = a * x + c**2 + x * w
coupled_coefficient = sp.expand(
    sp.diff(A_test, a) * bracket(F_coupled, G_coupled)
    - sp.diff(F_coupled, a) * bracket(A_test, G_coupled)
    + sp.diff(G_coupled, a) * bracket(A_test, F_coupled)
)
chart_coordinates = tuple(
    sp.cancel(item.subs(w, P / c))
    for item in (A_test, F_coupled, G_coupled)
)
chart_jacobian = sp.cancel(
    sp.Matrix(chart_coordinates).jacobian((a, c, x)).det() * c
)
assert sp.cancel(
    chart_jacobian - coupled_coefficient.subs(w, P / c)
) == 0
print("PASS: the coupled three-coordinate volume equation is exact")


# Exhaustive small-height regression.  Pure c and pure w are the only
# nonconstant affine-linear Hamiltonians without a critical point; the
# slice/fiber argument above excludes those two cases.
critical_cases = 0
critical_free_cases: set[tuple[int, int, int]] = set()
for alpha_value, beta_value, gamma_value in product(range(-2, 3), repeat=3):
    coefficients = (alpha_value, beta_value, gamma_value)
    if coefficients == (0, 0, 0):
        continue
    specialized = [
        component.subs(
            {
                alpha: alpha_value,
                beta: beta_value,
                gamma: gamma_value,
            }
        )
        for component in field
    ]
    ideal_basis = sp.groebner(
        [relation, *specialized],
        c,
        w,
        x,
        order="lex",
    )
    if ideal_basis.polys == [sp.Poly(1, c, w, x)]:
        critical_free_cases.add(coefficients)
    else:
        critical_cases += 1

expected_pure = {
    (value, 0, 0) for value in range(-2, 3) if value
} | {
    (0, 0, value) for value in range(-2, 3) if value
}
assert critical_free_cases == expected_pure
assert critical_cases == 116
print("PASS: all 116 small-height non-pure linear forms have a critical point")
print("PASS Danielewski Poisson-contraction frontier")

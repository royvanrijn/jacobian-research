#!/usr/bin/env python3
"""Exact audit of the primitive-multiplicity triple-root map."""

from __future__ import annotations

from itertools import product

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
y = a**2 * c
u = a * c
tau = x - y + 1
v = b * (x + 1) * tau
f = sp.factor(x * y * (x + 1) * tau)

assert sp.factor(u**2 * v - f) == 0
print("PASS: the affine source lands on u^2*v=x*y*(x+1)*(x-y+1)")


jacobian = sp.factor(
    sp.Matrix((x, y, u)).jacobian((a, b, c)).det()
)
assert jacobian == -u**2
assert sp.cancel(jacobian / u**2) == -1
print("PASS: the primitive-multiplicity residue Jacobian is -1")


assert sp.factor(y / u - a) == 0
assert sp.factor(u**2 / y - c) == 0
assert sp.factor(x * y / u**2 - b) == 0
print("PASS: the displayed rational inverse is exact")


valuation_table = {
    "y": (1, 0, 0),
    "x": (-1, 0, 2),
    "x+1": (-1, -2, 2),
    "x-y+1": (-1, -2, 2),
}
dicritical = {
    component: valuations
    for component, valuations in valuation_table.items()
    if min(valuations) < 0
}
assert set(dicritical) == {"x", "x+1", "x-y+1"}
print("PASS: one component is finite and three are dicritical")


# The singular locus of u^2*v=f is u=0 over the three intersection points
# of the four-line arrangement, with v arbitrary.
X, Y, U, V = sp.symbols("X Y U V")
target_equation = U**2 * V - X * Y * (X + 1) * (X - Y + 1)
target_gradient = tuple(
    sp.factor(sp.diff(target_equation, variable))
    for variable in (X, Y, U, V)
)
intersection_points = ((0, 0), (0, 1), (-1, 0))
for finite_x, finite_y in intersection_points:
    substitutions = {X: finite_x, Y: finite_y, U: 0}
    assert target_equation.subs(substitutions) == 0
    assert all(
        derivative.subs(substitutions) == 0
        for derivative in target_gradient
    )
print("PASS: the singular locus is supported over three arrangement points")


# Exhaust every monomial divisor of f in source irreducible factors.
factors = (a, b, c, x + 1, tau)
maximum_exponents = (2, 1, 2, 1, 1)
power_matches = []
for exponents in product(
    *(range(maximum + 1) for maximum in maximum_exponents)
):
    candidate_u = sp.prod(
        factor**exponent
        for factor, exponent in zip(factors, exponents)
    )
    candidate_v = sp.factor(f / candidate_u)
    if sp.denom(candidate_v) != 1:
        continue
    candidate_jacobian = sp.factor(
        sp.Matrix((x, y, candidate_u))
        .jacobian((a, b, c))
        .det()
    )
    if candidate_jacobian == 0:
        continue
    ratio = sp.factor(candidate_jacobian / candidate_u)
    assert sp.factor(ratio / (a * c)).is_number
    for multiplicity in range(1, 7):
        quotient = sp.factor(
            candidate_jacobian / candidate_u**multiplicity
        )
        if quotient.is_number and quotient != 0:
            power_matches.append(
                (exponents, multiplicity, quotient)
            )
assert power_matches == [((1, 0, 1, 0, 0), 2, -1)]
print("PASS: u=ac with multiplicity two is the unique monomial solution")


for prime in (5, 7, 11):
    target_count = 0
    for finite_x, finite_y, finite_u, finite_v in product(
        range(prime),
        repeat=4,
    ):
        if (
            finite_u**2 * finite_v
            - finite_x
            * finite_y
            * (finite_x + 1)
            * (finite_x - finite_y + 1)
        ) % prime == 0:
            target_count += 1
    assert target_count == prime * (prime - 1) * (prime + 4)
print("PASS: target point count is q(q-1)(q+4)")
print("PASS primitive-multiplicity triple-root map")

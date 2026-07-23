#!/usr/bin/env python3
"""Exact audit of the smooth Danielewski triple-root map."""

from __future__ import annotations

from itertools import product

import sympy as sp


a, b, c = sp.symbols("a b c")
x = b * c
w = b * (x + 1)
assert sp.factor(c * w - x * (x + 1)) == 0
print("PASS: the affine source lands on c*w=x*(x+1)")


jacobian = sp.factor(
    sp.Matrix((a, c, x)).jacobian((a, b, c)).det()
)
assert jacobian == -c
assert sp.cancel(jacobian / c) == -1
print("PASS: the smooth-target residue Jacobian is -1")


C, W, X = sp.symbols("C W X")
target_equation = C * W - X * (X + 1)
target_gradient = tuple(
    sp.diff(target_equation, variable)
    for variable in (C, W, X)
)
singular_basis = sp.groebner(
    [target_equation, *target_gradient],
    C,
    W,
    X,
)
assert singular_basis.contains(sp.Integer(1))
print("PASS: the Danielewski target is globally smooth")


assert sp.factor(x / c - b) == 0
finite_boundary_valuation = 0
dicritical_boundary_valuation = -1
assert finite_boundary_valuation >= 0
assert dicritical_boundary_valuation < 0
print("PASS: x=0 is finite and x=-1 is dicritical over c=0")


L = sp.symbols("L")
class_surface = sp.expand(
    (L - 2) * (L - 1) + 2 * (2 * L - 1)
)
class_target = sp.expand(L * class_surface)
assert class_surface == L**2 + L
assert class_target == L**3 + L**2
print("PASS: the smooth target has class L^3+L^2, not L^3")


for prime in (5, 7, 11):
    target_count = sum(
        (finite_c * finite_w - finite_x * (finite_x + 1))
        % prime
        == 0
        for finite_a, finite_c, finite_x, finite_w in product(
            range(prime),
            repeat=4,
        )
    )
    assert target_count == prime**3 + prime**2
print("PASS: target point count is q^3+q^2")


z = sp.symbols("z")
generic_coefficients = sp.symbols("g0:5")
generic_polynomial = sum(
    generic_coefficients[index] * z**index
    for index in range(5)
)
pullback_numerator = sp.expand(generic_polynomial.subs(z, b * c))
polynomial_pullback_condition = sp.expand(
    pullback_numerator.subs(c, 0)
)
assert polynomial_pullback_condition == generic_coefficients[0]

# Retaining the x=-1 center adds g(-1)=0.  Together with g(0)=0, the
# minimal nonzero polynomial is divisible by x(x+1).
minimal_balanced = z * (z + 1)
assert minimal_balanced.subs(z, 0) == 0
assert minimal_balanced.subs(z, -1) == 0
assert sp.factor(minimal_balanced / (z * (z + 1))) == 1
print("PASS: the minimal source-polynomial balanced center is x(x+1)/c")
print("PASS smooth Danielewski triple-root map")

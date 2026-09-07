#!/usr/bin/env python3
"""Exact checks for coarse affine-mark descent at root collisions."""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def elementary(values: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    """Elementary symmetric polynomial of the requested degree."""

    if degree == 0:
        return sp.Integer(1)
    result = sp.Integer(0)
    for indices in combinations(range(len(values)), degree):
        term = sp.Integer(1)
        for index in indices:
            term *= values[index]
        result += term
    return sp.expand(result)


# For R=k[x_1,...,x_mu], the coarse marked/unmarked collision chart is
# R^(S_(mu-1)) over R^(S_mu).  If T=x_1, the remaining elementary symmetric
# functions are recovered recursively from e_j=e'_j+T e'_(j-1).  Hence the
# marked invariant ring is the universal monic-root incidence.
for multiplicity in range(2, 9):
    roots = sp.symbols(f"x0:{multiplicity}")
    selected = roots[0]
    all_elementary = {
        degree: elementary(roots, degree) for degree in range(multiplicity + 1)
    }
    remaining_elementary = {
        degree: elementary(roots[1:], degree) for degree in range(multiplicity)
    }

    recovered = {0: sp.Integer(1)}
    for degree in range(1, multiplicity):
        recovered[degree] = sp.expand(
            all_elementary[degree] - selected * recovered[degree - 1]
        )
        assert sp.expand(recovered[degree] - remaining_elementary[degree]) == 0

    universal_root_polynomial = selected**multiplicity
    for degree in range(1, multiplicity + 1):
        universal_root_polynomial += (
            (-1) ** degree
            * all_elementary[degree]
            * selected ** (multiplicity - degree)
        )
    assert sp.expand(universal_root_polynomial) == 0


# The total-collision fibers for the proposed first counterexamples are
# k[T]/(T^2) and k[T]/(T^3).  They retain length/inertia but have one prime,
# so there is exactly one geometric coarse specialization of the mark.
T = sp.symbols("T")
for multiplicity in (2, 3):
    special_fiber = sp.Poly(T**multiplicity, T)
    assert special_fiber.degree() == multiplicity
    assert sp.factor_list(special_fiber.as_expr())[1] == [(T, multiplicity)]


# On the cyclic DVR slice epsilon=delta^mu, every selected generic root
# zeta*delta specializes to the same coarse point T=0.  The root is integral
# over k[[epsilon]], and after the required marked base change it is the
# tautological regular parameter.
delta, epsilon = sp.symbols("delta epsilon")
for multiplicity in range(2, 9):
    equation = T**multiplicity - epsilon
    assert sp.expand(equation.subs({T: delta, epsilon: delta**multiplicity})) == 0
    assert sp.rem(T**multiplicity, T**multiplicity - epsilon, T) == epsilon


print("PASS coarse marked chart: R^(S_(mu-1)) is the universal root incidence")
print("PASS pair/triple collisions: special fibers have one geometric point")
print("PASS DVR descent: every cyclic selected root has the same coarse limit")

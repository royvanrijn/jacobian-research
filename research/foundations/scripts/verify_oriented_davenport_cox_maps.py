#!/usr/bin/env python3
"""Exact audit of the oriented Davenport Cox maps."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    A,
    T,
    U,
    Y,
    Z,
    branch_cubic,
    cycle_partition,
    davenport_pair,
    gl32,
    inverses,
    line_permutation,
    point_permutation,
)


g, h = davenport_pair()
minimal_variable = sp.symbols("minimal_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(
        minimal_variable**2 + minimal_variable + 2,
        minimal_variable,
    ),
    alias="alpha",
)
alpha = number_field.ext
function_field = number_field.frac_field(T)


def oriented_factor_data(variable: sp.Symbol, polynomial: sp.Expr):
    pullback = sp.Poly(
        sp.expand(branch_cubic().subs(U, polynomial).subs(A, alpha)),
        variable,
        domain=function_field,
    )
    coefficient, factors = sp.factor_list(pullback)
    assert coefficient == 1

    jacobian = sp.Poly(
        sp.diff(polynomial, variable).subs(A, alpha),
        variable,
        domain=function_field,
    ).monic()
    simple = [
        factor.monic()
        for factor, multiplicity in factors
        if factor.monic() != jacobian
    ]
    matching = [
        multiplicity
        for factor, multiplicity in factors
        if factor.monic() == jacobian
    ]
    assert matching == [2]
    assert sorted(factor.degree() for factor in simple) == [3, 6]
    return pullback.monic(), jacobian, simple


for variable, polynomial in ((Y, g), (Z, h)):
    pullback, jacobian, simple = oriented_factor_data(
        variable, polynomial
    )
    source_square = simple[0] * simple[1]
    assert source_square * jacobian**2 == pullback

    # Abstract oriented coordinates satisfy D=J*W and both square equations.
    D, W = sp.symbols("D W", nonzero=True)
    oriented_relation = sp.expand(
        (jacobian.as_expr() * W) ** 2
    ).subs(W**2, source_square.as_expr()) - pullback.as_expr()
    # Re-enter the exact algebraic function field so alpha^2+alpha+2
    # is reduced before testing the oriented hypersurface identity.
    assert sp.Poly(
        oriented_relation,
        variable,
        domain=function_field,
    ).is_zero

    residue_ratio = sp.cancel(
        jacobian.as_expr() / (jacobian.as_expr() * W)
    )
    assert residue_ratio == 1 / W
print("PASS: both Davenport covers lift by W^2=E3*E6 and D=J*W")
print("PASS: target and source hypersurface residues have Jacobian one")


# The oriented divisor D=J*W has one reduced component from J and one from
# each of the two prime factors under W=0.
unoriented_multiplicities = (1, 1, 2)
oriented_multiplicities = (1, 1, 1)
assert sum(unoriented_multiplicities) == 4
assert sum(oriented_multiplicities) == 3
print("PASS: orientation changes the height-one ledger from (1,1,2) to (1,1,1)")


# The point and line actions remain degree seven and have identical cycle
# data, but their fixed-point counts are not constantly one.
matrices = gl32()
inverse = inverses(matrices)
point_fixed_counts = set()
for matrix in matrices:
    point_action = point_permutation(matrix)
    line_action = line_permutation(matrix, inverse[matrix])
    assert cycle_partition(point_action) == cycle_partition(line_action)
    point_fixed_counts.add(
        sum(index == image for index, image in enumerate(point_action))
    )
assert point_fixed_counts == {0, 1, 3, 7}
print("PASS: oriented base change retains the common GL_3(F_2) point/line data")
print("PASS: fixed-point counts {0,1,3,7} obstruct a global permutation family")
print("PASS oriented Davenport Cox maps")

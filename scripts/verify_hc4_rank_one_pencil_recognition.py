#!/usr/bin/env python3
"""Verify the constant-null-covector frontend for HC4 pencil recognition."""

from __future__ import annotations

import sympy as sp

from hc4_rank_one_pencil_recognition import (
    certify_constant_null_covector,
    constant_null_covector_system,
    has_rank_one_pencil_over_algebraic_closure,
    projective_null_covector_charts,
)


def symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    entries = {
        (row, column): sp.Symbol(f"{prefix}{row}{column}")
        for row in range(size)
        for column in range(row, size)
    }
    return sp.Matrix(
        size,
        size,
        lambda row, column: entries[min(row, column), max(row, column)],
    )


# Universal rank-one determinant and square-zero identities.
H = symmetric_matrix("h", 4)
ell = sp.Matrix(sp.symbols("l0:4"))
s = sp.symbols("s")
T = ell * ell.T
det_H = H.det(method="domain-ge")
adj_H = H.adjugate(method="domain-ge")
Q = sp.expand((ell.T * adj_H * ell)[0])
assert sp.expand((H + s * T).det(method="domain-ge") - det_H - s * Q) == 0

# If N=H^(-1)T, then det(H)*N=adj(H)*T.  Clearing denominators gives
# N^2=(Q/det(H))*N, so Q=0 makes the relative endomorphism square-zero.
N_numerator = adj_H * T
assert (N_numerator**2 - Q * N_numerator).applyfunc(sp.expand) == sp.zeros(4)


# A nonlinear cotangent control has a whole constant null source plane.
x, y, z, w = variables = sp.symbols("x y z w")
cotangent = z * (x + y**2) + w * y
cotangent_system = constant_null_covector_system(cotangent, variables)
assert cotangent_system.hessian_determinant == 1
assert certify_constant_null_covector(cotangent_system, (1, 0, 0, 0))
assert certify_constant_null_covector(cotangent_system, (0, 1, 0, 0))
cotangent_charts = projective_null_covector_charts(cotangent_system)
assert not cotangent_charts[0].is_empty
assert not cotangent_charts[1].is_empty

# The criterion is invariant under an oblique unimodular linear rechart.
linear_map = sp.Matrix(
    [
        [1, 1, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 1],
        [0, 0, 0, 1],
    ]
)
mapped_variables = linear_map * sp.Matrix(variables)
oblique = sp.expand(
    cotangent.subs(
        dict(zip(variables, mapped_variables, strict=True)),
        simultaneous=True,
    )
)
oblique_system = constant_null_covector_system(oblique, variables)
transported_covector = tuple(linear_map.T * sp.Matrix([1, 0, 0, 0]))
assert oblique_system.hessian_determinant == 1
assert certify_constant_null_covector(oblique_system, transported_covector)

# Even the nonsplit-looking constant quadratic control has a projective null
# quadric after scalar extension, as it must over an algebraically closed field.
quadratic = sum(variable**2 for variable in variables) / 2
quadratic_system = constant_null_covector_system(quadratic, variables)
assert quadratic_system.hessian_determinant == 1
assert has_rank_one_pencil_over_algebraic_closure(quadratic_system)

# Admission refuses to draw an HC4 conclusion before the constant-Hessian gate.
nonconstant = x**3 + y**2 / 2 + z**2 / 2 + w**2 / 2
nonconstant_system = constant_null_covector_system(nonconstant, variables)
assert not nonconstant_system.has_constant_nonzero_hessian_determinant
try:
    projective_null_covector_charts(nonconstant_system)
except ValueError:
    pass
else:
    raise AssertionError("nonconstant Hessian determinant passed admission")


print("PASS: rank-one determinant face is ell^T*adj(H)*ell")
print("PASS: a null metric covector gives a square-zero relative pencil")
print("PASS: nonlinear cotangent controls contain a constant null two-plane")
print("PASS: the recognition scheme survives an oblique linear rechart")
print("PASS: four exact projective charts decide algebraic-closure admission")
print("SCOPE: nonempty null-covector scheme is sufficient, not proved universal")


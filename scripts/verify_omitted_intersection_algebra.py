#!/usr/bin/env python3
"""Exact certificate for primitive merger faces and the first failure."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.omitted_intersection_algebra import (  # noqa: E402
    fixed_weight_branches,
    primitive_merger_hypergraph,
)


# Degree twelve: the two fixed-weight branches over two sixfold roots.
degree12 = primitive_merger_hypergraph(fixed_weight_branches(2, 1))
assert degree12.transfer_cluster_count == 2
assert degree12.cycle_rank == 0
assert degree12.predicted_transverse_length == 4
assert degree12.predicted_hilbert_vector == (1, 2, 1)

# Degree eighteen: singleton sheets form K_3.  Its one graph cycle does not
# remove any of the three independent Hensel-cluster variables.
degree18 = primitive_merger_hypergraph(fixed_weight_branches(3, 1))
assert len(degree18.branches) == 3
assert len(degree18.merger_edges) == 3
assert degree18.component_count == 1
assert degree18.cycle_rank == 1
assert degree18.transfer_cluster_count == 3
assert degree18.predicted_transverse_length == 8
assert degree18.predicted_hilbert_vector == (1, 3, 3, 1)

for left in range(3):
    for right in range(left + 1, 3):
        pair = primitive_merger_hypergraph(
            (degree18.branches[left], degree18.branches[right])
        )
        assert pair.transfer_cluster_count == 2
        assert pair.predicted_transverse_length == 4

# The next singleton face has merger graph K_4 and cycle rank three, but the
# strong Hensel product has four primitive factors and length sixteen.
degree24 = primitive_merger_hypergraph(fixed_weight_branches(4, 1))
assert len(degree24.branches) == 4
assert len(degree24.merger_edges) == 6
assert degree24.cycle_rank == 3
assert degree24.transfer_cluster_count == 4
assert degree24.predicted_transverse_length == 16
assert degree24.predicted_hilbert_vector == (1, 4, 6, 4, 1)


# Universal primitive block Q^2=T^3.
Z = sp.symbols("Z")
u, v, a, b = sp.symbols("u v a b")
Q = Z**3 + u * Z + v
T = Z**2 + a * Z + b
primitive_coefficients = sp.Poly(sp.expand(Q**2 - T**3), Z).all_coeffs()
primitive = sp.groebner(primitive_coefficients, v, u, b, a, order="lex")
primitive_expected = (a, v, 2 * u - 3 * b, b**2)
assert all(primitive.reduce(relation)[1] == 0 for relation in primitive_expected)
expected_primitive = sp.groebner(primitive_expected, v, u, b, a, order="lex")
assert all(
    expected_primitive.reduce(relation)[1] == 0
    for relation in primitive_coefficients
)


# First merged block.  Eliminate the six nonleading coefficients of U from
# U^2=V^3, then adapt V to V=S^2+XZ+Y.
v3, v2, v1, v0 = sp.symbols("v3 v2 v1 v0")
u5, u4, u3, u2, u1, u0 = sp.symbols("u5 u4 u3 u2 u1 u0")
V = Z**4 + v3 * Z**3 + v2 * Z**2 + v1 * Z + v0
U = Z**6 + u5 * Z**5 + u4 * Z**4 + u3 * Z**3 + u2 * Z**2 + u1 * Z + u0
difference = sp.Poly(sp.expand(U**2 - V**3), Z)

solution = {}
for degree, variable in zip(range(11, 5, -1), (u5, u4, u3, u2, u1, u0)):
    equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    solution[variable] = sp.factor(roots[0])

remaining = [
    sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
    for degree in range(5, -1, -1)
]
p, q, X, Y = sp.symbols("p q X Y")
adapted = {
    v3: 2 * p,
    v2: p**2 + 2 * q,
    v1: 2 * p * q + X,
    v0: q**2 + Y,
}
remaining = [sp.factor(equation.subs(adapted)) for equation in remaining]
candidate = (X**3, 2 * X * Y - p * X**2, Y**2 - q * X**2)

# Degrees >=2 are the actual affine-difference equations.  They already imply
# the discarded linear and constant equations in this local block.
affine = sp.groebner(remaining[:4], Y, X, q, p, order="lex", domain=sp.QQ)
candidate_basis = sp.groebner(candidate, Y, X, q, p, order="lex", domain=sp.QQ)
assert all(candidate_basis.reduce(relation)[1] == 0 for relation in remaining)
assert all(affine.reduce(relation)[1] == 0 for relation in candidate)

# Monic leading terms Y^2, XY, X^3 give the relative basis 1,X,Y,X^2.
relative_basis = (sp.Integer(1), X, Y, X**2)
assert all(candidate_basis.reduce(monomial)[1] == monomial for monomial in relative_basis)

# At the coincident root the multiplication differs from the Boolean B_2.
special = sp.groebner((X**3, X * Y, Y**2), Y, X, order="lex", domain=sp.QQ)
assert all(special.reduce(monomial)[1] == monomial for monomial in relative_basis)
assert special.reduce(X**2 * X)[1] == 0
assert special.reduce(X**2 * Y)[1] == 0
assert special.reduce(Y * X)[1] == 0
assert special.reduce(Y * Y)[1] == 0

# The special socle contains the independent classes X^2 and Y.  In the
# Boolean algebra the only socle monomial is e1*e2.
e1, e2 = sp.symbols("e1 e2")
boolean = sp.groebner((e1**2, e2**2), e1, e2, order="lex", domain=sp.QQ)
assert boolean.reduce(e1 * e2)[1] == e1 * e2
assert boolean.reduce(e1 * e2 * e1)[1] == 0
assert boolean.reduce(e1 * e2 * e2)[1] == 0

print("PASS: active sixfold hyperedges, not merger-cycle rank, give the Boolean exponent")
print("PASS: degree 12, degree 18, and the strong degree-24 face have lengths 4, 8, 16")
print("PASS: one primitive block is k[e]/(e^2)")
print("PASS: the first merged block is k[X,Y]/(X^3,XY,Y^2), not the Boolean B_2")

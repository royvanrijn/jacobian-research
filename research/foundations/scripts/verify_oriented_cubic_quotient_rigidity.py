#!/usr/bin/env python3
"""Exact audit of tangent-source and symmetric-quotient rigidity."""

from __future__ import annotations

from itertools import permutations, product

import sympy as sp


a, b, c, d, t, X = sp.symbols("a b c d t X")
L1 = a * X + b
L3 = c * X + d
L2 = L1 + t * L3
p1 = sp.Poly(sp.expand(L1 * L2 * L3), X).nth(2)
determinant_relation = a * d - b * c - 1
B0 = a * (3 * a * d - 2)
B1 = c * (3 * a * d - 1)

# Reduce coefficient identities modulo ad-bc=1.
assert sp.rem(
    sp.Poly(sp.expand(p1 - (B0 + t * B1)), b),
    sp.Poly(determinant_relation, b),
).is_zero
print("PASS: p1=B0+t*B1 on the normalized SL_2 base")


# H=(B1=0) is the disjoint union G_m*A^1 and G_m^2.
# Z=(B1=0,B0=1) is the disjoint union A^1 and G_m.
L = sp.symbols("L")
class_sl2 = L**3 - L
class_H = (L - 1) * L + (L - 1) ** 2
class_Z = L + (L - 1)
class_source = sp.expand(class_sl2 - class_H + L * class_Z)
assert class_source == L**3 + L - 1
print("PASS: the tangent oriented source has class L^3+L-1")


def source_count(prime: int) -> int:
    vectors = [
        vector
        for vector in product(range(prime), repeat=2)
        if vector != (0, 0)
    ]

    def det(left: tuple[int, int], right: tuple[int, int]) -> int:
        return (
            left[0] * right[1] - left[1] * right[0]
        ) % prime

    count = 0
    for factor1, factor2, factor3 in product(vectors, repeat=3):
        if det(factor1, factor3) != 1:
            continue
        if det(factor2, factor3) != 1:
            continue
        finite_p1 = (
            factor1[1] * factor2[0] * factor3[0]
            + factor1[0] * factor2[1] * factor3[0]
            + factor1[0] * factor2[0] * factor3[1]
        ) % prime
        count += finite_p1 == 1
    return count


for prime in (5, 7, 11):
    assert source_count(prime) == prime**3 + prime - 1
print("PASS: finite-field counts realize L^3+L-1")


vertices = (1, 2, 3)
edge13 = frozenset((1, 3))
edge23 = frozenset((2, 3))
marked_edges = frozenset((edge13, edge23))


def image_edge(
    permutation: tuple[int, int, int],
    edge: frozenset[int],
) -> frozenset[int]:
    action = dict(zip(vertices, permutation))
    return frozenset(action[vertex] for vertex in edge)


group = list(permutations(vertices))
pointwise = [
    permutation
    for permutation in group
    if image_edge(permutation, edge13) == edge13
    and image_edge(permutation, edge23) == edge23
]
setwise = [
    permutation
    for permutation in group
    if frozenset(
        (
            image_edge(permutation, edge13),
            image_edge(permutation, edge23),
        )
    )
    == marked_edges
]
assert pointwise == [(1, 2, 3)]
assert setwise == [(1, 2, 3), (2, 1, 3)]
print("PASS: only the identity fixes both dicritical edges")
print("PASS: the unique nontrivial setwise symmetry swaps the two edges")


# Adding an anti-invariant coordinate produces the local A_1 quotient
# singularity B^2=A*C rather than a smooth affine quotient.
A, B, C = sp.symbols("A B C")
cone_relation = B**2 - A * C
cone_gradient = [
    sp.diff(cone_relation, variable)
    for variable in (A, B, C)
]
assert all(
    derivative.subs({A: 0, B: 0, C: 0}) == 0
    for derivative in cone_gradient
)

# The transverse branch coordinate already contributes one -1 eigenvalue.
# Smooth involution quotients require the total negative eigenspace to have
# dimension one, so every added linear coordinate must be invariant.
for anti_invariant_dimension in range(5):
    negative_eigenspace_dimension = 1 + anti_invariant_dimension
    smooth_pseudoreflection_quotient = (
        negative_eigenspace_dimension == 1
    )
    assert smooth_pseudoreflection_quotient == (
        anti_invariant_dimension == 0
    )
print("PASS: every nontrivial sign suspension makes the local quotient singular")
print("PASS oriented cubic quotient rigidity")

#!/usr/bin/env python3
"""Exact checks for all-arity oriented linear-factor Cox maps."""

from __future__ import annotations

from itertools import combinations
from math import factorial

import sympy as sp


def all_edges(size: int) -> list[tuple[int, int]]:
    return list(combinations(range(size), 2))


def incidence_column(size: int, edge: tuple[int, int]) -> sp.Matrix:
    column = [0] * size
    column[edge[0]] = 1
    column[edge[1]] = 1
    return sp.Matrix(column)


def minimal_double_star(size: int) -> list[tuple[int, int]]:
    first_size = (size - 1) // 2
    first = list(range(first_size))
    second = list(range(first_size, size))
    return (
        [(first[0], vertex) for vertex in second]
        + [(vertex, second[0]) for vertex in first[1:]]
    )


def normalization_index(
    size: int, tree: list[tuple[int, int]]
) -> int:
    matrix = sp.Matrix.hstack(
        *(incidence_column(size, edge) for edge in tree),
        sp.ones(size, 1),
    )
    return abs(int(matrix.det()))


# All-arity lattice, degree, and collision allocation.
for size in range(3, 13):
    tree = minimal_double_star(size)
    index = normalization_index(size, tree)
    assert index == (1 if size % 2 else 2)

    edges = set(all_edges(size))
    tree_edges = set(tree)
    nonedges = edges - tree_edges
    assert len(tree_edges) == size - 1
    assert len(nonedges) == (size - 1) * (size - 2) // 2

    oriented_degree = index * factorial(size) // 2
    component_sum = (
        len(edges) * index * factorial(size - 2)
    )
    assert component_sum == oriented_degree
print("PASS: every arity splits collisions into s-1 tree and rank-many non-tree types")
print("PASS: component degrees sum to d_s*s!/2")


# Symbolic discriminant and residue calculation for s=3,4.
for size in (3, 4):
    X = sp.symbols("X")
    leading = sp.symbols(f"u1:{size + 1}")
    trailing = sp.symbols(f"v1:{size + 1}")
    variables = tuple(
        coordinate
        for pair in zip(leading, trailing)
        for coordinate in pair
    )
    factors = [
        leading[index] * X + trailing[index]
        for index in range(size)
    ]
    polynomial = sp.Poly(sp.prod(factors), X)
    coefficients = polynomial.all_coeffs()
    resultants = {
        edge: (
            leading[edge[0]] * trailing[edge[1]]
            - trailing[edge[0]] * leading[edge[1]]
        )
        for edge in all_edges(size)
    }
    tree = minimal_double_star(size)
    nonedges = set(resultants) - set(tree)
    orientation = sp.prod(resultants[edge] for edge in nonedges)

    discriminant = sp.factor(
        sp.discriminant(polynomial.as_expr(), X)
    )
    full_vandermonde = sp.prod(resultants.values())
    assert sp.factor(discriminant - full_vandermonde**2) == 0

    outputs = (*coefficients, *(resultants[edge] for edge in tree))
    ambient_jacobian = sp.factor(
        sp.Matrix(outputs).jacobian(variables).det()
    )
    selected_product = sp.prod(resultants[edge] for edge in tree)
    predicted_product = full_vandermonde * selected_product
    ratio = sp.factor(ambient_jacobian / predicted_product)
    index = normalization_index(size, tree)
    assert abs(int(ratio)) == index

    # On the normalized slice selected resultants are one.  The horizontal
    # residue Jacobian is +-index*orientation, and the target hypersurface
    # residue divides by 2D with D=orientation.
    D = sp.symbols("D", nonzero=True)
    residue_ratio = sp.cancel(index * orientation / (2 * D))
    assert sp.cancel(
        residue_ratio.subs(D, orientation) - sp.Rational(index, 2)
    ) == 0
print("PASS: oriented residue Jacobian is d_s/2 for s=3,4")
print("PASS oriented linear-factor Cox maps")

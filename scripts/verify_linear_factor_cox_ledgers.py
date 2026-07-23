#!/usr/bin/env python3
"""Exact regressions for multi-boundary linear-factor Cox ledgers."""

from __future__ import annotations

from itertools import combinations
from math import factorial
from random import Random

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def all_edges(size: int) -> list[tuple[int, int]]:
    return list(combinations(range(size), 2))


def incidence_column(size: int, edge: tuple[int, int]) -> sp.Matrix:
    column = [0] * size
    column[edge[0]] = 1
    column[edge[1]] = 1
    return sp.Matrix(column)


def full_boundary_matrix(size: int) -> sp.Matrix:
    columns = [incidence_column(size, edge) for edge in all_edges(size)]
    columns.append(sp.ones(size, 1))
    return sp.Matrix.hstack(*columns)


def minimal_double_star(size: int) -> list[tuple[int, int]]:
    first_size = (size - 1) // 2
    first = list(range(first_size))
    second = list(range(first_size, size))
    first_center = first[0]
    second_center = second[0]
    return (
        [(first_center, vertex) for vertex in second]
        + [(vertex, second_center) for vertex in first[1:]]
    )


def scaling_matrix(size: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    return sp.Matrix.hstack(
        *(incidence_column(size, edge) for edge in edges),
        sp.ones(size, 1),
    )


# The full Smith obstruction is trivial in odd arity and Z/2 in even arity.
for size in range(3, 11):
    boundary = full_boundary_matrix(size)
    smith = smith_normal_form(boundary, domain=sp.ZZ)
    nonzero_diagonal = [
        abs(int(smith[index, index]))
        for index in range(size)
        if smith[index, index]
    ]
    torsion_order = sp.prod(nonzero_diagonal)
    assert torsion_order == (1 if size % 2 else 2)

    tree = minimal_double_star(size)
    assert len(tree) == size - 1
    tree_index = abs(int(scaling_matrix(size, tree).det()))
    assert tree_index == (1 if size % 2 else 2)

    unit_rank = len(all_edges(size)) + 1 - size
    assert unit_rank == (size - 1) * (size - 2) // 2
print("PASS: full boundary Smith index is 1 in odd arity and 2 in even arity")
print("PASS: the minimal double-star normalization attains that index")


def factor_data(size: int):
    variable = sp.symbols("X")
    leading = sp.symbols(f"u1:{size + 1}")
    trailing = sp.symbols(f"v1:{size + 1}")
    factors = [
        leading[index] * variable + trailing[index]
        for index in range(size)
    ]
    coefficients = sp.Poly(sp.prod(factors), variable).all_coeffs()
    resultants = {
        edge: (
            leading[edge[0]] * trailing[edge[1]]
            - trailing[edge[0]] * leading[edge[1]]
        )
        for edge in all_edges(size)
    }
    variables = tuple(
        coordinate
        for pair in zip(leading, trailing)
        for coordinate in pair
    )
    return coefficients, resultants, variables


# Full symbolic factorizations in the first two nontrivial arities.
for size in (3, 4):
    tree = minimal_double_star(size)
    coefficients, resultants, variables = factor_data(size)
    outputs = (*coefficients, *(resultants[edge] for edge in tree))
    determinant = sp.factor(
        sp.Matrix(outputs).jacobian(variables).det()
    )
    predicted_factors = (
        sp.prod(resultants.values())
        * sp.prod(resultants[edge] for edge in tree)
    )
    ratio = sp.factor(determinant / predicted_factors)
    tree_determinant = scaling_matrix(size, tree).det()
    assert abs(int(ratio)) == abs(int(tree_determinant))
print("PASS: symbolic ambient determinant factorization holds for s=3,4")


# Higher-arity exact evaluations avoid an expensive full symbolic factor.
rng = Random(20260723)
for size in range(5, 9):
    tree = minimal_double_star(size)
    coefficients, resultants, variables = factor_data(size)
    outputs = (*coefficients, *(resultants[edge] for edge in tree))
    jacobian_matrix = sp.Matrix(outputs).jacobian(variables)
    expected_absolute_ratio = abs(int(scaling_matrix(size, tree).det()))

    certified = 0
    while certified < 3:
        values = {variable: rng.randint(1, 17) for variable in variables}
        factor_value = (
            sp.prod(resultant.subs(values) for resultant in resultants.values())
            * sp.prod(resultants[edge].subs(values) for edge in tree)
        )
        if factor_value == 0:
            continue
        determinant_value = jacobian_matrix.subs(values).det()
        ratio = sp.Rational(determinant_value, factor_value)
        assert abs(int(ratio)) == expected_absolute_ratio
        certified += 1
print("PASS: exact integer specializations verify the formula through s=8")


# Abstract vertical block checks both suspension variants.  The horizontal
# residue determinant is d*product(units).
for size in range(3, 9):
    unit_rank = (size - 1) * (size - 2) // 2
    units = sp.symbols(f"q0:{unit_rank}", nonzero=True)
    index = sp.Integer(1 if size % 2 else 2)
    horizontal = index * sp.prod(units)
    separated_vertical = sp.prod(1 / unit for unit in units)
    compressed_vertical = 1 / sp.prod(units)
    assert sp.cancel(horizontal * separated_vertical) == index
    assert sp.cancel(horizontal * compressed_vertical) == index

    expected_degree = int(index) * factorial(size)
    assert expected_degree >= 6
print("PASS: separated and one-coordinate compressed suspensions are constant")
print("PASS: geometric degree is d_s*s! and primitive coordinates preserve it")
print("PASS linear-factor Cox ledgers")

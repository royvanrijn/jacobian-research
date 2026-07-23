#!/usr/bin/env python3
"""Compute the first universal inverse-tree obstruction for ``(JH)^3=0``.

Trees encode vector covariants ``D^k H[child_1,...,child_k]`` with
``0 <= k <= 3``.  The formal inverse equation

    K = -H(X+K)

then generates exact rational linear combinations of non-plane rooted trees.
This script computes the pieces through degree eleven and reduces ``K_11``
by the universal differential/context closure of two consequences of
``N^3=0``:

* Euler homogeneity gives ``N^2H=0`` from ``N^3X=0``;
* the directional derivative of ``N^3 v=0``, followed by ``v=H`` and
  direction ``H``.

The reduction is a universal tensor identity, not a sampled calculation.
An exact triangular model independently checks the tree expansion.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path
from typing import TypeAlias

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_three_inverse_tree_obstruction.json"
)

Tree: TypeAlias = tuple["Tree", ...] | str
Vector = Counter[Tree]
H_TREE: Tree = ()


def tree_key(tree: Tree) -> str:
    if isinstance(tree, str):
        return tree
    return "(" + ",".join(tree_key(child) for child in tree) + ")"


def tree_formula(tree: Tree) -> str:
    if isinstance(tree, str):
        return tree
    if not tree:
        return "H"
    operator = {1: "N", 2: "B", 3: "C"}[len(tree)]
    return f"{operator}({','.join(tree_formula(child) for child in tree)})"


def node(children: tuple[Tree, ...]) -> Tree:
    if len(children) > 3:
        raise ValueError("a cubic H has no derivatives above order three")
    return tuple(sorted(children, key=tree_key))


def tree_degree(tree: Tree) -> int:
    if isinstance(tree, str):
        raise ValueError("marked leaves do not have a polynomial degree")
    return 3 + sum(tree_degree(child) - 1 for child in tree)


def add_scaled(target: Vector, source: Vector, scale: Fraction) -> None:
    for tree, coefficient in source.items():
        value = target[tree] + scale * coefficient
        if value:
            target[tree] = value
        else:
            target.pop(tree, None)


def inverse_parts(max_degree: int) -> dict[int, Vector]:
    parts: dict[int, Vector] = {}
    available: list[tuple[Tree, Fraction]] = []
    for degree in range(3, max_degree + 1, 2):
        answer: Vector = Counter()
        if degree == 3:
            answer[H_TREE] = Fraction(-1)
        for arity in range(1, 4):
            for choices in product(available, repeat=arity):
                children = tuple(choice[0] for choice in choices)
                root = node(children)
                if tree_degree(root) != degree:
                    continue
                coefficient = Fraction(-1, math.factorial(arity))
                for _, child_coefficient in choices:
                    coefficient *= child_coefficient
                answer[root] += coefficient
        parts[degree] = answer
        available.extend(answer.items())
    return parts


def derivative(tree: Tree, direction: Tree) -> Vector:
    """Differentiate a tree in an x-direction treated as constant."""

    if isinstance(tree, str):
        return Counter()
    answer: Vector = Counter()
    if len(tree) < 3:
        answer[node(tree + (direction,))] += 1
    for index, child in enumerate(tree):
        for replacement, coefficient in derivative(child, direction).items():
            children = list(tree)
            children[index] = replacement
            answer[node(tuple(children))] += coefficient
    return answer


def substitute(tree: Tree, replacements: dict[str, Tree]) -> Tree:
    if isinstance(tree, str):
        return replacements[tree]
    return node(tuple(substitute(child, replacements) for child in tree))


def substitute_vector(vector: Vector, replacements: dict[str, Tree]) -> Vector:
    answer: Vector = Counter()
    for tree, coefficient in vector.items():
        answer[substitute(tree, replacements)] += coefficient
    return answer


def vector_derivative(vector: Vector, direction: Tree) -> Vector:
    answer: Vector = Counter()
    for tree, coefficient in vector.items():
        add_scaled(answer, derivative(tree, direction), coefficient)
    return answer


def vector_context(vector: Vector, other_children: tuple[Tree, ...]) -> Vector:
    answer: Vector = Counter()
    for tree, coefficient in vector.items():
        answer[node((tree,) + other_children)] += coefficient
    return answer


def vector_degree(vector: Vector) -> int:
    degrees = {tree_degree(tree) for tree in vector}
    if len(degrees) != 1:
        raise AssertionError("relation is not homogeneous")
    return degrees.pop()


def normalized_signature(vector: Vector) -> tuple[tuple[str, Fraction], ...]:
    ordered = sorted(
        ((tree_key(tree), coefficient) for tree, coefficient in vector.items()),
        key=lambda item: item[0],
    )
    if not ordered:
        return ()
    scale = ordered[0][1]
    return tuple((key, coefficient / scale) for key, coefficient in ordered)


def relation_closure(
    seeds: list[Vector],
    tree_pool: list[Tree],
    max_degree: int,
) -> list[Vector]:
    """Close vector identities under differentiation and cubic contexts."""

    relations: list[Vector] = []
    seen: set[tuple[tuple[str, Fraction], ...]] = set()
    queue = list(seeds)
    while queue:
        relation = queue.pop(0)
        signature = normalized_signature(relation)
        if not signature or signature in seen:
            continue
        seen.add(signature)
        relations.append(relation)
        degree = vector_degree(relation)

        for direction in tree_pool:
            if degree + tree_degree(direction) - 1 <= max_degree:
                queue.append(vector_derivative(relation, direction))

        for arity in range(1, 4):
            for others in product(tree_pool, repeat=arity - 1):
                context_degree = (
                    3
                    + (degree - 1)
                    + sum(tree_degree(tree) - 1 for tree in others)
                )
                if context_degree <= max_degree:
                    queue.append(vector_context(relation, tuple(others)))
    return relations


def unary_chain(length: int, child: Tree) -> Tree:
    answer = child
    for _ in range(length):
        answer = node((answer,))
    return answer


def relation_span_contains(target: Vector, relations: list[Vector]) -> tuple[bool, list[Fraction]]:
    basis = sorted(
        set(target).union(*(set(relation) for relation in relations)),
        key=tree_key,
    )
    relation_matrix = sp.Matrix(
        [
            [sp.Rational(relation.get(tree, 0)) for relation in relations]
            for tree in basis
        ]
    )
    target_column = sp.Matrix(
        [sp.Rational(target.get(tree, 0)) for tree in basis]
    )
    augmented = relation_matrix.row_join(target_column)
    contained = relation_matrix.rank() == augmented.rank()
    if not contained:
        return False, []
    solution = sp.linsolve((relation_matrix, target_column))
    vector = next(iter(solution))
    substitutions = {
        symbol: 0
        for value in vector
        for symbol in value.free_symbols
    }
    return True, [Fraction(value.subs(substitutions)) for value in vector]


def quotient_normal_form(
    target: Vector, relations: list[Vector]
) -> dict[str, object]:
    basis = sorted(
        set(target).union(*(set(relation) for relation in relations)),
        key=tree_key,
    )
    rows = sp.Matrix(
        [
            [sp.Rational(relation.get(tree, 0)) for tree in basis]
            for relation in relations
        ]
    )
    reduced, pivots = rows.rref()
    nonzero_rows = [list(reduced.row(row)) for row in range(reduced.rows) if any(reduced.row(row))]
    remainder = [sp.Rational(target.get(tree, 0)) for tree in basis]
    for row in nonzero_rows:
        pivot = next(index for index, value in enumerate(row) if value)
        scale = remainder[pivot]
        if scale:
            remainder = [
                left - scale * right for left, right in zip(remainder, row)
            ]
    return {
        "tree_space_dimension": len(basis),
        "relation_rank": len(nonzero_rows),
        "quotient_dimension": len(basis) - len(nonzero_rows),
        "pivot_tree_indices": list(pivots),
        "target_normal_form": {
            tree_key(tree): str(value)
            for tree, value in zip(basis, remainder)
            if value
        },
        "target_normal_form_formula": {
            tree_formula(tree): str(value)
            for tree, value in zip(basis, remainder)
            if value
        },
    }


def apply_derivative_tensor(
    h: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    children: list[list[sp.Expr]],
) -> list[sp.Expr]:
    if not children:
        return h
    answer = []
    for component in h:
        value = sp.Integer(0)
        for indices in product(range(len(variables)), repeat=len(children)):
            derivative_value = component
            for index in indices:
                derivative_value = sp.diff(
                    derivative_value, variables[index]
                )
            value += derivative_value * sp.prod(
                children[slot][index] for slot, index in enumerate(indices)
            )
        answer.append(sp.expand(value))
    return answer


def evaluate_tree(
    tree: Tree,
    h: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    cache: dict[Tree, list[sp.Expr]],
) -> list[sp.Expr]:
    if tree in cache:
        return cache[tree]
    if isinstance(tree, str):
        raise ValueError("cannot evaluate a marked tree")
    children = [evaluate_tree(child, h, variables, cache) for child in tree]
    value = apply_derivative_tensor(h, variables, children)
    cache[tree] = value
    return value


def evaluate_vector(
    vector: Vector,
    h: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:
    cache: dict[Tree, list[sp.Expr]] = {}
    answer = [sp.Integer(0)] * len(h)
    for tree, coefficient in vector.items():
        value = evaluate_tree(tree, h, variables, cache)
        answer = [
            sp.expand(left + sp.Rational(coefficient) * right)
            for left, right in zip(answer, value)
        ]
    return answer


def main() -> None:
    parts = inverse_parts(11)
    assert [len(parts[degree]) for degree in parts] == [1, 1, 2, 4, 8]

    # Euler: N X=3H, hence N^3 X=3N^2 H=0.
    relation_euler = Counter({unary_chain(2, H_TREE): Fraction(1)})

    # D_w(N^3 v)=0, followed by v=H and w=H.
    marked_identity = unary_chain(3, "v")
    differentiated = derivative(marked_identity, "w")
    relation_derivative = substitute_vector(
        differentiated,
        {"v": H_TREE, "w": H_TREE},
    )
    tree_pool = sorted(
        {
            tree
            for degree, vector in parts.items()
            if degree <= 5
            for tree in vector
        },
        key=tree_key,
    )
    closure = relation_closure(
        [relation_euler, relation_derivative],
        tree_pool,
        11,
    )
    relations = [
        relation for relation in closure if vector_degree(relation) == 11
    ]
    assert relations

    contained, coefficients = relation_span_contains(parts[11], relations)
    quotient = quotient_normal_form(parts[11], relations)
    assert contained == (not quotient["target_normal_form"])

    # Independent exact check against the degree-nine triangular calibration.
    x = sp.symbols("x0:3")
    y = sp.symbols("y0:3")
    h = [x[1] ** 3, x[2] ** 3, sp.Integer(0)]
    exact_inverse_correction = [
        -(y[1] - y[2] ** 3) ** 3,
        -y[2] ** 3,
        sp.Integer(0),
    ]
    substitution = dict(zip(x, y))
    for degree, vector in parts.items():
        evaluated = [
            sp.expand(value.subs(substitution, simultaneous=True))
            for value in evaluate_vector(vector, h, x)
        ]
        expected = [
            sum(
                coefficient
                * sp.prod(
                    variable**power
                    for variable, power in zip(y, exponents)
                )
                for exponents, coefficient in sp.Poly(
                    component, *y, domain=sp.QQ
                ).terms()
                if sum(exponents) == degree
            )
            for component in exact_inverse_correction
        ]
        assert evaluated == expected

    payload = {
        "format": "cubic-index-three-inverse-tree-obstruction-v1",
        "field": "characteristic zero",
        "inverse_tree_counts_by_degree": {
            str(degree): len(vector) for degree, vector in parts.items()
        },
        "degree_11_tree_coefficients": {
            tree_key(tree): str(coefficient)
            for tree, coefficient in sorted(
                parts[11].items(), key=lambda item: tree_key(item[0])
            )
        },
        "universal_relation_counts_by_degree": {
            str(degree): sum(
                vector_degree(relation) == degree for relation in closure
            )
            for degree in (7, 9, 11)
        },
        "universal_degree_11_relations": [
            {
                tree_key(tree): str(coefficient)
                for tree, coefficient in sorted(
                    relation.items(), key=lambda item: tree_key(item[0])
                )
            }
            for relation in relations
        ],
        "degree_11_lies_in_relation_span": contained,
        "degree_11_quotient": quotient,
        "relation_combination_coefficients": [
            str(coefficient) for coefficient in coefficients
        ],
        "interpretation": (
            "if true, K_11 vanishes under the displayed universal closure; "
            "otherwise its displayed quotient class is the first unresolved "
            "inverse covariant and additional polarized identities are needed"
        ),
        "triangular_model_check": (
            "all tree pieces through degree 11 agree with the exact "
            "degree-nine inverse"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS index-three trees: exact inverse recursion through degree 11")
    print("PASS index-three trees: triangular calibration agrees degree by degree")
    print(
        "RESULT index-three trees: "
        f"K_11 in Euler/differential/context span = {contained}"
    )
    print(f"PASS index-three trees: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

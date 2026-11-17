#!/usr/bin/env python3
"""Verify the homological filling obstruction for the degree-six cubic E8 row."""

from __future__ import annotations

from itertools import permutations

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


Permutation = tuple[int, ...]


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(value: Permutation) -> Permutation:
    result = [0] * len(value)
    for i, image in enumerate(value):
        result[image] = i
    return tuple(result)


def power(value: Permutation, exponent: int) -> Permutation:
    result = tuple(range(len(value)))
    while exponent:
        if exponent & 1:
            result = compose(value, result)
        value = compose(value, value)
        exponent //= 2
    return result


def cycle_type(value: Permutation) -> tuple[int, ...]:
    seen: set[int] = set()
    result: list[int] = []
    for start in range(len(value)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = value[point]
        result.append(length)
    return tuple(sorted(result, reverse=True))


def lifted_path(
    start: int,
    word: tuple[tuple[str, int], ...],
    a: Permutation,
    b: Permutation,
) -> tuple[int, tuple[int, ...]]:
    """Return the endpoint and cellular one-chain of a lifted word."""

    inverses = {"a": inverse(a), "b": inverse(b)}
    generators = {"a": a, "b": b}
    offsets = {"a": 0, "b": len(a)}
    point = start
    chain = [0] * (2 * len(a))
    for generator, sign in word:
        if sign == 1:
            chain[offsets[generator] + point] += 1
            point = generators[generator][point]
        else:
            previous = inverses[generator][point]
            chain[offsets[generator] + previous] -= 1
            point = previous
    return point, tuple(chain)


def spanning_tree_chords(a: Permutation, b: Permutation) -> tuple[int, ...]:
    degree = len(a)
    parent = list(range(degree))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> bool:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return False
        parent[right_root] = left_root
        return True

    tree: list[int] = []
    edges = tuple((source, value[source]) for value in (a, b) for source in range(degree))
    for edge_index, (source, target) in enumerate(edges):
        if source != target and union(source, target):
            tree.append(edge_index)
    assert len(tree) == degree - 1
    return tuple(index for index in range(2 * degree) if index not in tree)


def filled_homology_for_meridian_word(
    meridian_word: tuple[tuple[str, int], ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    # The displayed unique cubic action.
    a = (1, 2, 0, 4, 5, 3)  # (1 2 3)(4 5 6)
    b = (0, 3, 4, 2, 5, 1)  # (2 4 3 5 6)
    identity = tuple(range(6))
    assert power(a, 3) == power(b, 5) == identity

    relator = (("a", 1),) * 3 + (("b", -1),) * 5
    relator_lifts = []
    for sheet in range(6):
        endpoint, chain = lifted_path(sheet, relator, a, b)
        assert endpoint == sheet
        relator_lifts.append(chain)

    fixed_lifts = []
    endpoints = []
    for sheet in range(6):
        endpoint, chain = lifted_path(sheet, meridian_word, a, b)
        endpoints.append(endpoint)
        if endpoint == sheet:
            fixed_lifts.append(chain)
    assert cycle_type(tuple(endpoints)) == (3, 1, 1, 1)
    assert len(fixed_lifts) == 3

    chords = spanning_tree_chords(a, b)
    assert len(chords) == 7

    unfilled_relations = sp.Matrix(
        [[chain[edge] for chain in relator_lifts] for edge in chords]
    )
    filled_relations = sp.Matrix(
        [
            [chain[edge] for chain in relator_lifts + fixed_lifts]
            for edge in chords
        ]
    )

    unfilled_smith = smith_normal_form(unfilled_relations, domain=ZZ)
    filled_smith = smith_normal_form(filled_relations, domain=ZZ)
    unfilled_diagonal = tuple(
        int(unfilled_smith[index, index])
        for index in range(min(unfilled_smith.shape))
    )
    filled_diagonal = tuple(
        int(filled_smith[index, index])
        for index in range(min(filled_smith.shape))
    )
    return unfilled_diagonal, filled_diagonal


def main() -> None:
    # Depending on whether fiber monodromy is encoded as a left or right
    # action, the geometric meridian is read as a^-1*b^2 or in reverse path
    # order.  Both cellular conventions give the same integral homology.
    direct = (("a", -1), ("b", 1), ("b", 1))
    reverse_path = (("b", 1), ("b", 1), ("a", -1))
    results = {
        filled_homology_for_meridian_word(direct),
        filled_homology_for_meridian_word(reverse_path),
    }
    assert results == {((1, 1, 1, 0, 0, 0), (1, 1, 1, 1, 1, 1, 0))}

    # Seven cycle generators modulo six primitive relations leave one free
    # class and no torsion after all three fixed-sheet meridians are filled.
    filled_diagonal = next(iter(results))[1]
    assert filled_diagonal.count(0) == 1
    assert all(value in (0, 1) for value in filled_diagonal)

    print(
        "PASS: the six-sheet E8 cusp cover has H_1=Z^4 before filling; "
        "after filling all three fixed-sheet affine meridians its Smith "
        "diagonal is (1,1,1,1,1,1,0), so H_1=Z and the filled space "
        "cannot be A^2"
    )


if __name__ == "__main__":
    main()

"""Permutation and forest helpers for polynomial admissible-cover bubbles."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import combinations, permutations, product


Permutation = tuple[int, ...]
Edge = tuple[int, int]


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return left after right."""

    if len(left) != len(right):
        raise ValueError("permutations must have the same degree")
    return tuple(left[right[index]] for index in range(len(left)))


def transposition(degree: int, left: int, right: int) -> Permutation:
    result = list(range(degree))
    result[left], result[right] = result[right], result[left]
    return tuple(result)


def permutation_product(
    permutations: Iterable[Permutation],
    degree: int,
) -> Permutation:
    result = identity(degree)
    for permutation in permutations:
        result = compose(permutation, result)
    return result


def cycles(permutation: Permutation) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    result = []
    while unseen:
        start = min(unseen)
        cycle = []
        current = start
        while current in unseen:
            unseen.remove(current)
            cycle.append(current)
            current = permutation[current]
        result.append(tuple(cycle))
    return tuple(result)


def cycle_lengths(permutation: Permutation) -> tuple[int, ...]:
    return tuple(sorted(len(cycle) for cycle in cycles(permutation)))


def forest_components(
    degree: int,
    edges: Sequence[Edge],
) -> tuple[tuple[int, ...], ...]:
    adjacency = {vertex: set() for vertex in range(degree)}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)

    unseen = set(range(degree))
    result = []
    while unseen:
        start = min(unseen)
        stack = [start]
        component = []
        unseen.remove(start)
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        result.append(tuple(sorted(component)))
    return tuple(sorted(result, key=lambda component: component[0]))


def is_tree(degree: int, edges: Sequence[Edge]) -> bool:
    return (
        len(edges) == degree - 1
        and len(forest_components(degree, edges)) == 1
    )


def standard_cycle(degree: int) -> Permutation:
    return tuple((index + 1) % degree for index in range(degree))


def reduced_cycle_factorizations(
    degree: int,
) -> tuple[tuple[Edge, ...], ...]:
    """Return all minimal transposition factorizations of a fixed d-cycle."""

    edges = tuple(combinations(range(degree), 2))
    edge_permutations = {
        edge: transposition(degree, *edge)
        for edge in edges
    }
    target = standard_cycle(degree)
    result = []
    for edge_sequence in product(edges, repeat=degree - 1):
        factors = tuple(edge_permutations[edge] for edge in edge_sequence)
        if permutation_product(factors, degree) == target:
            result.append(edge_sequence)
    return tuple(result)


def centralizer(generators: Sequence[Permutation]) -> tuple[Permutation, ...]:
    """Return the centralizer of a permutation subgroup in a symmetric group."""

    if not generators:
        raise ValueError("at least one generator is required")
    degree = len(generators[0])
    result = []
    for candidate in permutations(range(degree)):
        if all(
            compose(candidate, generator) == compose(generator, candidate)
            for generator in generators
        ):
            result.append(candidate)
    return tuple(result)


def fixes_labels(
    permutation: Permutation,
    labels: Iterable[int],
) -> bool:
    return all(permutation[label] == label for label in labels)

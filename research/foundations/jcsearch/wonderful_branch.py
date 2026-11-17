"""Nested-set combinatorics for labelled genus-zero target boundaries."""

from __future__ import annotations

from itertools import combinations
from typing import Iterable


Boundary = frozenset[int]


def boundary_subsets(mark_count: int) -> tuple[Boundary, ...]:
    """Represent Mbar_0,n boundary splits by the side not containing infinity."""

    if mark_count < 4:
        raise ValueError("at least four target marks are required")
    finite_marks = tuple(range(mark_count - 1))
    return tuple(
        frozenset(subset)
        for size in range(2, mark_count - 1)
        for subset in combinations(finite_marks, size)
    )


def compatible(left: Boundary, right: Boundary) -> bool:
    """Return whether two genus-zero boundary splits are compatible."""

    return (
        left <= right
        or right <= left
        or left.isdisjoint(right)
    )


def maximal_nested_sets(mark_count: int) -> tuple[tuple[Boundary, ...], ...]:
    """Enumerate zero-stratum nested sets of Mbar_0,n."""

    boundaries = boundary_subsets(mark_count)
    target_size = mark_count - 3
    return tuple(
        collection
        for collection in combinations(boundaries, target_size)
        if all(
            compatible(left, right)
            for left, right in combinations(collection, 2)
        )
    )


def act_on_boundary(
    boundary: Boundary,
    permutation: dict[int, int],
) -> Boundary:
    """Transport a boundary split by a permutation fixing infinity."""

    return frozenset(permutation.get(mark, mark) for mark in boundary)


def diagonal_generators(
    boundary: Boundary,
    value_names: dict[int, str],
    *,
    zero_mark: int = 0,
) -> tuple[str, ...]:
    """Return canonical generators for the corresponding affine diagonal."""

    ordered = sorted(boundary)
    if zero_mark in boundary:
        return tuple(value_names[mark] for mark in ordered if mark != zero_mark)
    base = ordered[0]
    return tuple(
        f"{value_names[mark]}-{value_names[base]}"
        for mark in ordered[1:]
    )


def is_nested(collection: Iterable[Boundary]) -> bool:
    """Return whether every pair in a collection is compatible."""

    values = tuple(collection)
    return all(
        compatible(left, right)
        for left, right in combinations(values, 2)
    )

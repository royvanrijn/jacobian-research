"""Finite phase algebra for normalized admissible-cover node charts."""

from __future__ import annotations

import math
from itertools import product
from typing import Iterable, Sequence


def node_lcm(indices: Sequence[int]) -> int:
    """Return the common Kummer index at one target node."""

    if not indices or any(index <= 0 for index in indices):
        raise ValueError("node indices must be a nonempty positive sequence")
    return math.lcm(*indices)


def phase_group(indices: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Return the product of source-node root-of-unity groups additively."""

    node_lcm(indices)
    return tuple(product(*(range(index) for index in indices)))


def diagonal_phase_group(
    indices: Sequence[int],
) -> frozenset[tuple[int, ...]]:
    """Return the simultaneous reparametrization subgroup mu_L."""

    common_index = node_lcm(indices)
    diagonal = frozenset(
        tuple(phase % index for index in indices)
        for phase in range(common_index)
    )
    assert len(diagonal) == common_index
    return diagonal


def normalization_branch_count(indices: Sequence[int]) -> int:
    """Return prod(e_j)/lcm(e_j) for the Harris--Mumford node ring."""

    return math.prod(indices) // node_lcm(indices)


def inertia_character_count(
    indices: Sequence[int],
    character_values: Iterable[Sequence[int]],
) -> int:
    """Count label-preserving automorphisms whose node phase is diagonal.

    ``character_values`` is the image, with multiplicity, of the
    label-preserving cover-automorphism group in ``prod_j Z/e_j``.
    """

    diagonal = diagonal_phase_group(indices)
    count = 0
    for value in character_values:
        reduced = tuple(
            phase % index for phase, index in zip(value, indices)
        )
        if len(reduced) != len(indices):
            raise ValueError("character value has the wrong length")
        if reduced in diagonal:
            count += 1
    return count


def compose_character_tables(
    moduli: Sequence[int],
    local_tables: Sequence[Sequence[Sequence[int]]],
) -> tuple[tuple[int, ...], ...]:
    """Compose independent local automorphism character tables.

    Every row of a local table is one automorphism character on the same
    flattened list of source-node phase coordinates.  Cartesian products of
    local automorphisms add their characters.  Duplicate rows are retained,
    because a character kernel contributes genuine inertia.
    """

    if not moduli or any(modulus <= 0 for modulus in moduli):
        raise ValueError("phase moduli must be a nonempty positive sequence")
    zero = (0,) * len(moduli)
    composed = (zero,)
    for table in local_tables:
        if not table:
            raise ValueError("local character tables must be nonempty")
        next_composed = []
        for accumulated in composed:
            for value in table:
                if len(value) != len(moduli):
                    raise ValueError("character value has the wrong length")
                next_composed.append(
                    tuple(
                        (left + right) % modulus
                        for left, right, modulus in zip(
                            accumulated,
                            value,
                            moduli,
                        )
                    )
                )
        composed = tuple(next_composed)
    return composed


def simultaneous_inertia_character_count(
    node_profiles: Sequence[Sequence[int]],
    character_values: Iterable[Sequence[int]],
) -> int:
    """Count automorphisms diagonal at every target node simultaneously."""

    return sum(
        simultaneous_phase_is_diagonal(node_profiles, value)
        for value in character_values
    )


def simultaneous_phase_is_diagonal(
    node_profiles: Sequence[Sequence[int]],
    value: Sequence[int],
) -> bool:
    """Return whether one phase row is diagonal at every target node."""

    if not node_profiles:
        raise ValueError("at least one target-node profile is required")
    profiles = tuple(tuple(profile) for profile in node_profiles)
    diagonals = tuple(diagonal_phase_group(profile) for profile in profiles)
    width = sum(len(profile) for profile in profiles)

    if len(value) != width:
        raise ValueError("simultaneous character has the wrong length")
    offset = 0
    for profile, diagonal in zip(profiles, diagonals):
        segment = tuple(
            phase % index
            for phase, index in zip(
                value[offset : offset + len(profile)],
                profile,
            )
        )
        if segment not in diagonal:
            return False
        offset += len(profile)
    return True


def gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    """Return the row rank of a binary matrix."""

    matrix = [[entry & 1 for entry in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows must have equal length")
    rank = 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column]:
                matrix[row] = [
                    left ^ right
                    for left, right in zip(matrix[row], matrix[rank])
                ]
        rank += 1
    return rank

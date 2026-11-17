"""Combinatorics of radial source expansions for polynomial root clusters."""

from __future__ import annotations

import math
from itertools import product
from typing import Hashable, Sequence, TypeVar


Label = TypeVar("Label", bound=Hashable)


def ordered_set_partitions(
    labels: Sequence[Label],
) -> tuple[tuple[tuple[Label, ...], ...], ...]:
    """Return all ordered set partitions, preserving label order in blocks."""

    labels = tuple(labels)
    result = []
    for level_count in range(1, len(labels) + 1):
        for assignment in product(range(level_count), repeat=len(labels)):
            if set(assignment) != set(range(level_count)):
                continue
            result.append(
                tuple(
                    tuple(
                        label
                        for label, level in zip(labels, assignment)
                        if level == block
                    )
                    for block in range(level_count)
                )
            )
    return tuple(result)


def radial_component_kinds(
    multiplicity: int,
    cluster_level: int,
    target_level: int,
) -> tuple[tuple[str, int, int], ...]:
    """Return (kind, degree, RH contribution) on one target bubble."""

    if multiplicity < 1:
        raise ValueError("cluster multiplicities must be positive")
    if cluster_level > target_level:
        return (("power_connector", multiplicity, 2 * multiplicity - 2),)
    if cluster_level == target_level:
        return (
            ("local_polynomial_tail", multiplicity, 2 * multiplicity - 2),
        )
    return tuple(("identity_strand", 1, 0) for _ in range(multiplicity))


def radial_node_indices(
    multiplicities: Sequence[int],
    cluster_levels: Sequence[int],
    after_level: int,
) -> tuple[int, ...]:
    """Return the source expansion indices over one radial target node."""

    if len(multiplicities) != len(cluster_levels):
        raise ValueError("multiplicities and cluster levels must have equal length")
    indices = []
    for multiplicity, cluster_level in zip(
        multiplicities,
        cluster_levels,
    ):
        if cluster_level <= after_level:
            indices.extend((1,) * multiplicity)
        else:
            indices.append(multiplicity)
    return tuple(sorted(indices))


def radial_inertia_factors(
    multiplicities: Sequence[int],
    ordered_blocks: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """Return the successive full-chain radial inertia factors.

    If ``B_0|...|B_k`` is the ordered scale partition, put

        M_j = lcm(mu_i : i in B_j),
        L_j = lcm(mu_i : i in B_j union ... union B_k).

    The diagonal phase at the node before level ``j`` is an element of
    ``mu_(L_j)``.  Rigidification at the active polynomial tails imposes its
    accumulated phase modulo ``M_j``.  Hence level ``j`` contributes
    ``L_j / M_j`` choices.  The last factor is always one.
    """

    multiplicities = tuple(multiplicities)
    blocks = tuple(tuple(block) for block in ordered_blocks)
    if not multiplicities:
        raise ValueError("at least one cluster multiplicity is required")
    if any(multiplicity < 1 for multiplicity in multiplicities):
        raise ValueError("cluster multiplicities must be positive")
    flattened = tuple(index for block in blocks for index in block)
    if (
        any(not block for block in blocks)
        or len(flattened) != len(multiplicities)
        or set(flattened) != set(range(len(multiplicities)))
    ):
        raise ValueError(
            "ordered blocks must partition the cluster indices"
        )

    result = []
    for level, block in enumerate(blocks):
        block_lcm = math.lcm(
            *(multiplicities[index] for index in block)
        )
        suffix_lcm = math.lcm(
            *(
                multiplicities[index]
                for suffix in blocks[level:]
                for index in suffix
            )
        )
        if suffix_lcm % block_lcm:
            raise AssertionError("a block lcm must divide its suffix lcm")
        result.append(suffix_lcm // block_lcm)
    return tuple(result)


def radial_inertia_order(
    multiplicities: Sequence[int],
    ordered_blocks: Sequence[Sequence[int]],
) -> int:
    """Return the normalized full-chain inertia order of a radial lift."""

    return math.prod(
        radial_inertia_factors(multiplicities, ordered_blocks)
    )

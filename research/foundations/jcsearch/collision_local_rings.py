"""Combinatorial invariants of tame admissible-cover local rings.

For one target node with source-node expansion indices ``e_j``, the raw
Harris--Mumford factor is

    k[[q,s_1,...,s_r]] / (s_j**e_j - q).

This module records its normalization branches, Fitting/different orders,
and conductor.  It also records the fan of a normalized two-generator
monomial blowup, the first contact-monoid test in degree five.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product
from typing import Iterable, Sequence


Phase = tuple[int, ...]


def _validate_profile(profile: Sequence[int]) -> tuple[int, ...]:
    indices = tuple(profile)
    if not indices or any(index <= 0 for index in indices):
        raise ValueError("an expansion profile must contain positive indices")
    return indices


def _diagonal_phase(profile: Sequence[int], phase: int) -> Phase:
    return tuple(phase % index for index in profile)


def phase_quotient_representatives(
    profile: Sequence[int],
) -> tuple[Phase, ...]:
    """Canonical representatives for prod(mu_ej)/mu_lcm.

    Phases are written additively.  Simultaneous reparametrization of the
    normalization parameter adds ``(a mod e_j)_j``.
    """

    indices = _validate_profile(profile)
    common_index = math.lcm(*indices)
    unvisited = set(product(*(range(index) for index in indices)))
    representatives: list[Phase] = []
    while unvisited:
        representative = min(unvisited)
        orbit = {
            tuple(
                (entry + shift) % index
                for entry, index in zip(representative, indices)
            )
            for shift in range(common_index)
        }
        representatives.append(min(orbit))
        unvisited.difference_update(orbit)
    return tuple(sorted(representatives))


@dataclass(frozen=True)
class NodeLocalAlgebra:
    """Numerical package for one Harris--Mumford node factor."""

    profile: tuple[int, ...]
    common_index: int
    source_orders: tuple[int, ...]
    phase_representatives: tuple[Phase, ...]
    branch_count: int
    raw_fitting_order: int
    normalized_different_order: int
    conductor_order: int


def node_local_algebra(profile: Sequence[int]) -> NodeLocalAlgebra:
    """Compile normalization, Fitting, different, and conductor data.

    On a normalization branch,

        q = t**L,  s_j = omega_j*t**(L/e_j).

    The raw relative Fitting order is the valuation of
    ``prod_j s_j**(e_j-1)``.  Since the node curve is a reduced complete
    intersection, adjunction identifies the conductor order with the raw
    Fitting order minus the different order ``L-1`` of ``k[[t]]/k[[q]]``.
    """

    indices = _validate_profile(profile)
    common_index = math.lcm(*indices)
    source_orders = tuple(common_index // index for index in indices)
    phases = phase_quotient_representatives(indices)
    branch_count = math.prod(indices) // common_index
    if len(phases) != branch_count:
        raise AssertionError("phase quotient has the wrong cardinality")
    raw_fitting_order = sum(
        (index - 1) * order
        for index, order in zip(indices, source_orders)
    )
    normalized_different_order = common_index - 1
    conductor_order = raw_fitting_order - normalized_different_order
    if conductor_order < 0:
        raise AssertionError("a conductor order cannot be negative")
    return NodeLocalAlgebra(
        profile=indices,
        common_index=common_index,
        source_orders=source_orders,
        phase_representatives=phases,
        branch_count=branch_count,
        raw_fitting_order=raw_fitting_order,
        normalized_different_order=normalized_different_order,
        conductor_order=conductor_order,
    )


@dataclass(frozen=True)
class NestedNodeLocalAlgebra:
    """Product package for independent target-node factors."""

    nodes: tuple[NodeLocalAlgebra, ...]
    branch_count: int
    conductor_exponent_vector: tuple[int, ...]
    raw_fitting_exponent_vector: tuple[int, ...]
    normalized_different_exponent_vector: tuple[int, ...]


def nested_node_local_algebra(
    profiles: Iterable[Sequence[int]],
) -> NestedNodeLocalAlgebra:
    """Compile the completed tensor product over several target nodes."""

    nodes = tuple(node_local_algebra(profile) for profile in profiles)
    if not nodes:
        raise ValueError("a nested package needs at least one target node")
    return NestedNodeLocalAlgebra(
        nodes=nodes,
        branch_count=math.prod(node.branch_count for node in nodes),
        conductor_exponent_vector=tuple(
            node.conductor_order for node in nodes
        ),
        raw_fitting_exponent_vector=tuple(
            node.raw_fitting_order for node in nodes
        ),
        normalized_different_exponent_vector=tuple(
            node.normalized_different_order for node in nodes
        ),
    )


@dataclass(frozen=True)
class MonomialBlowup:
    """Fan data for the normalized blowup of ``(x^a,y^b)``."""

    exponents: tuple[int, int]
    primitive_ray: tuple[int, int]
    chart_indices: tuple[int, int]


def monomial_blowup(exponent_x: int, exponent_y: int) -> MonomialBlowup:
    """Return the primitive subdivision ray and the two cone indices."""

    if exponent_x <= 0 or exponent_y <= 0:
        raise ValueError("monomial exponents must be positive")
    common = math.gcd(exponent_x, exponent_y)
    primitive_ray = (
        exponent_y // common,
        exponent_x // common,
    )
    return MonomialBlowup(
        exponents=(exponent_x, exponent_y),
        primitive_ray=primitive_ray,
        chart_indices=(primitive_ray[1], primitive_ray[0]),
    )

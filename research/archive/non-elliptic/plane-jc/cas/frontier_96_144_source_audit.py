#!/usr/bin/env python3
"""Exact source reconciliation for the repeated-tail (96,144) row.

The 2017 complete-chain table contains

    (8,40) -> (8,28) -> (11/4,7),  (m,n)=(3,2).

The 2016 lower-side paper says, without proving the transition, that this row
leads to the forbidden last lower corner (8,4).  The later 2022 proof for the
companion row beginning at (8,32) supplies the divisibility calculation that
can be reused here.  It does not kill the whole (8,40) row: it reduces its
vertical factor to the one remaining triple-root partition.

The remaining triple-root branch is then passed to the bounded
complete-chain audit.  Its translated edge has no complete-chain escape, so
the repeated-tail row is excluded before Laurent-band compilation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math

from complete_chain_no_escape import triple_root_no_escape_audit


Point = tuple[Fraction, Fraction]


def _primitive_direction(first: Point, second: Point) -> tuple[int, int]:
    """Return the primitive weight constant on the line through two points."""

    dx = first[0] - second[0]
    dy = first[1] - second[1]
    scale = math.lcm(dx.denominator, dy.denominator)
    integral_dx = int(dx * scale)
    integral_dy = int(dy * scale)
    divisor = math.gcd(abs(integral_dx), abs(integral_dy))
    rho, sigma = integral_dy // divisor, -integral_dx // divisor
    if rho < 0:
        rho, sigma = -rho, -sigma
    return rho, sigma


def _partitions(total: int, maximum: int | None = None) -> tuple[tuple[int, ...], ...]:
    """Integer partitions in nonincreasing order."""

    if total == 0:
        return ((),)
    maximum = total if maximum is None else min(total, maximum)
    result: list[tuple[int, ...]] = []
    for first in range(maximum, 0, -1):
        for tail in _partitions(total - first, first):
            result.append((first,) + tail)
    return tuple(result)


@dataclass(frozen=True)
class RootPartitionAudit:
    partition: tuple[int, ...]
    translated_corners: tuple[tuple[int, int], ...]
    has_simple_root: bool
    excluded_by_forbidden_corner: bool
    status: str


@dataclass(frozen=True)
class RepeatedTailAudit:
    degree_pair: tuple[int, int]
    chain: tuple[tuple[str, str], ...]
    final_edge_weight: tuple[int, int]
    forcing_fraction: tuple[int, int]
    q1: int
    d0_upper_bound: int
    d0: int
    primitive_vertical_factor_endpoints: tuple[tuple[int, int], ...]
    residual_factor_degree: int
    root_partitions: tuple[RootPartitionAudit, ...]
    surviving_partitions: tuple[tuple[int, ...], ...]
    surviving_factor: str
    triple_root_open_chain_counts: tuple[int, ...]
    triple_root_complete_chain_count: int
    repeated_tail_excluded: bool
    status: str
    source_boundary: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def repeated_tail_96_144_audit() -> RepeatedTailAudit:
    """Reduce the source conflict to the unique triple-root vertical factor."""

    a0: Point = (Fraction(8), Fraction(40))
    a1: Point = (Fraction(8), Fraction(28))
    a2: Point = (Fraction(11, 4), Fraction(7))

    weight = _primitive_direction(a1, a2)
    assert weight == (4, -1)
    assert weight[0] * a1[0] + weight[1] * a1[1] == 4
    assert weight[0] * a2[0] + weight[1] * a2[1] == 4

    # The 2022 calculation on the common tail gives
    # en(F_1)=(6,21)=(3/4)A_1, hence q_1=4.  The divisibility theorem gives
    # q_1 | d_0, while the primitive endpoint A_1/d_0 bounds
    # d_0 | gcd(8,28)=4.
    forcing_fraction = Fraction(3, 4)
    forced_end = (
        forcing_fraction * a1[0],
        forcing_fraction * a1[1],
    )
    assert forced_end == (6, 21)
    q1 = forcing_fraction.denominator
    d0_upper_bound = math.gcd(int(a1[0]), int(a1[1]))
    assert d0_upper_bound == 4 and d0_upper_bound % q1 == 0
    d0 = q1

    primitive_lower = (int(a1[0] / d0), int(a1[1] / d0))
    primitive_upper = (int(a0[0] / d0), int(a0[1] / d0))
    assert primitive_lower == (2, 7)
    assert primitive_upper == (2, 10)
    residual_degree = primitive_upper[1] - primitive_lower[1]
    assert residual_degree == 3

    partition_audits = []
    for partition in _partitions(residual_degree):
        translated = tuple((8, d0 * multiplicity) for multiplicity in partition)
        has_simple = 1 in partition
        partition_audits.append(
            RootPartitionAudit(
                partition=partition,
                translated_corners=translated,
                has_simple_root=has_simple,
                excluded_by_forbidden_corner=has_simple,
                status=(
                    "excluded: a simple nonzero root translates to (8,4)"
                    if has_simple
                    else "survives this argument: only a triple nonzero root remains"
                ),
            )
        )

    survivors = tuple(
        item.partition
        for item in partition_audits
        if not item.excluded_by_forbidden_corner
    )
    assert survivors == ((3,),)
    no_escape = triple_root_no_escape_audit()
    assert no_escape.excluded

    return RepeatedTailAudit(
        degree_pair=(96, 144),
        chain=(("8", "40"), ("8", "28"), ("11/4", "7")),
        final_edge_weight=weight,
        forcing_fraction=(3, 4),
        q1=q1,
        d0_upper_bound=d0_upper_bound,
        d0=d0,
        primitive_vertical_factor_endpoints=(primitive_lower, primitive_upper),
        residual_factor_degree=residual_degree,
        root_partitions=tuple(partition_audits),
        surviving_partitions=survivors,
        surviving_factor="R=kappa*x^2*y^7*(y-lambda)^3",
        triple_root_open_chain_counts=no_escape.open_chain_counts,
        triple_root_complete_chain_count=no_escape.complete_chain_count,
        repeated_tail_excluded=True,
        status=(
            "excluded: simple-root partitions give (8,4), while the sole "
            "triple-root partition has no complete-chain escape"
        ),
        source_boundary=(
            "The 2016 paper states the full (8,40) transition only in a remark.",
            "The 2017 table is an over-approximation; its Proposition-3.29 filter is commented out.",
            "The 2022 paper proves the analogous transition only for the degree-one residual factor in (8,32).",
            "The local complete-chain audit supplies the missing triple-root argument.",
        ),
    )


if __name__ == "__main__":
    print(json.dumps(repeated_tail_96_144_audit().as_dict(), indent=2, sort_keys=True))

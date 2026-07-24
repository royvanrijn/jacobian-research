#!/usr/bin/env python3
"""Bounded conductor-profile atlas for rational plane image curves.

A rational plane curve of degree c has total delta invariant

    delta_total = (c-1)(c-2)/2

and normalization-conductor degree 2*delta_total.  At a singular point with
r analytic branches, delta is at least binomial(r,2).  These constraints
bound the coarse singularity profiles at fixed c.

The atlas is deliberately numerical.  Passing it does not assert that a
plane singularity, curve, finite cover, or Keller map exists.  In particular,
large delta can be concentrated at a unibranch singularity, so conductor
degree alone does not force a multi-point collision packet.

The minimal profile is not merely formal: the degree-c curve
``y^(c-1) z = x^c`` is rational, smooth at infinity, and has one affine
unibranch singularity of delta ``(c-1)(c-2)/2``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from itertools import product

from finite_normalization_signatures import integer_partitions


@dataclass(frozen=True, order=True)
class TargetSingularityPacket:
    """Coarse local normalization data at one target singular point."""

    delta_invariant: int
    normalization_branches: int

    def __post_init__(self) -> None:
        if self.delta_invariant <= 0:
            raise ValueError("a singularity packet needs positive delta")
        if self.normalization_branches <= 0:
            raise ValueError("a singularity packet needs a positive branch count")
        if self.minimum_delta > self.delta_invariant:
            raise ValueError("the branch count exceeds the delta bound")

    @property
    def minimum_delta(self) -> int:
        branches = self.normalization_branches
        return branches * (branches - 1) // 2

    @property
    def conductor_degree(self) -> int:
        return 2 * self.delta_invariant

    @property
    def collision_excess(self) -> int:
        return self.normalization_branches - 1

    @property
    def is_unibranch(self) -> bool:
        return self.normalization_branches == 1

    @property
    def possibly_normalization_immersive(self) -> bool:
        """A singular unibranch normalization is never immersive."""

        return self.normalization_branches >= 2

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.update(
            minimum_delta=self.minimum_delta,
            conductor_degree=self.conductor_degree,
            collision_excess=self.collision_excess,
            is_unibranch=self.is_unibranch,
            possibly_normalization_immersive=(
                self.possibly_normalization_immersive
            ),
            conductor_weight_profiles=conductor_weight_profiles(
                self.delta_invariant,
                self.normalization_branches,
            ),
        )
        return result


@dataclass(frozen=True)
class TargetConductorProfile:
    """A coarse distribution of the genus defect among singular points."""

    curve_degree: int
    packets: tuple[TargetSingularityPacket, ...]

    def __post_init__(self) -> None:
        if self.curve_degree < 3:
            raise ValueError("a positive conductor profile needs curve degree >= 3")
        if not self.packets:
            raise ValueError("a conductor profile needs a singularity packet")
        if tuple(sorted(self.packets)) != self.packets:
            raise ValueError("singularity packets must be in canonical order")
        if self.total_delta != self.expected_total_delta:
            raise ValueError("the packets must exhaust the rational plane genus defect")

    @property
    def expected_total_delta(self) -> int:
        return (self.curve_degree - 1) * (self.curve_degree - 2) // 2

    @property
    def total_delta(self) -> int:
        return sum(packet.delta_invariant for packet in self.packets)

    @property
    def conductor_degree(self) -> int:
        return 2 * self.total_delta

    @property
    def singular_points(self) -> int:
        return len(self.packets)

    @property
    def normalization_points(self) -> int:
        return sum(packet.normalization_branches for packet in self.packets)

    @property
    def collision_excess(self) -> int:
        return sum(packet.collision_excess for packet in self.packets)

    @property
    def multibranch_points(self) -> int:
        return sum(not packet.is_unibranch for packet in self.packets)

    @property
    def maximum_branch_count(self) -> int:
        return max(packet.normalization_branches for packet in self.packets)

    @property
    def pareto_complexity(self) -> tuple[int, int, int, int, int]:
        return (
            self.curve_degree,
            self.singular_points,
            self.normalization_points,
            self.collision_excess,
            self.maximum_branch_count,
        )

    @property
    def conductor_forces_collision(self) -> bool:
        """Whether this selected profile, not merely its degree, is multibranch."""

        return self.collision_excess > 0

    @property
    def possibly_residue_immersive(self) -> bool:
        return all(
            packet.possibly_normalization_immersive
            for packet in self.packets
        )

    def lifted_conductor_degree(self, residue_degree: int) -> int:
        if residue_degree <= 0:
            raise ValueError("the residue degree must be positive")
        return residue_degree * self.conductor_degree

    def as_dict(self) -> dict[str, object]:
        return {
            "curve_degree": self.curve_degree,
            "packets": [packet.as_dict() for packet in self.packets],
            "total_delta": self.total_delta,
            "conductor_degree": self.conductor_degree,
            "singular_points": self.singular_points,
            "normalization_points": self.normalization_points,
            "collision_excess": self.collision_excess,
            "multibranch_points": self.multibranch_points,
            "maximum_branch_count": self.maximum_branch_count,
            "pareto_complexity": self.pareto_complexity,
            "conductor_forces_collision": self.conductor_forces_collision,
            "possibly_residue_immersive": self.possibly_residue_immersive,
        }


def conductor_weight_profiles(
    delta_invariant: int,
    normalization_branches: int,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate coarse branchwise conductor weights.

    A branch at an r-branch reduced plane singularity receives at least
    r-1 from intersections with the other branches.  The weights sum to
    2*delta.  These are necessary numerical profiles, not realizability
    assertions.
    """

    packet = TargetSingularityPacket(
        delta_invariant,
        normalization_branches,
    )
    minimum = max(1, normalization_branches - 1)
    return tuple(
        partition
        for partition in integer_partitions(packet.conductor_degree, minimum)
        if len(partition) == normalization_branches
    )


def _positive_partitions(total: int) -> tuple[tuple[int, ...], ...]:
    return integer_partitions(total, 1)


def _branch_counts(delta_invariant: int) -> tuple[int, ...]:
    return tuple(
        branches
        for branches in range(1, delta_invariant + 2)
        if branches * (branches - 1) // 2 <= delta_invariant
    )


def enumerate_target_conductor_profiles(
    curve_degree: int,
) -> tuple[TargetConductorProfile, ...]:
    """Enumerate all coarse delta/branch profiles at fixed curve degree."""

    if curve_degree < 3:
        return ()
    total_delta = (curve_degree - 1) * (curve_degree - 2) // 2
    profiles: set[TargetConductorProfile] = set()
    for delta_partition in _positive_partitions(total_delta):
        branch_options = tuple(
            _branch_counts(delta) for delta in delta_partition
        )
        for counts in product(*branch_options):
            packets = tuple(
                sorted(
                    TargetSingularityPacket(delta, branches)
                    for delta, branches in zip(delta_partition, counts)
                )
            )
            profiles.add(TargetConductorProfile(curve_degree, packets))
    return tuple(
        sorted(
            profiles,
            key=lambda profile: (
                profile.pareto_complexity,
                profile.packets,
            ),
        )
    )


def pareto_minimal_conductor_profiles(
    profiles: tuple[TargetConductorProfile, ...],
) -> tuple[TargetConductorProfile, ...]:
    """Keep profiles minimal in the four non-degree complexity coordinates."""

    def dominates(
        left: TargetConductorProfile,
        right: TargetConductorProfile,
    ) -> bool:
        left_complexity = left.pareto_complexity[1:]
        right_complexity = right.pareto_complexity[1:]
        return all(
            x <= y for x, y in zip(left_complexity, right_complexity)
        ) and any(
            x < y for x, y in zip(left_complexity, right_complexity)
        )

    return tuple(
        profile
        for profile in profiles
        if not any(
            dominates(other, profile)
            for other in profiles
            if other != profile
        )
    )


def concentrated_unibranch_profile(
    curve_degree: int,
) -> TargetConductorProfile:
    """Return the unique coordinatewise-minimal conductor profile.

    Every profile has at least one singular point, at least one normalization
    point, nonnegative collision excess, and maximum branch count at least
    one.  The curve ``y^(c-1)z=x^c`` concentrates the full genus defect at
    one unibranch point, attains ``(1,1,0,1)``, and hence dominates every
    other coarse profile.
    """

    if curve_degree < 3:
        raise ValueError("a positive conductor profile needs curve degree >= 3")
    total_delta = (curve_degree - 1) * (curve_degree - 2) // 2
    return TargetConductorProfile(
        curve_degree,
        (TargetSingularityPacket(total_delta, 1),),
    )


def concentrated_two_branch_profile(
    curve_degree: int,
) -> TargetConductorProfile:
    """Return the unique coarse Pareto minimum compatible with immersion.

    Once singular unibranch packets are forbidden, every singularity packet
    has at least two normalization branches.  Concentrating the full delta
    invariant at one two-branch point attains complexity ``(1,2,1,2)``.
    The two smooth branches may have arbitrarily high contact, so delta does
    not force more than these two normalization points.
    """

    if curve_degree < 3:
        raise ValueError("a positive conductor profile needs curve degree >= 3")
    total_delta = (curve_degree - 1) * (curve_degree - 2) // 2
    return TargetConductorProfile(
        curve_degree,
        (TargetSingularityPacket(total_delta, 2),),
    )


def residue_immersive_conductor_profiles(
    profiles: tuple[TargetConductorProfile, ...],
) -> tuple[TargetConductorProfile, ...]:
    """Discard profiles already impossible for an immersive normalization."""

    return tuple(
        profile for profile in profiles if profile.possibly_residue_immersive
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("degree", type=int)
    parser.add_argument("--pareto", action="store_true")
    args = parser.parse_args()
    profiles = enumerate_target_conductor_profiles(args.degree)
    if args.pareto:
        profiles = pareto_minimal_conductor_profiles(profiles)
    print(
        json.dumps(
            [profile.as_dict() for profile in profiles],
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

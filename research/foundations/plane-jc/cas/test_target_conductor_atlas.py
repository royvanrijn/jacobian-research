#!/usr/bin/env python3
"""Regression tests for rational-plane-curve conductor profiles."""

from target_conductor_atlas import (
    TargetSingularityPacket,
    concentrated_two_branch_profile,
    concentrated_unibranch_profile,
    conductor_weight_profiles,
    enumerate_target_conductor_profiles,
    pareto_minimal_conductor_profiles,
    residue_immersive_conductor_profiles,
)


assert enumerate_target_conductor_profiles(1) == ()
assert enumerate_target_conductor_profiles(2) == ()

# A rational cubic has delta one: its coarse alternatives are a unibranch
# cusp and a two-branch node.
cubic = enumerate_target_conductor_profiles(3)
assert len(cubic) == 2
assert tuple(
    tuple(
        (packet.delta_invariant, packet.normalization_branches)
        for packet in profile.packets
    )
    for profile in cubic
) == (((1, 1),), ((1, 2),))
assert conductor_weight_profiles(1, 1) == ((2,),)
assert conductor_weight_profiles(1, 2) == ((1, 1),)

# Degree four has total delta three.  The atlas includes concentrated,
# split, unibranch, and multibranch profiles without asserting existence.
quartic = enumerate_target_conductor_profiles(4)
assert len(quartic) == 11
assert all(profile.total_delta == 3 for profile in quartic)
assert all(profile.conductor_degree == 6 for profile in quartic)
assert any(profile.singular_points == 3 for profile in quartic)
assert any(profile.maximum_branch_count == 3 for profile in quartic)

# Low degrees can be enumerated exhaustively.  Conductor weight can always
# be concentrated at one unibranch point, so it does not by itself force a
# collision of distinct normalization points.
for curve_degree in range(3, 6):
    profiles = enumerate_target_conductor_profiles(curve_degree)
    pareto = pareto_minimal_conductor_profiles(profiles)
    total_delta = (curve_degree - 1) * (curve_degree - 2) // 2
    assert len(pareto) == 1
    assert pareto[0].packets == (
        TargetSingularityPacket(total_delta, 1),
    )
    assert not pareto[0].conductor_forces_collision
    assert pareto[0].lifted_conductor_degree(4) == (
        4 * (curve_degree - 1) * (curve_degree - 2)
    )
    immersive = residue_immersive_conductor_profiles(profiles)
    immersive_pareto = pareto_minimal_conductor_profiles(immersive)
    assert immersive_pareto == (
        concentrated_two_branch_profile(curve_degree),
    )

# The coordinatewise-minimal face has a direct construction and therefore
# scales without enumerating the rapidly growing dominated atlas.
for curve_degree in range(3, 126):
    concentrated = concentrated_unibranch_profile(curve_degree)
    assert concentrated.pareto_complexity[1:] == (1, 1, 0, 1)
    assert not concentrated.conductor_forces_collision
    assert not concentrated.possibly_residue_immersive
    immersive_minimum = concentrated_two_branch_profile(curve_degree)
    assert immersive_minimum.pareto_complexity[1:] == (1, 2, 1, 2)
    assert immersive_minimum.possibly_residue_immersive
    assert immersive_minimum.conductor_forces_collision

# At a three-branch delta-three point, pairwise intersections already use
# the entire genus defect and force the ordinary (2,2,2) weight profile.
assert conductor_weight_profiles(3, 3) == ((2, 2, 2),)

print("PASS: conductor profiles exhaust rational plane degrees three through five")
print("PASS: cubic profiles separate the cusp and node faces")
print("PASS: total conductor alone never forces a multibranch collision")
print("PASS: the concentrated unibranch face scales through degree 125")
print("PASS: residue immersion moves the Pareto minimum to two branches")

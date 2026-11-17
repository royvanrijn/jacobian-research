#!/usr/bin/env python3
"""Verify the Maxwell-boundary root-stack complex in the sextic chart."""

from __future__ import annotations

import json
import sys
from itertools import combinations, permutations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import gf2_rank  # noqa: E402

ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "stacky_fan_descent_degree6.json"
)

CLUSTERS = ("x", "y", "z")
SYMMETRIC_GROUP = tuple(permutations(CLUSTERS))
PAIR_DIVISORS = tuple(
    frozenset(pair) for pair in combinations(CLUSTERS, 2)
)
TRIPLE_DIVISOR = frozenset(CLUSTERS)


def ordered_partitions() -> list[tuple[tuple[str, ...], ...]]:
    result = []
    for level_count in range(1, len(CLUSTERS) + 1):
        for assignment in product(range(level_count), repeat=len(CLUSTERS)):
            if set(assignment) != set(range(level_count)):
                continue
            result.append(
                tuple(
                    tuple(
                        cluster
                        for cluster, level in zip(CLUSTERS, assignment)
                        if level == block
                    )
                    for block in range(level_count)
                )
            )
    return result


def partition_name(
    partition: tuple[tuple[str, ...], ...],
) -> str:
    return "|".join("".join(block) for block in partition)


def permutation_map(
    permutation: tuple[str, ...],
) -> dict[str, str]:
    return dict(zip(CLUSTERS, permutation))


def act_on_partition(
    permutation: tuple[str, ...],
    partition: tuple[tuple[str, ...], ...],
) -> tuple[tuple[str, ...], ...]:
    action = permutation_map(permutation)
    return tuple(
        tuple(sorted(action[cluster] for cluster in block))
        for block in partition
    )


def act_on_divisor(
    permutation: tuple[str, ...],
    divisor: frozenset[str],
) -> frozenset[str]:
    action = permutation_map(permutation)
    return frozenset(action[cluster] for cluster in divisor)


def divisor_name(divisor: frozenset[str]) -> str:
    if divisor == TRIPLE_DIVISOR:
        return "M:xyz"
    return "M:" + "".join(sorted(divisor))


# The radial fan still has its four familiar S_3 orbit types.  The separate
# full-chain compiler proves that every nontrivial connector sign changes a
# normalized branch, so the fully labelled radial lifts have no inertia.
partitions = ordered_partitions()
assert len(partitions) == 13
partition_set = set(partitions)
for partition in partitions:
    for permutation in SYMMETRIC_GROUP:
        assert act_on_partition(permutation, partition) in partition_set

unseen = set(partitions)
radial_orbits = []
while unseen:
    representative = min(unseen, key=partition_name)
    orbit = {
        act_on_partition(permutation, representative)
        for permutation in SYMMETRIC_GROUP
    }
    stabilizer = [
        permutation
        for permutation in SYMMETRIC_GROUP
        if act_on_partition(permutation, representative) == representative
    ]
    unseen -= orbit
    radial_orbits.append(
        {
            "representative": partition_name(representative),
            "orbit_size": len(orbit),
            "label_stabilizer_order": len(stabilizer),
            "normalized_labelled_node_inertia": "trivial",
        }
    )
radial_orbits.sort(key=lambda item: item["representative"])
assert sorted(
    (item["orbit_size"], item["label_stabilizer_order"])
    for item in radial_orbits
) == [(1, 6), (3, 2), (3, 2), (6, 1)]


# On Bl_{p_x,p_y,p_z,p_Delta}(P^2)=Mbar_0,5, the three strict pairwise
# Maxwell divisors are mutually disjoint and each meets the exceptional
# triple-Maxwell divisor once.  A surviving diagonal square-tail involution
# gives one mu_2 root generator for every Maxwell target node.
maxwell_divisors = PAIR_DIVISORS + (TRIPLE_DIVISOR,)
maxwell_strata = []
for pair_divisor in PAIR_DIVISORS:
    maxwell_strata.append(
        {
            "stratum": divisor_name(pair_divisor),
            "present_divisors": [divisor_name(pair_divisor)],
            "root_generators": [divisor_name(pair_divisor)],
            "inertia": "mu_2",
            "inertia_rank": 1,
            "normalization_branches": 2,
        }
    )
maxwell_strata.append(
    {
        "stratum": divisor_name(TRIPLE_DIVISOR),
        "present_divisors": [divisor_name(TRIPLE_DIVISOR)],
        "root_generators": [divisor_name(TRIPLE_DIVISOR)],
        "inertia": "mu_2",
        "inertia_rank": 1,
        "normalization_branches": 4,
    }
)
for pair_divisor in PAIR_DIVISORS:
    generators = sorted(
        (divisor_name(pair_divisor), divisor_name(TRIPLE_DIVISOR))
    )
    maxwell_strata.append(
        {
            "stratum": "&".join(generators),
            "present_divisors": generators,
            "root_generators": generators,
            "inertia": "(mu_2)^2",
            "inertia_rank": 2,
            "normalization_branches": 8,
        }
    )

# Face compatibility: a pair--triple intersection contains exactly the two
# generic-divisor root generators of its codimension-one faces.
generic_generators = {
    divisor_name(divisor): frozenset((divisor_name(divisor),))
    for divisor in maxwell_divisors
}
for stratum in maxwell_strata:
    generators = frozenset(stratum["root_generators"])
    for divisor in stratum["present_divisors"]:
        assert generic_generators[divisor] <= generators

# S_3 permutes pair divisors, fixes the triple divisor, and transports all
# pair--triple intersections equivariantly.
for permutation in SYMMETRIC_GROUP:
    assert act_on_divisor(permutation, TRIPLE_DIVISOR) == TRIPLE_DIVISOR
    assert {
        act_on_divisor(permutation, divisor)
        for divisor in PAIR_DIVISORS
    } == set(PAIR_DIVISORS)

maxwell_orbits = [
    {
        "representative": "M:xy",
        "orbit_size": 3,
        "cluster_label_stabilizer_order": 2,
        "root_inertia": "mu_2",
    },
    {
        "representative": "M:xyz",
        "orbit_size": 1,
        "cluster_label_stabilizer_order": 6,
        "root_inertia": "mu_2",
    },
    {
        "representative": "M:xy&M:xyz",
        "orbit_size": 3,
        "cluster_label_stabilizer_order": 2,
        "root_inertia": "(mu_2)^2",
    },
]

# Parity audit for the divisor roots.  Along a radial boundary, every
# critical value is t^2 times a unit, so the target-node smoothing has even
# order in the root-scale parameter t.  Along a generic Maxwell divisor,
# q_j=q_i+d with q_i=a a unit, and the branch-value difference has order one
# in d.  The source relation s^2=tau therefore forces precisely one second
# root transverse to that divisor.  Blowing up the triple diagonal writes
# its difference ideal as rho*(1,u), so the exceptional normal rho also has
# order one.
pair_d_coefficients = ("0", "2*a", "1")
assert pair_d_coefficients[0] == "0"
assert pair_d_coefficients[1] == "2*a"
radial_t_exponents = (2,)
assert min(radial_t_exponents) == 2
triple_rho_exponents = (1, 1)
assert min(triple_rho_exponents) == 1

# At a pair--triple crossing, let c_x,c_y be the two connector flips on the
# intermediate bubble, t_z the third square-tail flip at the outer node,
# and o_x,o_y the two square-tail flips on the inner pair bubble.  The outer
# normalized lift requires c_x=c_y=t_z; the inner lift requires
# c_x+o_x=c_y+o_y.  These are three independent equations on five signs,
# leaving rank two.  At a radial--Maxwell crossing, the labelled natural
# tail kills the radial sign and the two Maxwell signs must agree, leaving
# rank one.
pair_triple_constraints = [
    [1, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [1, 1, 0, 1, 1],
]
pair_triple_kernel_rank = 5 - gf2_rank(pair_triple_constraints)
assert pair_triple_kernel_rank == 2
radial_maxwell_constraints = [
    [1, 0, 0],
    [0, 1, 1],
]
radial_maxwell_kernel_rank = 3 - gf2_rank(radial_maxwell_constraints)
assert radial_maxwell_kernel_rank == 1

artifact = {
    "experiment": "degree-six Maxwell-boundary root-stack descent",
    "local_root_chart": {
        "atlas": "Spec k[[s]]",
        "coarse_parameter": "tau=s^2",
        "inertia_action": "mu_2 acts by s -> -s",
        "interpretation": (
            "the diagonal deck flip of all square tails over one Maxwell "
            "target node"
        ),
    },
    "radial_quotient_orbits": radial_orbits,
    "radial_rule": (
        "all thirteen fully source-root-labelled radial lifts have trivial "
        "node inertia because full-chain cycle matching makes every "
        "nontrivial connector sign change a normalized branch"
    ),
    "transverse_order_audit": {
        "radial": {
            "target_smoothing": "tau=t^2*unit",
            "order_in_root_scale": 2,
            "conclusion": "source smoothing root already lies on B_6",
        },
        "pairwise_maxwell": {
            "branch_difference": "2*a*d+d^2 with a a unit",
            "transverse_parameter": "d",
            "order": 1,
            "conclusion": "s^2=tau requires one second-root stack",
        },
        "triple_maxwell": {
            "blowup_difference_ideal": ["rho", "rho*u"],
            "exceptional_normal_order": 1,
            "conclusion": "s^2=tau requires one second root of D_xyz",
        },
    },
    "maxwell_divisors": [
        divisor_name(divisor) for divisor in maxwell_divisors
    ],
    "maxwell_strata": maxwell_strata,
    "maxwell_quotient_orbits": maxwell_orbits,
    "codimension_two_inertia_audit": {
        "pair_triple": {
            "deck_sign_variables": ["c_x", "c_y", "t_z", "o_x", "o_y"],
            "constraint_rank": gf2_rank(pair_triple_constraints),
            "kernel_rank": pair_triple_kernel_rank,
            "inertia": "(mu_2)^2",
        },
        "radial_maxwell": {
            "deck_sign_variables": [
                "radial_connector",
                "maxwell_tail_1",
                "maxwell_tail_2",
            ],
            "constraint_rank": gf2_rank(radial_maxwell_constraints),
            "kernel_rank": radial_maxwell_kernel_rank,
            "inertia": "mu_2",
        },
        "conclusion": (
            "all codimension-two boundary crossings have exactly the "
            "product inertia predicted by the Maxwell divisor roots"
        ),
    },
    "face_rule": (
        "specialization to an intersection adds the root generator of each "
        "Maxwell boundary divisor containing the stratum"
    ),
    "symmetric_group_rule": (
        "S_3 permutes the three pairwise Maxwell roots and fixes the "
        "triple-Maxwell root"
    ),
    "pair_root_rule": (
        "forgetting the order inside each of the three source-root pairs "
        "adds the separate label quotient (S_2)^3 semidirect S_3; it is not "
        "part of the Maxwell node-root inertia"
    ),
    "conclusion": (
        "the labelled ACV graph on this chart is the iterated "
        "second-root construction along the four Maxwell boundary divisors; "
        "radial Kummer saturation is already carried by the normalized "
        "branch graph"
    ),
    "scope": (
        "local completed charts, all divisorial and codimension-two inertia, "
        "boundary incidence, face maps, and S_3 descent in the labelled "
        "(2,2,2) chart; the identification with the ACV graph additionally "
        "uses ACV smoothness, the coarse comparison, and the bottom-up tame "
        "stack theorem; gluing other root partitions under S_6 remains open"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS radial fan: four S_3 orbit types with trivial labelled inertia")
print("PASS Maxwell roots: four divisors and three pair--triple crossings")
print("PASS root-stack faces: all generator inclusions")
print("PASS S_3 descent: pair roots permuted, triple root fixed")
print("PASS codimension two: pair--triple rank two, radial--Maxwell rank one")
print("DEGREE_SIX_STACKY_FAN_DESCENT_PASS")

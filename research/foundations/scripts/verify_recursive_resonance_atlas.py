#!/usr/bin/env python3
"""Verify recursive residue screens, contractions, and extracted inertia."""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.monodromy_forests import (  # noqa: E402
    compose,
    identity,
    reduced_cycle_factorizations,
    transposition,
)
from jcsearch.resonance_atlas import (  # noqa: E402
    affine_screen_composition,
    compile_nested_monodromy_characters,
    compile_resonance_characters,
    framed_screen_chart,
    normalized_polynomial_pullback,
    ramified_initial_form,
    screen_contraction_substitution,
    screen_frame_change_substitution,
)
from jcsearch.radial_sources import (  # noqa: E402
    ordered_set_partitions,
    radial_inertia_order,
)
from jcsearch.wonderful_branch import (  # noqa: E402
    is_nested,
    maximal_nested_sets,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "recursive_resonance_atlas.json"
)


def all_laminar_families(leaf_count: int):
    """Recover the branch-only nested complex from maximal Mbar strata."""

    branch_marks = frozenset(range(1, leaf_count + 1))
    families = {()}
    for maximal in maximal_nested_sets(leaf_count + 2):
        branch_clusters = tuple(
            sorted(
                (
                    frozenset(mark - 1 for mark in boundary)
                    for boundary in maximal
                    if boundary <= branch_marks
                ),
                key=lambda value: (len(value), tuple(sorted(value))),
            )
        )
        for mask in range(1 << len(branch_clusters)):
            family = tuple(
                branch_clusters[index]
                for index in range(len(branch_clusters))
                if mask & (1 << index)
            )
            families.add(family)
    return tuple(
        sorted(
            families,
            key=lambda family: (
                len(family),
                tuple(
                    (len(cluster), tuple(sorted(cluster)))
                    for cluster in family
                ),
            ),
        )
    )


def interval_laminar_families(edge_count: int):
    intervals = tuple(
        frozenset(range(left, right + 1))
        for left in range(edge_count)
        for right in range(left + 1, edge_count)
    )
    result = []
    for size in range(1, edge_count):
        result.extend(
            family
            for family in combinations(intervals, size)
            if is_nested(family)
        )
    return tuple(result)


def integer_partitions(total: int, maximum: int | None = None):
    if total == 0:
        return ((),)
    maximum = total if maximum is None else min(maximum, total)
    return tuple(
        (first,) + rest
        for first in range(maximum, 0, -1)
        for rest in integer_partitions(total - first, first)
    )


def permutation_with_cycle_profile(profile: tuple[int, ...]):
    degree = sum(profile)
    result = list(range(degree))
    offset = 0
    for length in profile:
        for index in range(length):
            result[offset + index] = offset + ((index + 1) % length)
        offset += length
    return tuple(result)


def verify_compiled_groups(atlas) -> int:
    """Check pointwise composition on the global and inertia subgroups."""

    for elements in (
        atlas.global_automorphisms,
        atlas.inertia_automorphisms,
    ):
        element_set = frozenset(elements)
        assert len(element_set) == len(elements)
        identity_element = tuple(
            identity(atlas.degree)
            for _ in elements[0]
        )
        assert identity_element in element_set
        for left in elements:
            for right in elements:
                composed = tuple(
                    compose(left_vertex, right_vertex)
                    for left_vertex, right_vertex in zip(left, right)
                )
                assert composed in element_set
    return 1


group_law_atlases_checked = 0


# Every framed wonderful chart has an exact polynomial blowdown, and deleting
# any screen is affine composition.  Check every branch-only nested family
# through five moving branch labels, including all two-step contractions.
screen_family_counts = {}
screen_frame_equivariance_checks = 0
screen_contractions = 0
two_step_contractions = 0
for leaf_count in range(2, 6):
    families = all_laminar_families(leaf_count)
    screen_family_counts[str(leaf_count)] = len(families)
    for family in families:
        chart = framed_screen_chart(range(leaf_count), family)
        screens = (None, *chart.family)
        translations = {
            screen: sp.Integer(index + 1)
            for index, screen in enumerate(screens)
        }
        dilations = {
            screen: sp.Integer(index + 2)
            for index, screen in enumerate(screens)
        }
        frame_change = screen_frame_change_substitution(
            chart,
            translations,
            dilations,
        )
        for leaf in chart.leaves:
            transformed = chart.leaf_expressions[leaf].xreplace(
                frame_change
            )
            assert sp.cancel(
                transformed
                - (
                    translations[None]
                    + dilations[None] * chart.leaf_expressions[leaf]
                )
            ) == 0
        screen_frame_equivariance_checks += 1

        for collision in family:
            coarser, substitution = screen_contraction_substitution(
                chart,
                collision,
            )
            for leaf in chart.leaves:
                contracted = coarser.leaf_expressions[leaf].xreplace(
                    substitution
                )
                assert sp.expand(
                    contracted - chart.leaf_expressions[leaf]
                ) == 0
            screen_contractions += 1

        for first, second in permutations(family, 2):
            once, first_substitution = screen_contraction_substitution(
                chart,
                first,
            )
            twice, second_substitution = screen_contraction_substitution(
                once,
                second,
            )
            for leaf in chart.leaves:
                contracted = (
                    twice.leaf_expressions[leaf]
                    .xreplace(second_substitution)
                    .xreplace(first_substitution)
                )
                assert sp.expand(
                    contracted - chart.leaf_expressions[leaf]
                ) == 0
            two_step_contractions += 1


# Normalized initial forms put every target-flag divisor on the source
# screen by one equation.  Exhaust ramification orders through degree seven.
W, U, V, source_scale, inner_scale = sp.symbols("W U V r s")
initial_form_checks = 0
deformation_initial_form_checks = 0
fiber_equation_checks = 0
source_frame_equivariance_checks = 0
transitivity_checks = 0
for degree in range(1, 8):
    for ramification_index in range(1, degree + 1):
        for point in (-2, 0, 3):
            shifted = W - point
            unit_polynomial = sum(
                (power + 2) * shifted**power
                for power in range(degree - ramification_index + 1)
            )
            polynomial = sp.expand(
                5 + shifted**ramification_index * unit_polynomial
            )
            family, initial = ramified_initial_form(
                polynomial,
                W,
                point,
                ramification_index,
                source_scale,
                U,
            )
            assert initial == 2 * U**ramification_index
            assert sp.Poly(family, U).degree() == degree
            initial_form_checks += 1

            source_translation = sp.Integer(3)
            source_dilation = sp.Integer(5)
            target_translation = sp.Integer(11)
            target_dilation = sp.Integer(7)
            target_flag = sp.Symbol("a")
            transformed_map = sp.expand(
                target_translation
                + target_dilation
                * family.subs(
                    U,
                    (V - source_translation) / source_dilation,
                )
            )
            transformed_flag = (
                target_translation + target_dilation * target_flag
            )
            assert sp.cancel(
                transformed_map
                - transformed_flag
                - target_dilation
                * (
                    family.subs(
                        U,
                        (V - source_translation) / source_dilation,
                    )
                    - target_flag
                )
            ) == 0
            source_frame_equivariance_checks += 1

            # A smoothing need not factor by (W-point)^e over the base.
            # Lower Taylor terms are allowed when their weighted orders make
            # the normalized quotient regular and make them vanish in the
            # special initial form.
            if ramification_index > 1:
                deformed_polynomial = sp.expand(
                    polynomial
                    + sum(
                        (power + 7)
                        * source_scale
                        ** (ramification_index - power + 1)
                        * shifted**power
                        for power in range(1, ramification_index)
                    )
                )
                deformed_family, deformed_initial = ramified_initial_form(
                    deformed_polynomial,
                    W,
                    point,
                    ramification_index,
                    source_scale,
                    U,
                )
                assert deformed_initial == initial
                assert sp.expand(deformed_family - family) != 0
                deformation_initial_form_checks += 1

            moving_target_value = (
                polynomial.subs(W, point)
                + source_scale**ramification_index * target_flag
            )
            fiber_identity = sp.expand(
                source_scale**ramification_index
                * (family - target_flag)
                - (
                    polynomial.subs(
                        W,
                        point + source_scale * U,
                    )
                    - moving_target_value
                )
            )
            assert fiber_identity == 0
            fiber_equation_checks += 1

            inner_center = sp.Integer(1)
            outer_target_center = polynomial.subs(W, point)
            outer_target_scale = source_scale**ramification_index
            outer_map = normalized_polynomial_pullback(
                polynomial,
                W,
                source_center=point,
                source_scale=source_scale,
                source_coordinate=U,
                target_center=outer_target_center,
                target_scale=outer_target_scale,
            )
            inner_target_center = outer_map.subs(U, inner_center)
            recursive_map = normalized_polynomial_pullback(
                outer_map,
                U,
                source_center=inner_center,
                source_scale=inner_scale,
                source_coordinate=V,
                target_center=inner_target_center,
                target_scale=inner_scale,
            )
            composed_center, composed_scale = affine_screen_composition(
                point,
                source_scale,
                inner_center,
                inner_scale,
            )
            direct_map = normalized_polynomial_pullback(
                polynomial,
                W,
                source_center=composed_center,
                source_scale=composed_scale,
                source_coordinate=V,
                target_center=(
                    outer_target_center
                    + outer_target_scale * inner_target_center
                ),
                target_scale=outer_target_scale * inner_scale,
            )
            assert sp.cancel(recursive_map - direct_map) == 0
            transitivity_checks += 1


# For a single collision, automatic full-centralizer extraction must recover
# the closed anchored/unanchored formula on every reduced polynomial
# factorization through degree five.
single_node_character_checks = 0
for degree in range(3, 6):
    edge_count = degree - 1
    intervals = tuple(
        frozenset(range(left, right + 1))
        for left in range(edge_count)
        for right in range(left + 1, edge_count)
    )
    for edge_sequence in reduced_cycle_factorizations(degree):
        for collision in intervals:
            atlas = compile_resonance_characters(
                degree,
                edge_sequence,
                (collision,),
            )
            common_index = math.lcm(*atlas.node_profiles[0])
            constrained = tuple(
                index
                for index in atlas.node_profiles[0]
                if index != 2
            )
            constrained_lcm = (
                math.lcm(*constrained)
                if constrained
                else 1
            )
            assert atlas.inertia_order == common_index // constrained_lcm

            anchored = compile_resonance_characters(
                degree,
                edge_sequence,
                (collision,),
                anchored_vertices=(collision,),
            )
            assert anchored.inertia_order == 1
            single_node_character_checks += 1


# The compiler itself is not restricted to transpositions.  Check every
# ordered pair of conjugacy-class representatives through degree six.
general_cycle_profile_checks = 0
for degree in range(2, 7):
    representatives = tuple(
        permutation_with_cycle_profile(profile)
        for profile in integer_partitions(degree)
    )
    collision = frozenset({0, 1})
    for left in representatives:
        for right in representatives:
            atlas = compile_nested_monodromy_characters(
                degree,
                (left, right),
                (collision,),
            )
            assert atlas.inertia_order <= atlas.global_automorphism_count
            anchored = compile_nested_monodromy_characters(
                degree,
                (left, right),
                (collision,),
                anchored_vertices=(collision,),
            )
            assert anchored.inertia_order == 1
            general_cycle_profile_checks += 1

zero_monodromy = identity(6)
for edge in ((0, 1), (2, 3), (4, 5)):
    zero_monodromy = compose(
        transposition(6, *edge),
        zero_monodromy,
    )
mixed_family = (
    frozenset({0, 1}),
    frozenset({0, 1, 2}),
)
mixed_atlas = compile_nested_monodromy_characters(
    6,
    (
        zero_monodromy,
        transposition(6, 0, 2),
        transposition(6, 1, 4),
    ),
    mixed_family,
    anchored_vertices=(mixed_family[0],),
    root_anchored=False,
)
assert mixed_atlas.node_profiles == ((4, 2), (6,))
assert mixed_atlas.inertia_order == 6
group_law_atlases_checked += verify_compiled_groups(mixed_atlas)


# Unequal cluster multiplicities reveal genuine radial stack structure.
# For an ordered scale partition B_0|...|B_k, put M_j=lcm(B_j) and
# L_j=lcm(B_j union ... union B_k).  Full-chain matching leaves
# product_j L_j/M_j normalized inertia.  Compare that closed formula with
# the automatic monodromy compiler on every ordered profile of total degree
# at most seven with two or three clusters and multiplicities two through
# four.
general_radial_profiles = tuple(
    multiplicities
    for cluster_count in (2, 3)
    for multiplicities in product(
        range(2, 5),
        repeat=cluster_count,
    )
    if sum(multiplicities) <= 7
)
general_radial_charts = 0
general_radial_inertia_distribution = Counter()
unequal_radial_examples = {}
for multiplicities in general_radial_profiles:
    degree = sum(multiplicities)
    sheets_by_cluster = []
    offset = 0
    moving_edges_by_cluster = []
    for multiplicity in multiplicities:
        sheets = tuple(range(offset, offset + multiplicity))
        sheets_by_cluster.append(sheets)
        moving_edges_by_cluster.append(
            tuple(
                (sheets[index], sheets[index + 1])
                for index in range(multiplicity - 1)
            )
        )
        offset += multiplicity
    stationary_edges = tuple(
        (
            sheets_by_cluster[index][-1],
            sheets_by_cluster[index + 1][0],
        )
        for index in range(len(multiplicities) - 1)
    )

    for blocks in ordered_set_partitions(
        tuple(range(len(multiplicities)))
    ):
        inner_to_outer_blocks = tuple(reversed(blocks))
        branch_permutations = [identity(degree)]
        cumulative = {0}
        position = 1
        radial_family = []
        for block in inner_to_outer_blocks:
            for cluster in block:
                for edge in moving_edges_by_cluster[cluster]:
                    branch_permutations.append(
                        transposition(degree, *edge)
                    )
                    cumulative.add(position)
                    position += 1
            radial_family.append(frozenset(cumulative))
        branch_permutations.extend(
            transposition(degree, *edge)
            for edge in stationary_edges
        )
        compiled = compile_nested_monodromy_characters(
            degree,
            tuple(branch_permutations),
            tuple(radial_family),
            anchored_vertices=(radial_family[0],),
            root_anchored=False,
        )
        predicted_order = radial_inertia_order(
            multiplicities,
            blocks,
        )
        assert compiled.inertia_order == predicted_order
        group_law_atlases_checked += verify_compiled_groups(compiled)
        general_radial_inertia_distribution[
            str(predicted_order)
        ] += 1
        general_radial_charts += 1

        example_key = (
            multiplicities,
            blocks,
        )
        if example_key in {
            ((2, 3), ((0,), (1,))),
            ((2, 3), ((1,), (0,))),
            ((2, 4), ((0,), (1,))),
        }:
            unequal_radial_examples[
                (
                    ",".join(map(str, multiplicities))
                    + ":"
                    + "|".join(
                        ",".join(map(str, block))
                        for block in blocks
                    )
                )
            ] = predicted_order

assert unequal_radial_examples == {
    "2,3:0|1": 3,
    "2,3:1|0": 2,
    "2,4:0|1": 2,
}


# The actual polynomial radial chart has an identity target-zero monodromy,
# three moving matching edges, and two stationary edges rigidifying the root
# vertex.  Put the zero mark on the deepest screen and exhaust all thirteen
# ordered scale types.  Special-fiber connector automorphisms can exist, but
# every nontrivial one fails at least one simultaneous diagonal condition.
moving_edges = ((0, 1), (2, 3), (4, 5))
stationary_edges = ((1, 2), (3, 4))
degree_six_radial_types = 0
degree_six_radial_global_automorphisms = Counter()
for blocks in ordered_set_partitions((0, 1, 2)):
    inner_to_outer_blocks = tuple(reversed(blocks))
    moving_order = tuple(
        label
        for block in inner_to_outer_blocks
        for label in block
    )
    branch_permutations = (
        identity(6),
        *(
            transposition(6, *moving_edges[label])
            for label in moving_order
        ),
        *(
            transposition(6, *edge)
            for edge in stationary_edges
        ),
    )
    cumulative = {0}
    position = 1
    radial_family = []
    for block in inner_to_outer_blocks:
        cumulative.update(range(position, position + len(block)))
        position += len(block)
        radial_family.append(frozenset(cumulative))
    radial = compile_nested_monodromy_characters(
        6,
        tuple(branch_permutations),
        tuple(radial_family),
        anchored_vertices=(radial_family[0],),
        root_anchored=False,
    )
    assert radial_inertia_order((2, 2, 2), blocks) == 1
    assert radial.inertia_order == 1
    group_law_atlases_checked += verify_compiled_groups(radial)
    degree_six_radial_global_automorphisms[
        str(radial.global_automorphism_count)
    ] += 1
    degree_six_radial_types += 1
assert degree_six_radial_types == 13


# A degree-six tree with its first three edges disjoint contains the complete
# nested pair/triple Maxwell test.  Exhaust every interval nested family.
matching_tree = (
    (0, 1),
    (2, 3),
    (4, 5),
    (1, 2),
    (3, 4),
)
degree_six_families = interval_laminar_families(5)
degree_six_inertia_distribution = Counter()
max_global_automorphisms = 0
for family in degree_six_families:
    atlas = compile_resonance_characters(
        6,
        matching_tree,
        family,
    )
    degree_six_inertia_distribution[str(atlas.inertia_order)] += 1
    group_law_atlases_checked += verify_compiled_groups(atlas)
    max_global_automorphisms = max(
        max_global_automorphisms,
        atlas.global_automorphism_count,
    )
    anchored = compile_resonance_characters(
        6,
        matching_tree,
        family,
        anchored_vertices=family,
    )
    assert anchored.inertia_order == 1

pair_triple_family = (
    frozenset({0, 1}),
    frozenset({0, 1, 2}),
)
pair_triple = compile_resonance_characters(
    6,
    matching_tree,
    pair_triple_family,
)
assert pair_triple.node_profiles == (
    (2, 2, 1, 1),
    (2, 2, 2),
)
assert pair_triple.inertia_order == 4
assert all(
    all(
        compose(vertex, vertex) == identity(pair_triple.degree)
        for vertex in automorphism
    )
    for automorphism in pair_triple.inertia_automorphisms
)
group_law_atlases_checked += verify_compiled_groups(pair_triple)


# Relabelling the six sheets conjugates every local monodromy group and
# leaves the extracted global character count unchanged.
label_equivariance_checks = 0
for relabelling in permutations(range(6)):
    relabelled_tree = tuple(
        (
            relabelling[left],
            relabelling[right],
        )
        for left, right in matching_tree
    )
    relabelled = compile_resonance_characters(
        6,
        relabelled_tree,
        pair_triple_family,
    )
    assert relabelled.inertia_order == pair_triple.inertia_order
    assert (
        relabelled.global_automorphism_count
        == pair_triple.global_automorphism_count
    )
    label_equivariance_checks += 1


artifact = {
    "experiment": (
        "recursive root-residue screens and automatic nested resonance "
        "character extraction"
    ),
    "screen_leaf_range": [2, 5],
    "screen_family_counts": screen_family_counts,
    "screen_frame_equivariance_checks": (
        screen_frame_equivariance_checks
    ),
    "screen_contractions_checked": screen_contractions,
    "two_step_contractions_checked": two_step_contractions,
    "initial_form_degree_range": [1, 7],
    "initial_form_checks": initial_form_checks,
    "deformation_initial_form_checks": deformation_initial_form_checks,
    "flag_fiber_equation_checks": fiber_equation_checks,
    "source_frame_equivariance_checks": (
        source_frame_equivariance_checks
    ),
    "initial_form_transitivity_checks": transitivity_checks,
    "single_node_character_checks": single_node_character_checks,
    "general_cycle_profile_checks": general_cycle_profile_checks,
    "unrigidified_tame_chain_node_profiles": [
        list(profile)
        for profile in mixed_atlas.node_profiles
    ],
    "unrigidified_tame_chain_inertia_order": mixed_atlas.inertia_order,
    "general_radial_profiles": len(general_radial_profiles),
    "general_radial_charts": general_radial_charts,
    "general_radial_inertia_distribution": dict(
        sorted(general_radial_inertia_distribution.items())
    ),
    "unequal_radial_examples": unequal_radial_examples,
    "degree_six_radial_types": degree_six_radial_types,
    "degree_six_radial_inertia_order": 1,
    "degree_six_radial_global_automorphism_distribution": dict(
        sorted(degree_six_radial_global_automorphisms.items())
    ),
    "degree_six_nested_families": len(degree_six_families),
    "degree_six_inertia_distribution": dict(
        sorted(degree_six_inertia_distribution.items())
    ),
    "maximum_global_automorphisms": max_global_automorphisms,
    "group_law_atlases_checked": group_law_atlases_checked,
    "pair_triple_maxwell": {
        "node_profiles": [
            list(profile)
            for profile in pair_triple.node_profiles
        ],
        "global_automorphisms": pair_triple.global_automorphism_count,
        "inertia_order": pair_triple.inertia_order,
        "inertia_structure": "mu_2 x mu_2",
    },
    "sheet_label_equivariance_checks": label_equivariance_checks,
    "rules": {
        "screen_blowdown": (
            "u_parent=a_cluster+q_cluster*u_cluster"
        ),
        "screen_contraction": (
            "delete a screen by composing its two adjacent affine maps"
        ),
        "screen_frame_change": (
            "independent affine frame changes act by "
            "q_C -> d_parent*q_C/d_C and "
            "a_(parent,C) -> c_parent+d_parent*a_(parent,C)-q'_C*c_C"
        ),
        "flag_divisor": (
            "on a normalized source screen R(U), the inverse image of a "
            "target flag a is the effective divisor R(U)-a=0"
        ),
        "character_extraction": (
            "take full vertex centralizers, match their permutations of "
            "source-node cycles across every target node, add the two "
            "half-node rotations, then impose every diagonal Kummer phase"
        ),
        "radial_inertia": (
            "for ordered blocks B_j, the full-chain order is the product "
            "of L_j/M_j, where M_j is the lcm on B_j and L_j is the lcm "
            "on the suffix B_j union ... union B_k"
        ),
    },
    "scope": (
        "all nested tame branch charts after choosing a collision-adapted "
        "Hurwitz system; exhaustive bounded checks cover simple polynomial "
        "trees, every pair of cycle profiles through degree six, and a "
        "mixed unrigidified tame chain; bounded unequal-multiplicity radial "
        "profiles and all thirteen equal-multiplicity sextic radial scale "
        "types are checked separately"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print(
    "PASS framed residue screens: "
    f"{sum(map(int, screen_family_counts.values()))} nested families and "
    f"{screen_frame_equivariance_checks} gauge changes"
)
print(
    "PASS screen contractions: "
    f"{screen_contractions} one-step and "
    f"{two_step_contractions} two-step"
)
print(
    "PASS normalized flag equations and transitivity: "
    f"{initial_form_checks} ramification profiles and "
    f"{deformation_initial_form_checks} nonfactorized smoothings"
)
print(
    "PASS source/target frame descent: "
    f"{source_frame_equivariance_checks} map and flag transitions"
)
print(
    "PASS automatic single-node characters: "
    f"{single_node_character_checks} collision charts"
)
print(
    "PASS arbitrary tame cycle profiles: "
    f"{general_cycle_profile_checks} pairs"
)
print("PASS unrigidified tame chain: surviving order-six ghost inertia")
print(
    "PASS general radial inertia formula: "
    f"{general_radial_charts} unequal/equal multiplicity charts"
)
print("PASS all thirteen degree-six polynomial radial inertias: trivial")
print(
    "PASS degree-six nested characters: "
    f"{len(degree_six_families)} families"
)
print(
    "PASS pair-triple Maxwell inertia: "
    f"order {pair_triple.inertia_order}"
)
print(
    "PASS concrete inertia subgroups: "
    f"{group_law_atlases_checked} compiled atlases"
)
print(
    "PASS sheet-label equivariance: "
    f"{label_equivariance_checks} permutations"
)
print("RECURSIVE_RESONANCE_ATLAS_PASS")

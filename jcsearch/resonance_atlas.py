"""Recursive screen charts and characters for nested branch resonances.

The screen-coordinate half is independent of admissible covers.  A laminar
family of collision subsets gives a rooted hierarchy.  Every cluster embeds
its affine screen into its parent by

    u_parent = a_cluster + q_cluster * u_cluster.

Deleting a screen is therefore literal composition of affine maps.

The character half takes a collision-adapted polynomial monodromy tree.
For every target screen it computes the deck centralizer of each source
component and records its rotations on all incident source-node cycles.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations, product
from typing import Iterable, Mapping, TypeAlias

import sympy as sp

from jcsearch.admissible_nodes import (
    simultaneous_phase_is_diagonal,
)
from jcsearch.monodromy_forests import (
    Edge,
    Permutation,
    compose,
    cycles,
    forest_components,
    identity,
    permutation_product,
    transposition,
)
from jcsearch.wonderful_branch import compatible


Collision: TypeAlias = frozenset[int]
Parent: TypeAlias = Collision | None
Child: TypeAlias = Collision | int
NodeId: TypeAlias = tuple[Collision, tuple[int, ...]]


def collision_key(collision: Collision) -> tuple[int, tuple[int, ...]]:
    return (len(collision), tuple(sorted(collision)))


def child_key(child: Child) -> tuple[int, tuple[int, ...]]:
    if isinstance(child, int):
        return (0, (child,))
    return (1, tuple(sorted(child)))


def validate_nested_family(
    leaves: Iterable[int],
    family: Iterable[Collision],
) -> tuple[tuple[int, ...], tuple[Collision, ...]]:
    """Validate and canonically order a laminar collision family."""

    leaf_tuple = tuple(sorted(set(leaves)))
    leaf_set = frozenset(leaf_tuple)
    collisions = tuple(
        sorted(set(family), key=collision_key)
    )
    if not leaf_tuple:
        raise ValueError("a screen chart needs at least one leaf")
    if any(len(collision) < 2 for collision in collisions):
        raise ValueError("collision subsets must contain at least two leaves")
    if any(not collision <= leaf_set for collision in collisions):
        raise ValueError("a collision contains an unknown leaf")
    if any(
        not compatible(left, right)
        for index, left in enumerate(collisions)
        for right in collisions[index + 1 :]
    ):
        raise ValueError("collision subsets must be nested or disjoint")
    return leaf_tuple, collisions


def collision_parents(
    family: Iterable[Collision],
) -> dict[Collision, Parent]:
    """Return the unique minimal strict containing cluster, if it exists."""

    collisions = tuple(sorted(set(family), key=collision_key))
    result: dict[Collision, Parent] = {}
    for collision in collisions:
        containers = tuple(
            candidate
            for candidate in collisions
            if collision < candidate
        )
        result[collision] = (
            min(containers, key=collision_key)
            if containers
            else None
        )
    return result


def screen_children(
    leaves: Iterable[int],
    family: Iterable[Collision],
) -> dict[Parent, tuple[Child, ...]]:
    """Return cluster and singleton children on every framed screen."""

    leaf_tuple, collisions = validate_nested_family(leaves, family)
    parents = collision_parents(collisions)
    result: dict[Parent, tuple[Child, ...]] = {}
    for parent in (None, *collisions):
        scope = set(leaf_tuple if parent is None else parent)
        cluster_children = tuple(
            collision
            for collision in collisions
            if parents[collision] == parent
        )
        covered = set().union(*cluster_children) if cluster_children else set()
        singleton_children = tuple(sorted(scope - covered))
        result[parent] = tuple(
            sorted(
                (*cluster_children, *singleton_children),
                key=child_key,
            )
        )
    return result


def _cluster_name(collision: Collision) -> str:
    return "_".join(map(str, sorted(collision)))


def _parent_name(parent: Parent) -> str:
    return "root" if parent is None else f"s_{_cluster_name(parent)}"


def _child_name(child: Child) -> str:
    if isinstance(child, int):
        return f"l_{child}"
    return f"s_{_cluster_name(child)}"


@dataclass
class FramedScreenChart:
    """A redundant affine-framed chart for one nested-set stratum.

    The two-dimensional affine gauge on each non-root screen may be removed
    by fixing any two child placements to zero and one.  Keeping the frames
    makes every contraction polynomial and choice-free.
    """

    leaves: tuple[int, ...]
    family: tuple[Collision, ...]
    children: dict[Parent, tuple[Child, ...]]
    placements: dict[tuple[Parent, Child], sp.Symbol]
    scales: dict[Collision, sp.Symbol]
    leaf_expressions: dict[int, sp.Expr]
    prefix: str


def framed_screen_chart(
    leaves: Iterable[int],
    family: Iterable[Collision],
    *,
    prefix: str = "z",
) -> FramedScreenChart:
    """Construct recursive affine screen coordinates and their blowdown."""

    leaf_tuple, collisions = validate_nested_family(leaves, family)
    children = screen_children(leaf_tuple, collisions)
    placements = {
        (parent, child): sp.Symbol(
            f"{prefix}_a_{_parent_name(parent)}_{_child_name(child)}"
        )
        for parent, child_list in children.items()
        for child in child_list
    }
    scales = {
        collision: sp.Symbol(
            f"{prefix}_q_{_cluster_name(collision)}"
        )
        for collision in collisions
    }

    def expression_on(parent: Parent, leaf: int) -> sp.Expr:
        containing_child = next(
            child
            for child in children[parent]
            if (
                child == leaf
                if isinstance(child, int)
                else leaf in child
            )
        )
        placement = placements[(parent, containing_child)]
        if isinstance(containing_child, int):
            return placement
        return sp.expand(
            placement
            + scales[containing_child]
            * expression_on(containing_child, leaf)
        )

    leaf_expressions = {
        leaf: expression_on(None, leaf)
        for leaf in leaf_tuple
    }
    return FramedScreenChart(
        leaves=leaf_tuple,
        family=collisions,
        children=children,
        placements=placements,
        scales=scales,
        leaf_expressions=leaf_expressions,
        prefix=prefix,
    )


def screen_frame_change_substitution(
    chart: FramedScreenChart,
    translations: Mapping[Parent, sp.Expr],
    dilations: Mapping[Parent, sp.Expr],
) -> dict[sp.Symbol, sp.Expr]:
    """Return the parameter change for independent affine screen frames.

    The coordinate on screen ``S`` changes by

        u'_S = translation[S] + dilation[S] * u_S.

    For a cluster child ``C`` of ``S`` this sends

        q_C -> dilation[S] * q_C / dilation[C],
        a_(S,C) -> translation[S] + dilation[S] * a_(S,C)
                   - q'_C * translation[C].

    A singleton placement has no child-frame correction.
    """

    screens = (None, *chart.family)
    screen_set = frozenset(screens)
    if (
        not set(translations) <= screen_set
        or not set(dilations) <= screen_set
    ):
        raise ValueError("a frame change names an unknown screen")
    translation = {
        screen: sp.sympify(translations.get(screen, 0))
        for screen in screens
    }
    dilation = {
        screen: sp.sympify(dilations.get(screen, 1))
        for screen in screens
    }
    if any(value.is_zero is True for value in dilation.values()):
        raise ValueError("screen dilations must be units")

    transformed_scales = {
        cluster: sp.cancel(
            dilation[parent]
            * chart.scales[cluster]
            / dilation[cluster]
        )
        for cluster, parent in collision_parents(chart.family).items()
    }
    substitution: dict[sp.Symbol, sp.Expr] = {
        chart.scales[cluster]: transformed_scale
        for cluster, transformed_scale in transformed_scales.items()
    }
    for (parent, child), placement in chart.placements.items():
        transformed = (
            translation[parent]
            + dilation[parent] * placement
        )
        if isinstance(child, frozenset):
            transformed -= (
                transformed_scales[child] * translation[child]
            )
        substitution[placement] = sp.cancel(transformed)
    return substitution


def screen_contraction_substitution(
    chart: FramedScreenChart,
    collision: Collision,
) -> tuple[FramedScreenChart, dict[sp.Symbol, sp.Expr]]:
    """Delete one screen and return the exact affine-composition map."""

    if collision not in chart.family:
        raise ValueError("the requested screen is not in the chart")
    parents = collision_parents(chart.family)
    parent = parents[collision]
    coarser = framed_screen_chart(
        chart.leaves,
        (
            candidate
            for candidate in chart.family
            if candidate != collision
        ),
        prefix=chart.prefix,
    )
    substitution: dict[sp.Symbol, sp.Expr] = {}

    for key, coarse_symbol in coarser.placements.items():
        if key in chart.placements:
            substitution[coarse_symbol] = chart.placements[key]
            continue
        coarse_parent, promoted_child = key
        if coarse_parent != parent:
            raise AssertionError("unexpected promoted screen child")
        substitution[coarse_symbol] = sp.expand(
            chart.placements[(parent, collision)]
            + chart.scales[collision]
            * chart.placements[(collision, promoted_child)]
        )

    old_parents = collision_parents(chart.family)
    for cluster, coarse_scale in coarser.scales.items():
        if old_parents[cluster] == collision:
            substitution[coarse_scale] = sp.expand(
                chart.scales[collision] * chart.scales[cluster]
            )
        else:
            substitution[coarse_scale] = chart.scales[cluster]

    return coarser, substitution


def normalized_polynomial_pullback(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    *,
    source_center: sp.Expr,
    source_scale: sp.Expr,
    source_coordinate: sp.Symbol,
    target_center: sp.Expr,
    target_scale: sp.Expr,
) -> sp.Expr:
    """Return the normalized map on one source/target screen.

    The caller works on the chart where the displayed quotient is regular.
    In the polynomial graph this regularity is exactly the initial-form
    divisibility supplied by the source and target scale construction.
    """

    numerator = sp.expand(
        polynomial.subs(
            variable,
            source_center + source_scale * source_coordinate,
        )
        - target_center
    )
    return sp.cancel(numerator / target_scale)


def ramified_initial_form(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    point: sp.Expr,
    ramification_index: int,
    scale: sp.Symbol,
    local_coordinate: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return a regular weighted Taylor family and its special fiber.

    The polynomial may itself depend on ``scale``.  Lower Taylor
    coefficients need not vanish identically: it is enough that, after
    substituting ``W=point+scale*U``, every coefficient is divisible by
    ``scale**ramification_index`` and that the special fiber has order
    exactly ``ramification_index`` at ``U=0``.
    """

    if ramification_index <= 0:
        raise ValueError("ramification index must be positive")
    target_center = sp.expand(polynomial.subs(variable, point))
    family = normalized_polynomial_pullback(
        polynomial,
        variable,
        source_center=point,
        source_scale=scale,
        source_coordinate=local_coordinate,
        target_center=target_center,
        target_scale=scale**ramification_index,
    )
    denominator = sp.denom(family)
    if sp.expand(denominator) != 1:
        raise ValueError(
            "the Taylor coefficients fail the asserted weighted divisibility"
        )
    initial = sp.expand(family.subs(scale, 0))
    if initial == 0:
        raise ValueError(
            "the chosen ramification index is smaller than the true order"
        )
    initial_polynomial = sp.Poly(initial, local_coordinate)
    initial_order = min(
        monomial[0]
        for monomial, coefficient in initial_polynomial.terms()
        if coefficient != 0
    )
    if initial_order != ramification_index:
        raise ValueError(
            "the special fiber does not have the asserted ramification order"
        )
    return sp.expand(family), initial


def affine_screen_composition(
    outer_center: sp.Expr,
    outer_scale: sp.Expr,
    inner_center: sp.Expr,
    inner_scale: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Compose u=c+rU and U=b+sV."""

    return (
        sp.expand(outer_center + outer_scale * inner_center),
        sp.expand(outer_scale * inner_scale),
    )


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def subset_monodromy(
    degree: int,
    edge_sequence: tuple[Edge, ...],
    positions: Collision,
) -> Permutation:
    branch_permutations = tuple(
        transposition(degree, *edge)
        for edge in edge_sequence
    )
    return subset_factor_monodromy(
        degree,
        branch_permutations,
        positions,
    )


@lru_cache(maxsize=None)
def subset_factor_monodromy(
    degree: int,
    branch_permutations: tuple[Permutation, ...],
    positions: Collision,
) -> Permutation:
    factors = tuple(
        branch_permutations[position]
        for position in sorted(positions)
    )
    return permutation_product(factors, degree)


def is_interval_collision(collision: Collision) -> bool:
    ordered = tuple(sorted(collision))
    return ordered == tuple(range(ordered[0], ordered[-1] + 1))


def is_collision_adapted(family: Iterable[Collision]) -> bool:
    """Return whether the current Hurwitz order makes every cluster an interval."""

    return all(is_interval_collision(collision) for collision in family)


@dataclass
class ResonanceCharacterAtlas:
    degree: int
    family: tuple[Collision, ...]
    node_ids: tuple[tuple[NodeId, ...], ...]
    node_profiles: tuple[tuple[int, ...], ...]
    vertex_automorphism_counts: tuple[tuple[str, int], ...]
    global_automorphism_count: int
    global_automorphisms: tuple[tuple[Permutation, ...], ...]
    character_values: tuple[tuple[int, ...], ...]
    inertia_order: int
    inertia_automorphisms: tuple[tuple[Permutation, ...], ...]


def _node_data(
    degree: int,
    branch_permutations: tuple[Permutation, ...],
    family: tuple[Collision, ...],
) -> tuple[
    tuple[tuple[NodeId, ...], ...],
    tuple[tuple[int, ...], ...],
    dict[NodeId, int],
]:
    node_ids = []
    node_profiles = []
    coordinate_by_id: dict[NodeId, int] = {}
    offset = 0
    for collision in family:
        boundary = subset_factor_monodromy(
            degree,
            branch_permutations,
            collision,
        )
        components = tuple(
            tuple(sorted(cycle))
            for cycle in cycles(boundary)
        )
        ids = tuple((collision, component) for component in components)
        node_ids.append(ids)
        profile = tuple(len(component) for component in components)
        node_profiles.append(profile)
        for index, node_id in enumerate(ids):
            coordinate_by_id[node_id] = offset + index
        offset += len(ids)
    return tuple(node_ids), tuple(node_profiles), coordinate_by_id


def _cycle_action(
    deck: Permutation,
    boundary_monodromy: Permutation,
    source_cycle: tuple[int, ...],
) -> tuple[tuple[int, ...], int]:
    """Return the image cycle and its oriented rotation exponent."""

    boundary_cycles = {
        frozenset(cycle): cycle
        for cycle in cycles(boundary_monodromy)
    }
    image_set = frozenset(deck[index] for index in source_cycle)
    target_cycle = boundary_cycles.get(image_set)
    if target_cycle is None:
        raise ValueError("the deck map does not preserve boundary cycles")
    image_of_base = deck[source_cycle[0]]
    exponent = target_cycle.index(image_of_base)
    if any(
        deck[source_cycle[index]]
        != target_cycle[(index + exponent) % len(source_cycle)]
        for index in range(len(source_cycle))
    ):
        raise ValueError("the deck map reverses a boundary-cycle orientation")
    return tuple(sorted(target_cycle)), exponent


@lru_cache(maxsize=None)
def _vertex_options(
    degree: int,
    branch_permutations: tuple[Permutation, ...],
    collision: Collision,
    child_clusters: tuple[Collision, ...],
    *,
    anchored: bool,
) -> tuple[
    tuple[
        Permutation,
        dict[NodeId, tuple[NodeId, int]],
    ],
    ...,
]:
    """Return deck permutations and half-node actions on one target screen."""

    covered_positions = (
        frozenset().union(*child_clusters)
        if child_clusters
        else frozenset()
    )
    direct_positions = tuple(sorted(collision - covered_positions))
    local_generators = tuple(
        branch_permutations[position]
        for position in direct_positions
    ) + tuple(
        subset_factor_monodromy(degree, branch_permutations, child)
        for child in child_clusters
    )
    candidates: Iterable[Permutation]
    if anchored:
        candidates = (identity(degree),)
    else:
        candidates = permutations(range(degree))

    incident_boundaries = [
        (
            collision,
            inverse(
                subset_factor_monodromy(
                    degree,
                    branch_permutations,
                    collision,
                )
            ),
        )
    ] + [
        (
            child,
            subset_factor_monodromy(
                degree,
                branch_permutations,
                child,
            ),
        )
        for child in child_clusters
    ]
    result = []
    for candidate in candidates:
        deck = tuple(candidate)
        if not all(
            compose(deck, generator) == compose(generator, deck)
            for generator in local_generators
        ):
            continue
        action: dict[NodeId, tuple[NodeId, int]] = {}
        for target_node, boundary in incident_boundaries:
            for cycle in cycles(boundary):
                node_id = (target_node, tuple(sorted(cycle)))
                image_cycle, exponent = _cycle_action(
                    deck,
                    boundary,
                    cycle,
                )
                action[node_id] = (
                    (target_node, image_cycle),
                    exponent,
                )
        result.append((deck, action))
    if not result:
        raise AssertionError("every target vertex has the identity deck map")
    return tuple(result)


@lru_cache(maxsize=None)
def _root_vertex_options(
    degree: int,
    branch_permutations: tuple[Permutation, ...],
    maximal_clusters: tuple[Collision, ...],
    *,
    anchored: bool,
) -> tuple[
    tuple[
        Permutation,
        dict[NodeId, tuple[NodeId, int]],
    ],
    ...,
]:
    """Return deck permutations and node actions on the root target screen."""

    covered_positions = (
        frozenset().union(*maximal_clusters)
        if maximal_clusters
        else frozenset()
    )
    direct_positions = tuple(
        sorted(
            frozenset(range(len(branch_permutations)))
            - covered_positions
        )
    )
    local_generators = tuple(
        branch_permutations[position]
        for position in direct_positions
    ) + tuple(
        subset_factor_monodromy(
            degree,
            branch_permutations,
            child,
        )
        for child in maximal_clusters
    )
    candidates: Iterable[Permutation]
    if anchored:
        candidates = (identity(degree),)
    else:
        candidates = permutations(range(degree))

    result = []
    for candidate in candidates:
        deck = tuple(candidate)
        if not all(
            compose(deck, generator) == compose(generator, deck)
            for generator in local_generators
        ):
            continue
        action: dict[NodeId, tuple[NodeId, int]] = {}
        for child in maximal_clusters:
            boundary = subset_factor_monodromy(
                degree,
                branch_permutations,
                child,
            )
            for cycle in cycles(boundary):
                node_id = (child, tuple(sorted(cycle)))
                image_cycle, exponent = _cycle_action(
                    deck,
                    boundary,
                    cycle,
                )
                action[node_id] = (
                    (child, image_cycle),
                    exponent,
                )
        result.append((deck, action))
    if not result:
        raise AssertionError("the root target vertex lost its identity map")
    return tuple(result)


def compile_resonance_characters(
    degree: int,
    edge_sequence: tuple[Edge, ...],
    family: Iterable[Collision],
    *,
    anchored_vertices: Iterable[Collision] = (),
) -> ResonanceCharacterAtlas:
    """Extract characters from a simple polynomial monodromy tree."""

    branch_permutations = tuple(
        transposition(degree, *edge)
        for edge in edge_sequence
    )
    return compile_nested_monodromy_characters(
        degree,
        branch_permutations,
        family,
        anchored_vertices=anchored_vertices,
        root_anchored=True,
    )


def compile_nested_monodromy_characters(
    degree: int,
    branch_permutations: tuple[Permutation, ...],
    family: Iterable[Collision],
    *,
    anchored_vertices: Iterable[Collision] = (),
    root_anchored: bool = True,
) -> ResonanceCharacterAtlas:
    """Extract all local node characters on a nested branch chart.

    The branch tuple may have arbitrary tame cycle types.  It must be written
    in a collision-adapted Hurwitz system: every collision subset is an
    interval in the tuple.  Any laminar family admits such a local branch-path
    system.  ``compile_resonance_characters`` is the transposition-tree
    wrapper used on the generic polynomial simple-branch locus.  Set
    ``root_anchored=False`` precisely when the labelled zero flag lies on a
    non-root collision screen named in ``anchored_vertices``.
    """

    positions, collisions = validate_nested_family(
        range(len(branch_permutations)),
        family,
    )
    if positions != tuple(range(len(branch_permutations))):
        raise AssertionError("edge positions must be consecutive")
    if any(
        len(permutation) != degree
        for permutation in branch_permutations
    ):
        raise ValueError("branch permutations have the wrong degree")
    if not is_collision_adapted(collisions):
        raise ValueError(
            "the transposition tuple is not collision-adapted"
        )
    anchored = frozenset(anchored_vertices)
    if not anchored <= frozenset(collisions):
        raise ValueError("an anchored vertex is not in the collision family")

    parents = collision_parents(collisions)
    children = screen_children(
        range(len(branch_permutations)),
        collisions,
    )
    cluster_children = {
        collision: tuple(
            child
            for child in children[collision]
            if isinstance(child, frozenset)
        )
        for collision in collisions
    }
    node_ids, node_profiles, coordinate_by_id = _node_data(
        degree,
        branch_permutations,
        collisions,
    )
    moduli = tuple(
        modulus
        for profile in node_profiles
        for modulus in profile
    )
    width = len(moduli)

    vertex_options = {}
    vertex_counts = []
    for collision in collisions:
        child_clusters = cluster_children[collision]
        options = _vertex_options(
            degree,
            branch_permutations,
            collision,
            child_clusters,
            anchored=collision in anchored,
        )
        vertex_options[collision] = options
        vertex_counts.append(
            (
                ",".join(map(str, sorted(collision))),
                len(options),
            )
        )

    maximal_clusters = tuple(
        collision
        for collision in collisions
        if parents[collision] is None
    )
    root_options = _root_vertex_options(
        degree,
        branch_permutations,
        maximal_clusters,
        anchored=root_anchored,
    )
    vertex_counts.insert(0, ("root", len(root_options)))

    global_automorphisms = []
    character_rows = []
    option_lists = tuple(vertex_options[collision] for collision in collisions)
    for root_choice, *choices in product(root_options, *option_lists):
        actions_by_vertex = {
            collision: choice[1]
            for collision, choice in zip(collisions, choices)
        }
        row = [0] * width
        compatible_actions = True
        for collision, ids in zip(collisions, node_ids):
            child_action = actions_by_vertex[collision]
            parent = parents[collision]
            parent_action = (
                root_choice[1]
                if parent is None
                else actions_by_vertex[parent]
            )
            for node_id in ids:
                child_image, child_phase = child_action[node_id]
                parent_image, parent_phase = parent_action[node_id]
                if child_image != parent_image:
                    compatible_actions = False
                    break
                coordinate = coordinate_by_id[node_id]
                row[coordinate] = (
                    child_phase + parent_phase
                ) % moduli[coordinate]
            if not compatible_actions:
                break
        if compatible_actions:
            global_automorphisms.append(
                (
                    root_choice[0],
                    *(choice[0] for choice in choices),
                )
            )
            character_rows.append(tuple(row))

    global_automorphisms = tuple(global_automorphisms)
    character_values = tuple(character_rows)
    if not character_values:
        raise AssertionError("the global identity automorphism was lost")
    inertia_automorphisms = tuple(
        automorphism
        for automorphism, value in zip(
            global_automorphisms,
            character_values,
        )
        if simultaneous_phase_is_diagonal(node_profiles, value)
    )
    inertia_order = len(inertia_automorphisms)
    return ResonanceCharacterAtlas(
        degree=degree,
        family=collisions,
        node_ids=node_ids,
        node_profiles=node_profiles,
        vertex_automorphism_counts=tuple(vertex_counts),
        global_automorphism_count=len(character_values),
        global_automorphisms=global_automorphisms,
        character_values=character_values,
        inertia_order=inertia_order,
        inertia_automorphisms=inertia_automorphisms,
    )

"""Bounded metric-assisted search for exact relation-complex embeddings.

Height angles are used only to order and prune proposals.  A returned map is
accepted solely after exact ternary replay and an integral lattice lift in
``relations.py``.  This separation is important for specializations, where
canonical-height forms are close only up to scale and specialization error.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations
from math import log, sqrt
from typing import Sequence

import numpy as np

from .codes import CandidateFiniteSignature, candidate_finite_signature
from .finite import FiniteQuotientBlock
from .integer import (
    canonical_unoriented,
    modular_rank,
    rational_nullspace,
    rational_rank,
    row_basis_coordinates,
)
from .local import ComponentBlock
from .pari import (
    primitive_column_closure,
    row_embedding_smith_invariant_factors,
)
from .relations import (
    RelationComplex,
    build_relation_complex,
    lift_relation_vertex_injection,
)


@dataclass(frozen=True)
class MetricRelationEmbedding:
    """One exact survivor of a bounded metric-assisted relation search."""

    vertex_map: tuple[int, ...]
    source_vertex_indices: tuple[int, ...]
    integral_matrices: tuple[tuple[tuple[int, ...], ...], ...]
    matched_seed_relations: int
    skipped_source_relations: int
    metric_rms_residual: float
    global_replay_ray_count: int
    global_replay_source_rank: int
    global_replay_relation_count: int
    metric_reseed_count: int


@dataclass(frozen=True)
class MetricRelationSearchLedger:
    """Auditable counts and exact survivors from one declared search box."""

    source_center_count: int
    target_center_count: int
    initial_state_count: int
    maximum_beam_population: int
    expanded_state_count: int
    exact_lift_attempt_count: int
    maximum_source_rank_reached: int
    maximum_mapped_vertex_count: int
    maximum_matched_relation_count: int
    maximum_global_replay_ray_count: int
    exact_lifts_below_global_support: int
    maximum_metric_reseed_count: int
    partial_replay_attempt_count: int
    maximum_partial_replay_ray_count: int
    maximum_partial_replay_relation_count: int
    embeddings: tuple[MetricRelationEmbedding, ...]


@dataclass(frozen=True)
class PartialRelationReplay:
    """Exact replay certificate for a mapped proper rational subspace.

    ``integral_matrix`` acts on coordinates in ``primitive_source_basis``.
    Thus a certificate remains meaningful when the mapped component has rank
    below the ambient source rank: no arbitrary completion to a square basis
    is introduced.
    """

    source_rank: int
    primitive_source_basis: tuple[tuple[int, ...], ...]
    mapped_source_vertex_indices: tuple[int, ...]
    intrinsic_relation_count: int
    integral_matrix: tuple[tuple[int, ...], ...]
    target_smith_invariant_factors: tuple[int, ...]
    primitive_target_image: bool
    source_subspace_ray_count: int
    replayed_source_vertex_indices: tuple[int, ...]
    replayed_source_rank: int
    replayed_relation_count: int


@dataclass(frozen=True)
class StarComponentCandidate:
    """One exact partial-replay survivor of a bounded star assignment."""

    vertex_map: tuple[int, ...]
    matched_star_relations: int
    skipped_star_relations: int
    metric_rms_residual: float
    replay: PartialRelationReplay


@dataclass(frozen=True)
class StarComponentSearchLedger:
    """Auditable box and survivors for center-star component matching."""

    source_center_count: int
    target_center_count: int
    center_pair_count: int
    processed_star_layer_count: int
    maximum_beam_population: int
    expanded_state_count: int
    finite_rank_rejection_count: int
    partial_replay_attempt_count: int
    maximum_source_rank_reached: int
    maximum_mapped_vertex_count: int
    maximum_partial_replay_ray_count: int
    maximum_partial_replay_relation_count: int
    candidates: tuple[StarComponentCandidate, ...]


@dataclass(frozen=True)
class _State:
    stratum: int
    mapping: tuple[int, ...]
    processed_edges: frozenset[int]
    matched_edges: int
    skipped_edges: int
    squared_residual: float
    residual_count: int
    scale: float
    stagnation: int
    reseeds: int
    partial_replay_rays: int = 0
    partial_replay_relations: int = 0
    partial_replay_rank: int = 0
    partial_audited_rank: int = 0
    partial_audited_mapped_count: int = 0

    @property
    def rms(self) -> float:
        return sqrt(self.squared_residual / max(1, self.residual_count))


def _quadratic_data(vectors, gram):
    matrix = np.asarray(vectors, dtype=float)
    form = np.asarray(gram, dtype=float)
    if form.shape != (matrix.shape[1], matrix.shape[1]):
        raise ValueError("Gram shape differs from relation-vector width")
    norms = np.einsum("ij,jk,ik->i", matrix, form, matrix)
    if np.any(norms <= 0):
        raise ValueError("metric-assisted search requires positive norms")
    pairings = matrix @ form @ matrix.T
    angles = np.abs(pairings / np.sqrt(norms[:, None] * norms[None, :]))
    return norms, angles


def _incidence(complex_: RelationComplex):
    by_vertex = [[] for _ in complex_.vertices]
    by_pair: dict[tuple[int, int], list[int]] = {}
    for edge_index, edge in enumerate(complex_.ternary_relations):
        for vertex in edge:
            by_vertex[vertex].append(edge_index)
        for left in range(3):
            for right in range(left):
                pair = tuple(sorted((edge[left], edge[right])))
                by_pair.setdefault(pair, []).append(edge_index)
    return tuple(tuple(items) for items in by_vertex), {
        pair: tuple(items) for pair, items in by_pair.items()
    }


def _star_features(complex_, by_vertex, norms, angles, center):
    answer = []
    for edge_index in by_vertex[center]:
        others = [
            vertex
            for vertex in complex_.ternary_relations[edge_index]
            if vertex != center
        ]
        pairs = sorted(
            (log(norms[vertex] / norms[center]), angles[vertex, center])
            for vertex in others
        )
        answer.append(tuple(value for pair in pairs for value in pair))
    return np.asarray(answer, dtype=float)


def _star_profile_distance(source_features, target_features, retained_edges=8):
    if not len(source_features) or not len(target_features):
        return float("inf")
    distances = np.sqrt(
        np.sum(
            (source_features[:, None, :] - target_features[None, :, :]) ** 2,
            axis=2,
        )
    )
    nearest = np.sort(np.min(distances, axis=1))
    return float(np.mean(nearest[: min(retained_edges, len(nearest))]))


def _metric_extension(
    state: _State,
    additions: Sequence[tuple[int, int]],
    source_norms,
    target_norms,
    source_angles,
    target_angles,
    *,
    norm_log_tolerance: float,
    angle_tolerance: float,
    angle_hard_tolerance: float,
) -> tuple[tuple[int, ...], float, int] | None:
    mapping = list(state.mapping)
    used = {target for target in mapping if target >= 0}
    if len({source for source, _target in additions}) != len(additions):
        return None
    if any(target in used for _source, target in additions):
        return None
    if len({target for _source, target in additions}) != len(additions):
        return None
    prior = [(source, target) for source, target in enumerate(mapping) if target >= 0]
    squared = 0.0
    count = 0
    accepted: list[tuple[int, int]] = []
    for source, target in additions:
        norm_error = abs(log((target_norms[target] / state.scale) / source_norms[source]))
        if norm_error > norm_log_tolerance:
            return None
        squared += norm_error * norm_error
        count += 1
        angle_errors = []
        for old_source, old_target in (*prior, *accepted):
            error = abs(
                target_angles[target, old_target]
                - source_angles[source, old_source]
            )
            if error > angle_hard_tolerance:
                return None
            angle_errors.append(error)
            squared += error * error
            count += 1
        if angle_errors and sqrt(
            sum(error * error for error in angle_errors) / len(angle_errors)
        ) > angle_tolerance:
            return None
        mapping[source] = target
        accepted.append((source, target))
    return tuple(mapping), squared, count


def _state_key(state: _State) -> tuple[float, ...]:
    mapped = sum(target >= 0 for target in state.mapping)
    # Rank is recomputed by the caller and placed first there.  This secondary
    # key rewards exact-relation support while retaining metric discrimination.
    return (
        float(state.partial_replay_rays),
        float(state.partial_replay_relations),
        float(state.matched_edges),
        float(mapped),
        -state.rms,
        -float(state.skipped_edges),
        -float(state.reseeds),
    )


def _diverse_truncate(states, width, key, minimum_per_stratum):
    """Retain bounded center-pair diversity before global score filling."""

    if len(states) <= width:
        return sorted(states, key=key, reverse=True)
    grouped: dict[int, list[_State]] = {}
    for state in states:
        grouped.setdefault(state.stratum, []).append(state)
    reserved = []
    for stratum in sorted(grouped):
        reserved.extend(
            sorted(grouped[stratum], key=key, reverse=True)[:minimum_per_stratum]
        )
    reserved = sorted(reserved, key=key, reverse=True)[:width]
    selected = {id(state) for state in reserved}
    if len(reserved) < width:
        for state in sorted(states, key=key, reverse=True):
            if id(state) in selected:
                continue
            reserved.append(state)
            if len(reserved) == width:
                break
    return reserved


def _grouped_truncate(states, width, key, group_key, minimum_per_group):
    """Retain a declared minimum from each caller-defined state group."""

    if len(states) <= width:
        return sorted(states, key=key, reverse=True)
    grouped = {}
    for state in states:
        grouped.setdefault(group_key(state), []).append(state)
    reserved = []
    for group in sorted(grouped, key=repr):
        reserved.extend(
            sorted(grouped[group], key=key, reverse=True)[:minimum_per_group]
        )
    reserved = sorted(reserved, key=key, reverse=True)[:width]
    selected = {id(state) for state in reserved}
    if len(reserved) < width:
        for state in sorted(states, key=key, reverse=True):
            if id(state) in selected:
                continue
            reserved.append(state)
            if len(reserved) == width:
                break
    return reserved


def _global_replay(source, target, matrix):
    target_index = {vertex: index for index, vertex in enumerate(target.vertices)}
    mapped = []
    source_indices = []
    width = len(matrix[0])
    for source_index, vector in enumerate(source.vertices):
        image = canonical_unoriented(
            tuple(
                sum(int(vector[row]) * int(matrix[row][column]) for row in range(len(matrix)))
                for column in range(width)
            )
        )
        target_vertex = target_index.get(image)
        mapped.append(target_vertex)
        if target_vertex is not None:
            source_indices.append(source_index)
    target_edges = set(target.ternary_relations)
    relation_count = sum(
        all(mapped[index] is not None for index in edge)
        and tuple(sorted(mapped[index] for index in edge)) in target_edges
        for edge in source.ternary_relations
    )
    return (
        len(source_indices),
        rational_rank(tuple(source.vertices[index] for index in source_indices)),
        int(relation_count),
    )


def _primitive_row_span_basis(
    rows: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    """Return a primitive integral basis of the rational row span."""

    independent = []
    for row in rows:
        if rational_rank((*independent, tuple(row))) > len(independent):
            independent.append(tuple(map(int, row)))
    if not independent:
        raise ValueError("partial replay needs at least one mapped source ray")
    columns = tuple(tuple(column) for column in zip(*independent))
    closure_columns = primitive_column_closure(columns)
    basis = tuple(tuple(column) for column in zip(*closure_columns))
    if rational_rank(basis) != len(independent):
        raise ArithmeticError("primitive closure changed partial source rank")
    return basis


def exact_partial_relation_replay(
    source: RelationComplex,
    target: RelationComplex,
    vertex_map: Sequence[int],
    *,
    maximum_sign_components: int = 16,
) -> tuple[PartialRelationReplay, ...]:
    """Lift and globally replay an exact injection on its primitive subspace.

    ``vertex_map`` has one entry per source vertex; negative entries are
    unmapped.  Mapped rays first define a rational subspace.  That subspace is
    saturated inside the source coordinate lattice, all mapped rays are
    rewritten in the resulting intrinsic integer coordinates, and the
    relation injection is lifted exactly.  Every supplied source ray in the
    same rational subspace is then replayed against the target cloud.

    The finite clouds and the sign-component bound are explicit boundaries.
    An empty return means that no integral lift survived within those bounds;
    it is not a statement about rays absent from either cloud.
    """

    mapping = tuple(map(int, vertex_map))
    if len(mapping) != len(source.vertices):
        raise ValueError("partial vertex_map length differs from source complex")
    mapped_indices = tuple(index for index, value in enumerate(mapping) if value >= 0)
    mapped_targets = tuple(mapping[index] for index in mapped_indices)
    if not mapped_indices:
        raise ValueError("partial vertex_map has no mapped source rays")
    if len(set(mapped_targets)) != len(mapped_targets):
        raise ValueError("partial vertex_map is not injective")
    if any(index >= len(target.vertices) for index in mapped_targets):
        raise ValueError("partial vertex_map target is outside the relation complex")

    primitive_basis = _primitive_row_span_basis(
        tuple(source.vertices[index] for index in mapped_indices)
    )
    rank = len(primitive_basis)
    normals = rational_nullspace(primitive_basis)
    subspace_indices = tuple(
        index
        for index, vector in enumerate(source.vertices)
        if all(
            sum(int(value) * int(normal[column]) for column, value in enumerate(vector))
            == 0
            for normal in normals
        )
    )
    subspace_coordinates = row_basis_coordinates(
        tuple(source.vertices[index] for index in subspace_indices),
        primitive_basis,
    )
    coordinates_by_source = dict(zip(subspace_indices, subspace_coordinates))
    mapped_coordinates = tuple(coordinates_by_source[index] for index in mapped_indices)
    intrinsic = build_relation_complex(mapped_coordinates)
    target_by_coordinate = {}
    source_by_coordinate = {}
    for source_index in mapped_indices:
        coordinate = canonical_unoriented(coordinates_by_source[source_index])
        previous = target_by_coordinate.get(coordinate)
        if previous is not None and previous != mapping[source_index]:
            raise ValueError("distinct mapped targets represent one intrinsic ray")
        target_by_coordinate[coordinate] = mapping[source_index]
        source_by_coordinate[coordinate] = source_index
    intrinsic_map = tuple(target_by_coordinate[vertex] for vertex in intrinsic.vertices)
    try:
        lifts = lift_relation_vertex_injection(
            intrinsic,
            target,
            intrinsic_map,
            maximum_sign_components=maximum_sign_components,
        )
    except ValueError:
        return ()

    target_index = {vertex: index for index, vertex in enumerate(target.vertices)}
    target_edges = set(target.ternary_relations)
    answers = []
    for matrix in lifts:
        smith = row_embedding_smith_invariant_factors(matrix)
        replay_map = [None] * len(source.vertices)
        replayed = []
        for source_index in subspace_indices:
            coordinate = coordinates_by_source[source_index]
            image = canonical_unoriented(
                tuple(
                    sum(
                        int(coordinate[row]) * int(matrix[row][column])
                        for row in range(rank)
                    )
                    for column in range(len(matrix[0]))
                )
            )
            target_vertex = target_index.get(image)
            replay_map[source_index] = target_vertex
            if target_vertex is not None:
                replayed.append(source_index)
        relation_count = sum(
            all(replay_map[index] is not None for index in edge)
            and tuple(sorted(replay_map[index] for index in edge)) in target_edges
            for edge in source.ternary_relations
        )
        answers.append(
            PartialRelationReplay(
                source_rank=rank,
                primitive_source_basis=primitive_basis,
                mapped_source_vertex_indices=tuple(
                    source_by_coordinate[vertex] for vertex in intrinsic.vertices
                ),
                intrinsic_relation_count=len(intrinsic.ternary_relations),
                integral_matrix=matrix,
                target_smith_invariant_factors=smith,
                primitive_target_image=all(value == 1 for value in smith),
                source_subspace_ray_count=len(subspace_indices),
                replayed_source_vertex_indices=tuple(replayed),
                replayed_source_rank=rational_rank(
                    tuple(source.vertices[index] for index in replayed)
                ),
                replayed_relation_count=int(relation_count),
            )
        )
    return tuple(answers)


def partial_replay_finite_signature(
    replay: PartialRelationReplay,
    target: RelationComplex,
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> CandidateFiniteSignature:
    """Attach source-free finite-code invariants to an exact partial replay.

    The signature is computed on the replay's integral target image and is
    invariant under integral rebasing of that image and quotient-code basis
    choices.  Comparing such signatures across fibres is a calibrated
    heuristic score; the underlying code restrictions are exact.
    """

    return candidate_finite_signature(
        replay.integral_matrix,
        target,
        finite_blocks=finite_blocks,
        component_blocks=component_blocks,
    )


def bounded_metric_star_component_search(
    source: RelationComplex,
    target: RelationComplex,
    source_gram: Sequence[Sequence[float]],
    target_gram: Sequence[Sequence[float]],
    *,
    source_center_indices: Sequence[int],
    target_center_indices: Sequence[int],
    source_edge_limit: int = 24,
    target_edges_per_source: int = 64,
    beam_width: int = 512,
    minimum_states_per_match_rank_group: int = 4,
    maximum_partial_replay_attempts: int = 500,
    minimum_partial_replay_rank: int = 10,
    partial_replay_rank_stride: int = 2,
    partial_replay_vertex_stride: int = 8,
    partial_replay_candidate_limit: int = 64,
    minimum_partial_replay_rays: int = 0,
    maximum_candidates: int = 32,
    norm_log_tolerance: float = 0.45,
    angle_tolerance: float = 0.20,
    angle_hard_tolerance: float | None = None,
    maximum_sign_components: int = 16,
    partial_replay_ray_weight: float = 0.25,
    partial_replay_relation_weight: float = 0.10,
    finite_matroid_primes: Sequence[int] = (2, 3),
    finite_matroid_subset_size: int = 1,
) -> StarComponentSearchLedger:
    """Match bounded center stars before proposing a proper subspace.

    Each layer assigns or skips one complete source relation incident to the
    chosen center.  An assignment adds both noncentral rays simultaneously,
    so it can increase source rank while preserving an exact ternary edge.
    Candidate target edges are preordered by center-relative height data, and
    all cross-angles to the already assigned star are checked.  Exact partial
    replay, not metric fit, determines the returned survivors.
    """

    if (
        not source_center_indices
        or not target_center_indices
        or source_edge_limit < 1
        or target_edges_per_source < 1
        or beam_width < 1
        or minimum_states_per_match_rank_group < 1
        or maximum_partial_replay_attempts < 1
        or minimum_partial_replay_rank < 1
        or partial_replay_rank_stride < 1
        or partial_replay_vertex_stride < 1
        or partial_replay_candidate_limit < 1
        or minimum_partial_replay_rays < 0
        or maximum_candidates < 1
        or not finite_matroid_primes
        or any(int(prime) < 2 for prime in finite_matroid_primes)
        or finite_matroid_subset_size not in (1, 2, 3)
    ):
        raise ValueError("star-component search bounds must be positive")
    if angle_hard_tolerance is None:
        angle_hard_tolerance = 2.25 * angle_tolerance
    if angle_hard_tolerance < angle_tolerance:
        raise ValueError("hard angle tolerance must not be smaller than RMS tolerance")
    if any(index < 0 or index >= len(source.vertices) for index in source_center_indices):
        raise ValueError("source center is outside the relation complex")
    if any(index < 0 or index >= len(target.vertices) for index in target_center_indices):
        raise ValueError("target center is outside the relation complex")

    source_rank = rational_rank(source.vertices)
    source_norms, source_angles = _quadratic_data(source.vertices, source_gram)
    target_norms, target_angles = _quadratic_data(target.vertices, target_gram)
    source_by_vertex, _source_by_pair = _incidence(source)
    target_by_vertex, _target_by_pair = _incidence(target)
    center_pairs = tuple(
        (int(source_center), int(target_center))
        for source_center in source_center_indices
        for target_center in target_center_indices
    )
    states = []
    star_orders = {}
    star_proposals = {}
    for stratum, (source_center, target_center) in enumerate(center_pairs):
        scale = target_norms[target_center] / source_norms[source_center]
        mapping = [-1] * len(source.vertices)
        mapping[source_center] = target_center
        base = _State(
            stratum,
            tuple(mapping),
            frozenset(),
            0,
            0,
            0.0,
            0,
            float(scale),
            0,
            0,
        )
        proposals_by_edge = {}
        for source_edge_index in source_by_vertex[source_center]:
            source_other = tuple(
                index
                for index in source.ternary_relations[source_edge_index]
                if index != source_center
            )
            proposals = []
            for target_edge_index in target_by_vertex[target_center]:
                target_other = tuple(
                    index
                    for index in target.ternary_relations[target_edge_index]
                    if index != target_center
                )
                for ordering in (target_other, tuple(reversed(target_other))):
                    extension = _metric_extension(
                        base,
                        tuple(zip(source_other, ordering)),
                        source_norms,
                        target_norms,
                        source_angles,
                        target_angles,
                        norm_log_tolerance=norm_log_tolerance,
                        angle_tolerance=angle_tolerance,
                        angle_hard_tolerance=angle_hard_tolerance,
                    )
                    if extension is None:
                        continue
                    _mapping, squared, count = extension
                    proposals.append((sqrt(squared / max(1, count)), ordering))
            proposals_by_edge[source_edge_index] = tuple(
                ordering
                for _score, ordering in sorted(proposals)[:target_edges_per_source]
            )
        ordered_edges = sorted(
            source_by_vertex[source_center],
            key=lambda edge: (
                len(proposals_by_edge[edge]) or 10**9,
                -sum(
                    source.additive_degrees[index]
                    for index in source.ternary_relations[edge]
                ),
                edge,
            ),
        )[:source_edge_limit]
        star_orders[stratum] = tuple(ordered_edges)
        star_proposals[stratum] = proposals_by_edge
        states.append(base)

    rank_cache = {}

    def mapped_rank(state: _State) -> int:
        indices = tuple(
            index for index, target_index in enumerate(state.mapping) if target_index >= 0
        )
        if indices not in rank_cache:
            rank_cache[indices] = rational_rank(
                tuple(source.vertices[index] for index in indices)
            )
        return rank_cache[indices]

    def beam_key(state: _State):
        rank = mapped_rank(state)
        mapped_count = sum(index >= 0 for index in state.mapping)
        objective = (
            state.matched_edges
            + 0.5 * mapped_count
            + 0.5 * rank
            + partial_replay_ray_weight * state.partial_replay_rays
            + partial_replay_relation_weight * state.partial_replay_relations
            - 4.0 * state.rms
        )
        return (objective, rank, *_state_key(state))

    maximum_population = len(states)
    expanded = 0
    partial_attempts = 0
    maximum_rank = 1
    maximum_mapped = 1
    maximum_partial_rays = 0
    maximum_partial_relations = 0
    finite_rank_rejections = 0
    answers = {}
    source_modular_subset_ranks = {}
    target_modular_subset_ranks = {}

    def cached_modular_subset_rank(complex_, indices, prime, cache):
        key = (int(prime), tuple(sorted(int(index) for index in indices)))
        if key not in cache:
            cache[key] = modular_rank(
                tuple(complex_.vertices[index] for index in key[1]), key[0]
            )
        return cache[key]

    layer_count = max((len(order) for order in star_orders.values()), default=0)
    for layer in range(layer_count):
        next_states = []
        for state in states:
            order = star_orders[state.stratum]
            if layer >= len(order):
                next_states.append(state)
                continue
            source_edge_index = order[layer]
            source_edge = source.ternary_relations[source_edge_index]
            source_center = center_pairs[state.stratum][0]
            source_other = tuple(index for index in source_edge if index != source_center)
            processed = state.processed_edges | {source_edge_index}
            next_states.append(
                replace(
                    state,
                    processed_edges=processed,
                    skipped_edges=state.skipped_edges + 1,
                )
            )
            for target_ordering in star_proposals[state.stratum][source_edge_index]:
                if any(
                    state.mapping[source_index] >= 0
                    and state.mapping[source_index] != target_index
                    for source_index, target_index in zip(
                        source_other, target_ordering
                    )
                ):
                    continue
                additions = tuple(
                    (source_index, target_index)
                    for source_index, target_index in zip(
                        source_other, target_ordering
                    )
                    if state.mapping[source_index] < 0
                )
                extension = _metric_extension(
                    state,
                    additions,
                    source_norms,
                    target_norms,
                    source_angles,
                    target_angles,
                    norm_log_tolerance=norm_log_tolerance,
                    angle_tolerance=angle_tolerance,
                    angle_hard_tolerance=angle_hard_tolerance,
                )
                if extension is None:
                    continue
                mapping, squared, count = extension
                mapped_pairs = tuple(
                    (source_index, target_index)
                    for source_index, target_index in enumerate(mapping)
                    if target_index >= 0
                )
                if any(
                    cached_modular_subset_rank(
                        source,
                        tuple(index for index, _target in mapped_pairs),
                        prime,
                        source_modular_subset_ranks,
                    )
                    != cached_modular_subset_rank(
                        target,
                        tuple(index for _source, index in mapped_pairs),
                        prime,
                        target_modular_subset_ranks,
                    )
                    for prime in finite_matroid_primes
                ):
                    finite_rank_rejections += 1
                    continue
                if finite_matroid_subset_size >= 2 and additions:
                    new_sources = {source_index for source_index, _target in additions}
                    finite_incompatible = False
                    for subset_size in range(2, finite_matroid_subset_size + 1):
                        for subset in combinations(mapped_pairs, subset_size):
                            if not any(source_index in new_sources for source_index, _ in subset):
                                continue
                            if any(
                                cached_modular_subset_rank(
                                    source,
                                    tuple(
                                        source_index
                                        for source_index, _target_index in subset
                                    ),
                                    prime,
                                    source_modular_subset_ranks,
                                )
                                != cached_modular_subset_rank(
                                    target,
                                    tuple(
                                        target_index
                                        for _source_index, target_index in subset
                                    ),
                                    prime,
                                    target_modular_subset_ranks,
                                )
                                for prime in finite_matroid_primes
                            ):
                                finite_incompatible = True
                                break
                        if finite_incompatible:
                            break
                    if finite_incompatible:
                        finite_rank_rejections += 1
                        continue
                next_states.append(
                    replace(
                        state,
                        mapping=mapping,
                        processed_edges=processed,
                        matched_edges=state.matched_edges + 1,
                        squared_residual=state.squared_residual + squared,
                        residual_count=state.residual_count + count,
                    )
                )
            expanded += 1
        deduplicated = {}
        for state in next_states:
            previous = deduplicated.get(state.mapping)
            if previous is None or beam_key(state) > beam_key(previous):
                deduplicated[state.mapping] = state

        remaining = maximum_partial_replay_attempts - partial_attempts
        eligible = []
        if remaining:
            for state in deduplicated.values():
                rank = mapped_rank(state)
                mapped_count = sum(index >= 0 for index in state.mapping)
                if rank < minimum_partial_replay_rank or rank >= source_rank:
                    continue
                if state.partial_audited_rank and not (
                    rank >= state.partial_audited_rank + partial_replay_rank_stride
                    or mapped_count
                    >= state.partial_audited_mapped_count + partial_replay_vertex_stride
                ):
                    continue
                eligible.append(state)
        audit_limit = min(remaining, partial_replay_candidate_limit, len(eligible))
        if audit_limit:
            audit_states = _grouped_truncate(
                eligible,
                audit_limit,
                beam_key,
                lambda state: (
                    state.stratum,
                    state.matched_edges,
                    mapped_rank(state),
                ),
                1,
            )
            for state in audit_states:
                partial_attempts += 1
                rank = mapped_rank(state)
                mapped_count = sum(index >= 0 for index in state.mapping)
                replays = exact_partial_relation_replay(
                    source,
                    target,
                    state.mapping,
                    maximum_sign_components=maximum_sign_components,
                )
                best = max(
                    replays,
                    key=lambda replay: (
                        len(replay.replayed_source_vertex_indices),
                        replay.replayed_relation_count,
                    ),
                    default=None,
                )
                updated = replace(
                    state,
                    partial_audited_rank=rank,
                    partial_audited_mapped_count=mapped_count,
                )
                if best is not None:
                    replay_rays = len(best.replayed_source_vertex_indices)
                    updated = replace(
                        updated,
                        partial_replay_rays=max(state.partial_replay_rays, replay_rays),
                        partial_replay_relations=max(
                            state.partial_replay_relations,
                            best.replayed_relation_count,
                        ),
                        partial_replay_rank=max(
                            state.partial_replay_rank, best.replayed_source_rank
                        ),
                    )
                    maximum_partial_rays = max(maximum_partial_rays, replay_rays)
                    maximum_partial_relations = max(
                        maximum_partial_relations, best.replayed_relation_count
                    )
                    if replay_rays >= minimum_partial_replay_rays:
                        key = (
                            best.source_rank,
                            rational_nullspace(best.integral_matrix),
                        )
                        candidate = StarComponentCandidate(
                            vertex_map=state.mapping,
                            matched_star_relations=state.matched_edges,
                            skipped_star_relations=state.skipped_edges,
                            metric_rms_residual=state.rms,
                            replay=best,
                        )
                        previous = answers.get(key)
                        if previous is None or (
                            replay_rays,
                            best.replayed_relation_count,
                            -state.rms,
                        ) > (
                            len(previous.replay.replayed_source_vertex_indices),
                            previous.replay.replayed_relation_count,
                            -previous.metric_rms_residual,
                        ):
                            answers[key] = candidate
                deduplicated[state.mapping] = updated

        states = _grouped_truncate(
            tuple(deduplicated.values()),
            beam_width,
            beam_key,
            lambda state: (
                state.stratum,
                state.matched_edges,
                mapped_rank(state),
            ),
            minimum_states_per_match_rank_group,
        )
        maximum_population = max(maximum_population, len(states))
        maximum_rank = max(maximum_rank, *(mapped_rank(state) for state in states))
        maximum_mapped = max(
            maximum_mapped,
            *(sum(index >= 0 for index in state.mapping) for state in states),
        )
    ordered_answers = sorted(
        answers.values(),
        key=lambda candidate: (
            len(candidate.replay.replayed_source_vertex_indices),
            candidate.replay.replayed_relation_count,
            candidate.matched_star_relations,
            -candidate.metric_rms_residual,
        ),
        reverse=True,
    )[:maximum_candidates]
    return StarComponentSearchLedger(
        source_center_count=len(source_center_indices),
        target_center_count=len(target_center_indices),
        center_pair_count=len(center_pairs),
        processed_star_layer_count=layer_count,
        maximum_beam_population=maximum_population,
        expanded_state_count=expanded,
        finite_rank_rejection_count=finite_rank_rejections,
        partial_replay_attempt_count=partial_attempts,
        maximum_source_rank_reached=maximum_rank,
        maximum_mapped_vertex_count=maximum_mapped,
        maximum_partial_replay_ray_count=maximum_partial_rays,
        maximum_partial_replay_relation_count=maximum_partial_relations,
        candidates=tuple(ordered_answers),
    )


def bounded_metric_relation_search(
    source: RelationComplex,
    target: RelationComplex,
    source_gram: Sequence[Sequence[float]],
    target_gram: Sequence[Sequence[float]],
    *,
    source_center_limit: int = 20,
    target_center_limit: int = 80,
    source_center_indices: Sequence[int] | None = None,
    target_center_indices: Sequence[int] | None = None,
    center_pair_limit: int | None = None,
    initial_states_per_center_pair: int = 128,
    minimum_states_per_center_pair: int = 4,
    seed_edges_per_center: int = 4,
    beam_width: int = 2_000,
    maximum_steps: int = 80,
    maximum_exact_lift_attempts: int = 500,
    maximum_embeddings: int = 1,
    minimum_global_replay_rays: int = 0,
    reseed_after_skips: int = 40,
    reseed_source_limit: int = 8,
    reseed_target_limit: int = 4,
    reseed_state_limit: int = 64,
    norm_log_tolerance: float = 0.45,
    angle_tolerance: float = 0.20,
    angle_hard_tolerance: float | None = None,
    maximum_sign_components: int = 16,
    maximum_partial_replay_attempts: int = 0,
    minimum_partial_replay_rank: int = 10,
    partial_replay_rank_stride: int = 2,
    partial_replay_vertex_stride: int = 12,
    partial_replay_candidate_limit: int = 64,
    partial_replay_ray_weight: float = 0.25,
    partial_replay_relation_weight: float = 0.10,
    preserve_initial_state_lineages: bool = False,
) -> MetricRelationSearchLedger:
    """Search a finite box for exact full-rank rectangular embeddings.

    The algorithm seeds a high-incidence source star in a high-incidence
    target star, grows one- and two-vertex relation frontiers with a bounded
    beam, and invokes the exact lift as soon as the mapped source rays span.
    A skipped-edge branch is always retained because height cutoffs make the
    source and target relation clouds unequal.
    """

    if min(source_center_limit, target_center_limit, seed_edges_per_center) < 1:
        raise ValueError("center and seed bounds must be positive")
    if (
        beam_width < 1
        or maximum_steps < 1
        or maximum_exact_lift_attempts < 1
        or maximum_embeddings < 1
        or minimum_global_replay_rays < 0
        or reseed_after_skips < 1
        or reseed_source_limit < 1
        or reseed_target_limit < 1
        or reseed_state_limit < 1
        or maximum_partial_replay_attempts < 0
        or minimum_partial_replay_rank < 1
        or partial_replay_rank_stride < 1
        or partial_replay_vertex_stride < 1
        or partial_replay_candidate_limit < 1
    ):
        raise ValueError("beam, step, and exact-lift bounds must be positive")
    if initial_states_per_center_pair < 1 or minimum_states_per_center_pair < 1:
        raise ValueError("center-pair diversity bounds must be positive")
    if angle_hard_tolerance is None:
        angle_hard_tolerance = 2.25 * angle_tolerance
    if angle_hard_tolerance < angle_tolerance:
        raise ValueError("hard angle tolerance must not be smaller than RMS tolerance")
    source_rank = rational_rank(source.vertices)
    if source_rank != len(source.vertices[0]):
        raise ValueError("source complex must span its ambient coordinate lattice")
    source_norms, source_angles = _quadratic_data(source.vertices, source_gram)
    target_norms, target_angles = _quadratic_data(target.vertices, target_gram)
    source_by_vertex, _source_by_pair = _incidence(source)
    target_by_vertex, target_by_pair = _incidence(target)
    source_centers = (
        list(map(int, source_center_indices))
        if source_center_indices is not None
        else sorted(
            range(len(source.vertices)),
            key=lambda index: (
                -source.additive_degrees[index],
                source_norms[index],
                index,
            ),
        )[:source_center_limit]
    )
    target_centers = (
        list(map(int, target_center_indices))
        if target_center_indices is not None
        else sorted(
            range(len(target.vertices)),
            key=lambda index: (
                -target.additive_degrees[index],
                target_norms[index],
                index,
            ),
        )[:target_center_limit]
    )
    if (
        not source_centers
        or not target_centers
        or len(set(source_centers)) != len(source_centers)
        or len(set(target_centers)) != len(target_centers)
        or any(index < 0 or index >= len(source.vertices) for index in source_centers)
        or any(index < 0 or index >= len(target.vertices) for index in target_centers)
    ):
        raise ValueError("explicit center indices are empty, duplicated, or out of range")

    if center_pair_limit is None:
        center_pairs = [
            (source_center, target_center)
            for source_center in source_centers
            for target_center in target_centers
        ]
    else:
        if center_pair_limit < 1:
            raise ValueError("center-pair bound must be positive")
        source_features = {
            center: _star_features(
                source, source_by_vertex, source_norms, source_angles, center
            )
            for center in source_centers
        }
        target_features = {
            center: _star_features(
                target, target_by_vertex, target_norms, target_angles, center
            )
            for center in target_centers
        }
        scored_pairs = sorted(
            (
                _star_profile_distance(
                    source_features[source_center], target_features[target_center]
                ),
                source_center,
                target_center,
            )
            for source_center in source_centers
            for target_center in target_centers
        )
        center_pairs = [
            (source_center, target_center)
            for _score, source_center, target_center in scored_pairs[:center_pair_limit]
        ]

    states: list[_State] = []
    for stratum, (source_center, target_center) in enumerate(center_pairs):
        stratum_states = []
        source_edges = sorted(
            source_by_vertex[source_center],
            key=lambda edge: (
                -sum(source.additive_degrees[index] for index in source.ternary_relations[edge]),
                edge,
            ),
        )[:seed_edges_per_center]
        scale = target_norms[target_center] / source_norms[source_center]
        base_mapping = [-1] * len(source.vertices)
        base_mapping[source_center] = target_center
        base = _State(
            stratum,
            tuple(base_mapping),
            frozenset(),
            0,
            0,
            0.0,
            0,
            float(scale),
            0,
            0,
        )
        for source_edge_index in source_edges:
            source_edge = source.ternary_relations[source_edge_index]
            source_other = [index for index in source_edge if index != source_center]
            for target_edge_index in target_by_vertex[target_center]:
                target_edge = target.ternary_relations[target_edge_index]
                target_other = [index for index in target_edge if index != target_center]
                for ordered_target in (target_other, tuple(reversed(target_other))):
                    extension = _metric_extension(
                        base,
                        tuple(zip(source_other, ordered_target)),
                        source_norms,
                        target_norms,
                        source_angles,
                        target_angles,
                        norm_log_tolerance=norm_log_tolerance,
                        angle_tolerance=angle_tolerance,
                        angle_hard_tolerance=angle_hard_tolerance,
                    )
                    if extension is None:
                        continue
                    mapping, squared, count = extension
                    stratum_states.append(
                        _State(
                            stratum,
                            mapping,
                            frozenset((source_edge_index,)),
                            1,
                            0,
                            squared,
                            count,
                            float(scale),
                            0,
                            0,
                        )
                    )
        states.extend(
            sorted(stratum_states, key=_state_key, reverse=True)[
                :initial_states_per_center_pair
            ]
        )
    if preserve_initial_state_lineages:
        # This is a declared search-cost tradeoff: when the beam is at least
        # this wide, every retained seed gets one descendant slot until exact
        # partial replay can distinguish the lineages.
        states = [replace(state, stratum=index) for index, state in enumerate(states)]
    initial_state_count = len(states)
    # Deduplicate identical partial maps before the first beam truncation.
    states_by_map = {}
    for state in states:
        previous = states_by_map.get(state.mapping)
        if previous is None or _state_key(state) > _state_key(previous):
            states_by_map[state.mapping] = state
    states = _diverse_truncate(
        tuple(states_by_map.values()),
        beam_width,
        _state_key,
        minimum_states_per_center_pair,
    )
    maximum_population = len(states)
    expanded = 0
    exact_attempts = 0
    maximum_rank_reached = 0
    maximum_mapped_vertices = max(
        (sum(target >= 0 for target in state.mapping) for state in states),
        default=0,
    )
    maximum_matched_relations = max(
        (state.matched_edges for state in states), default=0
    )
    maximum_global_replay = 0
    below_global_support = 0
    maximum_reseeds = 0
    partial_attempts = 0
    maximum_partial_replay = 0
    maximum_partial_relations = 0
    answers: list[MetricRelationEmbedding] = []
    seen_exact_maps = set()
    seen_embedding_spaces = set()
    source_rank_cache: dict[tuple[int, ...], int] = {}

    def mapped_source_rank(indices: Sequence[int]) -> int:
        key = tuple(indices)
        if key not in source_rank_cache:
            source_rank_cache[key] = rational_rank(
                tuple(source.vertices[index] for index in key)
            )
        return source_rank_cache[key]

    for _step in range(maximum_steps):
        next_states: list[_State] = []
        reseed_eligible = [
            state
            for state in states
            if state.stagnation == reseed_after_skips
            or (
                state.reseeds == 0
                and mapped_source_rank(
                    tuple(
                        index
                        for index, target_index in enumerate(state.mapping)
                        if target_index >= 0
                    )
                )
                == source_rank - 1
            )
        ]
        reseed_allowed = {
            id(state)
            for state in _diverse_truncate(
                reseed_eligible,
                reseed_state_limit,
                _state_key,
                1,
            )
        }
        for state in states:
            mapped_sources = [
                source_index
                for source_index, target_index in enumerate(state.mapping)
                if target_index >= 0
            ]
            current_rank = mapped_source_rank(mapped_sources)
            maximum_rank_reached = max(maximum_rank_reached, current_rank)
            if current_rank == source_rank and exact_attempts < maximum_exact_lift_attempts:
                key = tuple((index, state.mapping[index]) for index in mapped_sources)
                if key not in seen_exact_maps:
                    seen_exact_maps.add(key)
                    exact_attempts += 1
                    subcomplex = build_relation_complex(
                        tuple(source.vertices[index] for index in mapped_sources)
                    )
                    source_index = {
                        vertex: index for index, vertex in enumerate(source.vertices)
                    }
                    vertex_map = tuple(
                        state.mapping[source_index[vertex]]
                        for vertex in subcomplex.vertices
                    )
                    try:
                        lifts = lift_relation_vertex_injection(
                            subcomplex,
                            target,
                            vertex_map,
                            maximum_sign_components=maximum_sign_components,
                        )
                    except ValueError:
                        lifts = ()
                    if lifts:
                        space_key = rational_nullspace(lifts[0])
                        if space_key not in seen_embedding_spaces:
                            seen_embedding_spaces.add(space_key)
                            replay = _global_replay(source, target, lifts[0])
                            maximum_global_replay = max(
                                maximum_global_replay, replay[0]
                            )
                            if replay[0] < minimum_global_replay_rays:
                                below_global_support += 1
                                continue
                            answers.append(
                                MetricRelationEmbedding(
                                    vertex_map=vertex_map,
                                    source_vertex_indices=tuple(
                                        source_index[vertex]
                                        for vertex in subcomplex.vertices
                                    ),
                                    integral_matrices=lifts,
                                    matched_seed_relations=state.matched_edges,
                                    skipped_source_relations=state.skipped_edges,
                                    metric_rms_residual=state.rms,
                                    global_replay_ray_count=replay[0],
                                    global_replay_source_rank=replay[1],
                                    global_replay_relation_count=replay[2],
                                    metric_reseed_count=state.reseeds,
                                )
                            )
                        if len(answers) >= maximum_embeddings:
                            return MetricRelationSearchLedger(
                                len(source_centers),
                                len(target_centers),
                                initial_state_count,
                                maximum_population,
                                expanded,
                                exact_attempts,
                                maximum_rank_reached,
                                maximum_mapped_vertices,
                                maximum_matched_relations,
                                maximum_global_replay,
                                below_global_support,
                                maximum_reseeds,
                                partial_attempts,
                                maximum_partial_replay,
                                maximum_partial_relations,
                                tuple(answers),
                            )

            frontier = set()
            for vertex in mapped_sources:
                frontier.update(source_by_vertex[vertex])
            frontier.difference_update(state.processed_edges)
            if id(state) in reseed_allowed:
                used_targets = {target for target in state.mapping if target >= 0}
                span_normals = rational_nullspace(
                    tuple(source.vertices[index] for index in mapped_sources)
                )
                reseed_sources = []
                for source_index in sorted(
                    (
                        index
                        for index, target_index in enumerate(state.mapping)
                        if target_index < 0
                    ),
                    key=lambda index: (
                        -source.additive_degrees[index],
                        source_norms[index],
                        index,
                    ),
                ):
                    if not any(
                        sum(
                            int(value) * int(normal[column])
                            for column, value in enumerate(source.vertices[source_index])
                        )
                        for normal in span_normals
                    ):
                        continue
                    reseed_sources.append(source_index)
                    if len(reseed_sources) == reseed_source_limit:
                        break
                for source_index in reseed_sources:
                    old_sources = np.asarray(mapped_sources, dtype=int)
                    old_targets = np.asarray(
                        [state.mapping[index] for index in mapped_sources], dtype=int
                    )
                    norm_errors = np.abs(
                        np.log((target_norms / state.scale) / source_norms[source_index])
                    )
                    angle_errors = np.abs(
                        target_angles[:, old_targets]
                        - source_angles[source_index, old_sources][None, :]
                    )
                    mask = norm_errors <= norm_log_tolerance
                    if len(old_targets):
                        mask &= np.max(angle_errors, axis=1) <= angle_hard_tolerance
                        mask &= np.sqrt(np.mean(angle_errors**2, axis=1)) <= angle_tolerance
                    if used_targets:
                        mask[np.asarray(sorted(used_targets), dtype=int)] = False
                    metric_scores = norm_errors**2
                    if len(old_targets):
                        metric_scores += np.mean(angle_errors**2, axis=1)
                    target_candidates = np.flatnonzero(mask)
                    target_candidates = target_candidates[
                        np.argsort(metric_scores[target_candidates])[:reseed_target_limit]
                    ]
                    for target_index in target_candidates:
                        extension = _metric_extension(
                            state,
                            ((source_index, int(target_index)),),
                            source_norms,
                            target_norms,
                            source_angles,
                            target_angles,
                            norm_log_tolerance=norm_log_tolerance,
                            angle_tolerance=angle_tolerance,
                            angle_hard_tolerance=angle_hard_tolerance,
                        )
                        if extension is None:
                            continue
                        mapping, squared, count = extension
                        next_states.append(
                            replace(
                                state,
                                mapping=mapping,
                                squared_residual=state.squared_residual + squared,
                                residual_count=state.residual_count + count,
                                stagnation=0,
                                reseeds=state.reseeds + 1,
                            )
                        )
                        maximum_reseeds = max(maximum_reseeds, state.reseeds + 1)
            if not frontier:
                continue
            mapped_set = set(mapped_sources)
            edge_index = max(
                frontier,
                key=lambda edge: (
                    len(mapped_set & set(source.ternary_relations[edge])),
                    sum(
                        source.additive_degrees[index]
                        for index in source.ternary_relations[edge]
                    ),
                    -edge,
                ),
            )
            source_edge = source.ternary_relations[edge_index]
            known = [index for index in source_edge if state.mapping[index] >= 0]
            missing = [index for index in source_edge if state.mapping[index] < 0]
            processed = state.processed_edges | {edge_index}
            # The explicit skip branch is the cutoff-robust part of the search.
            next_states.append(
                replace(
                    state,
                    processed_edges=processed,
                    skipped_edges=state.skipped_edges + 1,
                    stagnation=state.stagnation + 1,
                )
            )
            if len(known) == 2:
                target_pair = tuple(sorted(state.mapping[index] for index in known))
                target_edges = target_by_pair.get(target_pair, ())
            elif len(known) == 1:
                target_edges = target_by_vertex[state.mapping[known[0]]]
            else:
                target_edges = ()
            for target_edge_index in target_edges:
                target_edge = target.ternary_relations[target_edge_index]
                if any(state.mapping[index] not in target_edge for index in known):
                    continue
                available = [
                    index
                    for index in target_edge
                    if index not in {state.mapping[source_index] for source_index in known}
                ]
                orderings = (available,) if len(available) < 2 else (
                    available,
                    tuple(reversed(available)),
                )
                for ordered_target in orderings:
                    extension = _metric_extension(
                        state,
                        tuple(zip(missing, ordered_target)),
                        source_norms,
                        target_norms,
                        source_angles,
                        target_angles,
                        norm_log_tolerance=norm_log_tolerance,
                        angle_tolerance=angle_tolerance,
                        angle_hard_tolerance=angle_hard_tolerance,
                    )
                    if extension is None:
                        continue
                    mapping, squared, count = extension
                    next_states.append(
                        replace(
                            state,
                            mapping=mapping,
                            processed_edges=processed,
                            matched_edges=state.matched_edges + 1,
                            squared_residual=state.squared_residual + squared,
                            residual_count=state.residual_count + count,
                            stagnation=0,
                        )
                    )
            expanded += 1

        if not next_states:
            break
        deduplicated = {}
        for state in next_states:
            previous = deduplicated.get(state.mapping)
            if previous is None or _state_key(state) > _state_key(previous):
                deduplicated[state.mapping] = state

        # Audit a bounded, stratum-diverse sample of proper subspaces before
        # beam truncation.  A successful replay is exact and remains a valid
        # lower bound along descendants of that state; height never enters the
        # certificate.  Re-audits require a declared rank or vertex increase.
        remaining_partial = maximum_partial_replay_attempts - partial_attempts
        if remaining_partial > 0:
            partial_eligible = []
            for state in deduplicated.values():
                mapped_indices = tuple(
                    index
                    for index, target_index in enumerate(state.mapping)
                    if target_index >= 0
                )
                rank = mapped_source_rank(mapped_indices)
                mapped_count = len(mapped_indices)
                if rank < minimum_partial_replay_rank or rank >= source_rank:
                    continue
                if state.partial_audited_rank and not (
                    rank >= state.partial_audited_rank + partial_replay_rank_stride
                    or mapped_count
                    >= state.partial_audited_mapped_count + partial_replay_vertex_stride
                ):
                    continue
                partial_eligible.append(state)
            audit_limit = min(
                remaining_partial,
                partial_replay_candidate_limit,
                len(partial_eligible),
            )
            if audit_limit:
                audit_states = _diverse_truncate(
                    partial_eligible,
                    audit_limit,
                    _state_key,
                    1,
                )
                for state in audit_states:
                    mapped_indices = tuple(
                        index
                        for index, target_index in enumerate(state.mapping)
                        if target_index >= 0
                    )
                    rank = mapped_source_rank(mapped_indices)
                    partial_attempts += 1
                    replays = exact_partial_relation_replay(
                        source,
                        target,
                        state.mapping,
                        maximum_sign_components=maximum_sign_components,
                    )
                    best = max(
                        replays,
                        key=lambda replay: (
                            len(replay.replayed_source_vertex_indices),
                            replay.replayed_relation_count,
                        ),
                        default=None,
                    )
                    updated = replace(
                        state,
                        partial_audited_rank=rank,
                        partial_audited_mapped_count=len(mapped_indices),
                    )
                    if best is not None:
                        replay_rays = len(best.replayed_source_vertex_indices)
                        updated = replace(
                            updated,
                            partial_replay_rays=max(
                                state.partial_replay_rays, replay_rays
                            ),
                            partial_replay_relations=max(
                                state.partial_replay_relations,
                                best.replayed_relation_count,
                            ),
                            partial_replay_rank=max(
                                state.partial_replay_rank, best.replayed_source_rank
                            ),
                        )
                        maximum_partial_replay = max(
                            maximum_partial_replay, replay_rays
                        )
                        maximum_partial_relations = max(
                            maximum_partial_relations,
                            best.replayed_relation_count,
                        )
                    deduplicated[state.mapping] = updated

        def beam_key(state):
            mapped_indices = tuple(
                index
                for index, target_index in enumerate(state.mapping)
                if target_index >= 0
            )
            rank = mapped_source_rank(mapped_indices)
            # Avoid lexicographically forcing a false rank gain: relation
            # support and metric fit remain in the same scalar objective.
            mapped_count = sum(target >= 0 for target in state.mapping)
            objective = (
                0.5 * rank
                + 1.0 * state.matched_edges
                + 0.5 * mapped_count
                + partial_replay_ray_weight * state.partial_replay_rays
                + partial_replay_relation_weight * state.partial_replay_relations
                - 4.0 * state.rms
                - 0.10 * state.reseeds
            )
            return (objective, rank, *_state_key(state))
        states = _diverse_truncate(
            tuple(deduplicated.values()),
            beam_width,
            beam_key,
            minimum_states_per_center_pair,
        )
        maximum_population = max(maximum_population, len(states))
        maximum_mapped_vertices = max(
            maximum_mapped_vertices,
            max(sum(target >= 0 for target in state.mapping) for state in states),
        )
        maximum_matched_relations = max(
            maximum_matched_relations,
            max(state.matched_edges for state in states),
        )
    return MetricRelationSearchLedger(
        len(source_centers),
        len(target_centers),
        initial_state_count,
        maximum_population,
        expanded,
        exact_attempts,
        maximum_rank_reached,
        maximum_mapped_vertices,
        maximum_matched_relations,
        maximum_global_replay,
        below_global_support,
        maximum_reseeds,
        partial_attempts,
        maximum_partial_replay,
        maximum_partial_relations,
        tuple(answers),
    )

"""Exact relation components and finite-aware proper-subspace merges.

The routines in this module deliberately separate proposal from acceptance.
Stars, overlapping stars, graph cores, biconnected blocks, and sampled dense
hyperplanes are finite proposal units.  A proposed merge is accepted only
after exact rational-rank checks and :func:`exact_partial_relation_replay`.
Height angles are rejection filters, never identities.
"""

from __future__ import annotations

from dataclasses import dataclass
from random import Random
from typing import Sequence

import numpy as np

from .integer import modular_rank, rational_nullspace, rational_rank
from .codes import CandidateFiniteSignature, finite_signature_distance
from .fingerprints import (
    CandidateRelationFingerprint,
    robust_standardize_fingerprint_families,
)
from .matching import PartialRelationReplay, exact_partial_relation_replay
from .relations import RelationComplex


@dataclass(frozen=True)
class RelationComponent:
    """A vertex-induced finite component of a ternary relation complex."""

    kind: str
    vertex_indices: tuple[int, ...]
    relation_indices: tuple[int, ...]
    rational_rank: int
    mod2_rank: int
    mod3_rank: int

    @property
    def relation_density(self) -> float:
        return len(self.relation_indices) / max(1, self.rational_rank)


@dataclass(frozen=True)
class DenseHyperplaneLedger:
    """Declared deterministic sampling box for dense proper subspaces."""

    ambient_rank: int
    random_seed: int
    sample_count: int
    independent_sample_count: int
    distinct_hyperplane_count: int
    components: tuple[RelationComponent, ...]


@dataclass(frozen=True)
class ExactComponentMerge:
    """One exact finite-aware merge, scored only after replay."""

    replay: PartialRelationReplay
    left_source_rank: int
    right_source_rank: int
    added_rank: int
    held_out_replayed_ray_count: int
    held_out_rays_per_added_rank: float
    height_angle_rms: float


@dataclass(frozen=True)
class ReplayedComponentCandidate:
    """One exact rational component after global replay and deduplication.

    ``origin_indices`` records every proposal occurrence which produced the
    same rational space.  Development and held-out sets belong to the ambient
    fibre, not to an enclosing parent; this is what makes replay counts
    comparable across different parent proposals.
    """

    basis_rows: tuple[tuple[int, ...], ...]
    origin_indices: tuple[int, ...]
    development_replayed_ray_indices: tuple[int, ...]
    held_out_replayed_ray_indices: tuple[int, ...]
    full_replayed_ray_indices: tuple[int, ...]
    full_replayed_relation_indices: tuple[int, ...]
    rational_rank: int
    modular_ranks: tuple[tuple[int, int], ...]

    @property
    def held_out_rays_per_rank(self) -> float:
        return len(self.held_out_replayed_ray_indices) / max(1, self.rational_rank)


@dataclass(frozen=True)
class JointComponentBundle:
    """One cross-fibre bundle under equal-weight structural channels."""

    candidate_indices: tuple[int, ...]
    relation_distance_percentile: float
    hermite_distance_percentile: float
    finite_distance_percentile: float
    mean_combined_distance_percentile: float
    maximum_combined_distance_percentile: float
    held_out_replayed_ray_count: int
    held_out_rays_per_rank: float


@dataclass(frozen=True)
class JointComponentBundleLedger:
    """Frozen structural pruning followed by exact replay optimization."""

    structural_quantile: float
    generated_bundle_count: int
    eligible_bundle_count: int
    structurally_retained_bundle_count: int
    selected: JointComponentBundle
    bundles: tuple[JointComponentBundle, ...]


def exact_rational_space_key(
    basis_rows: Sequence[Sequence[int]],
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    """Return a basis-independent exact key for a rational row space.

    The right nullspace is computed from an exact RREF with deterministic
    pivot order.  Its primitive oriented rows therefore canonically describe
    the rational space, including when the ambient coordinate width exceeds
    the candidate rank by more than one.
    """

    rank = rational_rank(basis_rows)
    if not basis_rows or rank != len(basis_rows):
        raise ValueError("space key needs a nonempty independent row basis")
    width = len(basis_rows[0])
    if any(len(row) != width for row in basis_rows):
        raise ValueError("space-key row widths differ")
    kernel = () if rank == width else rational_nullspace(basis_rows)
    return rank, kernel


def replay_and_deduplicate_components(
    vectors: Sequence[Sequence[int]],
    complex_: RelationComplex,
    candidate_basis_rows: Sequence[Sequence[Sequence[int]]],
    *,
    development_indices: Sequence[int],
    held_out_indices: Sequence[int],
    finite_primes: Sequence[int] = (2, 3),
) -> tuple[ReplayedComponentCandidate, ...]:
    """Replay proposed spaces globally and collapse exact duplicates.

    The input indices refer to ``vectors``.  They must form a disjoint split;
    vectors outside both sets are allowed but are counted only by full replay.
    Every membership and induced-relation count is exact.  The function does
    not saturate proposal bases: callers must do that immediately after each
    component construction or merge, before invoking this audit.
    """

    if not vectors or not candidate_basis_rows:
        return ()
    development = frozenset(map(int, development_indices))
    held_out = frozenset(map(int, held_out_indices))
    if development & held_out:
        raise ValueError("development and held-out ray sets overlap")
    if development | held_out and (
        min(development | held_out) < 0 or max(development | held_out) >= len(vectors)
    ):
        raise ValueError("development/held-out index outside the vector population")
    vertex_by_ray = {
        tuple(vector): index for index, vector in enumerate(complex_.vertices)
    }
    # RelationComplex vertices are canonical unoriented rays; normalize the
    # caller population through the same constructor-independent convention.
    from .integer import canonical_unoriented

    input_to_vertex = tuple(vertex_by_ray[canonical_unoriented(vector)] for vector in vectors)
    by_key: dict[
        tuple[int, tuple[tuple[int, ...], ...]],
        dict[str, object],
    ] = {}
    for origin, raw_basis in enumerate(candidate_basis_rows):
        basis = tuple(tuple(map(int, row)) for row in raw_basis)
        key = exact_rational_space_key(basis)
        rank = key[0]
        normals = key[1]
        masks = [
            all(
                sum(int(value) * int(coefficient) for value, coefficient in zip(vector, normal))
                == 0
                for normal in normals
            )
            for vector in vectors
        ]
        replayed = tuple(index for index, keep in enumerate(masks) if keep)
        replayed_vertices = frozenset(input_to_vertex[index] for index in replayed)
        relation_indices = tuple(
            index
            for index, relation in enumerate(complex_.ternary_relations)
            if set(relation) <= replayed_vertices
        )
        modular_ranks = tuple(
            (int(prime), modular_rank(basis, int(prime))) for prime in finite_primes
        )
        previous = by_key.get(key)
        if previous is None:
            by_key[key] = {
                "basis": basis,
                "origins": [origin],
                "development": tuple(index for index in replayed if index in development),
                "held_out": tuple(index for index in replayed if index in held_out),
                "full": replayed,
                "relations": relation_indices,
                "modular_ranks": modular_ranks,
            }
            continue
        previous["origins"].append(origin)
        # Equality of the exact key already proves equality of rational
        # spaces; replay equality is checked defensively against programming
        # or canonicalization mistakes.
        if previous["full"] != replayed or previous["relations"] != relation_indices:
            raise ArithmeticError("equal rational-space keys produced unequal replay")
    answer = [
        ReplayedComponentCandidate(
            basis_rows=record["basis"],
            origin_indices=tuple(record["origins"]),
            development_replayed_ray_indices=record["development"],
            held_out_replayed_ray_indices=record["held_out"],
            full_replayed_ray_indices=record["full"],
            full_replayed_relation_indices=record["relations"],
            rational_rank=key[0],
            modular_ranks=record["modular_ranks"],
        )
        for key, record in by_key.items()
    ]
    return tuple(
        sorted(
            answer,
            key=lambda item: (
                -item.held_out_rays_per_rank,
                -len(item.held_out_replayed_ray_indices),
                -len(item.full_replayed_ray_indices),
                -len(item.full_replayed_relation_indices),
                item.origin_indices,
            ),
        )
    )


def _distance_percentiles(matrix: np.ndarray) -> np.ndarray:
    """Replace finite distances by empirical lower-is-better percentiles."""

    values = np.asarray(matrix, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("distance channel contains no finite population")
    ordered = np.sort(np.unique(values.ravel()))
    if len(ordered) == 1:
        return np.zeros_like(values)
    return np.searchsorted(ordered, values, side="left") / (len(ordered) - 1)


def joint_component_bundle_ledger(
    relation_families: Sequence[Sequence[CandidateRelationFingerprint]],
    hermite_families: Sequence[Sequence[float | str]],
    finite_families: Sequence[Sequence[CandidateFiniteSignature]],
    held_out_replay_families: Sequence[Sequence[int]],
    candidate_ranks: Sequence[int],
    *,
    eligible_families: Sequence[Sequence[bool]] | None = None,
    structural_quantile: float = 0.10,
) -> JointComponentBundleLedger:
    """Match components jointly, then optimize held-out rays per rank.

    Relation/height, primitive-Hermite, and development finite-code distances
    are each converted to empirical percentiles for every pair of fibres and
    receive equal weight.  Candidate bundles are generated by nearest matches
    from every possible anchor.  The best declared structural quantile is an
    early-pruning set; only inside it does exact held-out replay decide.

    This routine deliberately accepts only *development* finite signatures.
    Supplying disjoint validation signatures is a caller responsibility, and
    they must not be passed here.
    """

    fibre_count = len(relation_families)
    if fibre_count < 3 or len(candidate_ranks) != fibre_count:
        raise ValueError("joint component selection needs at least three fibres")
    if not 0.0 < structural_quantile <= 1.0:
        raise ValueError("structural quantile must lie in (0,1]")
    if not (
        len(hermite_families)
        == len(finite_families)
        == len(held_out_replay_families)
        == fibre_count
    ):
        raise ValueError("joint component family counts differ")
    widths = tuple(len(family) for family in relation_families)
    if any(width == 0 for width in widths):
        raise ValueError("joint component families must be nonempty")
    if any(
        len(hermite_families[index]) != width
        or len(finite_families[index]) != width
        or len(held_out_replay_families[index]) != width
        for index, width in enumerate(widths)
    ):
        raise ValueError("component channel widths differ")
    if eligible_families is None:
        eligible_families = tuple(tuple(True for _ in range(width)) for width in widths)
    if len(eligible_families) != fibre_count or any(
        len(eligible_families[index]) != width for index, width in enumerate(widths)
    ):
        raise ValueError("eligibility channel widths differ")

    relation_matrices = robust_standardize_fingerprint_families(relation_families)
    hermite_values = tuple(np.asarray(family, dtype=float) for family in hermite_families)
    pair_channels: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for left in range(fibre_count):
        for right in range(left + 1, fibre_count):
            relation = np.sqrt(
                np.mean(
                    (
                        relation_matrices[left][:, None, :]
                        - relation_matrices[right][None, :, :]
                    )
                    ** 2,
                    axis=2,
                )
            )
            hermite = np.abs(
                hermite_values[left][:, None] - hermite_values[right][None, :]
            )
            finite = np.asarray(
                [
                    [
                        finite_signature_distance(
                            left_signature,
                            right_signature,
                            include_components=False,
                        )
                        for right_signature in finite_families[right]
                    ]
                    for left_signature in finite_families[left]
                ],
                dtype=float,
            )
            relation_percentile = _distance_percentiles(relation)
            hermite_percentile = _distance_percentiles(hermite)
            finite_percentile = _distance_percentiles(finite)
            combined = (
                relation_percentile + hermite_percentile + finite_percentile
            ) / 3.0
            pair_channels[(left, right)] = (
                relation_percentile,
                hermite_percentile,
                finite_percentile,
                combined,
            )

    keys = set()
    for anchor in range(fibre_count):
        for candidate in range(widths[anchor]):
            key = [None] * fibre_count
            key[anchor] = candidate
            for other in range(fibre_count):
                if other == anchor:
                    continue
                if anchor < other:
                    combined = pair_channels[(anchor, other)][3]
                    choice = int(np.argmin(combined[candidate]))
                else:
                    combined = pair_channels[(other, anchor)][3]
                    choice = int(np.argmin(combined[:, candidate]))
                key[other] = choice
            keys.add(tuple(map(int, key)))

    bundles = []
    for key in sorted(keys):
        if not all(eligible_families[fibre][candidate] for fibre, candidate in enumerate(key)):
            continue
        relation_pieces = []
        hermite_pieces = []
        finite_pieces = []
        combined_pieces = []
        for left in range(fibre_count):
            for right in range(left + 1, fibre_count):
                channels = pair_channels[(left, right)]
                indices = (key[left], key[right])
                relation_pieces.append(float(channels[0][indices]))
                hermite_pieces.append(float(channels[1][indices]))
                finite_pieces.append(float(channels[2][indices]))
                combined_pieces.append(float(channels[3][indices]))
        held_out = sum(
            int(held_out_replay_families[fibre][candidate])
            for fibre, candidate in enumerate(key)
        )
        total_rank = sum(map(int, candidate_ranks))
        bundles.append(
            JointComponentBundle(
                candidate_indices=key,
                relation_distance_percentile=float(np.mean(relation_pieces)),
                hermite_distance_percentile=float(np.mean(hermite_pieces)),
                finite_distance_percentile=float(np.mean(finite_pieces)),
                mean_combined_distance_percentile=float(np.mean(combined_pieces)),
                maximum_combined_distance_percentile=float(np.max(combined_pieces)),
                held_out_replayed_ray_count=held_out,
                held_out_rays_per_rank=held_out / max(1, total_rank),
            )
        )
    if not bundles:
        raise ArithmeticError("no joint component bundle survived eligibility pruning")
    structural = sorted(
        bundles,
        key=lambda item: (
            item.maximum_combined_distance_percentile,
            item.mean_combined_distance_percentile,
            -item.held_out_rays_per_rank,
            item.candidate_indices,
        ),
    )
    retained_count = max(1, int(np.ceil(structural_quantile * len(structural))))
    retained = structural[:retained_count]
    selected = max(
        retained,
        key=lambda item: (
            item.held_out_rays_per_rank,
            item.held_out_replayed_ray_count,
            -item.maximum_combined_distance_percentile,
            -item.mean_combined_distance_percentile,
            tuple(-index for index in item.candidate_indices),
        ),
    )
    ordered = tuple(
        sorted(
            bundles,
            key=lambda item: (
                item.maximum_combined_distance_percentile,
                item.mean_combined_distance_percentile,
                -item.held_out_rays_per_rank,
                item.candidate_indices,
            ),
        )
    )
    return JointComponentBundleLedger(
        structural_quantile=float(structural_quantile),
        generated_bundle_count=len(keys),
        eligible_bundle_count=len(bundles),
        structurally_retained_bundle_count=retained_count,
        selected=selected,
        bundles=ordered,
    )


def height_angle_profile(
    vectors: Sequence[Sequence[int]],
    gram: Sequence[Sequence[float | str]],
    *,
    quantile_count: int = 17,
) -> tuple[float, ...]:
    """Return scale-free quantiles of absolute pairwise height angles."""

    if len(vectors) < 2 or quantile_count < 2:
        raise ValueError("an angle profile needs two vectors and two quantiles")
    matrix = np.asarray(vectors, dtype=float)
    form = np.asarray(gram, dtype=float)
    if form.shape != (matrix.shape[1], matrix.shape[1]):
        raise ValueError("Gram shape differs from vector width")
    norms = np.einsum("ij,jk,ik->i", matrix, form, matrix)
    if np.any(norms <= 0):
        raise ValueError("height-angle profile requires a positive form")
    pairings = matrix @ form @ matrix.T
    angles = np.abs(pairings / np.sqrt(norms[:, None] * norms[None, :]))
    values = angles[np.triu_indices(len(matrix), 1)]
    return tuple(
        map(float, np.quantile(values, np.linspace(0.0, 1.0, quantile_count)))
    )


def height_angle_profile_distance(
    left: Sequence[float], right: Sequence[float]
) -> float:
    """Return RMS distance between two equal-length angle profiles."""

    if not left or len(left) != len(right):
        raise ValueError("height-angle profiles must have equal positive length")
    return float(
        np.sqrt(np.mean(np.square(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))))
    )


def _incidence(complex_: RelationComplex) -> tuple[tuple[int, ...], ...]:
    by_vertex = [[] for _ in complex_.vertices]
    for relation_index, relation in enumerate(complex_.ternary_relations):
        for vertex in relation:
            by_vertex[vertex].append(relation_index)
    return tuple(tuple(items) for items in by_vertex)


def _component(
    complex_: RelationComplex,
    vertices: Sequence[int],
    *,
    kind: str,
) -> RelationComponent:
    vertex_indices = tuple(sorted(set(map(int, vertices))))
    selected = set(vertex_indices)
    relations = tuple(
        index
        for index, relation in enumerate(complex_.ternary_relations)
        if set(relation) <= selected
    )
    rows = tuple(complex_.vertices[index] for index in vertex_indices)
    return RelationComponent(
        kind=kind,
        vertex_indices=vertex_indices,
        relation_indices=relations,
        rational_rank=rational_rank(rows),
        mod2_rank=modular_rank(rows, 2),
        mod3_rank=modular_rank(rows, 3),
    )


def maximal_star_components(
    complex_: RelationComplex,
    *,
    minimum_relations: int = 2,
) -> tuple[RelationComponent, ...]:
    """Return every sufficiently large center star, without naming centers."""

    if minimum_relations < 1:
        raise ValueError("minimum_relations must be positive")
    incidence = _incidence(complex_)
    answers = []
    for center, relation_indices in enumerate(incidence):
        if len(relation_indices) < minimum_relations:
            continue
        vertices = {center}
        for relation_index in relation_indices:
            vertices.update(complex_.ternary_relations[relation_index])
        answers.append(_component(complex_, vertices, kind="maximal_star"))
    return tuple(
        sorted(
            answers,
            key=lambda item: (
                -len(item.relation_indices),
                -len(item.vertex_indices),
                item.vertex_indices,
            ),
        )
    )


def overlapping_star_components(
    complex_: RelationComplex,
    *,
    minimum_shared_relation: bool = True,
    maximum_components: int | None = None,
) -> tuple[RelationComponent, ...]:
    """Return unions of two stars whose centers are relation-adjacent."""

    incidence = _incidence(complex_)
    adjacent = set()
    for relation in complex_.ternary_relations:
        for left in range(3):
            for right in range(left):
                adjacent.add(tuple(sorted((relation[left], relation[right]))))
    if not minimum_shared_relation:
        centers = [index for index, items in enumerate(incidence) if items]
        adjacent = {
            (centers[left], centers[right])
            for left in range(len(centers))
            for right in range(left)
        }
    answers = {}
    for left, right in sorted(adjacent):
        vertices = {left, right}
        for relation_index in (*incidence[left], *incidence[right]):
            vertices.update(complex_.ternary_relations[relation_index])
        component = _component(complex_, vertices, kind="overlapping_stars")
        answers.setdefault(component.vertex_indices, component)
    ordered = sorted(
        answers.values(),
        key=lambda item: (
            -item.relation_density,
            -len(item.relation_indices),
            item.vertex_indices,
        ),
    )
    return tuple(ordered if maximum_components is None else ordered[:maximum_components])


def dense_two_core_component(
    complex_: RelationComplex,
    *,
    minimum_relation_degree: int = 2,
) -> RelationComponent:
    """Peel vertices until every survivor meets the relation-degree bound."""

    if minimum_relation_degree < 1:
        raise ValueError("minimum_relation_degree must be positive")
    alive = set(range(len(complex_.vertices)))
    changed = True
    while changed:
        degrees = {vertex: 0 for vertex in alive}
        for relation in complex_.ternary_relations:
            if set(relation) <= alive:
                for vertex in relation:
                    degrees[vertex] += 1
        remove = {
            vertex for vertex, degree in degrees.items() if degree < minimum_relation_degree
        }
        changed = bool(remove)
        alive.difference_update(remove)
    return _component(
        complex_, alive, kind=f"relation_{minimum_relation_degree}_core"
    )


def biconnected_components(
    complex_: RelationComplex,
) -> tuple[RelationComponent, ...]:
    """Return exact vertex-biconnected blocks of the primal relation graph."""

    adjacency = [set() for _ in complex_.vertices]
    for relation in complex_.ternary_relations:
        for left in range(3):
            for right in range(left):
                adjacency[relation[left]].add(relation[right])
                adjacency[relation[right]].add(relation[left])
    discovery = [-1] * len(adjacency)
    low = [0] * len(adjacency)
    parent = [-1] * len(adjacency)
    edge_stack: list[tuple[int, int]] = []
    blocks: list[set[int]] = []
    time = 0

    def visit(vertex: int) -> None:
        nonlocal time
        discovery[vertex] = low[vertex] = time
        time += 1
        for neighbour in sorted(adjacency[vertex]):
            if discovery[neighbour] < 0:
                parent[neighbour] = vertex
                edge_stack.append((vertex, neighbour))
                visit(neighbour)
                low[vertex] = min(low[vertex], low[neighbour])
                if low[neighbour] >= discovery[vertex]:
                    block = set()
                    while edge_stack:
                        edge = edge_stack.pop()
                        block.update(edge)
                        if edge == (vertex, neighbour):
                            break
                    if block:
                        blocks.append(block)
            elif neighbour != parent[vertex] and discovery[neighbour] < discovery[vertex]:
                low[vertex] = min(low[vertex], discovery[neighbour])
                edge_stack.append((vertex, neighbour))

    for vertex in range(len(adjacency)):
        if discovery[vertex] >= 0:
            continue
        visit(vertex)
        if edge_stack:
            block = set()
            for edge in edge_stack:
                block.update(edge)
            edge_stack.clear()
            blocks.append(block)
    components = (_component(complex_, block, kind="biconnected") for block in blocks)
    return tuple(
        sorted(
            components,
            key=lambda item: (
                -len(item.vertex_indices),
                -len(item.relation_indices),
                item.vertex_indices,
            ),
        )
    )


def bounded_dense_hyperplanes(
    complex_: RelationComplex,
    *,
    sample_count: int,
    random_seed: int,
    maximum_components: int = 32,
) -> DenseHyperplaneLedger:
    """Sample exact rank-``r-1`` spans and rank them by full-cloud replay.

    Sampling is deterministic and coordinate-rebasing invariant in law.  The
    returned support and relation counts are exact, but exhaustion is only over
    the declared pseudorandom samples.
    """

    ambient_rank = rational_rank(complex_.vertices)
    if ambient_rank < 2 or sample_count < 1 or maximum_components < 1:
        raise ValueError("dense-hyperplane bounds must be positive")
    rng = Random(int(random_seed))
    seen: set[tuple[int, ...]] = set()
    independent = 0
    for _ in range(sample_count):
        indices = rng.sample(range(len(complex_.vertices)), ambient_rank - 1)
        rows = tuple(complex_.vertices[index] for index in indices)
        if rational_rank(rows) != ambient_rank - 1:
            continue
        normals = rational_nullspace(rows)
        if len(normals) != len(complex_.vertices[0]) - ambient_rank + 1:
            continue
        independent += 1
        # The ambient vectors may occupy a proper coordinate subspace.  Use
        # the full normal tuple as the exact rational-space key.
        vertices = tuple(
            index
            for index, vector in enumerate(complex_.vertices)
            if all(
                sum(int(value) * int(normal[column]) for column, value in enumerate(vector))
                == 0
                for normal in normals
            )
        )
        # The sampled rows already have rank r-1, and ``vertices`` is exactly
        # their common rational nullspace inside the supplied ambient cloud.
        seen.add(vertices)
    support_order = sorted(seen, key=lambda vertices: (-len(vertices), vertices))
    cutoff = (
        len(support_order[min(maximum_components, len(support_order)) - 1])
        if support_order
        else 0
    )
    finalists = (
        _component(complex_, vertices, kind="sampled_dense_hyperplane")
        for vertices in support_order
        if len(vertices) >= cutoff
    )
    ordered = sorted(
        finalists,
        key=lambda item: (
            -len(item.vertex_indices),
            -len(item.relation_indices),
            item.vertex_indices,
        ),
    )[:maximum_components]
    return DenseHyperplaneLedger(
        ambient_rank=ambient_rank,
        random_seed=int(random_seed),
        sample_count=int(sample_count),
        independent_sample_count=independent,
        distinct_hyperplane_count=len(seen),
        components=tuple(ordered),
    )


def _height_angle_rms(
    source: RelationComplex,
    target: RelationComplex,
    vertex_map: Sequence[int],
    source_gram: Sequence[Sequence[float]],
    target_gram: Sequence[Sequence[float]],
) -> float:
    pairs = tuple((index, target_index) for index, target_index in enumerate(vertex_map) if target_index >= 0)
    if len(pairs) < 2:
        return 0.0
    source_rows = np.asarray([source.vertices[index] for index, _ in pairs], dtype=float)
    target_rows = np.asarray([target.vertices[index] for _, index in pairs], dtype=float)
    source_form = np.asarray(source_gram, dtype=float)
    target_form = np.asarray(target_gram, dtype=float)
    source_norms = np.einsum("ij,jk,ik->i", source_rows, source_form, source_rows)
    target_norms = np.einsum("ij,jk,ik->i", target_rows, target_form, target_rows)
    source_pairings = source_rows @ source_form @ source_rows.T
    target_pairings = target_rows @ target_form @ target_rows.T
    source_angles = np.abs(source_pairings / np.sqrt(source_norms[:, None] * source_norms[None, :]))
    target_angles = np.abs(target_pairings / np.sqrt(target_norms[:, None] * target_norms[None, :]))
    errors = [target_angles[left, right] - source_angles[left, right] for left in range(len(pairs)) for right in range(left)]
    return float(np.sqrt(np.mean(np.square(errors))))


def merge_component_vertex_maps(
    source: RelationComplex,
    target: RelationComplex,
    left_map: Sequence[int],
    right_map: Sequence[int],
    *,
    held_out_source_indices: Sequence[int] = (),
    source_gram: Sequence[Sequence[float]] | None = None,
    target_gram: Sequence[Sequence[float]] | None = None,
    maximum_height_angle_rms: float | None = None,
    finite_primes: Sequence[int] = (2, 3),
    maximum_sign_components: int = 16,
) -> tuple[ExactComponentMerge, ...]:
    """Merge two partial maps and run exact proper-subspace replay immediately."""

    if len(left_map) != len(source.vertices) or len(right_map) != len(source.vertices):
        raise ValueError("component map length differs from source complex")
    merged = []
    for left, right in zip(left_map, right_map):
        if left >= 0 and right >= 0 and int(left) != int(right):
            return ()
        merged.append(int(left if left >= 0 else right))
    targets = [index for index in merged if index >= 0]
    if len(targets) != len(set(targets)):
        return ()
    source_indices = tuple(index for index, target_index in enumerate(merged) if target_index >= 0)
    left_indices = tuple(index for index, target_index in enumerate(left_map) if target_index >= 0)
    right_indices = tuple(index for index, target_index in enumerate(right_map) if target_index >= 0)
    for prime in finite_primes:
        if modular_rank(tuple(source.vertices[index] for index in source_indices), int(prime)) != modular_rank(tuple(target.vertices[index] for index in targets), int(prime)):
            return ()
    angle_rms = 0.0
    if source_gram is not None or target_gram is not None:
        if source_gram is None or target_gram is None:
            raise ValueError("both height Grams are required for angle pruning")
        angle_rms = _height_angle_rms(source, target, merged, source_gram, target_gram)
        if maximum_height_angle_rms is not None and angle_rms > maximum_height_angle_rms:
            return ()
    left_rank = rational_rank(tuple(source.vertices[index] for index in left_indices))
    right_rank = rational_rank(tuple(source.vertices[index] for index in right_indices))
    union_rank = rational_rank(tuple(source.vertices[index] for index in source_indices))
    added_rank = union_rank - max(left_rank, right_rank)
    held_out = set(map(int, held_out_source_indices))
    answers = []
    for replay in exact_partial_relation_replay(
        source,
        target,
        merged,
        maximum_sign_components=maximum_sign_components,
    ):
        held_out_count = len(held_out & set(replay.replayed_source_vertex_indices))
        answers.append(
            ExactComponentMerge(
                replay=replay,
                left_source_rank=left_rank,
                right_source_rank=right_rank,
                added_rank=added_rank,
                held_out_replayed_ray_count=held_out_count,
                held_out_rays_per_added_rank=(
                    held_out_count / max(1, added_rank)
                ),
                height_angle_rms=angle_rms,
            )
        )
    return tuple(
        sorted(
            answers,
            key=lambda item: (
                -item.held_out_rays_per_added_rank,
                -item.held_out_replayed_ray_count,
                -len(item.replay.replayed_source_vertex_indices),
                -item.replay.replayed_relation_count,
            ),
        )
    )

"""Exact candidate-intersection consensus after a numerical shape prefilter."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .pari import exact_rational_ranks


@dataclass(frozen=True)
class IntersectionConsensusCandidate:
    source_index: int
    shape_value: str
    normalized_shape_extremality: str
    codimension_one_intersection_count: int
    codimension_one_intersection_fraction: str
    combined_score: str

    def to_record(self) -> dict[str, object]:
        return {
            "source_index": self.source_index,
            "shape_value": self.shape_value,
            "normalized_shape_extremality": self.normalized_shape_extremality,
            "codimension_one_intersection_count": self.codimension_one_intersection_count,
            "codimension_one_intersection_fraction": self.codimension_one_intersection_fraction,
            "combined_score": self.combined_score,
        }


@dataclass(frozen=True)
class IntersectionConsensusLedger:
    dimension: int
    population_count: int
    pool_size: int
    exact_pair_count: int
    candidates: tuple[IntersectionConsensusCandidate, ...]

    @property
    def selected(self) -> IntersectionConsensusCandidate:
        return min(
            self.candidates,
            key=lambda item: (
                -float(item.combined_score),
                -float(item.shape_value),
                item.source_index,
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "population_count": self.population_count,
            "pool_size": self.pool_size,
            "exact_pair_count": self.exact_pair_count,
            "selected_source_index": self.selected.source_index,
            "candidates": [candidate.to_record() for candidate in self.candidates],
        }


def exact_intersection_consensus(
    basis_matrices: Sequence[Sequence[Sequence[int]]],
    shape_values: Sequence[float | str],
    *,
    pool_size: int = 64,
    timeout: float = 300.0,
) -> IntersectionConsensusLedger:
    """Rank candidates by shape extremality plus exact codimension-one support.

    Shape values are used only to choose and normalize a bounded pool.  Every
    intersection in that pool is then computed over ``Q``.  The equal-weight
    combined score lies in ``[0,2]`` and contains no withheld truth data.
    """

    matrices = tuple(
        tuple(tuple(map(int, row)) for row in matrix) for matrix in basis_matrices
    )
    values = tuple(float(value) for value in shape_values)
    if not matrices or len(matrices) != len(values):
        raise ValueError("candidate matrices and shape values must agree")
    dimension = len(matrices[0])
    ambient = len(matrices[0][0])
    if dimension < 1 or pool_size < 2:
        raise ValueError("consensus dimension and pool_size must be positive")
    if any(
        len(matrix) != dimension
        or any(len(row) != ambient for row in matrix)
        for matrix in matrices
    ):
        raise ValueError("candidate embedding shapes differ")
    retained_count = min(int(pool_size), len(matrices))
    retained = tuple(
        sorted(range(len(matrices)), key=lambda index: (-values[index], index))[
            :retained_count
        ]
    )
    pairs = []
    pair_indices = []
    for left in range(retained_count):
        for right in range(left):
            pairs.append(matrices[retained[left]] + matrices[retained[right]])
            pair_indices.append((left, right))
    ranks = exact_rational_ranks(pairs, batch_size=128, timeout=timeout)
    counts = [0] * retained_count
    for (left, right), rank in zip(pair_indices, ranks):
        intersection = 2 * dimension - rank
        if intersection >= dimension - 1:
            counts[left] += 1
            counts[right] += 1
    retained_values = [values[index] for index in retained]
    low, high = min(retained_values), max(retained_values)
    denominator = high - low
    candidates = []
    for local, source in enumerate(retained):
        normalized = (
            (values[source] - low) / denominator if denominator > 0 else 1.0
        )
        fraction = counts[local] / (retained_count - 1)
        candidates.append(
            IntersectionConsensusCandidate(
                source_index=source,
                shape_value=f"{values[source]:.17g}",
                normalized_shape_extremality=f"{normalized:.17g}",
                codimension_one_intersection_count=counts[local],
                codimension_one_intersection_fraction=f"{fraction:.17g}",
                combined_score=f"{normalized + fraction:.17g}",
            )
        )
    return IntersectionConsensusLedger(
        dimension=dimension,
        population_count=len(matrices),
        pool_size=retained_count,
        exact_pair_count=len(pairs),
        candidates=tuple(candidates),
    )


@dataclass(frozen=True)
class GraphWalkConsensusCandidate:
    """One candidate's exact graph counts and numerical rank score."""

    source_index: int
    shape_value: str
    shape_percentile: str
    triangle_count: int
    triangle_percentile: str
    length_four_walk_count: int
    length_four_walk_percentile: str
    combined_score: str

    def to_record(self) -> dict[str, object]:
        return {
            "source_index": self.source_index,
            "shape_value": self.shape_value,
            "shape_percentile": self.shape_percentile,
            "triangle_count": self.triangle_count,
            "triangle_percentile": self.triangle_percentile,
            "length_four_walk_count": self.length_four_walk_count,
            "length_four_walk_percentile": self.length_four_walk_percentile,
            "combined_score": self.combined_score,
        }


@dataclass(frozen=True)
class GraphWalkConsensusLedger:
    dimension: int
    population_count: int
    pool_size: int
    exact_pair_count: int
    codimension_one_edge_count: int
    shape_gap: str
    shape_gap_threshold: str
    graph_weight: str
    selector_mode: str
    candidates: tuple[GraphWalkConsensusCandidate, ...]

    @property
    def selected(self) -> GraphWalkConsensusCandidate:
        source = rescore_graph_walk_consensus(
            self,
            shape_gap_threshold=float(self.shape_gap_threshold),
            graph_weight=float(self.graph_weight),
        )
        return next(item for item in self.candidates if item.source_index == source)

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "population_count": self.population_count,
            "pool_size": self.pool_size,
            "exact_pair_count": self.exact_pair_count,
            "codimension_one_edge_count": self.codimension_one_edge_count,
            "shape_gap": self.shape_gap,
            "shape_gap_threshold": self.shape_gap_threshold,
            "graph_weight": self.graph_weight,
            "selector_mode": self.selector_mode,
            "selected_source_index": self.selected.source_index,
            "candidates": [candidate.to_record() for candidate in self.candidates],
        }


def _average_rank_percentiles(values: Sequence[int | float]) -> np.ndarray:
    """Return tie-averaged ascending percentiles in ``[0,1]``."""

    values = np.asarray(values)
    if not len(values):
        return np.empty(0, dtype=float)
    if len(values) == 1:
        return np.ones(1, dtype=float)
    answer = np.empty(len(values), dtype=float)
    for value in np.unique(values):
        indices = np.flatnonzero(values == value)
        answer[indices] = (
            int(np.sum(values < value)) + (len(indices) - 1) / 2
        ) / (len(values) - 1)
    return answer


def rescore_graph_walk_consensus(
    ledger: GraphWalkConsensusLedger,
    *,
    shape_gap_threshold: float,
    graph_weight: float,
) -> int:
    """Select from pinned graph counts under a declared score pair."""

    if shape_gap_threshold < 0 or graph_weight <= 0:
        raise ValueError("graph-walk score bounds must be nonnegative")
    if float(ledger.shape_gap) >= shape_gap_threshold:
        return max(
            ledger.candidates,
            key=lambda item: (float(item.shape_value), -item.source_index),
        ).source_index
    return max(
        ledger.candidates,
        key=lambda item: (
            graph_weight
            * (
                float(item.triangle_percentile)
                + float(item.length_four_walk_percentile)
            )
            + float(item.shape_percentile),
            float(item.shape_value),
            -item.source_index,
        ),
    ).source_index


def exact_graph_walk_consensus(
    basis_matrices: Sequence[Sequence[Sequence[int]]],
    shape_values: Sequence[float | str],
    *,
    pool_size: int = 64,
    shape_gap_threshold: float = 0.005,
    graph_weight: float = 1.5,
    timeout: float = 300.0,
) -> GraphWalkConsensusLedger:
    """Combine Hermite extremality with exact codimension-one graph walks.

    The graph has one vertex per shape-prefiltered primitive candidate and an
    edge when two ``k``-spaces meet in dimension at least ``k-1``.  Triangle
    counts and the entries of ``A^4 1`` are exact integers.  Only shape values,
    rank percentiles, and their final weighted score are numerical/search
    data.
    """

    matrices = tuple(
        tuple(tuple(map(int, row)) for row in matrix) for matrix in basis_matrices
    )
    values = tuple(float(value) for value in shape_values)
    if not matrices or len(matrices) != len(values):
        raise ValueError("candidate matrices and shape values must agree")
    if pool_size < 2 or shape_gap_threshold < 0 or graph_weight <= 0:
        raise ValueError("invalid graph-walk consensus bounds")
    dimension = len(matrices[0])
    ambient = len(matrices[0][0])
    if any(
        len(matrix) != dimension
        or any(len(row) != ambient for row in matrix)
        for matrix in matrices
    ):
        raise ValueError("candidate embedding shapes differ")
    retained_count = min(int(pool_size), len(matrices))
    retained = tuple(
        sorted(range(len(matrices)), key=lambda index: (-values[index], index))[
            :retained_count
        ]
    )
    pairs = []
    pair_indices = []
    for left in range(retained_count):
        for right in range(left):
            pairs.append(matrices[retained[left]] + matrices[retained[right]])
            pair_indices.append((left, right))
    ranks = exact_rational_ranks(pairs, batch_size=128, timeout=timeout)
    adjacency = np.zeros((retained_count, retained_count), dtype=np.int64)
    for (left, right), rank in zip(pair_indices, ranks):
        if 2 * dimension - rank >= dimension - 1:
            adjacency[left, right] = adjacency[right, left] = 1
    squared = adjacency @ adjacency
    triangles = np.diag(squared @ adjacency) // 2
    walks = np.ones(retained_count, dtype=np.int64)
    for _step in range(4):
        walks = adjacency @ walks
    if np.any(walks < 0):
        raise ArithmeticError("length-four graph walk count overflowed int64")
    retained_values = np.asarray([values[index] for index in retained])
    shape_percentiles = _average_rank_percentiles(retained_values)
    triangle_percentiles = _average_rank_percentiles(triangles)
    walk_percentiles = _average_rank_percentiles(walks)
    candidates = tuple(
        GraphWalkConsensusCandidate(
            source_index=source,
            shape_value=f"{retained_values[local]:.17g}",
            shape_percentile=f"{shape_percentiles[local]:.17g}",
            triangle_count=int(triangles[local]),
            triangle_percentile=f"{triangle_percentiles[local]:.17g}",
            length_four_walk_count=int(walks[local]),
            length_four_walk_percentile=f"{walk_percentiles[local]:.17g}",
            combined_score=f"{graph_weight * (triangle_percentiles[local] + walk_percentiles[local]) + shape_percentiles[local]:.17g}",
        )
        for local, source in enumerate(retained)
    )
    gap = retained_values[0] - retained_values[1]
    ledger = GraphWalkConsensusLedger(
        dimension=dimension,
        population_count=len(matrices),
        pool_size=retained_count,
        exact_pair_count=len(pairs),
        codimension_one_edge_count=int(np.sum(adjacency) // 2),
        shape_gap=f"{gap:.17g}",
        shape_gap_threshold=f"{float(shape_gap_threshold):.17g}",
        graph_weight=f"{float(graph_weight):.17g}",
        selector_mode=(
            "separated_shape_extremum"
            if gap >= shape_gap_threshold
            else "exact_graph_walk_consensus"
        ),
        candidates=candidates,
    )
    # Force selection now so malformed score records fail during construction.
    ledger.selected
    return ledger

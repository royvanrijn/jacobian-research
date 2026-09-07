"""Compact basis-independent fingerprints for candidate ray subspaces.

These summaries are designed for joint candidate selection across unequal
finite clouds.  Every exact combinatorial count is computed on the declared
cloud; normalized height and arithmetic profiles are numerical search
features, not identities or isometry certificates.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log1p
from typing import Mapping, Sequence

import numpy as np

from .integer import canonical_unoriented
from .relations import RelationComplex


@dataclass(frozen=True)
class CandidateRelationFingerprint:
    """Fixed-width candidate summary with exact audit counts."""

    dimension: int
    ray_count: int
    ternary_relation_count: int
    scaled_relation_count: int
    integral_ray_count: int
    feature_names: tuple[str, ...]
    feature_values: tuple[str, ...]

    def vector(self) -> np.ndarray:
        return np.asarray(self.feature_values, dtype=float)

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "ray_count": self.ray_count,
            "ternary_relation_count": self.ternary_relation_count,
            "scaled_relation_count": self.scaled_relation_count,
            "integral_ray_count": self.integral_ray_count,
            "feature_names": list(self.feature_names),
            "feature_values": list(self.feature_values),
        }


@dataclass(frozen=True)
class JointCandidateScore:
    """Nearest-neighbour consensus score for one fibre candidate."""

    fibre_index: int
    candidate_index: int
    other_fibre_neighbours: tuple[tuple[int, int, str], ...]
    mean_distance: str
    maximum_distance: str
    mutual_neighbour_count: int

    def to_record(self) -> dict[str, object]:
        return {
            "fibre_index": self.fibre_index,
            "candidate_index": self.candidate_index,
            "other_fibre_neighbours": [
                {
                    "fibre_index": fibre,
                    "candidate_index": candidate,
                    "distance": distance,
                }
                for fibre, candidate, distance in self.other_fibre_neighbours
            ],
            "mean_distance": self.mean_distance,
            "maximum_distance": self.maximum_distance,
            "mutual_neighbour_count": self.mutual_neighbour_count,
        }


def _quantiles(values, count: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.zeros(count, dtype=float)
    return np.quantile(values, np.linspace(0.0, 1.0, count))


def _projective_class(vector: Sequence[int], prime: int) -> tuple[int, ...]:
    values = tuple(int(value) % prime for value in vector)
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return ()
    inverse = pow(pivot, -1, prime)
    return tuple(value * inverse % prime for value in values)


def candidate_relation_fingerprint(
    vectors: Sequence[Sequence[int]],
    heights: Sequence[float | str],
    arithmetic: Sequence[Mapping[str, object]],
    inlier_indices: Sequence[int],
    complex_: RelationComplex,
    *,
    dimension: int,
    quantiles: int = 16,
    projective_multiplicities: int = 16,
    finite_primes: Sequence[int] = (2, 3),
) -> CandidateRelationFingerprint:
    """Summarize one candidate without reference to an ambient basis label.

    Coordinate changes are allowed to act simultaneously on every ray.  The
    induced relation complex and projective-class multiplicity partitions are
    invariant under unimodular ambient changes.  Height profiles are centered
    logarithmically, so a common specialization scale is forgotten.
    """

    if not vectors or len(vectors) != len(heights) or len(vectors) != len(arithmetic):
        raise ValueError("vectors, heights, and arithmetic populations must agree")
    if quantiles < 2 or projective_multiplicities < 1:
        raise ValueError("fingerprint profile bounds must be positive")
    input_by_ray = {
        canonical_unoriented(vector): index for index, vector in enumerate(vectors)
    }
    if len(input_by_ray) != len(vectors) or set(input_by_ray) != set(complex_.vertices):
        raise ValueError("vector population and relation complex differ")
    retained_input = tuple(sorted(set(map(int, inlier_indices))))
    if not retained_input or min(retained_input) < 0 or max(retained_input) >= len(vectors):
        raise ValueError("candidate inlier indices are empty or out of range")
    retained_input_set = set(retained_input)
    retained_complex = np.asarray(
        [input_by_ray[vertex] in retained_input_set for vertex in complex_.vertices],
        dtype=bool,
    )
    ternary = np.asarray(complex_.ternary_relations, dtype=np.int64)
    scaled = np.asarray(
        [edge[:3] for edge in complex_.scaled_relations], dtype=np.int64
    )
    kept_ternary = (
        ternary[np.all(retained_complex[ternary], axis=1)]
        if len(ternary)
        else np.empty((0, 3), dtype=np.int64)
    )
    kept_scaled = (
        scaled[np.all(retained_complex[scaled], axis=1)]
        if len(scaled)
        else np.empty((0, 3), dtype=np.int64)
    )
    degrees = np.zeros(len(complex_.vertices), dtype=float)
    if len(kept_ternary):
        degrees += np.bincount(kept_ternary.ravel(), minlength=len(degrees))
    scaled_degrees = np.zeros(len(complex_.vertices), dtype=float)
    if len(kept_scaled):
        scaled_degrees += np.bincount(kept_scaled.ravel(), minlength=len(degrees))
    retained_vertices = np.flatnonzero(retained_complex)

    height_values = np.asarray([float(value) for value in heights], dtype=float)
    if np.min(height_values) <= 0 or not np.all(np.isfinite(height_values)):
        raise ValueError("candidate heights must be finite and positive")
    retained_heights = np.asarray(
        [height_values[input_by_ray[complex_.vertices[index]]] for index in retained_vertices]
    )
    ray_logs = np.log(retained_heights)
    ray_logs -= float(np.mean(ray_logs))
    if len(kept_ternary):
        edge_logs = np.asarray(
            [
                sorted(
                    np.log(height_values[input_by_ray[complex_.vertices[index]]])
                    for index in edge
                )
                for edge in kept_ternary
            ],
            dtype=float,
        )
        edge_logs -= np.mean(edge_logs, axis=1)[:, None]
    else:
        edge_logs = np.zeros((1, 3), dtype=float)

    total_bits = np.asarray(
        [float(item.get("total_bits", 0)) for item in arithmetic], dtype=float
    )
    ambient_bit_center = float(np.median(np.log1p(total_bits)))
    retained_bits = np.asarray(
        [
            log1p(total_bits[input_by_ray[complex_.vertices[index]]])
            - ambient_bit_center
            for index in retained_vertices
        ]
    )
    integral_count = sum(
        bool(arithmetic[index].get("integral", False)) for index in retained_input
    )

    names = []
    values = []

    def add(name: str, value: float) -> None:
        names.append(name)
        values.append(float(value))

    add("log_ray_count", np.log(len(retained_vertices)))
    add("ternary_relations_per_ray", len(kept_ternary) / len(retained_vertices))
    add("scaled_relations_per_ray", len(kept_scaled) / len(retained_vertices))
    add("integral_ray_fraction", integral_count / len(retained_vertices))
    for index, value in enumerate(_quantiles(ray_logs, quantiles)):
        add(f"ray_log_height_q{index}", value)
    for column in range(3):
        for index, value in enumerate(_quantiles(edge_logs[:, column], quantiles)):
            add(f"edge_centered_log_height_{column}_q{index}", value)
    retained_degree_logs = np.log1p(degrees[retained_vertices])
    retained_degree_logs -= float(np.mean(retained_degree_logs))
    for index, value in enumerate(_quantiles(retained_degree_logs, quantiles)):
        add(f"ternary_log_degree_q{index}", value)
    retained_scaled_logs = np.log1p(scaled_degrees[retained_vertices])
    retained_scaled_logs -= float(np.mean(retained_scaled_logs))
    for index, value in enumerate(_quantiles(retained_scaled_logs, quantiles)):
        add(f"scaled_log_degree_q{index}", value)
    for index, value in enumerate(_quantiles(retained_bits, quantiles)):
        add(f"relative_log_complexity_q{index}", value)

    retained_vectors = [complex_.vertices[index] for index in retained_vertices]
    for prime in finite_primes:
        counts = {}
        for vector in retained_vectors:
            key = _projective_class(vector, int(prime))
            counts[key] = counts.get(key, 0) + 1
        multiplicities = sorted(counts.values(), reverse=True)
        normalized = [value / len(retained_vectors) for value in multiplicities]
        normalized.extend([0.0] * (projective_multiplicities - len(normalized)))
        add(f"mod{prime}_projective_class_count_per_ray", len(counts) / len(retained_vectors))
        for index, value in enumerate(normalized[:projective_multiplicities]):
            add(f"mod{prime}_projective_multiplicity_{index}", value)

    return CandidateRelationFingerprint(
        dimension=int(dimension),
        ray_count=len(retained_vertices),
        ternary_relation_count=len(kept_ternary),
        scaled_relation_count=len(kept_scaled),
        integral_ray_count=int(integral_count),
        feature_names=tuple(names),
        feature_values=tuple(f"{value:.17g}" for value in values),
    )


def robust_standardize_fingerprint_families(
    families: Sequence[Sequence[CandidateRelationFingerprint]],
) -> tuple[np.ndarray, ...]:
    """Apply one coordinatewise median/MAD scale to several fibre ledgers."""

    if len(families) < 2 or any(not family for family in families):
        raise ValueError("joint standardization needs at least two nonempty fibres")
    names = families[0][0].feature_names
    dimension = families[0][0].dimension
    if any(
        item.feature_names != names or item.dimension != dimension
        for family in families
        for item in family
    ):
        raise ValueError("fingerprint families have incompatible schemas")
    matrices = tuple(np.asarray([item.vector() for item in family]) for family in families)
    population = np.vstack(matrices)
    median = np.median(population, axis=0)
    absolute = np.abs(population - median)
    scale = np.median(absolute, axis=0)
    fallback = np.std(population, axis=0)
    scale = np.where(scale > 1e-12, scale, np.where(fallback > 1e-12, fallback, 1.0))
    return tuple((matrix - median) / scale for matrix in matrices)


def joint_nearest_candidate_scores(
    families: Sequence[Sequence[CandidateRelationFingerprint]],
    *,
    chunk_size: int = 256,
) -> tuple[tuple[JointCandidateScore, ...], ...]:
    """Score every candidate by nearest neighbours in all other fibres.

    Distances use the robustly standardized compact profiles.  Mutual-nearest
    counts are reported separately rather than folded into the numerical
    distance, so calibration can detect whether reciprocity helps.
    """

    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    matrices = robust_standardize_fingerprint_families(families)
    nearest = {}
    for left, left_matrix in enumerate(matrices):
        for right, right_matrix in enumerate(matrices):
            if left == right:
                continue
            indices = np.empty(len(left_matrix), dtype=np.int64)
            distances = np.empty(len(left_matrix), dtype=float)
            for start in range(0, len(left_matrix), chunk_size):
                stop = min(len(left_matrix), start + chunk_size)
                differences = left_matrix[start:stop, None, :] - right_matrix[None, :, :]
                squared = np.mean(differences * differences, axis=2)
                choice = np.argmin(squared, axis=1)
                indices[start:stop] = choice
                distances[start:stop] = np.sqrt(
                    squared[np.arange(stop - start), choice]
                )
            nearest[(left, right)] = (indices, distances)
    answer = []
    for fibre, family in enumerate(families):
        scores = []
        for candidate in range(len(family)):
            neighbours = []
            mutual = 0
            for other in range(len(families)):
                if other == fibre:
                    continue
                indices, distances = nearest[(fibre, other)]
                other_candidate = int(indices[candidate])
                neighbours.append(
                    (other, other_candidate, f"{float(distances[candidate]):.17g}")
                )
                reverse_indices, _reverse_distances = nearest[(other, fibre)]
                mutual += int(reverse_indices[other_candidate] == candidate)
            numeric = [float(item[2]) for item in neighbours]
            scores.append(
                JointCandidateScore(
                    fibre_index=fibre,
                    candidate_index=candidate,
                    other_fibre_neighbours=tuple(neighbours),
                    mean_distance=f"{float(np.mean(numeric)):.17g}",
                    maximum_distance=f"{float(np.max(numeric)):.17g}",
                    mutual_neighbour_count=mutual,
                )
            )
        answer.append(tuple(scores))
    return tuple(answer)

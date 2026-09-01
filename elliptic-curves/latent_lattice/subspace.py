"""Bounded exact/heuristic subspace operations for short-vector clouds."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import log
from random import Random
from typing import Any, Sequence

import numpy as np

from .codes import finite_rarity_scores
from .finite import FiniteQuotientBlock
from .integer import rational_nullspace, rational_rank
from .local import ComponentBlock
from .pari import ShortVectorRecord, gp_matrix, primitive_column_closure, run_gp
from .relations import RelationComplex, build_relation_complex


def modular_rank(rows: Sequence[Sequence[int]], prime: int = 1_000_003) -> int:
    matrix = np.array(
        [[int(value) % prime for value in row] for row in rows], dtype=np.int64
    )
    if matrix.size == 0:
        return 0
    rank = 0
    for column in range(matrix.shape[1]):
        pivot = next(
            (index for index in range(rank, matrix.shape[0]) if matrix[index, column]),
            None,
        )
        if pivot is None:
            continue
        matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inverse = pow(int(matrix[rank, column]), -1, prime)
        matrix[rank] = (matrix[rank] * inverse) % prime
        for index in range(matrix.shape[0]):
            if index != rank and matrix[index, column]:
                matrix[index] = (
                    matrix[index] - matrix[index, column] * matrix[rank]
                ) % prime
        rank += 1
        if rank == matrix.shape[0]:
            break
    return rank


def independent_row_basis(
    rows: Sequence[Sequence[int]], prime: int = 1_000_003
) -> tuple[tuple[int, ...], ...]:
    """Select original rows forming a basis modulo ``prime`` in one pass.

    The selected integer rows are automatically independent over ``Q``.  As
    with :func:`modular_rank`, a deficient modular rank can occur only when the
    chosen prime divides a relevant minor; final retained spaces are therefore
    checked by exact rational rank and saturation.
    """

    if not rows:
        return ()
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("row widths differ")
    echelon: list[tuple[int, list[int]]] = []
    selected: list[tuple[int, ...]] = []
    for raw_row in rows:
        row = [int(value) % prime for value in raw_row]
        for pivot, pivot_row in echelon:
            coefficient = row[pivot]
            if coefficient:
                row = [
                    (value - coefficient * basis_value) % prime
                    for value, basis_value in zip(row, pivot_row)
                ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], -1, prime)
        row = [(value * inverse) % prime for value in row]
        echelon.append((pivot, row))
        selected.append(tuple(map(int, raw_row)))
        if len(selected) == width:
            break
    return tuple(selected)


def modular_row_space_key(
    rows: Sequence[Sequence[int]], prime: int = 1_000_003
) -> tuple[tuple[int, ...], ...]:
    """Return a canonical reduced-row-echelon key over one finite field.

    This is a fast proposal/deduplication key, not an exact rational-space
    certificate.  A reported survivor must be replayed with exact rank,
    primitive closure, and exact span tests.
    """

    if not rows:
        return ()
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("row widths differ")
    matrix = [
        [int(value) % prime for value in row]
        for row in rows
        if any(int(value) % prime for value in row)
    ]
    rank = 0
    pivots: list[int] = []
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank:
                continue
            coefficient = matrix[index][column]
            if coefficient:
                matrix[index] = [
                    (value - coefficient * pivot_value) % prime
                    for value, pivot_value in zip(matrix[index], matrix[rank])
                ]
        pivots.append(column)
        rank += 1
        if rank == len(matrix):
            break
    return tuple(tuple(row) for row in matrix[:rank])


def modular_right_kernel_basis(
    rows: Sequence[Sequence[int]], prime: int = 1_000_003
) -> tuple[tuple[int, ...], ...]:
    """Return a row basis for the right kernel over ``F_prime``."""

    if not rows:
        raise ValueError("kernel input must be nonempty")
    rref = modular_row_space_key(rows, prime)
    width = len(rows[0])
    pivots = [next(index for index, value in enumerate(row) if value) for row in rref]
    free = [index for index in range(width) if index not in set(pivots)]
    answer = []
    for free_column in free:
        vector = [0] * width
        vector[free_column] = 1
        for row, pivot in zip(rref, pivots):
            vector[pivot] = (-row[free_column]) % prime
        answer.append(tuple(vector))
    return tuple(answer)


def exact_row_space_intersection(
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
    *,
    saturate: bool = True,
) -> tuple[tuple[int, ...], ...]:
    """Return an integer basis for two rational row-space intersection.

    ``saturate=False`` avoids the external exact-closure call and is useful for
    bulk scoring; it still returns the exact rational intersection, but not
    necessarily a primitive ambient embedding.
    """

    if not left or not right:
        return ()
    width = len(left[0])
    if any(len(row) != width for row in (*left, *right)):
        raise ValueError("row widths differ")
    if rational_rank(left) != len(left) or rational_rank(right) != len(right):
        raise ValueError("intersection inputs must be row bases")
    coefficient_rows = tuple(
        tuple(int(left[row][column]) for row in range(len(left)))
        + tuple(-int(right[row][column]) for row in range(len(right)))
        for column in range(width)
    )
    kernel_rows = rational_nullspace(coefficient_rows)
    nullity = len(kernel_rows)
    if nullity == 0:
        return ()
    intersection = tuple(
        tuple(
            sum(
                int(kernel_rows[basis_column][row]) * int(left[row][coordinate])
                for row in range(len(left))
            )
            for coordinate in range(width)
        )
        for basis_column in range(nullity)
    )
    if rational_rank(intersection) != nullity:
        raise ArithmeticError("intersection kernel lost rank")
    return primitive_span_basis(intersection) if saturate else intersection


def integer_right_kernel(rows: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    """Return rows of the matrix whose columns are ``ker_Z(rows)``."""

    if not rows:
        raise ValueError("kernel input must be nonempty")
    width = len(rows[0])
    program = f"""
B={gp_matrix(rows)};K=matkerint(B);print("SIZE|",matsize(K)[1],"|",matsize(K)[2]);
print("BEGIN");for(i=1,matsize(K)[1],for(j=1,matsize(K)[2],if(j>1,print1("|"));print1(K[i,j]));print());print("END");
"""
    lines = run_gp(program)
    height, nullity = map(
        int,
        next(line.split("|")[1:] for line in lines if line.startswith("SIZE|")),
    )
    if height != width:
        raise ArithmeticError("integer kernel has the wrong ambient dimension")
    start = lines.index("BEGIN") + 1
    stop = lines.index("END", start)
    if nullity == 0:
        return tuple(() for _ in range(width))
    answer = tuple(
        tuple(int(value) for value in line.split("|")) for line in lines[start:stop]
    )
    if len(answer) != width or any(len(row) != nullity for row in answer):
        raise ArithmeticError("integer kernel output has the wrong shape")
    return answer


def exact_span_mask(
    vectors: Sequence[Sequence[int]], basis_rows: Sequence[Sequence[int]]
) -> np.ndarray:
    if not basis_rows:
        return np.zeros(len(vectors), dtype=bool)
    if len(basis_rows) == len(basis_rows[0]) and rational_rank(basis_rows) == len(basis_rows):
        return np.ones(len(vectors), dtype=bool)
    kernel = rational_nullspace(basis_rows)
    if not kernel:
        return np.ones(len(vectors), dtype=bool)
    vectors_object = np.asarray(vectors, dtype=object)
    products = vectors_object @ np.asarray(kernel, dtype=object).T
    return np.all(products == 0, axis=1)


def cross_bound_intersection_proposals(
    records: Sequence[ShortVectorRecord],
    left_enclosures: Sequence[GrowthProposal],
    right_enclosures: Sequence[GrowthProposal],
    *,
    target_dimension: int,
    left_count: int = 200,
    right_count: int = 200,
) -> tuple[GrowthProposal, ...]:
    """Intersect independently generated enclosures from two finite clouds.

    The two ledgers must use the same displayed ambient basis.  Modular rank
    is only a rejection filter; every retained intersection is recomputed over
    ``Q``, primitively saturated, and replayed against ``records`` exactly.
    """

    if not records:
        raise ValueError("intersection proposals need an evaluation cloud")
    vectors = tuple(record.coordinates for record in records)
    integral = np.asarray(
        [bool(record.arithmetic.get("integral", False)) for record in records],
        dtype=bool,
    )
    features = arithmetic_feature_flags(records)
    retained: dict[tuple[tuple[int, ...], ...], GrowthProposal] = {}
    for left in left_enclosures[: int(left_count)]:
        for right in right_enclosures[: int(right_count)]:
            expected = left.dimension + right.dimension - target_dimension
            if modular_rank(left.basis_rows + right.basis_rows) != expected:
                continue
            basis = exact_row_space_intersection(left.basis_rows, right.basis_rows)
            if len(basis) != target_dimension:
                continue
            key = modular_row_space_key(basis)
            if key in retained:
                continue
            mask = exact_span_mask(vectors, basis)
            support = int(mask.sum())
            arithmetic_llr = arithmetic_enrichment_llr(mask, features)
            retained[key] = GrowthProposal(
                dimension=target_dimension,
                basis_rows=basis,
                inlier_indices=tuple(map(int, np.flatnonzero(mask))),
                support=support,
                integral_support=int(np.sum(mask & integral)),
                arithmetic_llr=arithmetic_llr,
                search_score=arithmetic_llr,
            )
    proposals = list(retained.values())
    proposals.sort(
        key=lambda proposal: (
            proposal.search_score,
            proposal.support,
            proposal.basis_rows,
        ),
        reverse=True,
    )
    return tuple(proposals)


def numerical_span_mask(vectors: np.ndarray, basis_rows: np.ndarray) -> np.ndarray:
    rank = len(basis_rows)
    _u, singular, vh = np.linalg.svd(basis_rows.astype(float), full_matrices=True)
    if len(singular) < rank or singular[-1] < 1e-8:
        return np.zeros(len(vectors), dtype=bool)
    if rank == vectors.shape[1]:
        return np.ones(len(vectors), dtype=bool)
    null = vh[rank:].T
    scale = np.maximum(1.0, np.linalg.norm(vectors, axis=1))
    residual = np.max(np.abs(vectors @ null), axis=1) / scale
    return residual < 1e-8


def primitive_span_basis(basis_rows: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    """Return rows of a primitive basis for the same rational row space."""

    columns = tuple(zip(*basis_rows))
    closure_rows = primitive_column_closure(columns)
    return tuple(tuple(column) for column in zip(*closure_rows))


def _bernoulli_log_likelihood(successes: int, total: int) -> float:
    if total == 0:
        return 0.0
    probability = min(1.0 - 1e-12, max(1e-12, successes / total))
    return successes * log(probability) + (total - successes) * log(1.0 - probability)


def integrality_llr(inlier: np.ndarray, integral: np.ndarray) -> float:
    inside_total = int(inlier.sum())
    outside_total = len(inlier) - inside_total
    inside_success = int(np.sum(inlier & integral))
    outside_success = int(np.sum(~inlier & integral))
    split = _bernoulli_log_likelihood(inside_success, inside_total) + _bernoulli_log_likelihood(
        outside_success, outside_total
    )
    null = _bernoulli_log_likelihood(int(integral.sum()), len(integral))
    return max(0.0, 2.0 * (split - null))


@dataclass(frozen=True)
class SubspaceCandidate:
    dimension: int
    primitive_basis_rows: tuple[tuple[int, ...], ...]
    inlier_indices: tuple[int, ...]
    support: int
    integral_support: int
    integral_rate: float
    outside_integral_rate: float
    integrality_llr: float
    relation_mass: int
    search_score: float
    arithmetic_llr: float = 0.0

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "primitive_basis_rows": [list(row) for row in self.primitive_basis_rows],
            "inlier_indices": list(self.inlier_indices),
            "support": self.support,
            "integral_support": self.integral_support,
            "integral_rate": f"{self.integral_rate:.17g}",
            "outside_integral_rate": f"{self.outside_integral_rate:.17g}",
            "integrality_likelihood_ratio": f"{self.integrality_llr:.17g}",
            "relation_mass": self.relation_mass,
            "search_score": f"{self.search_score:.17g}",
            "arithmetic_likelihood_ratio": f"{self.arithmetic_llr:.17g}",
        }


@dataclass(frozen=True)
class GrowthProposal:
    """An exact integral basis proposed by independent relation growth."""

    dimension: int
    basis_rows: tuple[tuple[int, ...], ...]
    inlier_indices: tuple[int, ...]
    support: int
    integral_support: int
    arithmetic_llr: float
    search_score: float

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": self.dimension,
            "basis_rows": [list(row) for row in self.basis_rows],
            "inlier_indices": list(self.inlier_indices),
            "support": self.support,
            "integral_support": self.integral_support,
            "arithmetic_likelihood_ratio": f"{self.arithmetic_llr:.17g}",
            "search_score": f"{self.search_score:.17g}",
        }


@dataclass(frozen=True)
class RecombinedSearchLedger:
    """Bounded high-recall proposal ledger for one chosen dimension."""

    dimension: int
    direct_proposals: tuple[GrowthProposal, ...]
    enclosure_proposals: tuple[GrowthProposal, ...]
    refined_proposals: tuple[GrowthProposal, ...]
    bounds: dict[str, int | float]

    def summary_record(self) -> dict[str, Any]:
        return {
            "dimension": self.dimension,
            "direct_proposal_count": len(self.direct_proposals),
            "enclosure_proposal_count": len(self.enclosure_proposals),
            "refined_proposal_count": len(self.refined_proposals),
            "bounds": dict(self.bounds),
        }


def arithmetic_feature_flags(
    records: Sequence[ShortVectorRecord],
) -> tuple[tuple[str, np.ndarray], ...]:
    """Return deterministic low-complexity indicator features.

    Thresholds are empirical quartiles of the complete retained cloud, so the
    construction is invariant under record ordering and uses no withheld
    subgroup labels.
    """

    if not records:
        return ()
    features: list[tuple[str, np.ndarray]] = [
        (
            "integral",
            np.asarray(
                [bool(record.arithmetic.get("integral", False)) for record in records],
                dtype=bool,
            ),
        )
    ]
    for key in (
        "total_bits",
        "x_numerator_bits",
        "x_denominator_bits",
        "y_numerator_bits",
        "y_denominator_bits",
    ):
        values = np.asarray(
            [int(record.arithmetic.get(key, 0)) for record in records], dtype=float
        )
        for quantile in (0.25, 0.5, 0.75):
            threshold = float(np.quantile(values, quantile))
            features.append((f"{key}<=q{int(100 * quantile)}", values <= threshold))
    # Identical flags (notably integrality and denominator q25/q50 on many
    # integral models) carry only one likelihood contribution.
    unique: dict[bytes, tuple[str, np.ndarray]] = {}
    for label, flag in features:
        unique.setdefault(flag.tobytes(), (label, flag))
    return tuple(unique.values())


def arithmetic_enrichment_llr(
    inlier: np.ndarray,
    features: Sequence[tuple[str, np.ndarray]],
) -> float:
    return sum(integrality_llr(inlier, flag) for _label, flag in features)


def independent_relation_growth_proposals(
    records: Sequence[ShortVectorRecord],
    complex_: RelationComplex,
    *,
    dimension: int,
    seed_edges: int = 3_000,
    maximum_proposals: int | None = None,
    priority_mode: str = "arithmetic",
    seed_strategy: str = "top",
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> tuple[GrowthProposal, ...]:
    """Grow one deterministic path from each high-arithmetic relation seed.

    Unlike a global beam, independent paths are never discarded because a
    different early seed has a better partial score.  This preserves candidate
    diversity for subsequent joint-fibre matching.  The result is bounded by
    ``seed_edges`` and is a proposal ledger, not an exhaustive subspace list.
    """

    if not 2 <= dimension <= len(records[0].coordinates):
        raise ValueError("dimension is outside the ambient range")
    if priority_mode not in {"arithmetic", "relations", "finite"}:
        raise ValueError("priority_mode must be 'arithmetic', 'relations', or 'finite'")
    if seed_strategy not in {"top", "stratified"}:
        raise ValueError("seed_strategy must be 'top' or 'stratified'")
    record_by_vector = {record.coordinates: record for record in records}
    if set(record_by_vector) != set(complex_.vertices):
        raise ValueError("records and relation-complex vertices differ")
    ordered = [record_by_vector[vector] for vector in complex_.vertices]
    input_index = {
        record.coordinates: index for index, record in enumerate(records)
    }
    vectors = np.asarray([record.coordinates for record in ordered], dtype=np.int64)
    edges = np.asarray(complex_.ternary_relations, dtype=np.int64)
    if not len(edges):
        return ()
    features = arithmetic_feature_flags(ordered)
    favorable = np.sum(
        np.asarray([flag for _label, flag in features], dtype=float), axis=0
    )
    integral = np.asarray(
        [bool(record.arithmetic.get("integral", False)) for record in ordered],
        dtype=bool,
    )
    total_bits = np.asarray(
        [int(record.arithmetic.get("total_bits", 0)) for record in ordered],
        dtype=float,
    )
    edge_priority = np.sum(
        5.0 * integral[edges] + favorable[edges] - 0.02 * total_bits[edges],
        axis=1,
    )
    relation_priority = np.log1p(
        np.asarray(complex_.additive_degrees, dtype=float)
        + np.asarray(complex_.divisibility_degrees, dtype=float)
    )
    finite_priority = np.asarray(
        finite_rarity_scores(
            [record.coordinates for record in ordered],
            finite_blocks=finite_blocks,
            component_blocks=component_blocks,
        ),
        dtype=float,
    )
    if priority_mode == "finite" and not (finite_blocks or component_blocks):
        raise ValueError("finite priority needs at least one finite code block")
    if priority_mode == "relations":
        edge_priority = np.sum(relation_priority[edges], axis=1)
    elif priority_mode == "finite":
        edge_priority = np.sum(
            finite_priority[edges] + 0.05 * relation_priority[edges], axis=1
        )
    ranked_edges = np.argsort(-edge_priority, kind="stable")
    retained_seed_count = min(seed_edges, len(edges))
    if seed_strategy == "top" or retained_seed_count == len(edges):
        seed_order = ranked_edges[:retained_seed_count]
    else:
        positions = np.linspace(
            0, len(ranked_edges) - 1, retained_seed_count, dtype=int
        )
        seed_order = ranked_edges[positions]
    proposals: list[GrowthProposal] = []
    seen: set[bytes] = set()
    for edge_index in seed_order:
        basis = vectors[edges[edge_index, :2]].copy()
        for current_dimension in range(3, dimension + 1):
            mask = numerical_span_mask(vectors, basis)
            occupancy = np.sum(mask[edges], axis=1)
            crossing = edges[occupancy == 1]
            outside = crossing[~mask[crossing]]
            counts = np.bincount(outside, minlength=len(vectors))
            if priority_mode == "arithmetic":
                scores = (
                    counts * (1.0 + 0.12 * favorable + 0.4 * integral)
                    - 0.005 * total_bits
                )
            elif priority_mode == "relations":
                scores = counts * (1.0 + 0.05 * relation_priority)
            else:
                scores = counts * (
                    1.0 + 0.05 * relation_priority + 0.05 * finite_priority
                )
            ranked = np.argsort(-scores, kind="stable")
            addition = next(
                (
                    int(index)
                    for index in ranked
                    if not mask[index]
                    and modular_rank(np.vstack((basis, vectors[index])))
                    == current_dimension
                ),
                None,
            )
            if addition is None:
                break
            basis = np.vstack((basis, vectors[addition]))
        if len(basis) != dimension:
            continue
        mask = numerical_span_mask(vectors, basis)
        key = np.packbits(mask).tobytes()
        if key in seen:
            continue
        seen.add(key)
        arithmetic_llr = arithmetic_enrichment_llr(mask, features)
        support = int(mask.sum())
        integral_support = int(np.sum(mask & integral))
        if priority_mode == "relations":
            base_score = support + 0.01 * float(np.sum(relation_priority[mask]))
        elif priority_mode == "finite":
            base_score = support + 0.01 * float(np.sum(finite_priority[mask]))
        else:
            base_score = arithmetic_llr
        # The population term prevents the Bernoulli likelihood from choosing
        # tiny, perfectly simple subspaces.  It does not assume a dimension-
        # specific shell size: the target is the median support of all paths
        # and is applied after proposal generation below.
        proposals.append(
            GrowthProposal(
                dimension=dimension,
                basis_rows=tuple(tuple(map(int, row)) for row in basis),
                # Public indices always refer to the caller's record order,
                # even though relation incidence is evaluated in the
                # complex's canonical vertex order.
                inlier_indices=tuple(
                    sorted(
                        input_index[ordered[index].coordinates]
                        for index in np.flatnonzero(mask)
                    )
                ),
                support=support,
                integral_support=integral_support,
                arithmetic_llr=arithmetic_llr,
                search_score=base_score,
            )
        )
    if not proposals:
        return ()
    median_support = float(np.median([proposal.support for proposal in proposals]))
    rescored = [
        GrowthProposal(
            dimension=proposal.dimension,
            basis_rows=proposal.basis_rows,
            inlier_indices=proposal.inlier_indices,
            support=proposal.support,
            integral_support=proposal.integral_support,
            arithmetic_llr=proposal.arithmetic_llr,
            search_score=(
                proposal.search_score
                - 0.01 * (proposal.support - median_support) ** 2
            ),
        )
        for proposal in proposals
    ]
    rescored.sort(key=lambda proposal: proposal.search_score, reverse=True)
    if maximum_proposals is not None:
        rescored = rescored[: int(maximum_proposals)]
    return tuple(rescored)


def extend_core_proposals(
    records: Sequence[ShortVectorRecord],
    core_bases: Sequence[Sequence[Sequence[int]]],
    *,
    target_dimension: int,
    target_support: float | None = None,
    support_penalty: float = 0.01,
) -> tuple[GrowthProposal, ...]:
    """Extend every rank-``k-1`` core by every retained ray, with deduplication."""

    if not records or not core_bases:
        return ()
    vectors = np.asarray([record.coordinates for record in records], dtype=np.int64)
    integral = np.asarray(
        [bool(record.arithmetic.get("integral", False)) for record in records],
        dtype=bool,
    )
    features = arithmetic_feature_flags(records)
    proposed: list[GrowthProposal] = []
    seen: set[bytes] = set()
    for raw_core in core_bases:
        core = np.asarray(raw_core, dtype=np.int64)
        if modular_rank(core) != target_dimension - 1:
            continue
        core_mask = numerical_span_mask(vectors, core)
        for index, vector in enumerate(vectors):
            if core_mask[index]:
                continue
            trial = np.vstack((core, vector))
            if modular_rank(trial) != target_dimension:
                continue
            mask = numerical_span_mask(vectors, trial)
            key = np.packbits(mask).tobytes()
            if key in seen:
                continue
            seen.add(key)
            support = int(mask.sum())
            arithmetic_llr = arithmetic_enrichment_llr(mask, features)
            proposed.append(
                GrowthProposal(
                    dimension=target_dimension,
                    basis_rows=tuple(tuple(map(int, row)) for row in trial),
                    inlier_indices=tuple(map(int, np.flatnonzero(mask))),
                    support=support,
                    integral_support=int(np.sum(mask & integral)),
                    arithmetic_llr=arithmetic_llr,
                    search_score=arithmetic_llr,
                )
            )
    if not proposed:
        return ()
    if target_support is None:
        target_support = float(np.median([proposal.support for proposal in proposed]))
    rescored = [
        GrowthProposal(
            dimension=proposal.dimension,
            basis_rows=proposal.basis_rows,
            inlier_indices=proposal.inlier_indices,
            support=proposal.support,
            integral_support=proposal.integral_support,
            arithmetic_llr=proposal.arithmetic_llr,
            search_score=(
                proposal.arithmetic_llr
                - support_penalty * (proposal.support - target_support) ** 2
            ),
        )
        for proposal in proposed
    ]
    rescored.sort(key=lambda proposal: proposal.search_score, reverse=True)
    return tuple(rescored)


def recombined_core_extension_search(
    records: Sequence[ShortVectorRecord],
    complex_: RelationComplex,
    *,
    dimension: int,
    seed_edges: int = 3_000,
    anchor_count: int = 20,
    enclosure_codimension: int = 3,
    enclosure_count: int = 20,
    inner_count: int = 2,
    enclosure_support_penalty: float = 0.002,
    refined_support_penalty: float = 0.002,
) -> RecombinedSearchLedger:
    """Build a diverse proposal ledger by bounded enclosure recombination.

    The routine deliberately optimizes recall, not single-fibre selection:
    independent rank-``k`` growth paths are paired with a small anchor set;
    rank-``k+codimension`` unions are ranked without labels; then rank-``k-1``
    cores inside the leading rank-``k`` paths are extended by every retained
    enclosure ray.  Cross-bound or cross-fibre matching must select from the
    returned ledger.
    """

    if not records:
        raise ValueError("proposal search needs nonempty records")
    ambient = len(records[0].coordinates)
    enclosure_dimension = dimension + enclosure_codimension
    if not 2 <= dimension < enclosure_dimension <= ambient:
        raise ValueError("invalid proposal/enclosure dimensions")
    direct = independent_relation_growth_proposals(
        records,
        complex_,
        dimension=dimension,
        seed_edges=seed_edges,
    )
    if not direct:
        return RecombinedSearchLedger(
            dimension, (), (), (), {
                "seed_edges": seed_edges,
                "anchor_count": anchor_count,
                "enclosure_codimension": enclosure_codimension,
                "enclosure_count": enclosure_count,
                "inner_count": inner_count,
            }
        )
    vectors = np.asarray([record.coordinates for record in records], dtype=np.int64)
    integral = np.asarray(
        [bool(record.arithmetic.get("integral", False)) for record in records],
        dtype=bool,
    )
    features = arithmetic_feature_flags(records)

    enclosures: list[GrowthProposal] = []
    seen_enclosures: set[tuple[tuple[int, ...], ...]] = set()
    for proposal in direct:
        for anchor in direct[:anchor_count]:
            basis = independent_row_basis(proposal.basis_rows + anchor.basis_rows)
            if len(basis) != enclosure_dimension:
                continue
            key = modular_row_space_key(basis)
            if key in seen_enclosures:
                continue
            seen_enclosures.add(key)
            mask = numerical_span_mask(vectors, np.asarray(basis, dtype=np.int64))
            support = int(mask.sum())
            arithmetic_llr = arithmetic_enrichment_llr(mask, features)
            enclosures.append(
                GrowthProposal(
                    dimension=enclosure_dimension,
                    basis_rows=basis,
                    inlier_indices=tuple(map(int, np.flatnonzero(mask))),
                    support=support,
                    integral_support=int(np.sum(mask & integral)),
                    arithmetic_llr=arithmetic_llr,
                    search_score=arithmetic_llr,
                )
            )
    if not enclosures:
        return RecombinedSearchLedger(
            dimension, direct, (), (), {
                "seed_edges": seed_edges,
                "anchor_count": anchor_count,
                "enclosure_codimension": enclosure_codimension,
                "enclosure_count": enclosure_count,
                "inner_count": inner_count,
            }
        )
    enclosure_target = float(np.median([item.support for item in enclosures]))
    enclosures = [
        GrowthProposal(
            dimension=item.dimension,
            basis_rows=item.basis_rows,
            inlier_indices=item.inlier_indices,
            support=item.support,
            integral_support=item.integral_support,
            arithmetic_llr=item.arithmetic_llr,
            search_score=(
                item.arithmetic_llr
                - enclosure_support_penalty
                * (item.support - enclosure_target) ** 2
            ),
        )
        for item in enclosures
    ]
    enclosures.sort(key=lambda item: item.search_score, reverse=True)

    refined_by_key: dict[tuple[tuple[int, ...], ...], GrowthProposal] = {}
    for enclosure in enclosures[:enclosure_count]:
        enclosure_records = tuple(records[index] for index in enclosure.inlier_indices)
        enclosure_complex = build_relation_complex(
            [record.coordinates for record in enclosure_records]
        )
        inner = independent_relation_growth_proposals(
            enclosure_records,
            enclosure_complex,
            dimension=dimension,
            seed_edges=seed_edges,
        )
        for proposal in inner[:inner_count]:
            proposal_records = tuple(
                enclosure_records[index] for index in proposal.inlier_indices
            )
            proposal_complex = build_relation_complex(
                [record.coordinates for record in proposal_records]
            )
            cores = independent_relation_growth_proposals(
                proposal_records,
                proposal_complex,
                dimension=dimension - 1,
                seed_edges=seed_edges,
            )
            extensions = extend_core_proposals(
                enclosure_records,
                [core.basis_rows for core in cores],
                target_dimension=dimension,
                target_support=proposal.support,
            )
            for extension in extensions:
                key = modular_row_space_key(extension.basis_rows)
                if key in refined_by_key:
                    continue
                mask = numerical_span_mask(
                    vectors, np.asarray(extension.basis_rows, dtype=np.int64)
                )
                arithmetic_llr = arithmetic_enrichment_llr(mask, features)
                refined_by_key[key] = GrowthProposal(
                    dimension=dimension,
                    basis_rows=extension.basis_rows,
                    inlier_indices=tuple(map(int, np.flatnonzero(mask))),
                    support=int(mask.sum()),
                    integral_support=int(np.sum(mask & integral)),
                    arithmetic_llr=arithmetic_llr,
                    search_score=arithmetic_llr,
                )
    refined = list(refined_by_key.values())
    if refined:
        refined_target = float(np.median([item.support for item in refined]))
        refined = [
            GrowthProposal(
                dimension=item.dimension,
                basis_rows=item.basis_rows,
                inlier_indices=item.inlier_indices,
                support=item.support,
                integral_support=item.integral_support,
                arithmetic_llr=item.arithmetic_llr,
                search_score=(
                    item.arithmetic_llr
                    - refined_support_penalty
                    * (item.support - refined_target) ** 2
                ),
            )
            for item in refined
        ]
        refined.sort(key=lambda item: item.search_score, reverse=True)
    return RecombinedSearchLedger(
        dimension=dimension,
        direct_proposals=direct,
        enclosure_proposals=tuple(enclosures),
        refined_proposals=tuple(refined),
        bounds={
            "seed_edges": seed_edges,
            "anchor_count": anchor_count,
            "enclosure_codimension": enclosure_codimension,
            "enclosure_count": enclosure_count,
            "inner_count": inner_count,
            "enclosure_support_penalty": enclosure_support_penalty,
            "refined_support_penalty": refined_support_penalty,
        },
    )


def beam_subspace_scan(
    records: Sequence[ShortVectorRecord],
    additive_degrees: Sequence[int] | RelationComplex,
    *,
    dimensions: Sequence[int] = tuple(range(10, 21)),
    pool: int = 300,
    beam_width: int = 8,
    branch_width: int = 80,
    candidates_per_dimension: int = 1,
    arithmetic_weight: float = 0.0,
    seed: int = 1729,
) -> tuple[SubspaceCandidate, ...]:
    """Greedily grow enriched exact subspaces without fixing their dimension.

    Search uses floating nullspaces only as a proposal score.  Every retained
    dimension winner is saturated and rescored by its exact integer kernel.
    The finite vector cloud, pool, beam and branching widths are explicit
    completeness boundaries; this routine is not an exhaustive Grassmannian
    search.
    """

    if isinstance(additive_degrees, RelationComplex):
        degree_by_vector = dict(
            zip(additive_degrees.vertices, additive_degrees.additive_degrees)
        )
        additive_degrees = tuple(
            degree_by_vector[record.coordinates] for record in records
        )
    if len(records) != len(additive_degrees):
        raise ValueError("record and degree counts differ")
    dimensions = tuple(sorted(set(map(int, dimensions))))
    if not dimensions or dimensions[0] < 1:
        raise ValueError("dimensions must be positive")
    if candidates_per_dimension < 1 or candidates_per_dimension > beam_width:
        raise ValueError("candidates_per_dimension must lie between one and beam_width")
    vectors = np.array([record.coordinates for record in records], dtype=np.int64)
    integral = np.array(
        [bool(record.arithmetic.get("integral", False)) for record in records],
        dtype=bool,
    )
    complexities = np.array(
        [max(1, int(record.arithmetic.get("total_bits", 1))) for record in records],
        dtype=float,
    )
    degrees = np.asarray(additive_degrees, dtype=float)
    arithmetic_features = arithmetic_feature_flags(records)
    complexity_scale = max(1.0, float(np.median(complexities)))
    favorable_count = np.sum(
        np.asarray([flag for _label, flag in arithmetic_features], dtype=float), axis=0
    )
    priority = (
        4.0 * integral.astype(float)
        + arithmetic_weight * favorable_count
        + np.log1p(degrees)
        - 0.35 * np.log1p(complexities / complexity_scale)
    )
    ordered = np.argsort(-priority, kind="stable")
    pool_indices = ordered[: min(pool, len(records))]
    rng = Random(seed)

    def state_score(basis: np.ndarray) -> tuple[float, np.ndarray]:
        mask = numerical_span_mask(vectors, basis)
        support = max(1, int(mask.sum()))
        llr = integrality_llr(mask, integral)
        arithmetic_llr = arithmetic_enrichment_llr(mask, arithmetic_features)
        relation_mass = float(np.sum(degrees[mask]))
        integral_support = float(np.sum(integral[mask]))
        if arithmetic_weight:
            score = (
                0.5 * support
                + 0.05 * relation_mass
                + arithmetic_weight * arithmetic_llr
            )
        else:
            score = (
                1.5 * support
                + 2.0 * integral_support
                + 0.35 * relation_mass
                + 3.0 * llr
            )
        return score, mask

    initial = list(pool_indices[: min(beam_width * 5, len(pool_indices))])
    rng.shuffle(initial)
    beam: list[tuple[float, np.ndarray]] = []
    for index in initial:
        basis = vectors[[index]]
        score, _mask = state_score(basis)
        beam.append((score, basis))
    beam.sort(key=lambda item: item[0], reverse=True)
    beam = beam[:beam_width]
    winners: dict[int, tuple[np.ndarray, ...]] = {}
    maximum = max(dimensions)
    for dimension in range(2, maximum + 1):
        proposed: dict[bytes, tuple[float, np.ndarray]] = {}
        candidate_indices = list(pool_indices[: min(branch_width, len(pool_indices))])
        # A small deterministic tail prevents an early priority mistake from
        # making all later dimensions identical.
        tail = list(pool_indices[min(branch_width, len(pool_indices)) :])
        rng.shuffle(tail)
        candidate_indices.extend(tail[: max(4, branch_width // 8)])
        for _old_score, basis in beam:
            for index in candidate_indices:
                trial = np.vstack((basis, vectors[index]))
                if modular_rank(trial) != dimension:
                    continue
                score, _mask = state_score(trial)
                _u, _s, vh = np.linalg.svd(trial.astype(float), full_matrices=False)
                # The orthogonal projector is invariant under changing the
                # proposal basis.  Using abs(vh) here would not be a row-space
                # key and could silently retain duplicate spans.
                key = np.round(vh.T @ vh, 8).tobytes()
                current = proposed.get(key)
                if current is None or score > current[0]:
                    proposed[key] = (score, trial)
        beam = sorted(proposed.values(), key=lambda item: item[0], reverse=True)[:beam_width]
        if not beam:
            break
        if dimension in dimensions:
            winners[dimension] = tuple(
                item[1] for item in beam[:candidates_per_dimension]
            )

    answers = []
    vector_rows = [tuple(map(int, row)) for row in vectors]
    degree_int = np.asarray(additive_degrees, dtype=np.int64)
    for dimension in dimensions:
        proposals = winners.get(dimension)
        if proposals is None:
            continue
        seen: set[tuple[tuple[int, ...], ...]] = set()
        for proposal in proposals:
            primitive = primitive_span_basis(
                [tuple(map(int, row)) for row in proposal]
            )
            if rational_rank(primitive) != dimension:
                raise ArithmeticError("primitive candidate changed dimension")
            if primitive in seen:
                continue
            seen.add(primitive)
            mask = exact_span_mask(vector_rows, primitive)
            support = int(mask.sum())
            integral_support = int(np.sum(mask & integral))
            outside_total = len(mask) - support
            outside_integral = int(np.sum(~mask & integral))
            llr = integrality_llr(mask, integral)
            arithmetic_llr = arithmetic_enrichment_llr(mask, arithmetic_features)
            relation_mass = int(np.sum(degree_int[mask]))
            if arithmetic_weight:
                score = (
                    0.5 * support
                    + 0.05 * relation_mass
                    + arithmetic_weight * arithmetic_llr
                )
            else:
                score = 1.5 * support + 2.0 * integral_support + 0.35 * relation_mass + 3.0 * llr
            answers.append(
                SubspaceCandidate(
                    dimension=dimension,
                    primitive_basis_rows=primitive,
                    inlier_indices=tuple(map(int, np.flatnonzero(mask))),
                    support=support,
                    integral_support=integral_support,
                    integral_rate=integral_support / support if support else 0.0,
                    outside_integral_rate=(
                        outside_integral / outside_total if outside_total else 0.0
                    ),
                    integrality_llr=llr,
                    relation_mass=relation_mass,
                    search_score=score,
                    arithmetic_llr=arithmetic_llr,
                )
            )
    return tuple(answers)


def relation_seeded_subspace_scan(
    records: Sequence[ShortVectorRecord],
    complex_: RelationComplex,
    *,
    dimension: int,
    seed_edges: int = 2_000,
    beam_width: int = 32,
    branch_width: int = 80,
    candidates: int = 8,
) -> tuple[SubspaceCandidate, ...]:
    """Grow subspaces along exact additive hyperedges.

    If one ray of ``a +/- b = c`` is already in a subspace, adjoining either
    of the other two adjoins both.  This bounded beam search exploits that
    basis-independent fact.  The declared edge, beam, and branch limits are
    part of the result's search boundary.
    """

    if not 2 <= dimension <= len(records[0].coordinates):
        raise ValueError("dimension is outside the ambient range")
    if candidates < 1 or candidates > beam_width:
        raise ValueError("candidates must lie between one and beam_width")
    record_by_vector = {record.coordinates: record for record in records}
    if set(record_by_vector) != set(complex_.vertices):
        raise ValueError("records and relation-complex vertices differ")
    ordered_records = [record_by_vector[vector] for vector in complex_.vertices]
    input_index = {
        record.coordinates: index for index, record in enumerate(records)
    }
    vectors = np.asarray([record.coordinates for record in ordered_records], dtype=np.int64)
    integral = np.asarray(
        [bool(record.arithmetic.get("integral", False)) for record in ordered_records],
        dtype=bool,
    )
    degrees = np.asarray(complex_.additive_degrees, dtype=np.int64)
    complexities = np.asarray(
        [max(1, int(record.arithmetic.get("total_bits", 1))) for record in ordered_records],
        dtype=float,
    )
    complexity_scale = max(1.0, float(np.median(complexities)))
    vertex_priority = (
        4.0 * integral.astype(float)
        + np.log1p(degrees)
        - 0.35 * np.log1p(complexities / complexity_scale)
    )
    edges = np.asarray(complex_.ternary_relations, dtype=np.int64)

    def state_score(basis: np.ndarray) -> tuple[float, np.ndarray]:
        mask = numerical_span_mask(vectors, basis)
        support = max(1, int(mask.sum()))
        llr = integrality_llr(mask, integral)
        relation_mass = float(np.sum(degrees[mask]))
        integral_support = float(np.sum(integral[mask]))
        return (
            1.5 * support
            + 2.0 * integral_support
            + 0.35 * relation_mass
            + 3.0 * llr,
            mask,
        )

    edge_priority = np.sum(vertex_priority[edges], axis=1)
    seed_order = np.argsort(-edge_priority, kind="stable")[: min(seed_edges, len(edges))]
    proposed: dict[bytes, tuple[float, np.ndarray, np.ndarray]] = {}
    for edge_index in seed_order:
        edge = edges[edge_index]
        basis = vectors[edge[:2]]
        if modular_rank(basis) != 2:
            continue
        score, mask = state_score(basis)
        _u, _s, vh = np.linalg.svd(basis.astype(float), full_matrices=False)
        key = np.round(vh.T @ vh, 8).tobytes()
        old = proposed.get(key)
        if old is None or score > old[0]:
            proposed[key] = score, basis, mask
    beam = sorted(proposed.values(), key=lambda item: item[0], reverse=True)[:beam_width]

    for current_dimension in range(3, dimension + 1):
        proposed = {}
        for _score, basis, mask in beam:
            occupancy = np.sum(mask[edges], axis=1)
            crossing = edges[occupancy == 1]
            counts = Counter(
                int(vertex)
                for edge in crossing
                for vertex in edge
                if not mask[vertex]
            )
            ranked = sorted(
                counts,
                key=lambda index: (counts[index], vertex_priority[index]),
                reverse=True,
            )[:branch_width]
            for index in ranked:
                trial = np.vstack((basis, vectors[index]))
                if modular_rank(trial) != current_dimension:
                    continue
                score, trial_mask = state_score(trial)
                _u, _s, vh = np.linalg.svd(trial.astype(float), full_matrices=False)
                key = np.round(vh.T @ vh, 8).tobytes()
                old = proposed.get(key)
                if old is None or score > old[0]:
                    proposed[key] = score, trial, trial_mask
        beam = sorted(proposed.values(), key=lambda item: item[0], reverse=True)[:beam_width]
        if not beam:
            return ()

    answers: list[SubspaceCandidate] = []
    vector_rows = [tuple(map(int, row)) for row in vectors]
    for _score, proposal, _numerical_mask in beam[:candidates]:
        primitive = primitive_span_basis([tuple(map(int, row)) for row in proposal])
        mask = exact_span_mask(vector_rows, primitive)
        support = int(mask.sum())
        integral_support = int(np.sum(mask & integral))
        outside_total = len(mask) - support
        outside_integral = int(np.sum(~mask & integral))
        llr = integrality_llr(mask, integral)
        relation_mass = int(np.sum(degrees[mask]))
        score = 1.5 * support + 2.0 * integral_support + 0.35 * relation_mass + 3.0 * llr
        answers.append(
            SubspaceCandidate(
                dimension=dimension,
                primitive_basis_rows=primitive,
                inlier_indices=tuple(
                    sorted(
                        input_index[ordered_records[index].coordinates]
                        for index in np.flatnonzero(mask)
                    )
                ),
                support=support,
                integral_support=integral_support,
                integral_rate=integral_support / support if support else 0.0,
                outside_integral_rate=outside_integral / outside_total if outside_total else 0.0,
                integrality_llr=llr,
                relation_mass=relation_mass,
                search_score=score,
            )
        )
    return tuple(answers)

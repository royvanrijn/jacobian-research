"""Intrinsic numerical height-lattice signatures with explicit enumeration bounds."""

from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import Sequence

import numpy as np

from .integer import canonical_unoriented, content
from .pari import run_gp
from .relations import RelationComplex, build_relation_complex


@dataclass(frozen=True)
class CloudHeightSignature:
    """Scale-free height quantiles of a retained candidate short-vector cloud.

    This is inexpensive and basis-independent, but unlike :class:`ThetaSignature`
    it is only complete through the ambient fibre's declared height cutoff.
    """

    vector_count: int
    quantile_count: int
    normalized_log_quantiles: tuple[str, ...]

    def profile(self) -> np.ndarray:
        return np.asarray(self.normalized_log_quantiles, dtype=float)

    def to_record(self) -> dict[str, object]:
        return {
            "retained_vector_count": self.vector_count,
            "quantile_count": self.quantile_count,
            "normalized_log_height_quantiles": list(self.normalized_log_quantiles),
        }


def cloud_height_signature(
    heights: Sequence[float | str],
    inlier_indices: Sequence[int],
    *,
    quantiles: int = 64,
) -> CloudHeightSignature:
    """Return a scale-free quantile sketch for selected unoriented heights."""

    if quantiles < 2:
        raise ValueError("at least two height quantiles are required")
    indices = tuple(map(int, inlier_indices))
    if not indices:
        raise ValueError("a cloud height signature needs at least one inlier")
    if min(indices) < 0 or max(indices) >= len(heights):
        raise ValueError("an inlier index is outside the height population")
    values = np.asarray([float(heights[index]) for index in indices], dtype=float)
    if np.min(values) <= 0 or not np.all(np.isfinite(values)):
        raise ValueError("canonical heights must be finite and positive")
    logs = np.log(values)
    logs -= float(np.mean(logs))
    profile = np.quantile(logs, np.linspace(0.0, 1.0, int(quantiles)))
    return CloudHeightSignature(
        vector_count=len(values),
        quantile_count=int(quantiles),
        normalized_log_quantiles=tuple(f"{value:.17g}" for value in profile),
    )


@dataclass(frozen=True)
class HermiteSignature:
    """Scale- and basis-invariant determinant/minimum shape statistic."""

    rank: int
    minimum: str
    determinant: str
    log_hermite_invariant: str

    def to_record(self) -> dict[str, object]:
        return {
            "rank": self.rank,
            "minimum": self.minimum,
            "determinant": self.determinant,
            "log_hermite_invariant": self.log_hermite_invariant,
        }


def hermite_signature(
    gram: Sequence[Sequence[float | str]],
    *,
    digits: int = 80,
    timeout: float = 300.0,
) -> HermiteSignature:
    """Compute ``log(min)-log(det)/rank`` after exact-to-precision SVP.

    The determinant and minimum scale compatibly, so the final statistic is
    unchanged by a scalar specialization factor.  It is numerical at the
    supplied Gram precision.
    """

    literal = _gp_real_matrix(gram)
    program = f"""
default(realprecision,{int(digits)});
G={literal};n=matsize(G)[1];U=qflllgram(G);R=U~*G*U;
B=vecmin(vector(n,i,R[i,i]));Q=qfminim(R,B,100000,2);V=Q[3];
L=vector(matsize(V)[2],j,V[,j]~*R*V[,j]);mn=vecmin(L);D=matdet(G);
print("HERMITE|",n,"|",mn,"|",D,"|",log(mn)-log(D)/n);
"""
    lines = run_gp(program, timeout=timeout)
    rank, minimum, determinant, invariant = next(
        line.split("|")[1:] for line in lines if line.startswith("HERMITE|")
    )
    return HermiteSignature(int(rank), minimum, determinant, invariant)


def hermite_signature_distance(left: HermiteSignature, right: HermiteSignature) -> float:
    if left.rank != right.rank:
        return float("inf")
    return abs(float(left.log_hermite_invariant) - float(right.log_hermite_invariant))


@dataclass(frozen=True)
class RelationMetricSignature:
    """Scale-free length geometry on an induced additive relation complex."""

    relation_count: int
    quantile_count: int
    normalized_sorted_log_length_quantiles: tuple[str, ...]

    def profile(self) -> np.ndarray:
        return np.asarray(self.normalized_sorted_log_length_quantiles, dtype=float)

    def to_record(self) -> dict[str, object]:
        return {
            "retained_ternary_relation_count": self.relation_count,
            "quantile_count": self.quantile_count,
            "normalized_sorted_log_length_quantiles": list(
                self.normalized_sorted_log_length_quantiles
            ),
        }


def relation_metric_signature(
    vectors: Sequence[Sequence[int]],
    heights: Sequence[float | str],
    inlier_indices: Sequence[int],
    complex_: RelationComplex,
    *,
    quantiles: int = 64,
) -> RelationMetricSignature:
    """Sketch normalized edge-length triples on the candidate subcomplex."""

    if len(vectors) != len(heights):
        raise ValueError("vector and height populations differ")
    if quantiles < 2:
        raise ValueError("at least two relation-metric quantiles are required")
    input_index = {
        canonical_unoriented(vector): index for index, vector in enumerate(vectors)
    }
    if len(input_index) != len(vectors):
        raise ValueError("vector population repeats an unoriented ray")
    if set(input_index) != set(complex_.vertices):
        raise ValueError("vectors and relation-complex vertices differ")
    retained_input = set(map(int, inlier_indices))
    retained = tuple(input_index[vertex] in retained_input for vertex in complex_.vertices)
    log_heights = np.log(np.asarray([float(value) for value in heights], dtype=float))
    triples = []
    for edge in complex_.ternary_relations:
        if not all(retained[index] for index in edge):
            continue
        values = np.sort(
            np.asarray(
                [log_heights[input_index[complex_.vertices[index]]] for index in edge],
                dtype=float,
            )
        )
        triples.append(values - float(np.mean(values)))
    if not triples:
        raise ValueError("candidate has no retained ternary relation")
    matrix = np.asarray(triples, dtype=float)
    positions = np.linspace(0.0, 1.0, int(quantiles))
    profile = np.concatenate(
        [np.quantile(matrix[:, column], positions) for column in range(3)]
    )
    return RelationMetricSignature(
        relation_count=len(triples),
        quantile_count=int(quantiles),
        normalized_sorted_log_length_quantiles=tuple(
            f"{value:.17g}" for value in profile
        ),
    )


def relation_metric_profile_distance(
    left: RelationMetricSignature, right: RelationMetricSignature
) -> float:
    if left.quantile_count != right.quantile_count:
        return float("inf")
    return float(np.sqrt(np.mean((left.profile() - right.profile()) ** 2)))


def cloud_height_profile_distance(
    left: CloudHeightSignature, right: CloudHeightSignature
) -> float:
    """RMS distance between two scale-free bounded-cloud profiles."""

    if left.quantile_count != right.quantile_count:
        return float("inf")
    return float(np.sqrt(np.mean((left.profile() - right.profile()) ** 2)))


def cloud_height_signature_from_record(
    record: dict[str, object],
) -> CloudHeightSignature:
    return CloudHeightSignature(
        vector_count=int(record["retained_vector_count"]),
        quantile_count=int(record["quantile_count"]),
        normalized_log_quantiles=tuple(
            map(str, record["normalized_log_height_quantiles"])
        ),
    )


def restricted_height_gram(
    ambient_gram: Sequence[Sequence[float | str]],
    basis_rows: Sequence[Sequence[int]],
) -> np.ndarray:
    """Restrict an ambient numerical Gram to an integral row basis."""

    ambient = np.asarray(ambient_gram, dtype=float)
    basis = np.asarray(basis_rows, dtype=float)
    if ambient.ndim != 2 or ambient.shape[0] != ambient.shape[1]:
        raise ValueError("ambient Gram must be square")
    if basis.ndim != 2 or basis.shape[1] != ambient.shape[0]:
        raise ValueError("basis and ambient Gram dimensions differ")
    answer = basis @ ambient @ basis.T
    if np.min(np.linalg.eigvalsh(answer)) <= 0:
        raise ValueError("restricted Gram is not positive definite")
    return answer


def _gp_real_matrix(matrix: Sequence[Sequence[float | str]]) -> str:
    rows = [tuple(row) for row in matrix]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("Gram must be nonempty and square")
    return "[" + ";".join(
        ",".join(str(value) for value in row) for row in rows
    ) + "]"


@dataclass(frozen=True)
class ThetaSignature:
    """Complete unoriented-vector lengths through one adaptive norm bound."""

    rank: int
    minimum: str
    enumeration_bound: str
    growth_steps: int
    vector_count: int
    lengths: tuple[str, ...]

    def normalized_log_lengths(self, count: int | None = None) -> np.ndarray:
        values = np.asarray([float(value) for value in self.lengths], dtype=float)
        if count is not None:
            values = values[: int(count)]
        if not len(values):
            return values
        logs = np.log(values)
        return logs - float(np.mean(logs))

    def quantile_profile(self, count: int = 64) -> np.ndarray:
        values = self.normalized_log_lengths()
        if not len(values):
            return values
        positions = np.linspace(0.0, 1.0, int(count))
        return np.quantile(values, positions)

    def to_record(self) -> dict[str, object]:
        return {
            "rank": self.rank,
            "minimum": self.minimum,
            "enumeration_bound": self.enumeration_bound,
            "growth_steps": self.growth_steps,
            "unoriented_vector_count": self.vector_count,
            "lengths": list(self.lengths),
        }


@dataclass(frozen=True)
class IntrinsicShellSignature:
    """Exact additive complex on a complete numerically bounded lattice shell.

    The selected vectors and their relations are exact integer data.  Only the
    decision that a vector lies below the declared canonical-height boundary
    is numerical.
    """

    rank: int
    minimum: str
    enumeration_bound: str
    growth_steps: int
    raw_vector_count: int
    primitive_vector_count: int
    normalized_log_height_quantiles: tuple[str, ...]
    relation_complex: RelationComplex

    def to_record(self, *, include_relations: bool = False) -> dict[str, object]:
        return {
            "rank": self.rank,
            "minimum": self.minimum,
            "enumeration_bound": self.enumeration_bound,
            "growth_steps": self.growth_steps,
            "raw_unoriented_vector_count": self.raw_vector_count,
            "primitive_unoriented_vector_count": self.primitive_vector_count,
            "normalized_log_height_quantiles": list(
                self.normalized_log_height_quantiles
            ),
            "relation_complex": self.relation_complex.to_record(
                include_relations=include_relations
            ),
        }


@dataclass(frozen=True)
class IntrinsicShellDistance:
    """Decomposed soft distance between two exact intrinsic shell complexes."""

    height_quantile_rms: float
    centered_degree_quantile_rms: float
    log_population_difference: float
    log_edge_density_difference: float
    total: float


def intrinsic_shell_profile_distance(
    left: IntrinsicShellSignature,
    right: IntrinsicShellSignature,
) -> IntrinsicShellDistance:
    """Compare deformed shells without requiring literal graph isomorphism."""

    if left.rank != right.rank:
        return IntrinsicShellDistance(*(float("inf"),) * 5)
    left_height = np.asarray(left.normalized_log_height_quantiles, dtype=float)
    right_height = np.asarray(right.normalized_log_height_quantiles, dtype=float)
    if len(left_height) != len(right_height):
        return IntrinsicShellDistance(*(float("inf"),) * 5)
    quantiles = len(left_height)

    def degree_profile(signature: IntrinsicShellSignature) -> np.ndarray:
        values = np.log1p(
            np.asarray(signature.relation_complex.additive_degrees, dtype=float)
        )
        values -= float(np.mean(values))
        return np.quantile(values, np.linspace(0.0, 1.0, quantiles))

    height_distance = float(np.sqrt(np.mean((left_height - right_height) ** 2)))
    degree_distance = float(
        np.sqrt(np.mean((degree_profile(left) - degree_profile(right)) ** 2))
    )
    population_distance = abs(
        log(left.primitive_vector_count / right.primitive_vector_count)
    )
    left_density = (len(left.relation_complex.ternary_relations) + 1) / (
        left.primitive_vector_count + 1
    )
    right_density = (len(right.relation_complex.ternary_relations) + 1) / (
        right.primitive_vector_count + 1
    )
    density_distance = abs(log(left_density / right_density))
    total = (
        height_distance
        + degree_distance
        + 0.25 * population_distance
        + 0.25 * density_distance
    )
    return IntrinsicShellDistance(
        height_quantile_rms=height_distance,
        centered_degree_quantile_rms=degree_distance,
        log_population_difference=population_distance,
        log_edge_density_difference=density_distance,
        total=total,
    )


def intrinsic_shell_signature(
    gram: Sequence[Sequence[float | str]],
    *,
    minimum_vectors: int = 128,
    quantiles: int = 64,
    growth: float = 1.12,
    maximum_steps: int = 40,
    maximum_vectors: int = 100_000,
    digits: int = 80,
    timeout: float = 300.0,
) -> IntrinsicShellSignature:
    """Enumerate a complete intrinsic shell and retain its exact relations.

    PARI chooses a complete norm boundary on an LLL-reduced Gram.  Its vectors
    are mapped back to the supplied integral basis, filtered to primitive
    unoriented rays, and canonicalized before relation construction.
    """

    if minimum_vectors < 1:
        raise ValueError("minimum_vectors must be positive")
    if quantiles < 2:
        raise ValueError("at least two height quantiles are required")
    if not 1.0 < growth < 2.0:
        raise ValueError("growth must lie strictly between one and two")
    if maximum_vectors < minimum_vectors:
        raise ValueError("maximum_vectors must be at least minimum_vectors")
    literal = _gp_real_matrix(gram)
    program = f"""
default(realprecision,{int(digits)});
G={literal};U=qflllgram(G);R=U~*G*U;
B=vecmin(vector(matsize(R)[1],i,R[i,i]));
Q=qfminim(R,B,{int(maximum_vectors)},2);s=0;
while(matsize(Q[3])[2] < {int(minimum_vectors)} && s < {int(maximum_steps)},{{B=B*{float(growth):.17g};Q=qfminim(R,B,{int(maximum_vectors)},2);s++}});
V=U*Q[3];L=vector(matsize(V)[2],j,V[,j]~*G*V[,j]);mn=vecmin(L);
print("META|",matsize(G)[1],"|",mn,"|",B,"|",s,"|",#L);
print("BEGIN");for(j=1,matsize(V)[2],{{for(i=1,matsize(V)[1],if(i>1,print1(","));print1(V[i,j]));print("|",L[j])}});print("END");
"""
    lines = run_gp(program, timeout=timeout)
    rank, minimum, bound, steps, raw_count = next(
        line.split("|")[1:] for line in lines if line.startswith("META|")
    )
    vectors: dict[tuple[int, ...], float] = {}
    for line in lines[lines.index("BEGIN") + 1 : lines.index("END")]:
        raw_text, height_text = line.split("|", 1)
        raw = tuple(map(int, raw_text.split(",")))
        if content(raw) != 1:
            continue
        vector = canonical_unoriented(raw)
        height = float(height_text)
        old = vectors.setdefault(vector, height)
        if abs(old - height) > 1e-8 * max(1.0, abs(height)):
            raise ArithmeticError("opposite shell vectors have inconsistent heights")
    if len(vectors) < minimum_vectors:
        raise ArithmeticError(
            "complete raw shell did not contain the requested primitive population"
        )
    if int(raw_count) >= maximum_vectors:
        raise ArithmeticError("intrinsic shell reached maximum_vectors cap")
    ordered = tuple(sorted(vectors))
    log_heights = np.log(np.asarray([vectors[vector] for vector in ordered]))
    log_heights -= float(np.mean(log_heights))
    profile = np.quantile(log_heights, np.linspace(0.0, 1.0, int(quantiles)))
    return IntrinsicShellSignature(
        rank=int(rank),
        minimum=minimum,
        enumeration_bound=bound,
        growth_steps=int(steps),
        raw_vector_count=int(raw_count),
        primitive_vector_count=len(ordered),
        normalized_log_height_quantiles=tuple(
            f"{value:.17g}" for value in profile
        ),
        relation_complex=build_relation_complex(ordered),
    )


def theta_signature(
    gram: Sequence[Sequence[float | str]],
    *,
    minimum_vectors: int = 256,
    growth: float = 1.12,
    maximum_steps: int = 40,
    maximum_vectors: int = 100_000,
    digits: int = 80,
    timeout: float = 300.0,
) -> ThetaSignature:
    """Enumerate a complete adaptive initial theta segment.

    The norm bound starts at the exact-to-current-precision lattice minimum
    returned by PARI and grows geometrically until at least
    ``minimum_vectors`` unoriented vectors have been enumerated.  No ``m`` cap
    is supplied to ``qfminim``: every vector through the final bound is kept.
    Canonical heights remain numerical, but enumeration completeness is exact
    for the supplied real Gram and bound.
    """

    if minimum_vectors < 1:
        raise ValueError("minimum_vectors must be positive")
    if not 1.0 < growth < 2.0:
        raise ValueError("growth must lie strictly between one and two")
    if maximum_steps < 0:
        raise ValueError("maximum_steps must be nonnegative")
    if maximum_vectors < minimum_vectors:
        raise ValueError("maximum_vectors must be at least minimum_vectors")
    literal = _gp_real_matrix(gram)
    program = f"""
default(realprecision,{int(digits)});
G={literal};U=qflllgram(G);R=U~*G*U;
B=vecmin(vector(matsize(R)[1],i,R[i,i]));
Q=qfminim(R,B,{int(maximum_vectors)},2);V=Q[3];
L=vector(matsize(V)[2],j,V[,j]~*R*V[,j]);mn=vecmin(L);
B=mn;Q=qfminim(R,B,{int(maximum_vectors)},2);s=0;
while(matsize(Q[3])[2] < {int(minimum_vectors)} && s < {int(maximum_steps)},{{B=B*{float(growth):.17g};Q=qfminim(R,B,{int(maximum_vectors)},2);s++}});
V=Q[3];L=vector(matsize(V)[2],j,V[,j]~*R*V[,j]);L=vecsort(L);
print("META|",matsize(G)[1],"|",mn,"|",B,"|",s,"|",#L);
print("BEGIN");for(i=1,#L,print(L[i]));print("END");
"""
    lines = run_gp(program, timeout=timeout)
    meta = next(line.split("|")[1:] for line in lines if line.startswith("META|"))
    rank, minimum, bound, steps, count = meta
    start = lines.index("BEGIN") + 1
    stop = lines.index("END", start)
    lengths = tuple(lines[start:stop])
    if len(lengths) != int(count):
        raise ArithmeticError("PARI theta length count changed during parsing")
    if len(lengths) < minimum_vectors and int(steps) >= maximum_steps:
        raise ArithmeticError("adaptive theta enumeration exhausted maximum_steps")
    if len(lengths) >= maximum_vectors:
        raise ArithmeticError("theta enumeration reached maximum_vectors cap")
    return ThetaSignature(
        rank=int(rank),
        minimum=minimum,
        enumeration_bound=bound,
        growth_steps=int(steps),
        vector_count=int(count),
        lengths=lengths,
    )


def theta_profile_distance(
    left: ThetaSignature,
    right: ThetaSignature,
    *,
    quantiles: int = 64,
) -> float:
    """Scale-invariant RMS distance between normalized theta quantiles."""

    if left.rank != right.rank:
        return float("inf")
    left_profile = left.quantile_profile(quantiles)
    right_profile = right.quantile_profile(quantiles)
    if len(left_profile) != len(right_profile):
        return float("inf")
    return float(np.sqrt(np.mean((left_profile - right_profile) ** 2)))

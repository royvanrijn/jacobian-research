"""Intrinsic numerical height-lattice signatures with explicit enumeration bounds."""

from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import Sequence

import numpy as np

from .pari import run_gp


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

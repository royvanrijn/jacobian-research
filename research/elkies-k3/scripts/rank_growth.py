from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import os
import tempfile
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ExtensionMetrics:
    base_rank: int
    candidate_index: int
    candidate_height: float
    projection_height: float
    orthogonal_height: float
    orthogonal_ratio: float
    pairing_vector: np.ndarray
    projection_coefficients: np.ndarray


@dataclass(frozen=True)
class CascadeMetrics:
    base_rank: int
    original_rank: int
    candidate_index: int
    new_block_dimension: int
    new_block_projection_height: float
    new_block_projection_share: float
    last_increment_corr: float
    last_increment_energy_share: float
    new_block_coefficients: np.ndarray


def load_gram(path: str | Path) -> np.ndarray:
    gram = np.loadtxt(path, dtype=np.float64)
    if gram.ndim != 2 or gram.shape[0] != gram.shape[1]:
        raise ValueError(f"Gram matrix must be square, got {gram.shape}: {path}")
    if not np.allclose(gram, gram.T, rtol=1e-10, atol=1e-10):
        raise ValueError(f"Gram matrix is not symmetric: {path}")
    return gram


def _solve_psd(gram: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve in a positive-definite height Gram with a useful error."""
    try:
        return np.linalg.solve(gram, rhs)
    except np.linalg.LinAlgError as exc:
        eig = np.linalg.eigvalsh(gram)
        raise ValueError(
            "Base height Gram is singular/ill-conditioned; "
            f"smallest eigenvalue={float(eig[0]):.12g}"
        ) from exc


def extension_metrics(
    gram: np.ndarray,
    base_rank: int,
    candidate_index: int,
    *,
    roundoff: float = 1e-8,
) -> ExtensionMetrics:
    """Schur-complement geometry of one candidate above a known base span.

    ``gram`` contains the base vectors first. ``candidate_index`` selects a
    later vector from the same height Gram.  The invariant quantity

        delta = h(P) - b^T G^-1 b

    is the squared height of the component perpendicular to the base span.
    Positive delta is the numerical rank-growth/ignition signal.
    """
    n = gram.shape[0]
    if not 0 < base_rank < n:
        raise ValueError(f"base_rank={base_rank} incompatible with {n}x{n} Gram")
    if not base_rank <= candidate_index < n:
        raise ValueError(
            f"candidate_index={candidate_index} must be in [{base_rank},{n})"
        )

    base = gram[:base_rank, :base_rank]
    pairing = np.asarray(gram[:base_rank, candidate_index], dtype=np.float64)
    height = float(gram[candidate_index, candidate_index])
    coeff = _solve_psd(base, pairing)
    projection = float(pairing @ coeff)
    delta = float(height - projection)
    if delta < 0.0 and delta > -roundoff:
        delta = 0.0
    ratio = delta / height if height > 0.0 else math.nan

    return ExtensionMetrics(
        base_rank=base_rank,
        candidate_index=candidate_index,
        candidate_height=height,
        projection_height=projection,
        orthogonal_height=delta,
        orthogonal_ratio=ratio,
        pairing_vector=pairing,
        projection_coefficients=coeff,
    )


def cascade_metrics(
    gram: np.ndarray,
    original_rank: int,
    base_rank: int,
    candidate_index: int,
    *,
    roundoff: float = 1e-8,
) -> CascadeMetrics:
    """Measure coupling to the sequentially-added rank-growth block.

    The first ``original_rank`` vectors are the original core.  Vectors
    ``original_rank:base_rank`` are prior ignition/cascade winners.  We
    quotient the latter block by the original core via a Schur complement.

    ``last_increment_corr`` is the absolute correlation of the candidate's
    projected new-block component with the newest *orthogonal increment* in
    that sequential block.  This is the quantity that was ~0.95 median in the
    controlled rank-21 experiment for the 19->20 and 20->21 steps.
    """
    if not 0 < original_rank < base_rank:
        raise ValueError("cascade metrics require at least one added direction")
    if candidate_index < base_rank or candidate_index >= gram.shape[0]:
        raise ValueError("candidate_index must follow the current base")

    current = gram[:base_rank, :base_rank]
    b = np.asarray(gram[:base_rank, candidate_index], dtype=np.float64)
    candidate_height = float(gram[candidate_index, candidate_index])

    c0 = original_rank
    g00 = current[:c0, :c0]
    g01 = current[:c0, c0:]
    g10 = current[c0:, :c0]
    g11 = current[c0:, c0:]
    b0 = b[:c0]
    b1 = b[c0:]

    s = g11 - g10 @ _solve_psd(g00, g01)
    residual_pair = b1 - g10 @ _solve_psd(g00, b0)
    alpha = _solve_psd(s, residual_pair)
    new_block_height = float(residual_pair @ alpha)
    if new_block_height < 0.0 and new_block_height > -roundoff:
        new_block_height = 0.0

    full_coeff = _solve_psd(current, b)
    full_projection = float(b @ full_coeff)
    share = (
        new_block_height / full_projection
        if full_projection > 0.0
        else 0.0
    )

    try:
        chol = np.linalg.cholesky(s)
        z = chol.T @ alpha
        energy = z * z
        total = float(energy.sum())
        if total > 0.0:
            last_share = float(energy[-1] / total)
            last_share = min(1.0, max(0.0, last_share))
            last_corr = math.sqrt(last_share)
        else:
            last_share = 0.0
            last_corr = 0.0
    except np.linalg.LinAlgError:
        last_share = math.nan
        last_corr = math.nan

    return CascadeMetrics(
        base_rank=base_rank,
        original_rank=original_rank,
        candidate_index=candidate_index,
        new_block_dimension=base_rank - original_rank,
        new_block_projection_height=new_block_height,
        new_block_projection_share=share,
        last_increment_corr=last_corr,
        last_increment_energy_share=last_share,
        new_block_coefficients=alpha,
    )


def matrix_numerical_rank(gram: np.ndarray, relative_tol: float = 1e-10) -> int:
    eig = np.linalg.eigvalsh(gram)
    scale = max(1.0, float(np.max(np.abs(eig))))
    return int(np.count_nonzero(eig > relative_tol * scale))


def jsonable_metrics(metrics: ExtensionMetrics | CascadeMetrics) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in vars(metrics).items():
        if isinstance(value, np.ndarray):
            result[key] = [float(x) for x in value]
        elif isinstance(value, (np.floating, float)):
            result[key] = float(value)
        elif isinstance(value, (np.integer, int)):
            result[key] = int(value)
        else:
            result[key] = value
    return result


def atomic_write_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise

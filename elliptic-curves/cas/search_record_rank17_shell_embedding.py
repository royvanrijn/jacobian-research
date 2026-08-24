#!/usr/bin/env python3
"""Search for a GL(17,Z) embedding with a coherent full R17 minimal shell.

Unlike the basis-entry calibration, this objective evaluates all 1,311
unoriented norm-four vectors after every proposed basis change.  It remains a
bounded numerical search, and the known curve-245 generic-rank-12 control must
be used to calibrate any claimed separation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import random

import numpy as np

from calibrate_record_rank17_fingerprint import (
    curve245_control,
    exact_r17_reduced_form_and_shell,
    integral_inverse,
    numerical_short_vectors,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_core_candidates_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_shell_embedding_search_v1.json"
)


def shell_statistics(
    gram: np.ndarray,
    basis: np.ndarray,
    shell: np.ndarray,
) -> tuple[float, dict[str, object]]:
    images = shell @ basis.T
    heights = np.einsum("ij,jk,ik->i", images, gram, images)
    mean = float(np.mean(heights))
    cv = float(np.std(heights) / mean)
    normalized = 4.0 * heights / mean
    quantiles = np.quantile(normalized, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])
    return cv * cv, {
        "raw_mean_height": f"{mean:.17g}",
        "coefficient_of_variation": f"{cv:.17g}",
        "normalized_rms_deviation_from_4": f"{4.0 * cv:.17g}",
        "normalized_quantiles_min_10_25_50_75_90_max": [
            f"{value:.17g}" for value in quantiles
        ],
        "fraction_with_normalized_height_between_3_and_5": f"{float(np.mean((normalized >= 3.0) & (normalized <= 5.0))):.17g}",
    }


def shell_descent(
    gram: np.ndarray,
    shell: np.ndarray,
    signed_candidates: np.ndarray,
    *,
    restarts: int,
    seed: int,
    random_steps: int,
    maximum_sweeps: int,
) -> dict[str, object]:
    """Coordinate-descent through unimodular short-vector bases."""

    rng = random.Random(seed)
    order_rng = np.random.default_rng(seed)

    def valid_candidates(basis: np.ndarray, column: int) -> np.ndarray:
        coordinates = signed_candidates @ integral_inverse(basis).T
        return np.flatnonzero(np.abs(coordinates[:, column]) == 1)

    def descend(basis: np.ndarray) -> tuple[np.ndarray, float, int]:
        score, _ = shell_statistics(gram, basis, shell)
        accepted = 0
        for _sweep in range(maximum_sweeps):
            changed = False
            for column_value in order_rng.permutation(17):
                column = int(column_value)
                indices = valid_candidates(basis, column)
                candidates = signed_candidates[indices]
                images = shell @ basis.T
                heights = np.einsum("ij,jk,ik->i", images, gram, images)
                differences = candidates - basis[:, column]
                mixed = (images @ gram) @ differences.T
                norms = np.einsum("ij,jk,ik->i", differences, gram, differences)
                coefficient = shell[:, column]
                trial_heights = (
                    heights[:, None]
                    + 2.0 * coefficient[:, None] * mixed
                    + (coefficient * coefficient)[:, None] * norms[None, :]
                )
                means = np.mean(trial_heights, axis=0)
                scores = np.mean(trial_heights * trial_heights, axis=0) / (means * means) - 1.0
                best_offset = int(np.argmin(scores))
                best_score = max(0.0, float(scores[best_offset]))
                if best_score < score - 1e-14:
                    basis[:, column] = candidates[best_offset]
                    score = best_score
                    accepted += 1
                    changed = True
            if not changed:
                break
        return basis, score, accepted

    best_basis = None
    best_score = float("inf")
    best_accepted = 0
    restart_scores = []
    for restart in range(restarts):
        basis = np.eye(17, dtype=np.int64)
        if restart:
            for _step in range(random_steps):
                column = rng.randrange(17)
                indices = valid_candidates(basis, column)
                basis[:, column] = signed_candidates[rng.choice(indices)]
        basis, score, accepted = descend(basis)
        restart_scores.append(math.sqrt(score))
        if score < best_score:
            best_basis = basis.copy()
            best_score = score
            best_accepted = accepted

    assert best_basis is not None
    if round(np.linalg.det(best_basis)) not in (-1, 1):
        raise RuntimeError("reported shell-aware target basis is not unimodular")
    _, statistics = shell_statistics(gram, best_basis, shell)
    return {
        "target_basis_columns": best_basis.tolist(),
        "accepted_descent_moves_in_best_restart": best_accepted,
        "restart_cv_min_median_max": [
            f"{min(restart_scores):.17g}",
            f"{float(np.median(restart_scores)):.17g}",
            f"{max(restart_scores):.17g}",
        ],
        **statistics,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("labels", nargs="*", default=["rank29", "curve245-negative-control"])
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--restarts", type=int, default=16)
    parser.add_argument("--random-steps", type=int, default=80)
    parser.add_argument("--maximum-sweeps", type=int, default=25)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    source = json.loads(args.input.read_text())
    grams = {
        record["label"]: np.array(record["core_lll_gram"], dtype=float)
        for record in source["curves"]
    }
    _control, control_gram = curve245_control(args.digits)
    grams["curve245-negative-control"] = control_gram
    _r17, signed_shell = exact_r17_reduced_form_and_shell()
    shell = signed_shell[: len(signed_shell) // 2]
    bounds = {
        "rank28": 75.0,
        "rank29": 75.0,
        "curve273": 75.0,
        "curve302": 82.0,
        "curve245-negative-control": 35.0,
    }
    results = []
    for index, label in enumerate(args.labels):
        if label not in grams:
            raise SystemExit(f"unknown label: {label}")
        candidates = numerical_short_vectors(grams[label], bounds[label])
        results.append(
            {
                "label": label,
                "target_short_vector_bound": bounds[label],
                "signed_target_candidates": len(candidates),
                **shell_descent(
                    grams[label],
                    shell,
                    candidates,
                    restarts=args.restarts,
                    seed=1_311_000 + 101 * index,
                    random_steps=args.random_steps,
                    maximum_sweeps=args.maximum_sweeps,
                ),
            }
        )

    result_by_label = {record["label"]: record for record in results}
    calibration = None
    if "rank29" in result_by_label and "curve245-negative-control" in result_by_label:
        positive = float(result_by_label["rank29"]["coefficient_of_variation"])
        control = float(
            result_by_label["curve245-negative-control"]["coefficient_of_variation"]
        )
        ratio = control / positive
        calibration = {
            "known_r17_positive": "rank29",
            "known_non_r17_parent_control": "curve245-negative-control",
            "control_cv_over_rank29_cv": f"{ratio:.17g}",
            "declared_separation_ratio": 1.20,
            "discrimination_result": (
                "FAILS_NEGATIVE_CONTROL"
                if ratio < 1.20
                else "SEPARATES_CONTROL_AT_DECLARED_THRESHOLD"
            ),
        }

    payload = {
        "schema": "record-rank17-shell-embedding-search-v1",
        "status": "bounded numerical shell-aware search",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "unoriented_r17_minimal_lines": len(shell),
        "objective": (
            "Minimize the coefficient of variation of all 1,311 mapped R17 "
            "minimal-vector heights over unimodular bases drawn from a bounded "
            "target short-vector cloud."
        ),
        "warning": (
            "A low value is evidence only after separation from the exact curve-245 "
            "negative control; bounded coordinate descent supplies no optimum proof."
        ),
        "parameters": {
            "restarts": args.restarts,
            "random_steps": args.random_steps,
            "maximum_sweeps": args.maximum_sweeps,
        },
        "calibration": calibration,
        "results": results,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(f"PASS|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"WROTE|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()

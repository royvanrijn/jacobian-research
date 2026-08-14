#!/usr/bin/env python3
"""Replay Nagao's four unpublished-rank finalists in the rank-21 family.

Section 6 of Nagao's 1994 paper lists four specializations that survived his
large-score search before the published rank-21 curve.  The paper uses half
the parameter of this repository's symmetric constructor, so the exact
constructor parameters are

``1393/108, 1649/6, 6629/174, 8057/438``.

This standalone pass recomputes conductor and root number, records the
selection-biased Nagao score only as provenance, and performs staged exact
quartic searches at heights 50k, 250k, and 1m.  Numerical height ranks select
the next stage but never certify rank.  Every stable rank at least 18 is sent
to saturation and the finite-reduction independence engine; only that exact
engine can create a rank lower-bound checkpoint or target hit.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import PRIMARY_SOURCE, RANK21_CONSTRUCTION, short_jacobian_coefficients
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank21_unbiased import (
    ExactScoredCandidate,
    PrefilterCandidate,
    build_residue_tables,
    exact_score_candidates,
    finite_reduction_certificate,
    point_pool_record,
    pool_priority,
    rank_one_pool,
    residue_score,
    search_one_pool,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
PAPER_PARAMETERS = (Q(1393, 216), Q(1649, 12), Q(6629, 348), Q(8057, 876))
CONSTRUCTOR_PARAMETERS = tuple(2 * value for value in PAPER_PARAMETERS)
EXPECTED_CONSTRUCTOR_PARAMETERS = (Q(1393, 108), Q(1649, 6), Q(6629, 174), Q(8057, 438))
STAGE_HEIGHTS = (50_000, 250_000, 1_000_000)
STAGE_KEEPS = (3, 2)
STAGE_TIMEOUTS = (10.0, 30.0, 90.0)
CHECKPOINT_RANK = 18
TARGET_RANK = 21
REPOSITORY = Path(__file__).resolve().parents[2]
CALIBRATION_ARTIFACT = (
    REPOSITORY / "artifacts/generated-results/elliptic_nagao_rank21_neighborhood.json"
)
CALIBRATION_ARTIFACT_SHA256 = (
    "7d59fe9a91c0f3e46604794e8931ae27e26eeea1ebf252176dffd6be8d6010fe"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_historical_finalists.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_source_convention() -> dict[str, Any]:
    if CONSTRUCTOR_PARAMETERS != EXPECTED_CONSTRUCTOR_PARAMETERS:
        raise AssertionError("the historical factor-two convention changed")
    if sha256_file(CALIBRATION_ARTIFACT) != CALIBRATION_ARTIFACT_SHA256:
        raise AssertionError("the published-rank21 calibration artifact changed")
    data = json.loads(CALIBRATION_ARTIFACT.read_text(encoding="utf-8"))
    calibration = data["published_record_calibration"]
    if (
        calibration["constructor_parameter"] != "14721/188"
        or not calibration["factor_two_identity_checked"]
        or not calibration["pari_replay_matches_printed_model_and_conductor"]
    ):
        raise AssertionError("the published parameter convention lost calibration")
    return data


def scored_candidates(args: argparse.Namespace) -> tuple[ExactScoredCandidate, ...]:
    tables = build_residue_tables(200)
    prefilter = []
    for parameter in CONSTRUCTOR_PARAMETERS:
        score, good, bad = residue_score(
            parameter.numerator, parameter.denominator, tables
        )
        prefilter.append(
            PrefilterCandidate(
                parameter.numerator,
                parameter.denominator,
                score,
                good,
                bad,
            )
        )
    return exact_score_candidates(
        tuple(prefilter),
        cutoff=2_000,
        batch_size=len(prefilter),
        timeout=args.score_timeout,
        stack_bytes=args.stack_bytes,
    )


def conductor_records(
    candidates: Sequence[ExactScoredCandidate], args: argparse.Namespace
) -> tuple[dict[Fraction, dict[str, Any]], list[dict[str, Any]]]:
    completed: dict[Fraction, dict[str, Any]] = {}
    errors = []
    for candidate in candidates:
        try:
            record = minimal_curve_data(
                short_jacobian_coefficients(
                    RANK21_CONSTRUCTION, candidate.parameter
                ),
                timeout=args.conductor_timeout,
                stack_bytes=args.stack_bytes,
                local_primes=(2, 3, 5, 7, 13, 17, 23),
            )
            record["below_strict_log_conductor_target"] = (
                Decimal(record["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
            completed[candidate.parameter] = record
        except (subprocess.TimeoutExpired, RuntimeError, ValueError) as error:
            errors.append(
                {
                    "constructor_parameter": rational_to_string(candidate.parameter),
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "error": str(error)[:500],
                }
            )
    return completed, errors


@dataclass(frozen=True)
class StageCheckpoint:
    candidate: ExactScoredCandidate
    pool: Any
    rank: dict[str, Any]


def staged_point_search(
    candidates: Sequence[ExactScoredCandidate], args: argparse.Namespace
) -> tuple[list[dict[str, Any]], dict[str, StageCheckpoint]]:
    retained = tuple(candidates)
    best: dict[str, StageCheckpoint] = {}
    stages = []
    for stage_index, (height, timeout) in enumerate(
        zip(STAGE_HEIGHTS, STAGE_TIMEOUTS), start=1
    ):
        pools = []
        ranked = []
        for candidate in retained:
            pool = search_one_pool(
                candidate,
                height_bound=height,
                timeout=timeout,
                stack_bytes=args.stack_bytes,
            )
            rank = rank_one_pool(
                pool,
                precisions=(72, 120),
                timeout=args.height_timeout,
                stack_bytes=args.stack_bytes,
            )
            pools.append(pool)
            ranked.append((pool, rank))
            if rank.get("status") == "completed":
                prior = best.get(candidate.identifier)
                if prior is None or int(rank["stable_numerical_rank"]) >= int(
                    prior.rank["stable_numerical_rank"]
                ):
                    best[candidate.identifier] = StageCheckpoint(candidate, pool, rank)
        ranked.sort(key=lambda item: pool_priority(item[0], item[1]))
        keep = len(ranked) if stage_index == len(STAGE_HEIGHTS) else STAGE_KEEPS[stage_index - 1]
        retained = tuple(pool.candidate for pool, _ in ranked[:keep])
        stages.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height,
                "population_searched": len(ranked),
                "completed_point_searches": sum(
                    pool.status == "completed" for pool, _ in ranked
                ),
                "point_search_timeouts": sum(
                    pool.status == "timeout" for pool, _ in ranked
                ),
                "point_search_errors": sum(
                    pool.status == "error" for pool, _ in ranked
                ),
                "ranked_population": [
                    point_pool_record(
                        pool,
                        rank,
                        include_points=(stage_index == len(STAGE_HEIGHTS)),
                    )
                    for pool, rank in ranked
                ],
                "retained_constructor_parameters": [
                    rational_to_string(candidate.parameter) for candidate in retained
                ],
            }
        )
    return stages, best


def exact_checkpoints(
    best: dict[str, StageCheckpoint],
    conductors: dict[Fraction, dict[str, Any]],
    args: argparse.Namespace,
) -> list[dict[str, Any]]:
    records = []
    for checkpoint in sorted(
        best.values(),
        key=lambda value: (
            -int(value.rank["stable_numerical_rank"]),
            value.candidate.parameter,
        ),
    ):
        numerical_rank = int(checkpoint.rank["stable_numerical_rank"])
        if numerical_rank < CHECKPOINT_RANK:
            continue
        conductor = conductors.get(checkpoint.candidate.parameter)
        if conductor is None:
            continue
        certificate = finite_reduction_certificate(
            checkpoint.pool,
            checkpoint.rank,
            saturation_timeout=args.saturation_timeout,
            certificate_prime_bound=args.certificate_prime_bound,
            stack_bytes=args.stack_bytes,
        )
        exact_rank = certificate["certified_algebraic_rank_lower_bound"]
        target_hit = bool(
            exact_rank is not None
            and exact_rank >= TARGET_RANK
            and conductor["below_strict_log_conductor_target"]
        )
        records.append(
            {
                "constructor_parameter": rational_to_string(
                    checkpoint.candidate.parameter
                ),
                "deepest_completed_height": checkpoint.pool.height_bound,
                "stable_numerical_rank": numerical_rank,
                "conductor": conductor,
                "exact_rank_certificate": certificate,
                "target_rank21_under_log_conductor_hit": target_hit,
            }
        )
    return records


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    validate_source_convention()
    candidates = scored_candidates(args)
    by_parameter = {candidate.parameter: candidate for candidate in candidates}
    if set(by_parameter) != set(CONSTRUCTOR_PARAMETERS):
        raise AssertionError("the historical finalist population changed")
    conductors, conductor_errors = conductor_records(candidates, args)
    stages, best = staged_point_search(candidates, args)
    checkpoints = exact_checkpoints(best, conductors, args)
    hits = [
        record
        for record in checkpoints
        if record["target_rank21_under_log_conductor_hit"]
    ]
    candidate_records = []
    for paper_parameter, constructor_parameter in zip(
        PAPER_PARAMETERS, CONSTRUCTOR_PARAMETERS
    ):
        candidate = by_parameter[constructor_parameter]
        candidate_records.append(
            {
                "paper_parameter_t": rational_to_string(paper_parameter),
                "constructor_parameter_T": rational_to_string(constructor_parameter),
                "factor_two_checked": constructor_parameter == 2 * paper_parameter,
                "selection_status": (
                    "listed by Nagao as a high-score survivor; no rank claimed"
                ),
                "residue_table_b200_score": candidate.prefilter.residue_score_b200,
                "exact_pari_b2000_score": candidate.exact_score_b2000,
                "conductor": conductors.get(constructor_parameter),
            }
        )
    return {
        "schema_version": 1,
        "status": (
            "bounded historical-finalist replay; numerical height ranks are "
            "triage only and exact finite reductions certify checkpoints"
        ),
        "primary_source": PRIMARY_SOURCE,
        "source_location": "Nagao 1994, section 6, page 217",
        "source_parameter_convention": {
            "paper": "E_t",
            "constructor": "q(X-T)q(X+T)",
            "identity": "T=2t",
            "published_rank21_calibration_sha256": CALIBRATION_ARTIFACT_SHA256,
        },
        "target": {
            "rank_at_least": TARGET_RANK,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "certified_hits": hits,
        },
        "historical_finalists": candidate_records,
        "conductor_errors": conductor_errors,
        "point_search_stages": stages,
        "exact_checkpoints_stable_numerical_rank_at_least_18": checkpoints,
        "bounds_and_caveats": {
            "stage_heights": list(STAGE_HEIGHTS),
            "stage_keeps": list(STAGE_KEEPS),
            "stage_timeouts_seconds": list(STAGE_TIMEOUTS),
            "conductor_timeout_seconds": args.conductor_timeout,
            "height_timeout_seconds": args.height_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "stack_bytes": args.stack_bytes,
            "historical_score_selection_is_biased": True,
            "bounded_search_is_not_a_rank_upper_bound": True,
            "all_subprocesses_synchronous_with_finite_timeouts": True,
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=90.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_rank21_historical_finalists.json"
        ),
    )
    args = parser.parse_args()
    if any(
        not 0 < value <= 120
        for value in (
            args.score_timeout,
            args.conductor_timeout,
            args.height_timeout,
            args.saturation_timeout,
        )
    ):
        raise SystemExit("all subprocess timeouts must lie in (0,120]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}", flush=True)
    for record in artifact["exact_checkpoints_stable_numerical_rank_at_least_18"]:
        certificate = record["exact_rank_certificate"]
        print(
            f"T={record['constructor_parameter']} numerical="
            f"{record['stable_numerical_rank']} exact="
            f"{certificate['certified_algebraic_rank_lower_bound']} target="
            f"{record['target_rank21_under_log_conductor_hit']}",
            flush=True,
        )


if __name__ == "__main__":
    main()

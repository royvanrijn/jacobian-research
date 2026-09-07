#!/usr/bin/env python3
"""Terminal direct H=10^6 Fermigier published-direction pair screen.

All 220 cross-direction genus-3 products from the exact published quartic
preimages P6 and P13--P22 are searched once by PARI ``hyperellratpoints`` at
naive rational height one million.  The exact record fiber remains a positive
calibration in every completed search.  Every parameter already present in
the Fermigier artifacts or in the terminal H=50000 pair stage is excluded.

A product square is only a necessary condition.  Candidates survive only if
both quartic factors are individually rational squares, both forced points
are checked exactly on the specialized quartic, neither point is generic,
and the two based Jacobian classes remain distinct.  Conductor comes before
the H=50000 specialized point/rank screen.  Stable numerical rank at least 21
triggers exact saturation and finite-reduction certification.

This script performs one direct pass, checkpoints after every pair, never
retries a subprocess, and stops at H=10^6.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any

import sympy as sp

from ek_k3 import rational_to_string
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_LABELS,
    EXPECTED_PAIR_COUNT,
    EXPECTED_PRIOR_PARAMETER_COUNT,
    EXPECTED_PRIOR_PARAMETER_SHA256,
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    SPECIALIZATION_HEIGHT,
    T0,
    TARGET_LOG_CONDUCTOR,
    aggregate_candidates,
    exact_slices,
    finalized_candidate_record,
    load_primary,
    pair_identifier,
    pair_population,
    polynomial_coprime,
    prior_fermigier_parameters,
    published_accidentals,
    published_preimage_digest,
    rational_digest,
    search_pair,
    sha256_bytes,
    sha256_file,
    triage_specialization,
)
from search_fermigier_published_pair_fiber_products_h50000 import (
    EXPECTED_H50000_RESULT_SHA256,
    h50000_result_digest,
)
from search_fermigier_rank22_accidental_slices import conductor_probe


HEIGHT = 1_000_000
H50000_ARTIFACT = (
    "artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products_h50000.json"
)
EXPECTED_TERMINAL_PRIOR_COUNT = 593
EXPECTED_TERMINAL_PRIOR_SHA256 = (
    "a4d06e4662d2e30c1a0f8873f91d8d348dae10f2abaffce88dcc0f480cfeede0"
)


def parse_precisions(value: str) -> tuple[int, ...]:
    values = tuple(int(part) for part in value.split(",") if part)
    if len(values) < 2 or tuple(sorted(set(values))) != values:
        raise argparse.ArgumentTypeError("provide increasing distinct precisions")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=HEIGHT)
    parser.add_argument("--pair-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--finalize-disproportionate-partial",
        action="store_true",
        help="finalize the checkpoint after the declared direct lane is stopped",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_published_pair_fiber_products_h1000000.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def terminal_result_digest(pair_rows: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in pair_rows:
        search = row["search"]
        digest.update(
            (
                f"{row['pair_id']}|{search['product_polynomial_sha256']}|"
                f"{search['search']['status']}|"
                f"{search['search'].get('signed_point_count')}|"
                f"{search['qualifying_new_parameter_count']}|"
                f"{row['record_T0_positive_calibration_count']}\n"
            ).encode()
        )
    return digest.hexdigest()


def main() -> None:
    args = build_parser().parse_args()
    if args.height != HEIGHT:
        raise SystemExit("this terminal direct tranche is pinned at H=1000000")
    if min(
        args.pair_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.pair_timeout > 60:
        raise SystemExit("pair timeout may not exceed 60 seconds")

    root = Path(__file__).resolve().parents[2]
    if args.finalize_disproportionate_partial:
        artifact = json.loads(args.output.read_bytes())
        rows = artifact["pair_searches"]
        if len(rows) != 2 or any(
            row["search"]["search"]["status"] != "timeout" for row in rows
        ):
            raise AssertionError("the disproportionate direct checkpoint changed")
        artifact["status"] = (
            "stopped as computationally disproportionate after two consecutive "
            "20-second direct H=1000000 pair timeouts; no retry"
        )
        artifact["execution"].update(
            {
                "phase": "stopped-computationally-disproportionate",
                "stopping_rule": "two consecutive first-pair-order timeouts",
                "owned_processes_remaining": 0,
            }
        )
        artifact["outcome"] = {
            "declared_pair_count": EXPECTED_PAIR_COUNT,
            "pairs_attempted": 2,
            "pairs_completed": 0,
            "pairs_timed_out": 2,
            "stopped_as_computationally_disproportionate": True,
            "conductor_calls": 0,
            "rank_triage_calls": 0,
        }
        artifact["target"]["reason"] = (
            "direct H=1000000 was disproportionate; deterministic transformed "
            "charts are the declared fallback"
        )
        artifact["script_sha256"] = sha256_file(Path(__file__).resolve())
        artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_artifact(args.output, artifact)
        return
    artifact_dir = root / "artifacts" / "generated-results"
    primary_path = artifact_dir / PRIMARY_ARTIFACT
    primary, primary_raw = load_primary(primary_path)
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published-preimage population changed")
    accidentals = published_accidentals(primary)
    slices = exact_slices(accidentals)
    pairs = pair_population(slices)

    h50000_path = root / H50000_ARTIFACT
    h50000 = json.loads(h50000_path.read_bytes())
    h50000_digest = h50000_result_digest(h50000["pair_searches"])
    if h50000_digest != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 pair-result digest changed")
    base_prior, base_prior_record = prior_fermigier_parameters(
        artifact_dir, args.output
    )
    if (
        len(base_prior) != EXPECTED_PRIOR_PARAMETER_COUNT
        or base_prior_record["prior_parameter_sha256"]
        != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact base prior-Fermigier population changed")
    h50000_seen = {
        abs(Fraction(incidence["canonical_parameter_t"]))
        for row in h50000["pair_searches"]
        for incidence in row["search"]["incidences"]
    }
    prior_parameters = base_prior | h50000_seen
    terminal_prior_sha256 = rational_digest(sorted(prior_parameters))
    if (
        len(prior_parameters) != EXPECTED_TERMINAL_PRIOR_COUNT
        or terminal_prior_sha256 != EXPECTED_TERMINAL_PRIOR_SHA256
    ):
        raise AssertionError("the terminal exact prior population changed")

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "in-progress bounded direct H=1000000 pair screen; exact factor-square "
            "and quartic checks, numerical height ranks only unless certified"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "hit": False,
        },
        "source": {
            "published_preimage_sha256": published_preimage_digest(primary),
            "published_accidental_labels": list(EXPECTED_LABELS),
            "primary_artifact_sha256_observed": sha256_bytes(primary_raw),
            "H50000_artifact": H50000_ARTIFACT,
            "H50000_artifact_sha256_observed": sha256_file(h50000_path),
            "H50000_exact_pair_result_sha256": h50000_digest,
            "record_parameter_t": rational_to_string(T0),
        },
        "prior_decontamination": {
            "base_prior_parameter_count": len(base_prior),
            "base_prior_parameter_sha256": base_prior_record[
                "prior_parameter_sha256"
            ],
            "H50000_seen_parameters": [
                rational_to_string(value) for value in sorted(h50000_seen)
            ],
            "terminal_prior_parameter_count": len(prior_parameters),
            "terminal_prior_parameter_sha256": terminal_prior_sha256,
        },
        "parameters": {
            "pair_height": args.height,
            "pair_timeout_seconds": args.pair_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "specialization_height": SPECIALIZATION_HEIGHT,
            "specialization_timeout_seconds": args.specialization_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": list(args.precisions),
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
            "terminal_bound": True,
            "checkpoint_after_every_pair": True,
        },
        "execution": {
            "phase": "pair-search-in-progress",
            "pairs_completed_or_attempted": 0,
        },
        "pair_searches": [],
        "candidates": [],
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    write_artifact(args.output, artifact)

    qualifying_all: list[tuple[str, dict[str, Any]]] = []
    started = time.monotonic()
    for left, right in pairs:
        if not polynomial_coprime(left.coefficients, right.coefficients):
            raise AssertionError("two distinct source slices shared a polynomial factor")
        pair_id = pair_identifier(left, right)
        search, qualifying = search_pair(
            left,
            right,
            height=args.height,
            timeout=args.pair_timeout,
            stack_bytes=args.stack_bytes,
            prior_parameters=prior_parameters,
        )
        calibrations = [
            incidence
            for incidence in search["incidences"]
            if incidence["canonical_parameter_t"] == str(abs(T0))
            and incidence["classification"] == "record-fiber-excluded"
            and incidence["left_factor_is_square"]
            and incidence["right_factor_is_square"]
            and len(incidence.get("exact_forced_quartic_points", [])) == 2
        ]
        if search["search"]["status"] == "completed" and len(calibrations) != 1:
            raise AssertionError(f"{pair_id}: completed search lost T0 calibration")
        artifact["pair_searches"].append(
            {
                "pair_id": pair_id,
                "left_source_label": left.accidental_label,
                "left_slice_id": left.identifier,
                "right_source_label": right.accidental_label,
                "right_slice_id": right.identifier,
                "distinct_source_labels": True,
                "factor_polynomials_coprime": True,
                "search": search,
                "record_T0_positive_calibration_count": len(calibrations),
            }
        )
        qualifying_all.extend((pair_id, incidence) for incidence in qualifying)
        artifact["execution"]["pairs_completed_or_attempted"] = len(
            artifact["pair_searches"]
        )
        artifact["execution"]["last_pair_id"] = pair_id
        artifact["execution"]["wall_seconds_so_far"] = time.monotonic() - started
        write_artifact(args.output, artifact)

    aggregated = aggregate_candidates(qualifying_all)
    candidates = [
        finalized_candidate_record(candidate)
        for _, candidate in sorted(aggregated.items())
    ]
    candidates = [
        candidate
        for candidate in candidates
        if candidate["distinct_published_source_direction_count"] >= 2
        and candidate["distinct_group_pullback_classes_modulo_inverse"] >= 2
    ]
    artifact["candidates"] = candidates
    artifact["execution"]["phase"] = "conductor-first"
    write_artifact(args.output, artifact)
    for candidate in candidates:
        parameter = Fraction(candidate["parameter_t"])
        candidate["conductor_probe"] = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        write_artifact(args.output, artifact)
        if candidate["conductor_probe"].get("below_strict_log_conductor_target"):
            try:
                candidate["rank_triage"] = triage_specialization(
                    candidate,
                    search_timeout=args.specialization_timeout,
                    height_timeout=args.height_timeout,
                    precisions=args.precisions,
                    stack_bytes=args.stack_bytes,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except subprocess.TimeoutExpired as error:
                candidate["rank_triage"] = {
                    "status": "timeout-no-retry",
                    "error": str(error)[:1000],
                }
            except (RuntimeError, AssertionError, ValueError) as error:
                candidate["rank_triage"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            write_artifact(args.output, artifact)

    classifications = Counter(
        incidence["classification"]
        for row in artifact["pair_searches"]
        for incidence in row["search"]["incidences"]
    )
    rank_records = [
        candidate["rank_triage"]
        for candidate in candidates
        if "rank_triage" in candidate
        and "full_pool_stable_numerical_rank" in candidate["rank_triage"]
    ]
    artifact["outcome"] = {
        "declared_pair_count": EXPECTED_PAIR_COUNT,
        "pairs_attempted": len(artifact["pair_searches"]),
        "pairs_completed": sum(
            row["search"]["search"]["status"] == "completed"
            for row in artifact["pair_searches"]
        ),
        "pairs_timed_out_or_errored": sum(
            row["search"]["search"]["status"] != "completed"
            for row in artifact["pair_searches"]
        ),
        "pair_wall_seconds": time.monotonic() - started,
        "record_T0_calibrated_pairs": sum(
            row["record_T0_positive_calibration_count"] == 1
            for row in artifact["pair_searches"]
        ),
        "incidence_classification_counts": dict(sorted(classifications.items())),
        "genuinely_new_double_forced_fibers": len(candidates),
        "completed_conductors": sum(
            candidate.get("conductor_probe", {}).get("status") == "completed"
            for candidate in candidates
        ),
        "subtarget_conductors": sum(
            candidate.get("conductor_probe", {}).get(
                "below_strict_log_conductor_target"
            )
            is True
            for candidate in candidates
        ),
        "rank_triage_count": len(rank_records),
        "maximum_stable_numerical_rank": max(
            (
                record["full_pool_stable_numerical_rank"]
                for record in rank_records
            ),
            default=None,
        ),
        "exact_terminal_pair_result_sha256": terminal_result_digest(
            artifact["pair_searches"]
        ),
    }
    artifact["target"]["hit"] = any(
        candidate.get("rank_triage", {})
        .get("finite_reduction_attempt", {})
        .get("certified_algebraic_rank_lower_bound", 0)
        >= 21
        and candidate.get("conductor_probe", {}).get(
            "below_strict_log_conductor_target"
        )
        for candidate in candidates
    )
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no new H=1000000 double-forced subtarget fiber received an exact rank-21 certificate"
        )
    artifact["status"] = (
        "completed bounded direct H=1000000 pair screen; exact factor-square and "
        "quartic checks, conductor-first triage, numerical height ranks unless certified"
    )
    artifact["execution"]["phase"] = "complete"
    artifact["execution"]["wall_seconds"] = time.monotonic() - started
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()

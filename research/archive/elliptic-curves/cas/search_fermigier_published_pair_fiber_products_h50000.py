#!/usr/bin/env python3
"""One-pass H=50000 extension of the Fermigier pairwise fiber products.

This is deliberately separate from the pinned H=5000 pilot.  It searches all
220 cross-direction products exactly once at H=50000, which includes the
record parameter T0=39508/39.  Every pair must therefore recover T0 with both
quartic factors individually square; that record fiber is an exact positive
calibration and is decontaminated before candidate aggregation.

Only genuinely new fibers forcing at least two distinct published source
directions receive conductor attempts.  H=50000 specialized quartic/rank
triage is reserved for completed conductors with log(N)<182.72.  No pair is
retried and there is no search beyond H=50000 in this extension.
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
    pair_identifier,
    pair_population,
    pair_result_digest,
    polynomial_coprime,
    prior_fermigier_parameters,
    published_accidentals,
    published_preimage_digest,
    search_pair,
    sha256_bytes,
    sha256_file,
    triage_specialization,
)
from search_fermigier_rank22_accidental_slices import conductor_probe


HEIGHT = 50_000
PILOT_SCRIPT = "elliptic-curves/cas/search_fermigier_published_pair_fiber_products.py"
PILOT_ARTIFACT = (
    "artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products.json"
)
EXPECTED_H5000_RESULT_SHA256 = (
    "80413701447b6468a826fa2185528da74057faa02619bdf02f237e7efb8b1b8b"
)
EXPECTED_H50000_RESULT_SHA256 = (
    "dea8b716c5aec56817a172afd6e894e7748aaddc482a2d29c0a3360abe55bf4b"
)
BOUNDED_EXECUTION_SCRIPT_SHA256 = (
    "61e2435221d523a9c7b3af28b3aad8be8477c57b44dd9709b7f86ddf709030d7"
)
PRE_NORMALIZATION_ARTIFACT_SHA256 = (
    "592c6f22abd5335d412e551d6397272f66272a628d9490d21ab1041153f1e56e"
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
    parser.add_argument("--pair-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--normalize-existing-artifact",
        action="store_true",
        help="refresh only stable dependency metadata; perform no bounded searches",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_published_pair_fiber_products_h50000.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def h50000_result_digest(pair_rows: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in pair_rows:
        search = row["search"]
        digest.update(
            (
                f"{row['pair_id']}|{search['product_polynomial_sha256']}|"
                f"{search['search'].get('signed_point_count')}|"
                f"{search['qualifying_new_parameter_count']}|"
                f"{row['record_T0_positive_calibration_count']}\n"
            ).encode()
        )
    return digest.hexdigest()


def normalize_existing_artifact(args: argparse.Namespace, root: Path) -> None:
    """Repin stable inputs without rerunning any pair/conductor/rank process."""

    if sha256_file(args.output) != PRE_NORMALIZATION_ARTIFACT_SHA256:
        raise AssertionError("the pre-normalization H=50000 artifact changed")
    artifact = json.loads(args.output.read_bytes())
    if artifact.get("script_sha256") != BOUNDED_EXECUTION_SCRIPT_SHA256:
        raise AssertionError("the bounded execution script provenance changed")
    exact_result_sha256 = h50000_result_digest(artifact["pair_searches"])
    if exact_result_sha256 != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 result digest changed")

    artifact_dir = root / "artifacts" / "generated-results"
    primary_path = artifact_dir / PRIMARY_ARTIFACT
    primary_raw = primary_path.read_bytes()
    primary = json.loads(primary_raw)
    preimage_sha256 = published_preimage_digest(primary)
    if preimage_sha256 != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published-preimage population changed")
    pilot_script_path = root / PILOT_SCRIPT
    pilot_artifact_path = root / PILOT_ARTIFACT
    pilot = json.loads(pilot_artifact_path.read_bytes())
    pilot_result_sha256 = pair_result_digest(pilot["pair_searches"])
    if pilot_result_sha256 != EXPECTED_H5000_RESULT_SHA256:
        raise AssertionError("the exact H=5000 result digest changed")
    prior_parameters, prior_record = prior_fermigier_parameters(
        artifact_dir, args.output
    )
    if (
        len(prior_parameters) != EXPECTED_PRIOR_PARAMETER_COUNT
        or prior_record["prior_parameter_sha256"]
        != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact prior-Fermigier parameter population changed")

    source = artifact["source"]
    source.update(
        {
            "finalized_accidental_artifact_sha256_observed": sha256_bytes(
                primary_raw
            ),
            "published_accidental_preimage_sha256": preimage_sha256,
            "H5000_pilot_script_sha256_observed": sha256_file(
                pilot_script_path
            ),
            "H5000_pilot_artifact_sha256_observed": sha256_file(
                pilot_artifact_path
            ),
            "H5000_exact_pair_result_sha256": pilot_result_sha256,
        }
    )
    # Remove whole-file fields which had previously acted as replay guards.
    source.pop("finalized_accidental_artifact_sha256", None)
    source.pop("finalized_accidental_script_sha256", None)
    source.pop("H5000_pilot_script_sha256", None)
    source.pop("H5000_pilot_artifact_sha256", None)
    artifact["prior_decontamination"] = prior_record
    artifact["outcome"]["exact_H50000_pair_result_sha256"] = exact_result_sha256
    current_script_sha256 = sha256_file(Path(__file__).resolve())
    artifact["execution_provenance"] = {
        "bounded_execution_script_sha256": BOUNDED_EXECUTION_SCRIPT_SHA256,
        "pre_normalization_artifact_sha256": PRE_NORMALIZATION_ARTIFACT_SHA256,
        "stable_replay_script_sha256": current_script_sha256,
        "metadata_normalization_only": True,
        "bounded_searches_rerun_during_normalization": 0,
    }
    artifact["script_sha256"] = current_script_sha256
    artifact["metadata_normalized_at_utc"] = datetime.now(timezone.utc).isoformat()
    artifact["reproducing_command"] = (
        f"PYTHONPATH=elliptic-curves/cas .venv/bin/python "
        f"elliptic-curves/cas/{Path(__file__).name}"
    )
    write_artifact(args.output, artifact)


def main() -> None:
    args = build_parser().parse_args()
    if args.height != HEIGHT:
        raise SystemExit("this extension is pinned at H=50000")
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
    if args.normalize_existing_artifact:
        normalize_existing_artifact(args, root)
        return
    artifact_dir = root / "artifacts" / "generated-results"
    pilot_script_path = root / PILOT_SCRIPT
    pilot_artifact_path = root / PILOT_ARTIFACT
    pilot = json.loads(pilot_artifact_path.read_bytes())
    if (
        pilot["outcome"]["pilot_pairs_completed"] != EXPECTED_PAIR_COUNT
        or pilot["outcome"]["productive_pilot_pairs"] != 0
    ):
        raise AssertionError("the H=5000 pilot checkpoint changed")
    pilot_result_sha256 = pair_result_digest(pilot["pair_searches"])
    if pilot_result_sha256 != EXPECTED_H5000_RESULT_SHA256:
        raise AssertionError("the exact H=5000 pair-result digest changed")

    primary_path = artifact_dir / PRIMARY_ARTIFACT
    primary_raw = primary_path.read_bytes()
    primary = json.loads(primary_raw)
    primary_script_path = (
        root / "elliptic-curves" / "cas" / "search_fermigier_rank22_accidental_slices.py"
    )
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published-preimage population changed")
    accidentals = published_accidentals(primary)
    slices = exact_slices(accidentals)
    pairs = pair_population(slices)
    prior_parameters, prior_record = prior_fermigier_parameters(
        artifact_dir, args.output
    )
    if (
        prior_record["unique_prior_parameter_count"] != EXPECTED_PRIOR_PARAMETER_COUNT
        or prior_record["prior_parameter_sha256"]
        != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact prior-Fermigier parameter population changed")

    pair_rows = []
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
        record_calibrations = [
            incidence
            for incidence in search["incidences"]
            if incidence["canonical_parameter_t"] == str(abs(T0))
            and incidence["classification"] == "record-fiber-excluded"
            and incidence["left_factor_is_square"]
            and incidence["right_factor_is_square"]
            and len(incidence.get("exact_forced_quartic_points", [])) == 2
            and all(
                point["exact_membership_checked"]
                for point in incidence["exact_forced_quartic_points"]
            )
        ]
        if search["search"]["status"] == "completed" and len(record_calibrations) != 1:
            raise AssertionError(
                f"{pair_id}: completed H=50000 search did not recover T0 exactly once"
            )
        pair_rows.append(
            {
                "pair_id": pair_id,
                "left_source_label": left.accidental_label,
                "left_slice_id": left.identifier,
                "right_source_label": right.accidental_label,
                "right_slice_id": right.identifier,
                "distinct_source_labels": True,
                "factor_polynomials_coprime": True,
                "search": search,
                "record_T0_positive_calibration_count": len(record_calibrations),
            }
        )
        qualifying_all.extend((pair_id, incidence) for incidence in qualifying)
    pair_wall_seconds = time.monotonic() - started
    if len(pair_rows) != EXPECTED_PAIR_COUNT:
        raise AssertionError("the H=50000 extension did not attempt all 220 pairs")

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
    for candidate in candidates:
        parameter = Fraction(candidate["parameter_t"])
        candidate["conductor_probe"] = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
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

    rank_records = [
        candidate["rank_triage"]
        for candidate in candidates
        if "rank_triage" in candidate
        and "full_pool_stable_numerical_rank" in candidate["rank_triage"]
    ]
    classifications = Counter(
        incidence["classification"]
        for row in pair_rows
        for incidence in row["search"]["incidences"]
    )
    exact_H50000_result_sha256 = h50000_result_digest(pair_rows)
    if exact_H50000_result_sha256 != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 pair-result population changed")
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded all-pairs H=50000 extension; both factor-square and forced-point "
            "checks are exact, T0 is an exact positive calibration, conductors are "
            "PARI computations, and height ranks remain numerical without a certificate"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "hit": any(
                candidate.get("rank_triage", {})
                .get("finite_reduction_attempt", {})
                .get("certified_algebraic_rank_lower_bound", 0)
                >= 21
                and candidate["conductor_probe"].get(
                    "below_strict_log_conductor_target"
                )
                for candidate in candidates
            ),
        },
        "source": {
            "finalized_accidental_artifact": str(primary_path.relative_to(root)),
            "finalized_accidental_artifact_sha256": sha256_bytes(primary_raw),
            "finalized_accidental_script_sha256": sha256_file(primary_script_path),
            "published_accidental_preimage_sha256": published_preimage_digest(
                primary
            ),
            "H5000_pilot_script": PILOT_SCRIPT,
            "H5000_pilot_script_sha256": sha256_file(pilot_script_path),
            "H5000_pilot_artifact": PILOT_ARTIFACT,
            "H5000_pilot_artifact_sha256": sha256_file(pilot_artifact_path),
            "H5000_exact_pair_result_sha256": pilot_result_sha256,
            "record_parameter_t": rational_to_string(T0),
            "record_parameter_projective_height": max(
                abs(T0.numerator), T0.denominator
            ),
            "published_accidental_labels": list(EXPECTED_LABELS),
        },
        "prior_decontamination": prior_record,
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
        },
        "pair_searches": pair_rows,
        "candidates": candidates,
        "outcome": {
            "declared_pair_count": EXPECTED_PAIR_COUNT,
            "pairs_attempted": len(pair_rows),
            "pairs_completed": sum(
                row["search"]["search"]["status"] == "completed"
                for row in pair_rows
            ),
            "pairs_timed_out_or_errored": sum(
                row["search"]["search"]["status"] != "completed"
                for row in pair_rows
            ),
            "pair_wall_seconds": pair_wall_seconds,
            "record_T0_calibrated_pairs": sum(
                row["record_T0_positive_calibration_count"] == 1
                for row in pair_rows
            ),
            "incidence_classification_counts": dict(sorted(classifications.items())),
            "genuinely_new_double_forced_fibers": len(candidates),
            "completed_conductors": sum(
                candidate["conductor_probe"]["status"] == "completed"
                for candidate in candidates
            ),
            "subtarget_conductors": sum(
                candidate["conductor_probe"].get("below_strict_log_conductor_target")
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
            "exact_H50000_pair_result_sha256": exact_H50000_result_sha256,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no new H=50000 double-forced subtarget fiber received an exact rank-21 certificate"
        )
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()

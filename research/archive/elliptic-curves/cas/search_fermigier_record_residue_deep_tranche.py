#!/usr/bin/env python3
"""Leakage-controlled deep tranche from Fermigier's H=50000 CRT population.

The earlier record-residue experiment retained candidates by cumulative
``B=500`` or staged cumulative scores.  This script regenerates the complete
23,769-member exact population and selects a disjoint tranche without using
that cumulative score.  It splits the primes into a discovery band and two
held-forward bands, selects independently within four height strata, and
adds a small height-first tail which must lie in the upper quartile of both
held-forward bands.

Every selected specialization receives a PARI conductor attempt before any
point search.  Only completed conductors below the strict target proceed to
the exact-quartic ``H=50000`` search.  A specialization advances to
``H=250000`` and then ``H=1000000`` only after the preceding stage has found
a stable numerical height rank above the exact generic rank-12 seed.  These
height ranks are numerical triage evidence, not rank certificates.  A stable
rank at least 20 triggers a bounded saturation and exact finite-reduction
independence attempt.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
from dataclasses import dataclass
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
from typing import Any, Iterable, Sequence

from batch_rank_triage import TriageCandidate, run_stage, stable_rank
from ek_k3 import rational_to_string
from exhaustive_multiple_root_height import (
    DEFAULT_HEIGHT,
    ProjectiveCandidate,
    enumerate_projective_height,
    verify_population_locally,
)
from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from pari_bridge import minimal_curve_data, pari_version
from search_fermigier_rank22_accidental_slices import finite_reduction_attempt
from search_multiple_root_crt import DEFAULT_GROUP_ORDER, choose_groups, crt_classes
from search_record_residue_class import PrimeScoreTable, build_score_tables, score_rational
from search_extra_points import height_replay


TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXPECTED_POPULATION_COUNT = 23_769
EXPECTED_POPULATION_SHA256 = (
    "b3c6751977ca35febe9d4f1d974183226bd0d77ec34321ed9d59244cffe4f086"
)
EXPECTED_EXCLUSION_COUNT = 111
EXPECTED_SELECTED_COUNT = 48
EXPECTED_SELECTED_SHA256 = (
    "22452d49f569aee1bc2b0cf6a3206445271d82c704f1b997f5f19cd7c03b3c44"
)
EXPECTED_FERMIGIER_ACCIDENTAL_PARAMETER_SHA256 = (
    "e8410fbcba4491165fd114e86cff11a1eabc1373f74cebda8dbc856bbcf0045f"
)
BENCHMARK = NORMALIZED_RECORD_PARAMETER
SCORE_BOUND = 500
DISCOVERY_MAX_PRIME = 199
HELD_ONE_MIN_PRIME = 211
HELD_ONE_MAX_PRIME = 349
HELD_TWO_MIN_PRIME = 353
HELD_TWO_MAX_PRIME = 499
HEIGHT_STRATA = ((1, 5_000), (5_001, 12_500), (12_501, 25_000), (25_001, 50_000))
SOURCE_FILENAMES = (
    "elliptic_fermigier_multiple_root_height_h50000.json",
    "elliptic_fermigier_batch_rank_triage.json",
    "elliptic_fermigier_rank22_accidental_slices.json",
    "elliptic_nagao_rank21_accidental_slices.json",
)
EXPECTED_SOURCE_SHA256 = {
    "elliptic_fermigier_multiple_root_height_h50000.json": (
        "34769e4f16833717fc20da13ab737f324d4cb0086b1e0f97f0830940dcc4a08d"
    ),
    "elliptic_fermigier_batch_rank_triage.json": (
        "d1fb509ecab0f1dcdb97bb9db852e807f7fc7b7c001211c6fcf11b2fe52f33f6"
    ),
    "elliptic_nagao_rank21_accidental_slices.json": (
        "125a6b0df7941099547039302b6f1878b5009dcde774328527952699877b1670"
    ),
}
LEGACY_FERMIGIER_ACCIDENTAL_PARAMETERS = tuple(
    Fraction(value)
    for value in ("19033/135", "22253/114", "31331/104", "38633/138")
)
LEGACY_FERMIGIER_ACCIDENTAL_ARTIFACT_SHA256 = (
    "d6670272bf88e225b44f5a7b0b6e45023652c62960248cb8674d7798d8e3af30"
)


@dataclass(frozen=True)
class ScoredCandidate:
    candidate: ProjectiveCandidate
    discovery: dict[str, Any]
    held_one: dict[str, Any]
    held_two: dict[str, Any]

    @property
    def t(self) -> Fraction:
        return self.candidate.parameter


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter_digest(parameters: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for parameter in parameters:
        digest.update((rational_to_string(parameter) + "\n").encode())
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def exact_prior_exclusions(source_dir: Path) -> tuple[set[Fraction], dict[str, Any]]:
    """Load only previously searched fiber parameters, never score-derived rows."""

    batch_path = source_dir / "elliptic_fermigier_batch_rank_triage.json"
    fermigier_path = source_dir / "elliptic_fermigier_rank22_accidental_slices.json"
    nagao_path = source_dir / "elliptic_nagao_rank21_accidental_slices.json"
    batch = read_json(batch_path)
    fermigier = read_json(fermigier_path)
    nagao = read_json(nagao_path)

    batch_parameters = {abs(Fraction(row["t"])) for row in batch["results"]}
    fermigier_parameters = {
        abs(Fraction(row["parameter_t"]))
        for row in fermigier["candidate_conductor_screen"]
    }
    fermigier_parameters.add(abs(Fraction(fermigier["record_parameter_normalized_T"])))
    fermigier_parameter_sha256 = parameter_digest(sorted(fermigier_parameters))
    if fermigier_parameter_sha256 != EXPECTED_FERMIGIER_ACCIDENTAL_PARAMETER_SHA256:
        raise AssertionError("the pinned Fermigier accidental parameter set changed")
    legacy_fermigier_parameters = {
        abs(parameter) for parameter in LEGACY_FERMIGIER_ACCIDENTAL_PARAMETERS
    }
    nagao_parameters = {
        abs(Fraction(row["T"]))
        for row in nagao["auxiliary_search"]["pinned_identity_pilot"][
            "association_records"
        ]
    }
    exclusions = (
        batch_parameters
        | fermigier_parameters
        | legacy_fermigier_parameters
        | nagao_parameters
        | {abs(BENCHMARK)}
    )
    if len(batch_parameters) != 11:
        raise AssertionError("the pinned batch-rank-triage exclusion set changed")
    if len(exclusions) != EXPECTED_EXCLUSION_COUNT:
        raise AssertionError(
            f"expected {EXPECTED_EXCLUSION_COUNT} exact prior parameters, got {len(exclusions)}"
        )
    return exclusions, {
        "canonicalization": "T -> abs(T), because the Fermigier family is exactly even",
        "batch_rank_triage_parameters": len(batch_parameters),
        "fermigier_accidental_slice_parameters_including_benchmark": len(
            fermigier_parameters
        ),
        "fermigier_accidental_parameter_sha256": fermigier_parameter_sha256,
        "legacy_fermigier_accidental_slice_parameters": len(
            legacy_fermigier_parameters
        ),
        "nagao_accidental_slice_H200000_parameters": len(nagao_parameters),
        "unique_prior_parameters_including_benchmark": len(exclusions),
        "benchmark_t": rational_to_string(abs(BENCHMARK)),
        "source_sha256": {
            batch_path.name: sha256_file(batch_path),
            fermigier_path.name: sha256_file(fermigier_path),
            nagao_path.name: sha256_file(nagao_path),
        },
        "legacy_fermigier_accidental_artifact_sha256": (
            LEGACY_FERMIGIER_ACCIDENTAL_ARTIFACT_SHA256
        ),
        "legacy_fermigier_parameters": sorted(
            rational_to_string(parameter) for parameter in legacy_fermigier_parameters
        ),
    }


def split_score_tables(
    tables: Sequence[PrimeScoreTable],
) -> tuple[tuple[PrimeScoreTable, ...], tuple[PrimeScoreTable, ...], tuple[PrimeScoreTable, ...]]:
    discovery = tuple(table for table in tables if table.prime <= DISCOVERY_MAX_PRIME)
    held_one = tuple(
        table
        for table in tables
        if HELD_ONE_MIN_PRIME <= table.prime <= HELD_ONE_MAX_PRIME
    )
    held_two = tuple(
        table
        for table in tables
        if HELD_TWO_MIN_PRIME <= table.prime <= HELD_TWO_MAX_PRIME
    )
    prime_sets = [
        {table.prime for table in band}
        for band in (discovery, held_one, held_two)
    ]
    if not all(prime_sets) or any(
        left & right
        for index, left in enumerate(prime_sets)
        for right in prime_sets[index + 1 :]
    ):
        raise AssertionError("score bands must be nonempty and pairwise disjoint")
    if set.union(*prime_sets) != {table.prime for table in tables}:
        raise AssertionError("the declared prime bands did not partition B=500")
    return discovery, held_one, held_two


def score_population(
    population: Sequence[ProjectiveCandidate],
    exclusions: set[Fraction],
) -> tuple[list[ScoredCandidate], dict[str, Any]]:
    tables = build_score_tables(SCORE_BOUND, "fermigier-good")
    discovery, held_one, held_two = split_score_tables(tables)
    scored: list[ScoredCandidate] = []
    excluded_population = []
    for candidate in population:
        if candidate.parameter in exclusions:
            excluded_population.append(candidate.identifier)
            continue
        score = lambda band: score_rational(  # noqa: E731 - compact exact replay
            candidate.numerator, candidate.denominator, band
        )
        scored.append(
            ScoredCandidate(
                candidate=candidate,
                discovery=score(discovery),
                held_one=score(held_one),
                held_two=score(held_two),
            )
        )
    return scored, {
        "score_definition": "sum (2-a_p)/#E(F_p)*log(p), good primes only",
        "selection_does_not_use": "the cumulative B=500 score or any point/rank result",
        "discovery_prime_interval": [5, DISCOVERY_MAX_PRIME],
        "held_forward_prime_intervals": [
            [HELD_ONE_MIN_PRIME, HELD_ONE_MAX_PRIME],
            [HELD_TWO_MIN_PRIME, HELD_TWO_MAX_PRIME],
        ],
        "band_prime_counts": {
            "discovery": len(discovery),
            "held_one": len(held_one),
            "held_two": len(held_two),
        },
        "prior_parameters_intersecting_population": sorted(excluded_population),
        "prior_parameter_intersection_count": len(excluded_population),
        "eligible_population_count": len(scored),
    }


def percentile(value: float, ordered_values: Sequence[float]) -> float:
    return bisect_right(ordered_values, value) / len(ordered_values)


def deterministic_key(row: ScoredCandidate) -> tuple[int, int, int]:
    candidate = row.candidate
    return candidate.height, candidate.numerator, candidate.denominator


def select_stratified_tranche(
    rows: Sequence[ScoredCandidate],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Select 10 per height stratum plus 8 globally low held-tail fibers."""

    selected: list[dict[str, Any]] = []
    selected_parameters: set[Fraction] = set()
    stratum_counts: dict[str, int] = {}
    for low, high in HEIGHT_STRATA:
        stratum = [row for row in rows if low <= row.candidate.height <= high]
        if not stratum:
            raise AssertionError("an announced height stratum was empty")
        held_one_values = sorted(float(row.held_one["value"]) for row in stratum)
        held_two_values = sorted(float(row.held_two["value"]) for row in stratum)
        percentile_rows = []
        for row in stratum:
            held_one_percentile = percentile(
                float(row.held_one["value"]), held_one_values
            )
            held_two_percentile = percentile(
                float(row.held_two["value"]), held_two_values
            )
            percentile_rows.append(
                (
                    row,
                    held_one_percentile,
                    held_two_percentile,
                    min(held_one_percentile, held_two_percentile),
                )
            )

        def take(reason: str, field_index: int, count: int) -> None:
            available = [
                item for item in percentile_rows if item[0].t not in selected_parameters
            ]
            available.sort(
                key=lambda item: (-item[field_index], *deterministic_key(item[0]))
            )
            for row, p_one, p_two, robust in available[:count]:
                selected_parameters.add(row.t)
                selected.append(
                    selection_record(
                        row,
                        reason=f"height-{low}-{high}:{reason}",
                        held_one_percentile=p_one,
                        held_two_percentile=p_two,
                        robust_percentile=robust,
                    )
                )

        take("robust-held-maximin", 3, 4)
        take("held-one", 1, 3)
        take("held-two", 2, 3)
        stratum_counts[f"{low}-{high}"] = len(stratum)

    all_held_one = sorted(float(row.held_one["value"]) for row in rows)
    all_held_two = sorted(float(row.held_two["value"]) for row in rows)
    height_tail = []
    for row in rows:
        p_one = percentile(float(row.held_one["value"]), all_held_one)
        p_two = percentile(float(row.held_two["value"]), all_held_two)
        if p_one >= 0.75 and p_two >= 0.75 and row.t not in selected_parameters:
            height_tail.append((row, p_one, p_two))
    height_tail.sort(key=lambda item: deterministic_key(item[0]))
    for row, p_one, p_two in height_tail[:8]:
        selected_parameters.add(row.t)
        selected.append(
            selection_record(
                row,
                reason="global-low-height-among-upper-quartile-of-both-held-bands",
                held_one_percentile=p_one,
                held_two_percentile=p_two,
                robust_percentile=min(p_one, p_two),
            )
        )

    if len(selected) != EXPECTED_SELECTED_COUNT or len(selected_parameters) != len(selected):
        raise AssertionError("the declared disjoint tranche did not have 48 fibers")
    selected.sort(key=lambda record: (record["height"], record["numerator"], record["denominator"]))
    selected_sha256 = parameter_digest(Fraction(record["t"]) for record in selected)
    if selected_sha256 != EXPECTED_SELECTED_SHA256:
        raise AssertionError("the held-forward 48-fiber tranche changed")
    return selected, {
        "height_strata": [list(interval) for interval in HEIGHT_STRATA],
        "eligible_counts_by_height_stratum": stratum_counts,
        "per_stratum": (
            "top 4 by min(held-one percentile, held-two percentile), then top 3 "
            "unused by each held-forward score"
        ),
        "height_tail": (
            "8 globally lowest-height unused fibers in the upper quartile of both "
            "held-forward bands"
        ),
        "selected_count": len(selected),
        "selected_parameter_sha256": selected_sha256,
    }


def selection_record(
    row: ScoredCandidate,
    *,
    reason: str,
    held_one_percentile: float,
    held_two_percentile: float,
    robust_percentile: float,
) -> dict[str, Any]:
    candidate = row.candidate
    return {
        "t": candidate.identifier,
        "numerator": candidate.numerator,
        "denominator": candidate.denominator,
        "height": candidate.height,
        "class_index": candidate.class_index,
        "crt_residue": candidate.crt_residue,
        "selection_reason": reason,
        "scores": {
            "discovery": score_record(row.discovery),
            "held_one": score_record(row.held_one),
            "held_two": score_record(row.held_two),
        },
        "held_forward_percentiles": {
            "held_one": held_one_percentile,
            "held_two": held_two_percentile,
            "maximin": robust_percentile,
        },
    }


def score_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "value": format(float(record["value"]), ".17g"),
        "primes_used": int(record["primes_used"]),
        "skipped_denominator_primes": int(record["skipped_denominator_primes"]),
        "skipped_bad_primes": int(record["skipped_bad_primes"]),
    }


def conductor_probe(parameter: Fraction, *, timeout: float, stack_bytes: int) -> dict[str, Any]:
    started = time.monotonic()
    try:
        data = minimal_curve_data(
            FermigierMestreFamily.coefficients(parameter),
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        return {
            "status": "completed",
            "wall_seconds": time.monotonic() - started,
            **data,
            "below_strict_log_conductor_target": (
                Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
        }
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
        }
    except (RuntimeError, AssertionError, ValueError, FileNotFoundError) as error:
        return {
            "status": "error",
            "error": str(error)[:1000],
            "wall_seconds": time.monotonic() - started,
        }


def exact_basis_points(stage: dict[str, Any]) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(
        (Fraction(record["jacobian_x"]), Fraction(record["jacobian_y"]))
        for record in stage["explicit_numerical_basis"]
    )


def write_checkpoint(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--point-timeout", type=float, default=45.0)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--precisions", default="72,120")
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_record_residue_deep_tranche.json",
    )
    return parser


def parse_precisions(value: str) -> tuple[int, ...]:
    values = tuple(int(part) for part in value.split(",") if part)
    if len(values) < 2 or any(value < 38 for value in values):
        raise argparse.ArgumentTypeError("provide at least two precisions >=38")
    return values


def main() -> None:
    args = build_parser().parse_args()
    args.precisions = parse_precisions(args.precisions)
    if args.height != DEFAULT_HEIGHT:
        # This experiment is deliberately pinned, not a tunable population search.
        raise SystemExit(f"--height must remain pinned at {DEFAULT_HEIGHT}")
    if min(args.conductor_timeout, args.point_timeout, args.saturation_timeout) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("PARI stack must be at least 8,000,000 bytes")

    root = Path(__file__).resolve().parents[2]
    source_dir = root / "artifacts" / "generated-results"
    script_path = Path(__file__).resolve()
    experiment_started = time.monotonic()

    source_hashes = {
        filename: sha256_file(source_dir / filename) for filename in SOURCE_FILENAMES
    }
    if any(
        source_hashes[filename] != expected
        for filename, expected in EXPECTED_SOURCE_SHA256.items()
    ):
        raise AssertionError("one or more pinned input artifacts changed")
    prior_exclusions, exclusion_record = exact_prior_exclusions(source_dir)
    groups = choose_groups(DEFAULT_GROUP_ORDER, ())
    classes = crt_classes(groups)
    population = enumerate_projective_height(classes, args.height)
    local_records, population_summary = verify_population_locally(population, groups)
    if len(population) != EXPECTED_POPULATION_COUNT:
        raise AssertionError("the exact H=50000 population count changed")
    if population_summary["population_and_valuation_sha256"] != EXPECTED_POPULATION_SHA256:
        raise AssertionError("the exact H=50000 population digest changed")
    if len(local_records) != len(population):
        raise AssertionError("the pinned population unexpectedly contains a singular fiber")

    scored, scoring_record = score_population(population, prior_exclusions)
    selected, selection_record_data = select_stratified_tranche(scored)
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "bounded leakage-controlled experiment; CRT membership and point membership "
            "are exact, conductor is a PARI computation, height ranks are numerical "
            "triage evidence, and only a successful finite-reduction block is a rank proof"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "hit": False,
        },
        "population_replay": population_summary,
        "prior_exclusions": exclusion_record,
        "scoring": scoring_record,
        "selection": selection_record_data,
        "candidates": selected,
        "execution": {
            "phase": "selection-complete",
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "point_timeout_seconds_per_subprocess": args.point_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "stack_bytes": args.stack_bytes,
            "precisions": list(args.precisions),
            "point_height_stages": [50_000, 250_000, 1_000_000],
            "no_retries": True,
        },
        "sources": {
            "sha256": source_hashes,
            "full_population_note": (
                "the stored H=50000 artifact retained only staged prefixes, so this "
                "script reran its exact scanner and verified the pinned digest"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    write_checkpoint(args.output, artifact)

    # Conductor first for every selected fiber, without score-dependent stopping.
    for record in artifact["candidates"]:
        record["conductor_probe"] = conductor_probe(
            Fraction(record["t"]),
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        write_checkpoint(args.output, artifact)
    artifact["execution"]["phase"] = "conductors-complete"
    write_checkpoint(args.output, artifact)

    subtarget = [
        record
        for record in artifact["candidates"]
        if record["conductor_probe"].get("below_strict_log_conductor_target") is True
    ]
    # Root number is parity triage only.  Searching root -1 first makes an
    # interruption most useful but does not alter membership in the point tranche.
    subtarget.sort(
        key=lambda record: (
            record["conductor_probe"].get("root_number") != -1,
            Decimal(record["conductor_probe"]["log_conductor"]),
            record["height"],
        )
    )
    by_t = {record["t"]: record for record in artifact["candidates"]}
    for selected_record in subtarget:
        parameter = Fraction(selected_record["t"])
        triage = TriageCandidate(
            identifier=f"deep-{parameter.numerator}-{parameter.denominator}",
            parameter=parameter,
            provenance="held-forward stratified H=50000 record-residue tranche",
            selection_metadata=selected_record["selection_reason"],
            known_log_conductor=selected_record["conductor_probe"]["log_conductor"],
            known_global_root_number=selected_record["conductor_probe"]["root_number"],
        )
        result = by_t[selected_record["t"]]
        result["point_search"] = {"status": "started", "stages": []}
        try:
            seeds = FermigierMestreFamily.known_jacobian_points(parameter)[1:]
            seed_height_runs = height_replay(
                parameter,
                seeds,
                precisions=args.precisions,
                timeout=args.point_timeout,
                stack_bytes=args.stack_bytes,
            )
            if stable_rank(seed_height_runs) != 12:
                raise AssertionError("the generic seed did not have stable numerical rank 12")
            result["point_search"]["seed_height_matrix_runs"] = list(seed_height_runs)
            for height_bound in (50_000, 250_000, 1_000_000):
                if height_bound > 50_000:
                    previous = result["point_search"]["stages"][-1]
                    if previous["stable_pool_numerical_rank"] <= 12:
                        break
                stage = run_stage(
                    triage,
                    height_bound=height_bound,
                    precisions=args.precisions,
                    timeout=args.point_timeout,
                    stack_bytes=args.stack_bytes,
                    seed_height_runs=tuple(seed_height_runs),
                )
                result["point_search"]["stages"].append(stage)
                write_checkpoint(args.output, artifact)
                if stage["stable_pool_numerical_rank"] >= 20:
                    result["point_search"]["finite_reduction_attempt"] = (
                        finite_reduction_attempt(
                            FermigierMestreFamily.coefficients(parameter),
                            exact_basis_points(stage),
                            saturation_timeout=args.saturation_timeout,
                            stack_bytes=args.stack_bytes,
                            certificate_prime_bound=args.certificate_prime_bound,
                        )
                    )
                    certificate = result["point_search"]["finite_reduction_attempt"]
                    if (
                        certificate.get("certified_algebraic_rank_lower_bound", 0) >= 21
                        and selected_record["conductor_probe"][
                            "below_strict_log_conductor_target"
                        ]
                    ):
                        artifact["target"] = {
                            "rank_at_least": certificate[
                                "certified_algebraic_rank_lower_bound"
                            ],
                            "strict_log_conductor_upper_bound": str(
                                TARGET_LOG_CONDUCTOR
                            ),
                            "hit": True,
                            "t": selected_record["t"],
                            "log_conductor": selected_record["conductor_probe"][
                                "log_conductor"
                            ],
                        }
                    write_checkpoint(args.output, artifact)
            result["point_search"]["status"] = "completed-declared-stages"
        except subprocess.TimeoutExpired as error:
            result["point_search"]["status"] = "timeout-no-retry"
            result["point_search"]["error"] = str(error)[:1000]
        except (RuntimeError, AssertionError, ValueError, FileNotFoundError) as error:
            result["point_search"]["status"] = "error-no-retry"
            result["point_search"]["error"] = str(error)[:1000]
        write_checkpoint(args.output, artifact)

    completed_conductors = [
        record
        for record in artifact["candidates"]
        if record["conductor_probe"]["status"] == "completed"
    ]
    ranked_stages = [
        stage
        for record in artifact["candidates"]
        for stage in record.get("point_search", {}).get("stages", [])
    ]
    artifact["outcome"] = {
        "selected_fibers": len(artifact["candidates"]),
        "completed_conductors": len(completed_conductors),
        "conductor_timeouts_or_errors": len(artifact["candidates"])
        - len(completed_conductors),
        "subtarget_completed_conductors": len(subtarget),
        "point_searched_fibers": sum(
            "point_search" in record for record in artifact["candidates"]
        ),
        "maximum_stable_numerical_rank": max(
            (stage["stable_pool_numerical_rank"] for stage in ranked_stages),
            default=12,
        ),
        "certified_target_hits": [
            record["t"]
            for record in artifact["candidates"]
            if record.get("point_search", {})
            .get("finite_reduction_attempt", {})
            .get("certified_algebraic_rank_lower_bound", 0)
            >= 21
            and record["conductor_probe"].get("below_strict_log_conductor_target")
        ],
    }
    artifact["execution"]["phase"] = "complete"
    artifact["execution"]["wall_seconds"] = time.monotonic() - experiment_started
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no stable numerical rank >=21 received a successful exact "
            "finite-reduction independence certificate"
        )
    write_checkpoint(args.output, artifact)


if __name__ == "__main__":
    main()

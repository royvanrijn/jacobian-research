#!/usr/bin/env python3
"""Exhaust the five-group multiple-root locus by projective height.

This is the height-box counterpart to ``search_multiple_root_crt.py``.  It
does not enumerate a box in a reduced lattice basis.  Instead it visits every
primitive rational ``T=a/b`` with ``b>0`` and

    max(abs(a), b) <= H

whose residue belongs to one of the 144 CRT classes modulo 6,441,589.  The
Fermigier family is even in ``T`` and the class union is stable under
negation, so the identical specializations ``T`` and ``-T`` are represented
once, with positive numerator.

Every candidate is checked by exact integer arithmetic against all five
local conditions before scoring.  The full population enters the first
good-prime score stage; later stages see only the declared retained prefix
from the preceding stage.  A small deterministic final prefix is sent to
PARI/GP for conductor and local-reduction computations.  Scores and root
numbers are triage features only, and no Mordell--Weil rank is computed or
claimed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any, Iterable, Sequence

from compare_score_cutoffs import parse_cutoffs
from ek_k3 import rational_to_string, valuation
from fermigier_mestre import (
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
)
from pari_bridge import minimal_curve_data, pari_version
from search_multiple_root_crt import (
    DEFAULT_GROUP_ORDER,
    LocalConstraintGroup,
    choose_groups,
    crt_classes,
    exact_group_certificate,
)
from search_record_residue_class import build_score_tables, score_rational
from staged_record_rescore import score_candidates_with_pari


TARGET_LOG_CONDUCTOR = Decimal("182.72")
DEFAULT_HEIGHT = 50_000
DEFAULT_STAGE_CUTOFFS = (200, 2_000, 10_000)
DEFAULT_KEEP_COUNTS = (256, 32, 12)


@dataclass(frozen=True)
class ProjectiveCandidate:
    """One sign-normalized rational in an allowed CRT class."""

    numerator: int
    denominator: int
    class_index: int
    crt_residue: int

    @property
    def parameter(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)

    @property
    def identifier(self) -> str:
        return rational_to_string(self.parameter)

    @property
    def height(self) -> int:
        return max(self.numerator, self.denominator)


def ceiling_division(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def parse_keep_counts(value: str) -> tuple[int, ...]:
    """Parse a nonincreasing comma-separated survivor schedule."""

    try:
        counts = tuple(int(part) for part in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("keep counts must be integers") from error
    if not counts or any(count < 1 for count in counts):
        raise argparse.ArgumentTypeError("every keep count must be positive")
    if any(left < right for left, right in zip(counts, counts[1:])):
        raise argparse.ArgumentTypeError("keep counts must be nonincreasing")
    return counts


def negation_orbits(
    residues: Iterable[int], modulus: int
) -> tuple[tuple[int, int], ...]:
    """Partition a negation-stable residue set into unordered sign orbits."""

    if modulus < 2:
        raise ValueError("the CRT modulus must be at least two")
    residue_set = {residue % modulus for residue in residues}
    if not residue_set:
        raise ValueError("at least one CRT residue is required")
    missing = sorted((-residue) % modulus for residue in residue_set if (-residue) % modulus not in residue_set)
    if missing:
        raise ValueError(f"the residue set is not stable under negation: {missing}")

    unseen = set(residue_set)
    orbits: list[tuple[int, int]] = []
    while unseen:
        residue = min(unseen)
        negative = (-residue) % modulus
        orbit = (residue, negative) if residue <= negative else (negative, residue)
        orbits.append(orbit)
        unseen.discard(residue)
        unseen.discard(negative)
    return tuple(sorted(orbits))


def enumerate_projective_height(
    classes: Sequence[dict[str, Any]], height: int
) -> tuple[ProjectiveCandidate, ...]:
    """Enumerate the exact sign-quotiented projective-height population.

    The implementation works for arbitrary height relative to the CRT
    modulus: every integer translate in ``[-height,height]`` is visited.
    Only one member of each residue-negation orbit is traversed, with negative
    numerators reflected to the canonical positive representative.
    """

    if height < 1:
        raise ValueError("the projective-height bound must be positive")
    if not classes:
        raise ValueError("at least one CRT class is required")
    moduli = {int(item["crt_modulus"]) for item in classes}
    if len(moduli) != 1:
        raise ValueError("all CRT classes must use one common modulus")
    modulus = moduli.pop()
    by_residue: dict[int, dict[str, Any]] = {}
    for item in classes:
        residue = int(item["crt_residue"]) % modulus
        if residue in by_residue:
            raise ValueError("CRT residues must be unique")
        by_residue[residue] = item

    orbits = negation_orbits(by_residue, modulus)
    records: dict[tuple[int, int], ProjectiveCandidate] = {}
    for denominator in range(1, height + 1):
        if gcd(denominator, modulus) != 1:
            continue
        for source_residue, _ in orbits:
            residue = source_residue * denominator % modulus
            minimum_multiplier = ceiling_division(-height - residue, modulus)
            maximum_multiplier = (height - residue) // modulus
            for multiplier in range(minimum_multiplier, maximum_multiplier + 1):
                signed_numerator = residue + multiplier * modulus
                if gcd(signed_numerator, denominator) != 1:
                    continue
                numerator = abs(signed_numerator)
                if numerator == 0:
                    canonical_residue = 0
                elif signed_numerator > 0:
                    canonical_residue = source_residue
                else:
                    canonical_residue = (-source_residue) % modulus
                source_class = by_residue[canonical_residue]
                candidate = ProjectiveCandidate(
                    numerator=numerator,
                    denominator=denominator,
                    class_index=int(source_class["class_index"]),
                    crt_residue=canonical_residue,
                )
                key = (numerator, denominator)
                previous = records.get(key)
                if previous is not None and previous != candidate:
                    raise AssertionError("a rational belongs to two distinct CRT classes")
                records[key] = candidate

    ordered = tuple(
        sorted(
            records.values(),
            key=lambda candidate: (
                candidate.height,
                candidate.numerator,
                candidate.denominator,
            ),
        )
    )
    assert_projective_population(ordered, classes, height)
    return ordered


def brute_force_projective_height(
    classes: Sequence[dict[str, Any]], height: int
) -> tuple[tuple[int, int], ...]:
    """Slow reference enumerator used to audit the optimized traversal."""

    if height < 1:
        raise ValueError("the projective-height bound must be positive")
    if not classes:
        raise ValueError("at least one CRT class is required")
    modulus = int(classes[0]["crt_modulus"])
    if any(int(item["crt_modulus"]) != modulus for item in classes):
        raise ValueError("all CRT classes must use one common modulus")
    residues = {int(item["crt_residue"]) % modulus for item in classes}
    negation_orbits(residues, modulus)
    records: list[tuple[int, int]] = []
    for denominator in range(1, height + 1):
        if gcd(denominator, modulus) != 1:
            continue
        inverse = pow(denominator, -1, modulus)
        for numerator in range(0, height + 1):
            if gcd(numerator, denominator) != 1:
                continue
            if numerator * inverse % modulus in residues:
                records.append((numerator, denominator))
    return tuple(
        sorted(records, key=lambda item: (max(item[0], item[1]), item[0], item[1]))
    )


def assert_projective_population(
    candidates: Sequence[ProjectiveCandidate],
    classes: Sequence[dict[str, Any]],
    height: int,
) -> None:
    """Check normalization, bounds, congruences, uniqueness, and class labels."""

    modulus = int(classes[0]["crt_modulus"])
    class_lookup = {
        int(item["crt_residue"]) % modulus: int(item["class_index"])
        for item in classes
    }
    seen: set[tuple[int, int]] = set()
    for candidate in candidates:
        key = (candidate.numerator, candidate.denominator)
        if key in seen:
            raise AssertionError("the projective population contains a duplicate")
        seen.add(key)
        if candidate.numerator < 0 or candidate.denominator <= 0:
            raise AssertionError("the sign-normalized coordinates are invalid")
        if candidate.height > height:
            raise AssertionError("a candidate exceeds the declared height")
        if gcd(candidate.numerator, candidate.denominator) != 1:
            raise AssertionError("a candidate is not primitive")
        if gcd(candidate.denominator, modulus) != 1:
            raise AssertionError("a constrained denominator is not a CRT unit")
        residue = candidate.numerator * pow(candidate.denominator, -1, modulus) % modulus
        if residue != candidate.crt_residue:
            raise AssertionError("a candidate lost its CRT residue")
        if class_lookup.get(residue) != candidate.class_index:
            raise AssertionError("a candidate has the wrong CRT class index")


def homogeneous_discriminant_factor(numerator: int, denominator: int) -> int:
    """Return ``denominator^20 * H(numerator/denominator)`` exactly."""

    degree = len(DISCRIMINANT_FACTOR_COEFFICIENTS) - 1
    value = DISCRIMINANT_FACTOR_COEFFICIENTS[degree]
    denominator_power = denominator
    for index in range(degree - 1, -1, -1):
        value = value * numerator + DISCRIMINANT_FACTOR_COEFFICIENTS[index] * denominator_power
        denominator_power *= denominator
    return value


def verify_local_constraints(
    candidate: ProjectiveCandidate,
    groups: Sequence[LocalConstraintGroup],
) -> dict[str, Any]:
    """Verify all allowed residues and actual forced valuations exactly."""

    homogeneous_h = homogeneous_discriminant_factor(
        candidate.numerator, candidate.denominator
    )
    if homogeneous_h == 0:
        return {
            "singular": True,
            "residues": {},
            "h_valuations": {},
        }
    residues: dict[str, int] = {}
    valuations: dict[str, int] = {}
    for group in groups:
        if gcd(candidate.denominator, group.modulus) != 1:
            raise AssertionError("a constrained denominator is not a local unit")
        residue = (
            candidate.numerator
            * pow(candidate.denominator, -1, group.modulus)
            % group.modulus
        )
        if residue not in group.residues:
            raise AssertionError(
                f"T={candidate.identifier} is outside the p={group.prime} union"
            )
        actual = valuation(homogeneous_h, group.prime)
        if actual < group.forced_h_valuation:
            raise AssertionError(
                f"T={candidate.identifier} lost the v_{group.prime}(H) guarantee"
            )
        residues[str(group.prime)] = residue
        valuations[str(group.prime)] = actual
    return {
        "singular": False,
        "residues": residues,
        "h_valuations": valuations,
    }


def verify_population_locally(
    candidates: Sequence[ProjectiveCandidate],
    groups: Sequence[LocalConstraintGroup],
) -> tuple[dict[tuple[int, int], dict[str, Any]], dict[str, Any]]:
    """Run exact local verification over the complete enumerated population."""

    records: dict[tuple[int, int], dict[str, Any]] = {}
    singular: list[str] = []
    minima = {str(group.prime): None for group in groups}
    maxima = {str(group.prime): None for group in groups}
    residue_counts = {str(group.prime): Counter() for group in groups}
    class_counts: Counter[int] = Counter()
    digest = hashlib.sha256()
    for candidate in candidates:
        local = verify_local_constraints(candidate, groups)
        digest.update(f"{candidate.numerator}/{candidate.denominator}\n".encode())
        if local["singular"]:
            singular.append(candidate.identifier)
            continue
        base = {
            "t": candidate.identifier,
            "numerator": candidate.numerator,
            "denominator": candidate.denominator,
            "height": candidate.height,
            "class_index": candidate.class_index,
            "crt_residue": candidate.crt_residue,
            "residues": local["residues"],
            "h_valuations": local["h_valuations"],
        }
        records[(candidate.numerator, candidate.denominator)] = base
        class_counts[candidate.class_index] += 1
        for prime, actual in local["h_valuations"].items():
            residue_counts[prime][local["residues"][prime]] += 1
            minima[prime] = actual if minima[prime] is None else min(minima[prime], actual)
            maxima[prime] = actual if maxima[prime] is None else max(maxima[prime], actual)
        digest.update(
            (
                "|"
                + ",".join(
                    f"{prime}:{local['h_valuations'][prime]}"
                    for prime in sorted(local["h_valuations"], key=int)
                )
                + "\n"
            ).encode()
        )
    summary = {
        "candidates_checked_exactly": len(candidates),
        "nonsingular_candidates": len(records),
        "singular_candidates": singular,
        "minimum_observed_h_valuations": minima,
        "maximum_observed_h_valuations": maxima,
        "class_population_counts": {
            str(index): class_counts[index]
            for index in sorted(class_counts)
        },
        "local_residue_population_counts": {
            prime: {str(residue): count for residue, count in sorted(counts.items())}
            for prime, counts in residue_counts.items()
        },
        "population_and_valuation_sha256": digest.hexdigest(),
    }
    return records, summary


def score_sort_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -Decimal(str(record["score"])),
        int(record["height"]),
        int(record["numerator"]),
        int(record["denominator"]),
    )


def score_with_tables(
    candidates: Sequence[ProjectiveCandidate], cutoff: int
) -> list[dict[str, Any]]:
    """Score a finite population from dependency-light finite-field tables."""

    tables = build_score_tables(cutoff, "fermigier-good")
    records: list[dict[str, Any]] = []
    for candidate in candidates:
        score = score_rational(candidate.numerator, candidate.denominator, tables)
        records.append(
            {
                "t": candidate.identifier,
                "numerator": candidate.numerator,
                "denominator": candidate.denominator,
                "height": candidate.height,
                "score": format(score["value"], ".17g"),
                "primes_used": score["primes_used"],
                "skipped_denominator_primes": score["skipped_denominator_primes"],
                "skipped_bad_primes": score["skipped_bad_primes"],
            }
        )
    records.sort(key=score_sort_key)
    return records


def candidate_lookup(
    candidates: Sequence[ProjectiveCandidate],
) -> dict[tuple[int, int], ProjectiveCandidate]:
    return {
        (candidate.numerator, candidate.denominator): candidate
        for candidate in candidates
    }


def verify_pari_local_constraints(
    record: dict[str, Any],
    pari: dict[str, Any],
    groups: Sequence[LocalConstraintGroup],
) -> dict[str, bool]:
    """Replay the exact engineered local prediction on PARI's minimal model."""

    checks: dict[str, bool] = {}
    for group in groups:
        prime = str(group.prime)
        local = pari["local_reduction"][prime]
        expected_delta = (
            int(record["h_valuations"][prime])
            - 12 * group.presented_model_scaling
        )
        valid = (
            int(local["conductor_exponent"]) == 1
            and int(local["minimal_c4_valuation"]) == 0
            and int(local["minimal_discriminant_valuation"]) == expected_delta
            and int(local["ellap"]) == 1
        )
        checks[prime] = valid
        if not valid:
            raise AssertionError(
                f"PARI contradicted the p={prime} certificate for T={record['t']}"
            )
    return checks


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument(
        "--stage-cutoffs", type=parse_cutoffs, default=DEFAULT_STAGE_CUTOFFS
    )
    parser.add_argument(
        "--keep-counts", type=parse_keep_counts, default=DEFAULT_KEEP_COUNTS
    )
    parser.add_argument("--conductor-count", type=int, default=6)
    parser.add_argument("--score-timeout", type=float, default=300.0)
    parser.add_argument("--conductor-timeout", type=float, default=60.0)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_multiple_root_height_h50000.json"
        ),
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.height < 1:
        raise SystemExit("--height must be positive")
    if len(args.stage_cutoffs) != len(args.keep_counts):
        raise SystemExit("--keep-counts must have one entry per score cutoff")
    if args.conductor_count < 0:
        raise SystemExit("--conductor-count must be nonnegative")
    if args.conductor_count > args.keep_counts[-1]:
        raise SystemExit("--conductor-count cannot exceed the final keep count")
    if args.score_timeout <= 0 or args.conductor_timeout <= 0:
        raise SystemExit("PARI timeouts must be positive")
    if args.pari_stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack must be at least 8,000,000 bytes")


def main() -> None:
    args = build_parser().parse_args()
    validate_arguments(args)
    groups = choose_groups(DEFAULT_GROUP_ORDER, ())
    classes = crt_classes(groups)
    modulus = int(classes[0]["crt_modulus"])
    orbits = negation_orbits(
        (int(item["crt_residue"]) for item in classes), modulus
    )
    if len(classes) != 144 or modulus != 6_441_589 or len(orbits) != 72:
        raise AssertionError("the pinned five-group CRT geometry changed")
    if any(
        coefficient
        for index, coefficient in enumerate(DISCRIMINANT_FACTOR_COEFFICIENTS)
        if index % 2
    ):
        raise AssertionError("H(T) is no longer even, so sign quotienting is invalid")

    enumeration_start = time.monotonic()
    population = enumerate_projective_height(classes, args.height)
    local_records, local_summary = verify_population_locally(population, groups)
    nonsingular = tuple(
        candidate
        for candidate in population
        if (candidate.numerator, candidate.denominator) in local_records
    )
    enumeration_seconds = time.monotonic() - enumeration_start
    if not nonsingular:
        raise SystemExit("the exhaustive height box contains no nonsingular candidate")

    stage_input = nonsingular
    stages: list[dict[str, Any]] = []
    histories: dict[tuple[int, int], dict[str, Any]] = {}
    benchmark_scores: dict[str, dict[str, Any]] = {}
    for stage_index, (cutoff, keep_count) in enumerate(
        zip(args.stage_cutoffs, args.keep_counts, strict=True)
    ):
        stage_start = time.monotonic()
        engine = (
            "dependency-light exact finite-field tables"
            if stage_index == 0
            else "PARI/GP ellap batch on exact minimal models"
        )
        if stage_index == 0:
            records = score_with_tables(stage_input, cutoff)
            benchmark_record = score_with_tables(
                (
                    ProjectiveCandidate(
                        NORMALIZED_RECORD_PARAMETER.numerator,
                        NORMALIZED_RECORD_PARAMETER.denominator,
                        0,
                        0,
                    ),
                ),
                cutoff,
            )[0]
        else:
            records = score_candidates_with_pari(
                stage_input,
                cutoff,
                timeout=args.score_timeout,
                stack_bytes=args.pari_stack_bytes,
            )
            benchmark_record = score_candidates_with_pari(
                (
                    ProjectiveCandidate(
                        NORMALIZED_RECORD_PARAMETER.numerator,
                        NORMALIZED_RECORD_PARAMETER.denominator,
                        0,
                        0,
                    ),
                ),
                cutoff,
                timeout=args.score_timeout,
                stack_bytes=args.pari_stack_bytes,
            )[0]
        records.sort(key=score_sort_key)
        retained = records[: min(keep_count, len(records))]
        lookup = candidate_lookup(stage_input)
        retained_candidates: list[ProjectiveCandidate] = []
        retained_rows: list[dict[str, Any]] = []
        for position, score_record in enumerate(retained, 1):
            key = (
                int(score_record["numerator"]),
                int(score_record["denominator"]),
            )
            candidate = lookup[key]
            retained_candidates.append(candidate)
            history = histories.setdefault(key, dict(local_records[key]))
            history.setdefault("scores", {})[str(cutoff)] = {
                "score": str(score_record["score"]),
                "position_among_stage_input": position,
                "primes_used": int(score_record["primes_used"]),
                "skipped_denominator_primes": int(
                    score_record["skipped_denominator_primes"]
                ),
                "skipped_bad_primes": int(score_record["skipped_bad_primes"]),
            }
            retained_rows.append(
                {
                    "t": score_record["t"],
                    "score": str(score_record["score"]),
                    "height": int(score_record["height"]),
                    "position": position,
                }
            )
        elapsed = time.monotonic() - stage_start
        benchmark_scores[str(cutoff)] = {
            "score": str(benchmark_record["score"]),
            "source": "separate comparison; never part of survivor selection",
        }
        stages.append(
            {
                "numeric_prime_cutoff": cutoff,
                "engine": engine,
                "input_count": len(records),
                "requested_keep_count": keep_count,
                "retained_count": len(retained),
                "elapsed_seconds": round(elapsed, 6),
                "best": retained_rows[0] if retained_rows else None,
                "worst_retained": retained_rows[-1] if retained_rows else None,
                "retained": retained_rows,
                "benchmark": benchmark_scores[str(cutoff)],
            }
        )
        stage_input = tuple(retained_candidates)

    final_cutoff = str(args.stage_cutoffs[-1])
    finalists = [
        histories[(candidate.numerator, candidate.denominator)]
        for candidate in stage_input
    ]
    finalists.sort(
        key=lambda record: (
            -Decimal(record["scores"][final_cutoff]["score"]),
            int(record["height"]),
            int(record["numerator"]),
            int(record["denominator"]),
        )
    )
    conductor_errors: list[dict[str, str]] = []
    for record in finalists[: args.conductor_count]:
        parameter = Fraction(record["numerator"], record["denominator"])
        try:
            pari = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.conductor_timeout,
                local_primes=tuple(group.prime for group in groups),
                stack_bytes=args.pari_stack_bytes,
            )
            record["pari"] = pari
            record["pari_local_constraints_verified"] = verify_pari_local_constraints(
                record, pari, groups
            )
            record["below_strict_log_conductor_target"] = (
                Decimal(pari["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
            record["root_number_role"] = "parity triage only; no rank inference"
            record["rank_status"] = "not computed or inferred"
        except Exception as error:
            conductor_errors.append({"t": record["t"], "error": str(error)})

    completed_conductors = [record for record in finalists if "pari" in record]
    below_target = [
        record
        for record in completed_conductors
        if record["below_strict_log_conductor_target"]
    ]
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive projective-height experiment; all local constraints "
            "are exact, score and root number are triage only, conductor/local "
            "reduction are PARI computations, and no Mordell--Weil rank is claimed"
        ),
        "family": "normalized Fermigier--Mestre family, exactly even in T",
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
            "reason": "no rank computation or independence certificate was made",
        },
        "enumeration": {
            "definition": (
                "every primitive T=a/b with b>0 and max(abs(a),b)<=height "
                "in the 144-class union; T and -T are one identical curve"
            ),
            "height": args.height,
            "crt_classes": len(classes),
            "crt_modulus": modulus,
            "negation_orbits_enumerated": len(orbits),
            "sign_normalization": "a>0; no a=0 member occurs in this locus",
            "primitive_sign_quotiented_candidates": len(population),
            "elapsed_seconds_including_exact_local_checks": round(
                enumeration_seconds, 6
            ),
            **local_summary,
        },
        "constraint_groups": [exact_group_certificate(group) for group in groups],
        "score_definition": {
            "formula": (
                "sum_{5<=p<=B, p not dividing b, good reduction} "
                "((2-a_p)/(p+1-a_p))*log(p)"
            ),
            "prime_bound_semantics": "numerical prime cutoff p<=B",
            "selection_protocol": (
                "the complete nonsingular height-box population enters stage one; "
                "each later stage receives only the retained prefix of its immediate "
                "predecessor"
            ),
            "rank_inference": "none",
        },
        "parameters": {
            "height": args.height,
            "stage_cutoffs": list(args.stage_cutoffs),
            "keep_counts": list(args.keep_counts),
            "conductor_count": args.conductor_count,
            "score_timeout_seconds_per_batch": args.score_timeout,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "pari_stack_bytes": args.pari_stack_bytes,
            "output": str(args.output),
        },
        "benchmark": {
            "parameter": rational_to_string(NORMALIZED_RECORD_PARAMETER),
            "role": "separate score comparison only; absent from selection",
            "scores": benchmark_scores,
        },
        "stages": stages,
        "finalists": finalists,
        "conductor_summary": {
            "requested": args.conductor_count,
            "completed": len(completed_conductors),
            "errors": len(conductor_errors),
            "below_strict_log_conductor_target": len(below_target),
            "best_log_conductor": (
                str(
                    min(
                        Decimal(record["pari"]["log_conductor"])
                        for record in completed_conductors
                    )
                )
                if completed_conductors
                else None
            ),
            "root_number_use": "parity triage only; no rank inference",
        },
        "conductor_errors": conductor_errors,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"height={args.height} population={len(population)} "
        f"nonsingular={len(nonsingular)} finalists={len(finalists)}"
    )
    for stage in stages:
        print(
            f"B={stage['numeric_prime_cutoff']} input={stage['input_count']} "
            f"retained={stage['retained_count']} best={stage['best']['t']} "
            f"score={stage['best']['score']}"
        )
    for record in completed_conductors:
        print(
            f"T={record['t']} logN={record['pari']['log_conductor']} "
            f"root={record['pari']['root_number']}"
        )


if __name__ == "__main__":
    main()

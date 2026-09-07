#!/usr/bin/env python3
"""Rank-blind bounded slice search from exact weight-three record orbits.

This extends ``search_fermigier_rank22_record_group_directions.py`` by exactly
one coefficient level.  It enumerates the global-sign quotient of the 6,160
vectors in ``{-1,0,1}^22`` having exact l1 norm three, transports every vector
through the pointed record-quartic group law, and excludes the complete
known/l1<=2 abscissa population.

Searching all resulting slices would require 12,320 genus-one calls.  The
declared rank-blind pilot instead retains the union of four deterministic
strata: lowest record-abscissa height, highest exact small-prime square yield,
best square yield for each first support index, and best square yield for each
coefficient-sign pattern.  No specialized rank or conductor enters selection.
Both slopes of every retained direction are searched once at H=50,000.

Every second parameter is exactly checked, decontaminated against the pinned
1,239-parameter Fermigier prior, and sent to exact conductor before rank work.
Only subtarget conductors receive H=50,000 rank triage; stable rank at least 21
triggers saturation and finite-reduction certification.  Unselected triples
are recorded by exact population and digest but are not claimed negative.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
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

import sympy as sp

from ek_k3 import fraction_mod, legendre_symbol, rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    published_preimage_digest,
    rational_digest,
    sha256_file,
    triage_specialization,
)
from search_fermigier_rank22_accidental_slices import T0, conductor_probe, slice_polynomial
from search_fermigier_rank22_record_group_directions import (
    EXPECTED_PUBLISHED_BASIS_POINT_SHA256,
    EXPECTED_TRANSPORT_SOURCE_SHA256,
    RecordQuarticAuxiliary,
    aggregate_candidates,
    auxiliary_combination,
    direction_digest,
    known_record_abscissas,
    load_prior_parameters,
    load_transport_source,
    projective_height,
    search_direction,
    short_combination,
    slice_result_digest,
    transport_source_digest,
)
from search_fermigier_rank22_accidental_slices import quartic_group_pullback
from triage_nagao_rank13_finalists import point_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

EXPECTED_TRIPLE_VECTOR_COUNT = 6_160
PREVIOUS_ARTIFACT = "elliptic_fermigier_rank22_record_group_directions.json"
EXPECTED_PREVIOUS_ARTIFACT_SHA256 = (
    "4928c44e27cada74b7a558dd97edfba20a554a508ac0e98ca051df7dea66a3c1"
)
EXPECTED_PREVIOUS_SLICE_RESULT_SHA256 = (
    "867ebd2e8c2b1a3608ea5577dbaa356de6c1b8e2574c45290fcbbbbf2ea4baf5"
)
EXPECTED_WEIGHT2_X_SHA256 = (
    "d1a922966f538dd14c3063d913488563908706ab218e81106881728cb9d36d62"
)
EXPECTED_PRIOR_RECORD_X_COUNT = 494
EXPECTED_PRIOR_RECORD_X_SHA256 = (
    "62fcbbda7493f940b8001582929fedf27e822fc6d40f6ddf15adf8b4d820dc82"
)
MODULAR_PRIMES = (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
LOW_HEIGHT_KEEP = 192
MODULAR_KEEP = 192
PER_FIRST_INDEX_KEEP = 4
PER_SIGN_PATTERN_KEEP = 16
SLICE_HEIGHT = 50_000
BOUNDED_EXECUTION_SCRIPT_SHA256 = (
    "1147c06cecec8cccd7439357ccffc98acda0574c89031bb4d80236bb40107d4f"
)


def triple_vectors() -> tuple[tuple[int, ...], ...]:
    vectors = []
    for first in range(22):
        for second in range(first + 1, 22):
            for third in range(second + 1, 22):
                for second_sign in (-1, 1):
                    for third_sign in (-1, 1):
                        vector = [0] * 22
                        vector[first] = 1
                        vector[second] = second_sign
                        vector[third] = third_sign
                        vectors.append(tuple(vector))
    if len(vectors) != EXPECTED_TRIPLE_VECTOR_COUNT:
        raise AssertionError("the exact sign-quotiented triple population changed")
    return tuple(vectors)


def triple_id(vector: Sequence[int]) -> str:
    return "_".join(
        f"{'p' if value > 0 else 'm'}{index:02d}"
        for index, value in enumerate(vector, start=1)
        if value
    )


def full_orbit_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: item["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{record['quartic_z']}|"
                f"{','.join(map(str, record['coefficient_vector']))}|"
                f"{record['classification']}\n"
            ).encode()
        )
    return digest.hexdigest()


def generate_triple_directions(
    auxiliary: RecordQuarticAuxiliary,
    auxiliary_basis: Sequence[tuple[Fraction, Fraction]],
    short_basis: Sequence[tuple[Fraction, Fraction]],
    prior_x: set[Fraction],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    short_coefficients = FermigierMestreFamily.coefficients(T0)
    raw_records = []
    by_x: dict[Fraction, dict[str, Any]] = {}
    duplicate_abscissas = 0
    prior_exclusions = 0
    exceptional = 0
    for vector in triple_vectors():
        auxiliary_point = auxiliary_combination(
            auxiliary.weierstrass_coefficients, auxiliary_basis, vector
        )
        inverse = auxiliary.inverse(auxiliary_point)
        if inverse is None:
            exceptional += 1
            continue
        x_value, ordinate = inverse
        expected_short = short_combination(short_coefficients, short_basis, vector)
        if quartic_group_pullback(T0, inverse) != expected_short:
            raise AssertionError("a triple orbit inverse lost its short group coordinate")
        classification = "new-weight3-direction"
        if x_value in prior_x:
            classification = "prior-quartic-abscissa-excluded"
            prior_exclusions += 1
        record = {
            "direction_id": triple_id(vector),
            "coefficient_vector": list(vector),
            "quartic_x": rational_to_string(x_value),
            "quartic_z": rational_to_string(ordinate),
            "projective_height": projective_height(x_value),
            "classification": classification,
            "exact_auxiliary_inverse_checked": True,
            "exact_short_group_combination_checked": True,
        }
        raw_records.append(record)
        if classification != "new-weight3-direction":
            continue
        existing = by_x.get(x_value)
        if existing is None:
            by_x[x_value] = record
        else:
            duplicate_abscissas += 1
            existing.setdefault("additional_coefficient_vectors", []).append(list(vector))
    directions = sorted(
        by_x.values(),
        key=lambda record: (
            record["projective_height"],
            Q(record["quartic_x"]),
            record["direction_id"],
        ),
    )
    return directions, {
        "global_sign_quotient": True,
        "coefficient_alphabet": [-1, 0, 1],
        "exact_l1_norm": 3,
        "full_vector_count": EXPECTED_TRIPLE_VECTOR_COUNT,
        "full_vector_direction_sha256": full_orbit_digest(raw_records),
        "exceptional_inverse_count": exceptional,
        "prior_quartic_abscissa_exclusions": prior_exclusions,
        "duplicate_new_abscissas": duplicate_abscissas,
        "genuinely_new_unique_abscissa_count": len(directions),
        "genuinely_new_direction_sha256": direction_digest(directions),
        "minimum_new_abscissa_projective_height": min(
            record["projective_height"] for record in directions
        ),
        "maximum_new_abscissa_projective_height": max(
            record["projective_height"] for record in directions
        ),
    }


def quartic_point_count_mod_prime(
    coefficients: Sequence[Fraction], prime: int
) -> int | None:
    if any(Q(value).denominator % prime == 0 for value in coefficients):
        return None
    reduced = tuple(fraction_mod(Q(value), prime) for value in coefficients)
    if reduced[-1] == 0:
        return None
    affine = 0
    for parameter in range(prime):
        value = 0
        for coefficient in reversed(reduced):
            value = (value * parameter + coefficient) % prime
        character = legendre_symbol(value, prime)
        affine += 1 if character == 0 else 2 if character == 1 else 0
    points_at_infinity = 1 + legendre_symbol(reduced[-1], prime)
    return affine + points_at_infinity


def modular_square_yield(direction: dict[str, Any]) -> dict[str, Any]:
    source_x = Q(direction["quartic_x"])
    slope_records = []
    for slope in (-1, 1):
        intercept = source_x - slope * T0
        coefficients = slice_polynomial(slope, intercept)
        counts = []
        score = Q(0)
        for prime in MODULAR_PRIMES:
            point_count = quartic_point_count_mod_prime(coefficients, prime)
            if point_count is None:
                continue
            excess = point_count - (prime + 1)
            score += Q(excess, prime + 1)
            counts.append(
                {
                    "prime": prime,
                    "projective_point_count": point_count,
                    "point_count_excess": excess,
                }
            )
        if not counts:
            raise AssertionError("a triple direction lost every modular pilot prime")
        score /= len(counts)
        slope_records.append(
            {
                "slope": slope,
                "usable_prime_count": len(counts),
                "average_normalized_point_count_excess": rational_to_string(score),
                "average_normalized_point_count_excess_float": float(score),
                "prime_counts": counts,
            }
        )
    scores = [Q(record["average_normalized_point_count_excess"]) for record in slope_records]
    return {
        "slope_records": slope_records,
        "maximum_slope_score": rational_to_string(max(scores)),
        "sum_slope_score": rational_to_string(sum(scores, Q(0))),
        "maximum_slope_score_float": float(max(scores)),
        "sum_slope_score_float": float(sum(scores, Q(0))),
    }


def modular_sort_key(record: dict[str, Any]) -> tuple[Any, ...]:
    modular = record["modular_square_yield"]
    return (
        -Q(modular["maximum_slope_score"]),
        -Q(modular["sum_slope_score"]),
        record["projective_height"],
        Q(record["quartic_x"]),
        record["direction_id"],
    )


def score_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: item["direction_id"]):
        modular = record["modular_square_yield"]
        digest.update(
            (
                f"{record['direction_id']}|{modular['maximum_slope_score']}|"
                f"{modular['sum_slope_score']}|"
                + ";".join(
                    f"{slope['slope']}:"
                    + ",".join(
                        f"{item['prime']}={item['projective_point_count']}"
                        for item in slope["prime_counts"]
                    )
                    for slope in modular["slope_records"]
                )
                + "\n"
            ).encode()
        )
    return digest.hexdigest()


def selection_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: item["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{','.join(record['selection_strata'])}\n"
            ).encode()
        )
    return digest.hexdigest()


def selection_audit_summary(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    coverages = [
        slope["usable_prime_count"]
        for record in records
        for slope in record["modular_square_yield"]["slope_records"]
    ]
    basis_participation = Counter()
    sign_patterns = Counter()
    for record in records:
        support = [
            (index, value)
            for index, value in enumerate(record["coefficient_vector"], start=1)
            if value
        ]
        for index, _ in support:
            basis_participation[index] += 1
        sign_patterns[f"{support[1][1]:+d},{support[2][1]:+d}"] += 1
    return {
        "selected_slope_usable_prime_count_minimum": min(coverages),
        "selected_slope_usable_prime_count_maximum": max(coverages),
        "selected_slope_usable_prime_count_distribution": {
            str(key): value for key, value in sorted(Counter(coverages).items())
        },
        "selected_basis_participation_minimum": min(basis_participation.values()),
        "selected_basis_participation_maximum": max(basis_participation.values()),
        "selected_relative_sign_pattern_counts": dict(sorted(sign_patterns.items())),
        "audit_caveats": [
            (
                "usable-prime coverage varies by denominator, so averaging over "
                "usable primes is exact but not denominator-robust"
            ),
            (
                "singular reductions are not excluded; the finite-field point "
                "counts remain exact but are not uniformly smooth genus-one signals"
            ),
            (
                "the score is normalized projective point-count excess, counting "
                "two ordinates over a nonzero square and one over zero; it is not "
                "a literal count of square x-residues"
            ),
            (
                "the deterministic quota union preserves some support/sign diversity "
                "but is not claimed balanced"
            ),
        ],
    }


def select_directions(
    directions: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    for index, direction in enumerate(directions, start=1):
        direction["modular_square_yield"] = modular_square_yield(direction)
        if index % 500 == 0:
            print(f"modular pilot {index}/{len(directions)}", flush=True)

    reasons: dict[str, set[str]] = defaultdict(set)
    by_id = {record["direction_id"]: record for record in directions}
    for record in sorted(
        directions,
        key=lambda item: (
            item["projective_height"], Q(item["quartic_x"]), item["direction_id"]
        ),
    )[:LOW_HEIGHT_KEEP]:
        reasons[record["direction_id"]].add("lowest-abscissa-height")
    for record in sorted(directions, key=modular_sort_key)[:MODULAR_KEEP]:
        reasons[record["direction_id"]].add("highest-modular-square-yield")

    by_first: dict[int, list[dict[str, Any]]] = defaultdict(list)
    by_pattern: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for record in directions:
        support = [
            (index, value)
            for index, value in enumerate(record["coefficient_vector"])
            if value
        ]
        by_first[support[0][0]].append(record)
        by_pattern[(support[1][1], support[2][1])].append(record)
    for first_index, population in sorted(by_first.items()):
        for record in sorted(population, key=modular_sort_key)[:PER_FIRST_INDEX_KEEP]:
            reasons[record["direction_id"]].add(
                f"first-support-index-{first_index + 1:02d}-diversity"
            )
    for pattern, population in sorted(by_pattern.items()):
        label = f"relative-sign-pattern-{pattern[0]:+d}-{pattern[1]:+d}"
        for record in sorted(population, key=modular_sort_key)[:PER_SIGN_PATTERN_KEEP]:
            reasons[record["direction_id"]].add(label)

    selected = []
    for identifier, strata in reasons.items():
        record = dict(by_id[identifier])
        record["selection_strata"] = sorted(strata)
        selected.append(record)
    selected.sort(
        key=lambda record: (
            record["projective_height"], Q(record["quartic_x"]), record["direction_id"]
        )
    )
    return selected, {
        "selection_uses_specialized_rank": False,
        "selection_uses_conductor": False,
        "modular_primes": list(MODULAR_PRIMES),
        "modular_score": (
            "for each slope, exact average of (#C(F_p)-(p+1))/(p+1) over "
            "usable primes; direction ordered by max slope then sum of slopes"
        ),
        "lowest_height_keep": LOW_HEIGHT_KEEP,
        "highest_modular_keep": MODULAR_KEEP,
        "per_first_support_index_keep": PER_FIRST_INDEX_KEEP,
        "per_relative_sign_pattern_keep": PER_SIGN_PATTERN_KEEP,
        "full_modular_score_sha256": score_digest(directions),
        "selected_direction_count": len(selected),
        "selected_direction_sha256": direction_digest(selected),
        "selected_population_and_strata_sha256": selection_digest(selected),
        "unselected_direction_count": len(directions) - len(selected),
        "scope_warning": (
            "only the explicitly selected union receives H=50000 searches; "
            "unselected exact triples are not negative search results"
        ),
        "postrun_independent_selector_audit": selection_audit_summary(selected),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument(
        "--finalize-selector-audit-metadata",
        action="store_true",
        help="add the independent selector caveats without rerunning bounded searches",
    )
    parser.add_argument("--slice-height", type=int, default=SLICE_HEIGHT)
    parser.add_argument("--slice-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_rank22_record_group_triples.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != SLICE_HEIGHT:
        raise SystemExit("this bounded triple tranche is pinned at slice H=50000")
    if min(
        args.slice_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0 or max(
        args.slice_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")

    root = Path(__file__).resolve().parents[2]
    artifact_directory = root / "artifacts" / "generated-results"
    if args.finalize_selector_audit_metadata:
        artifact = json.loads(args.output.read_text())
        if (
            artifact.get("status")
            != "completed exact full weight3 orbit, rank-blind selected H=50000 slice tranche"
            or artifact.get("outcome", {}).get("slice_calls_completed") != 798
            or len(artifact.get("selected_directions", [])) != 399
        ):
            raise AssertionError("the completed bounded triple artifact changed")
        artifact["rank_blind_selection"]["postrun_independent_selector_audit"] = (
            selection_audit_summary(artifact["selected_directions"])
        )
        artifact["execution_provenance"] = {
            "bounded_execution_script_sha256": BOUNDED_EXECUTION_SCRIPT_SHA256,
            "stable_replay_script_sha256": sha256_file(Path(__file__).resolve()),
            "metadata_normalization_only": True,
            "bounded_searches_rerun_during_normalization": 0,
            "normalization_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
        }
        artifact["script_sha256"] = sha256_file(Path(__file__).resolve())
        artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_artifact(args.output, artifact)
        return
    primary_path = artifact_directory / PRIMARY_ARTIFACT
    primary = json.loads(primary_path.read_text())
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published accidental preimages changed")
    previous_path = artifact_directory / PREVIOUS_ARTIFACT
    if sha256_file(previous_path) != EXPECTED_PREVIOUS_ARTIFACT_SHA256:
        raise AssertionError("the exact l1<=2 direction artifact changed")
    previous = json.loads(previous_path.read_text())
    if (
        previous["outcome"]["exact_slice_result_sha256"]
        != EXPECTED_PREVIOUS_SLICE_RESULT_SHA256
    ):
        raise AssertionError("the exact l1<=2 slice result changed")

    auxiliary = RecordQuarticAuxiliary.construct()
    auxiliary_basis, short_basis, _ = load_transport_source(primary, auxiliary)
    if (
        transport_source_digest(primary["published_point_preimages"])
        != EXPECTED_TRANSPORT_SOURCE_SHA256
        or point_digest(short_basis) != EXPECTED_PUBLISHED_BASIS_POINT_SHA256
    ):
        raise AssertionError("the exact transported rank22 basis changed")
    known_x, known_record = known_record_abscissas(primary)
    weight2_x = {Q(row["quartic_x"]) for row in previous["direction_searches"]}
    prior_x = known_x | weight2_x
    if (
        len(weight2_x) != 462
        or rational_digest(sorted(weight2_x)) != EXPECTED_WEIGHT2_X_SHA256
        or len(prior_x) != EXPECTED_PRIOR_RECORD_X_COUNT
        or rational_digest(sorted(prior_x)) != EXPECTED_PRIOR_RECORD_X_SHA256
    ):
        raise AssertionError("the exact prior record-quartic x population changed")

    started_preflight = time.monotonic()
    directions, full_population = generate_triple_directions(
        auxiliary, auxiliary_basis, short_basis, prior_x
    )
    selected, selection = select_directions(directions)
    prior_parameters, prior_record = load_prior_parameters(artifact_directory)
    preflight_seconds = time.monotonic() - started_preflight

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "selection-only exact weight3 orbit prepared"
            if args.selection_only
            else "in-progress selected weight3 manufactured-slice search"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "hit": False,
        },
        "source": {
            "primary_artifact": PRIMARY_ARTIFACT,
            "primary_artifact_sha256_observed": sha256_file(primary_path),
            "published_accidental_preimage_sha256": published_preimage_digest(primary),
            "transport_source_sha256": transport_source_digest(
                primary["published_point_preimages"]
            ),
            "certified_published_rank22_basis_point_sha256": point_digest(short_basis),
            "previous_l1_le_2_artifact": PREVIOUS_ARTIFACT,
            "previous_l1_le_2_artifact_sha256": EXPECTED_PREVIOUS_ARTIFACT_SHA256,
            "previous_exact_slice_result_sha256": EXPECTED_PREVIOUS_SLICE_RESULT_SHA256,
        },
        "prior_record_quartic_abscissas": {
            "known_source_union": known_record,
            "weight2_direction_count": len(weight2_x),
            "weight2_direction_x_sha256": rational_digest(sorted(weight2_x)),
            "full_prior_x_count": len(prior_x),
            "full_prior_x_sha256": rational_digest(sorted(prior_x)),
        },
        "full_triple_population": full_population,
        "rank_blind_selection": selection,
        "selected_directions": selected,
        "prior_parameter_decontamination": prior_record,
        "parameters": {
            "exact_l1_norm": 3,
            "coefficient_alphabet": [-1, 0, 1],
            "global_sign_quotient": True,
            "selected_direction_count": len(selected),
            "slopes_per_direction": [-1, 1],
            "declared_slice_call_count": 2 * len(selected),
            "slice_height": args.slice_height,
            "slice_timeout_seconds": args.slice_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "specialization_timeout_seconds": args.specialization_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": [72, 120],
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
            "checkpoint_after_each_selected_direction": True,
        },
        "direction_searches": [],
        "candidates": [],
        "execution": {
            "phase": "selection-only-complete" if args.selection_only else "slice-search-in-progress",
            "preflight_wall_seconds": preflight_seconds,
            "directions_completed": 0,
            "slice_calls_completed_or_attempted": 0,
            "owned_processes_remaining": 0 if args.selection_only else None,
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
    }
    write_artifact(args.output, artifact)
    if args.selection_only:
        artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_artifact(args.output, artifact)
        print(
            f"selection only: full={len(directions)} selected={len(selected)} "
            f"calls={2 * len(selected)}",
            flush=True,
        )
        return

    qualifying_all: list[dict[str, Any]] = []
    started = time.monotonic()
    for index, direction in enumerate(selected, start=1):
        direction_row, qualifying = search_direction(
            direction,
            height_bound=args.slice_height,
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
            prior_parameters=prior_parameters,
        )
        artifact["direction_searches"].append(direction_row)
        qualifying_all.extend(qualifying)
        artifact["execution"].update(
            {
                "directions_completed": index,
                "slice_calls_completed_or_attempted": index * 2,
                "last_direction_id": direction["direction_id"],
                "wall_seconds_so_far": time.monotonic() - started,
            }
        )
        write_artifact(args.output, artifact)
        if index % 25 == 0 or index == len(selected):
            print(
                f"selected directions {index}/{len(selected)} "
                f"qualifying_incidences={len(qualifying_all)}",
                flush=True,
            )

    candidates = aggregate_candidates(qualifying_all)
    artifact["candidates"] = candidates
    artifact["execution"]["phase"] = "conductor-first"
    write_artifact(args.output, artifact)
    for candidate in candidates:
        parameter = Q(candidate["parameter_t"])
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
                    precisions=(72, 120),
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

    slice_rows = [
        row
        for direction in artifact["direction_searches"]
        for row in direction["slice_searches"]
    ]
    classifications = Counter(
        incidence["classification"]
        for row in slice_rows
        for incidence in row["incidences"]
    )
    rank_records = [
        candidate["rank_triage"]
        for candidate in candidates
        if "rank_triage" in candidate
        and "full_pool_stable_numerical_rank" in candidate["rank_triage"]
    ]
    artifact["outcome"] = {
        "full_genuinely_new_triple_direction_count": len(directions),
        "searched_selected_direction_count": len(selected),
        "unsearched_triple_direction_count": len(directions) - len(selected),
        "unsearched_triples_claimed_negative": False,
        "declared_slice_call_count": 2 * len(selected),
        "slice_calls_attempted": len(slice_rows),
        "slice_calls_completed": sum(
            row["search"]["status"] == "completed" for row in slice_rows
        ),
        "slice_calls_timed_out_or_errored": sum(
            row["search"]["status"] != "completed" for row in slice_rows
        ),
        "record_T0_calibrated_slices": sum(
            row["record_T0_calibration_count"] == 1 for row in slice_rows
        ),
        "incidence_classification_counts": dict(sorted(classifications.items())),
        "genuinely_new_forced_fibres": len(candidates),
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
        "exact_selected_slice_result_sha256": slice_result_digest(
            artifact["direction_searches"]
        ),
        "wall_seconds": time.monotonic() - started,
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
            "no searched selected weight3 direction produced a certified subtarget rank21 fibre"
        )
    artifact["status"] = (
        "completed exact full weight3 orbit, rank-blind selected H=50000 slice tranche"
    )
    artifact["execution"].update(
        {
            "phase": "complete",
            "wall_seconds": time.monotonic() - started,
            "owned_processes_remaining": 0,
        }
    )
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()

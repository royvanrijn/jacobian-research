#!/usr/bin/env python3
"""Search a rank/conductor-blind tranche of exact weight-four record orbits.

The Fermigier rank-22 record fibre at ``T0=39508/39`` supplies a certified,
ordered basis of 22 rational points.  This script enumerates every one of the

    binomial(22,4) * 2^3 = 58,520

global-sign-quotiented vectors in ``{-1,0,1}^22`` of exact l1 norm four.  It
transports them through the exact pointed-quartic group law and removes the
complete l1<=3 record-abscissa population pinned by predecessor artifacts.

The bounded slice tranche is selected without looking at any weight-four
specialized rank or conductor.  The selector is the union of the lowest 7%
record-abscissa heights, both 1% tails of an exact small-prime square-density
score, and small support/sign diversity quotas.  The score is computed from
``#C(F_p)-(p+1)`` for both x=+/-T+n charts using a precomputed residue lookup.
The exact same selector is replayed on the fully searched l1=2 and l1=3
populations.  Its retention of the sole l1=3 second fibre is reported as an
in-sample calibration, not as an independent validation: the quota design
was chosen after that lower-weight outcome was known.

Every selected direction receives one H=50,000 search in both charts.  Work
is performed in deterministic parallel batches; a complete batch is appended
to an fsync'd JSONL stream before the summary checkpoint advances.  Each GP
subprocess has its own foreground process-group timeout, is never retried,
and is reaped by the predecessor search primitive.  An outer wall cap is
checked only between complete batches.

Every genuinely new fibre gets exact conductor first.  Subtarget conductors
then receive forced-point, H=50,000, and H=250,000 rank triage (with the
predeclared conditional H=1,000,000 escalation).  Stable numerical rank at
least 21 immediately enters saturation and exact finite-reduction
certification.  Unselected weight-four directions are explicitly left open.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Iterator, Sequence

import sympy as sp

from ek_k3 import fraction_mod, rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    prior_fermigier_parameters,
    published_preimage_digest,
    rational_digest,
    sha256_file,
)
from search_fermigier_rank22_accidental_slices import (
    T0,
    conductor_probe,
    quartic_group_pullback,
    short_add,
    short_negate,
    slice_polynomial,
)
from search_fermigier_rank22_record_group_directions import (
    EXPECTED_PUBLISHED_BASIS_POINT_SHA256,
    EXPECTED_TRANSPORT_SOURCE_SHA256,
    RecordQuarticAuxiliary,
    aggregate_candidates,
    direction_digest,
    known_record_abscissas,
    load_transport_source,
    projective_height,
    search_direction,
    transport_source_digest,
)
from search_fermigier_rank22_record_group_triples import (
    quartic_point_count_mod_prime,
)
from search_fermigier_rank22_record_group_triples_remainder import (
    staged_rank_triage,
)
from search_nagao_section7_auxiliary_jacobians import (
    weierstrass_add,
    weierstrass_multiply,
)
from triage_nagao_rank13_finalists import point_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

EXPECTED_QUAD_VECTOR_COUNT = 58_520
SLICE_HEIGHT = 50_000
MODULAR_PRIMES = (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
LOW_HEIGHT_NUMERATOR = 7
LOW_HEIGHT_DENOMINATOR = 100
MODULAR_TAIL_NUMERATOR = 1
MODULAR_TAIL_DENOMINATOR = 100
PER_FIRST_SUPPORT_HIGH_KEEP = 2
PER_RELATIVE_SIGN_TAIL_KEEP = 4

WEIGHT2_ARTIFACT = "elliptic_fermigier_rank22_record_group_directions.json"
WEIGHT3_ARTIFACT = "elliptic_fermigier_rank22_record_group_triples.json"
WEIGHT3_REMAINDER_ARTIFACT = (
    "elliptic_fermigier_rank22_record_group_triples_remainder.json"
)
WEIGHT3_REMAINDER_STREAM = (
    "elliptic_fermigier_rank22_record_group_triples_remainder_stream.jsonl"
)
EXPECTED_WEIGHT2_ARTIFACT_SHA256 = (
    "4928c44e27cada74b7a558dd97edfba20a554a508ac0e98ca051df7dea66a3c1"
)
EXPECTED_WEIGHT3_ARTIFACT_SHA256 = (
    "2803b1fa276c80eccceac5ce83215f8678d9fb771abdccb8e5043a9962b1ed36"
)
EXPECTED_WEIGHT3_REMAINDER_ARTIFACT_SHA256 = (
    "829bf44e50fb8d3583592190732ce4c9057c1bc54c9ababbea4e6cc002e6e028"
)
EXPECTED_WEIGHT3_REMAINDER_STREAM_SHA256 = (
    "66e69019bb53b28310bc1a4fa0d40989fe2fd95e677b53bca55feaafa2f3b5de"
)
EXPECTED_WEIGHT2_DIRECTION_COUNT = 462
EXPECTED_WEIGHT2_DIRECTION_SHA256 = (
    "2819df851b876fe430465e2d7fbb838c34b60d597aea2ee48c8d912b07c19939"
)
EXPECTED_WEIGHT3_DIRECTION_COUNT = 6_160
EXPECTED_WEIGHT3_DIRECTION_SHA256 = (
    "4673b556e0a60943d86b23c7f293bfc3a9952f6acd25ef5e152297e84a106455"
)
EXPECTED_WEIGHT3_SELECTED_COUNT = 399
EXPECTED_WEIGHT3_REMAINDER_COUNT = 5_761
EXPECTED_LOWER_WEIGHT_X_COUNT: int | None = 6_654
EXPECTED_LOWER_WEIGHT_X_SHA256: str | None = (
    "741591c009f0e3c8a09232e2e73a66c7155b38a0d43bc44d44e1a020fcac4796"
)
EXPECTED_FULL_QUAD_DIRECTION_SHA256: str | None = None
EXPECTED_SELECTED_QUAD_COUNT: int | None = None
EXPECTED_SELECTED_QUAD_SHA256: str | None = None
EXPECTED_SELECTION_AND_STRATA_SHA256: str | None = None
EXPECTED_PRIOR_PARAMETER_COUNT = 1_686
EXPECTED_PRIOR_PARAMETER_SHA256 = (
    "f2e060657e16e4b0b57f3ae210afba8c5147b69ed23c496416f2899339da3548"
)
KNOWN_WEIGHT3_SECOND_FIBRE_DIRECTION = "p02_p06_m12"
DEFAULT_SLICE_WALL_CAP_SECONDS = 1_800.0


def quad_vectors() -> Iterator[tuple[int, ...]]:
    """Yield the exact global-sign quotient in stable lexicographic order."""

    for support in itertools.combinations(range(22), 4):
        for relative_signs in itertools.product((-1, 1), repeat=3):
            vector = [0] * 22
            vector[support[0]] = 1
            for index, sign in zip(support[1:], relative_signs, strict=True):
                vector[index] = sign
            yield tuple(vector)


def vector_id(vector: Sequence[int]) -> str:
    return "_".join(
        f"{'p' if value > 0 else 'm'}{index:02d}"
        for index, value in enumerate(vector, start=1)
        if value
    )


def signed_auxiliary_basis(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
) -> dict[tuple[int, int], tuple[Fraction, Fraction]]:
    return {
        (index, sign): weierstrass_multiply(coefficients, point, sign)
        for index, point in enumerate(basis)
        for sign in (-1, 1)
    }


def signed_short_basis(
    basis: Sequence[tuple[Fraction, Fraction]],
) -> dict[tuple[int, int], tuple[Fraction, Fraction]]:
    return {
        (index, sign): point if sign == 1 else short_negate(point)
        for index, point in enumerate(basis)
        for sign in (-1, 1)
    }


def pair_cache(
    coefficients: Sequence[Fraction],
    signed_basis: dict[tuple[int, int], tuple[Fraction, Fraction]],
    add: Any,
) -> dict[tuple[int, int, int, int], tuple[Fraction, Fraction] | None]:
    return {
        (first, first_sign, second, second_sign): add(
            coefficients,
            signed_basis[first, first_sign],
            signed_basis[second, second_sign],
        )
        for first in range(22)
        for second in range(first + 1, 22)
        for first_sign in (-1, 1)
        for second_sign in (-1, 1)
    }


def generate_quad_directions(
    auxiliary: RecordQuarticAuxiliary,
    auxiliary_basis: Sequence[tuple[Fraction, Fraction]],
    short_basis: Sequence[tuple[Fraction, Fraction]],
    prior_x: set[Fraction],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Transport all weight-four vectors and classify every abscissa exactly."""

    auxiliary_signed = signed_auxiliary_basis(
        auxiliary.weierstrass_coefficients, auxiliary_basis
    )
    short_coefficients = FermigierMestreFamily.coefficients(T0)
    short_signed = signed_short_basis(short_basis)
    auxiliary_pairs = pair_cache(
        auxiliary.weierstrass_coefficients, auxiliary_signed, weierstrass_add
    )
    short_pairs = pair_cache(short_coefficients, short_signed, short_add)

    directions: list[dict[str, Any]] = []
    seen_x: dict[Fraction, str] = {}
    prior_exclusion_records: list[dict[str, str]] = []
    full_orbit_digest = hashlib.sha256()
    exceptional = 0
    prior_exclusions = 0
    duplicate_abscissas = 0
    vector_count = 0
    for vector in quad_vectors():
        vector_count += 1
        support = [index for index, value in enumerate(vector) if value]
        first, second, third, fourth = support
        signs = [vector[index] for index in support]
        auxiliary_point = weierstrass_add(
            auxiliary.weierstrass_coefficients,
            auxiliary_pairs[first, signs[0], second, signs[1]],
            auxiliary_pairs[third, signs[2], fourth, signs[3]],
        )
        inverse = auxiliary.inverse(auxiliary_point)
        if inverse is None:
            exceptional += 1
            continue
        expected_short = short_add(
            short_coefficients,
            short_pairs[first, signs[0], second, signs[1]],
            short_pairs[third, signs[2], fourth, signs[3]],
        )
        if quartic_group_pullback(T0, inverse) != expected_short:
            raise AssertionError("a weight-four inverse lost its short group coordinate")
        x_value, ordinate = inverse
        identifier = vector_id(vector)
        if x_value in prior_x:
            prior_exclusions += 1
            classification = "prior-l1-le-3-or-record-search-abscissa-excluded"
            prior_exclusion_records.append(
                {
                    "direction_id": identifier,
                    "quartic_x": rational_to_string(x_value),
                    "quartic_z": rational_to_string(ordinate),
                }
            )
            full_orbit_digest.update(
                (
                    f"{identifier}|{x_value}|{ordinate}|"
                    f"{','.join(map(str, vector))}|{classification}\n"
                ).encode()
            )
            continue
        if x_value in seen_x:
            duplicate_abscissas += 1
            classification = "duplicate-weight4-abscissa-excluded"
            full_orbit_digest.update(
                (
                    f"{identifier}|{x_value}|{ordinate}|"
                    f"{','.join(map(str, vector))}|{classification}\n"
                ).encode()
            )
            continue
        seen_x[x_value] = identifier
        classification = "genuinely-new-weight4-direction"
        full_orbit_digest.update(
            (
                f"{identifier}|{x_value}|{ordinate}|"
                f"{','.join(map(str, vector))}|{classification}\n"
            ).encode()
        )
        directions.append(
            {
                "direction_id": identifier,
                "coefficient_vector": list(vector),
                "quartic_x": rational_to_string(x_value),
                "quartic_z": rational_to_string(ordinate),
                "projective_height": projective_height(x_value),
                "exact_auxiliary_inverse_checked": True,
                "exact_short_group_combination_checked": True,
            }
        )
        if vector_count % 5_000 == 0:
            print(f"exact weight4 transport {vector_count}/{EXPECTED_QUAD_VECTOR_COUNT}", flush=True)

    if vector_count != EXPECTED_QUAD_VECTOR_COUNT:
        raise AssertionError("the exact weight-four vector population changed")
    directions.sort(
        key=lambda row: (
            row["projective_height"],
            Q(row["quartic_x"]),
            row["direction_id"],
        )
    )
    return directions, {
        "global_sign_quotient": True,
        "coefficient_alphabet": [-1, 0, 1],
        "exact_l1_norm": 4,
        "full_vector_count": vector_count,
        "full_vector_classified_direction_sha256": full_orbit_digest.hexdigest(),
        "exceptional_inverse_count": exceptional,
        "prior_quartic_abscissa_exclusions": prior_exclusions,
        "prior_quartic_abscissa_exclusion_records": prior_exclusion_records,
        "duplicate_new_abscissas": duplicate_abscissas,
        "genuinely_new_unique_abscissa_count": len(directions),
        "genuinely_new_direction_sha256": direction_digest(directions),
        "minimum_new_abscissa_projective_height": min(
            row["projective_height"] for row in directions
        ),
        "maximum_new_abscissa_projective_height": max(
            row["projective_height"] for row in directions
        ),
    }


def modular_lookup() -> dict[int, dict[int, tuple[int | None, ...]]]:
    """Precompute exact projective point counts by intercept residue."""

    answer: dict[int, dict[int, tuple[int | None, ...]]] = {}
    for prime in MODULAR_PRIMES:
        answer[prime] = {}
        for slope in (-1, 1):
            counts = []
            for intercept in range(prime):
                count = quartic_point_count_mod_prime(
                    slice_polynomial(slope, Q(intercept)), prime
                )
                counts.append(count)
            answer[prime][slope] = tuple(counts)
    return answer


def lookup_digest(
    lookup: dict[int, dict[int, tuple[int | None, ...]]]
) -> str:
    digest = hashlib.sha256()
    for prime in MODULAR_PRIMES:
        for slope in (-1, 1):
            digest.update(
                (
                    f"{prime}|{slope}|"
                    + ",".join(
                        "NA" if value is None else str(value)
                        for value in lookup[prime][slope]
                    )
                    + "\n"
                ).encode()
            )
    return digest.hexdigest()


def modular_profile(
    direction: dict[str, Any],
    lookup: dict[int, dict[int, tuple[int | None, ...]]],
) -> dict[str, Any]:
    source_x = Q(direction["quartic_x"])
    slope_scores = []
    coverages = []
    for slope in (-1, 1):
        intercept = source_x - slope * T0
        score = Q(0)
        usable = 0
        for prime in MODULAR_PRIMES:
            try:
                residue = fraction_mod(intercept, prime)
            except ValueError:
                continue
            point_count = lookup[prime][slope][residue]
            if point_count is None:
                continue
            score += Q(point_count - (prime + 1), prime + 1)
            usable += 1
        if usable == 0:
            raise AssertionError("a direction lost every modular pilot prime")
        slope_scores.append(score / usable)
        coverages.append(usable)
    return {
        "maximum_slope_score": rational_to_string(max(slope_scores)),
        "sum_slope_score": rational_to_string(sum(slope_scores, Q(0))),
        "usable_prime_counts_by_slope_m1_p1": coverages,
    }


def modular_sort_key(record: dict[str, Any]) -> tuple[Any, ...]:
    profile = record["modular_square_profile"]
    return (
        -Q(profile["maximum_slope_score"]),
        -Q(profile["sum_slope_score"]),
        record["projective_height"],
        Q(record["quartic_x"]),
        record["direction_id"],
    )


def modular_profile_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda row: row["direction_id"]):
        profile = record["modular_square_profile"]
        digest.update(
            (
                f"{record['direction_id']}|{profile['maximum_slope_score']}|"
                f"{profile['sum_slope_score']}|"
                f"{','.join(map(str, profile['usable_prime_counts_by_slope_m1_p1']))}\n"
            ).encode()
        )
    return digest.hexdigest()


def selection_and_strata_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda row: row["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{','.join(record['selection_strata'])}\n"
            ).encode()
        )
    return digest.hexdigest()


def select_directions(
    directions: Sequence[dict[str, Any]],
    lookup: dict[int, dict[int, tuple[int | None, ...]]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    scored = [dict(row) for row in directions]
    for row in scored:
        row["modular_square_profile"] = modular_profile(row, lookup)

    count = len(scored)
    low_height_keep = math.ceil(
        count * LOW_HEIGHT_NUMERATOR / LOW_HEIGHT_DENOMINATOR
    )
    modular_tail_keep = math.ceil(
        count * MODULAR_TAIL_NUMERATOR / MODULAR_TAIL_DENOMINATOR
    )
    reasons: dict[str, set[str]] = defaultdict(set)
    by_id = {row["direction_id"]: row for row in scored}

    height_order = sorted(
        scored,
        key=lambda row: (
            row["projective_height"], Q(row["quartic_x"]), row["direction_id"]
        ),
    )
    modular_order = sorted(scored, key=modular_sort_key)
    for row in height_order[:low_height_keep]:
        reasons[row["direction_id"]].add("lowest-seven-percent-abscissa-height")
    for row in modular_order[:modular_tail_keep]:
        reasons[row["direction_id"]].add("highest-one-percent-modular-square-yield")
    for row in modular_order[-modular_tail_keep:]:
        reasons[row["direction_id"]].add("lowest-one-percent-modular-square-yield")

    by_first: dict[int, list[dict[str, Any]]] = defaultdict(list)
    by_pattern: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in scored:
        support = [
            (index, value)
            for index, value in enumerate(row["coefficient_vector"])
            if value
        ]
        by_first[support[0][0]].append(row)
        by_pattern[tuple(value for _, value in support[1:])].append(row)
    for first, population in sorted(by_first.items()):
        for row in sorted(population, key=modular_sort_key)[
            :PER_FIRST_SUPPORT_HIGH_KEEP
        ]:
            reasons[row["direction_id"]].add(
                f"first-support-{first + 1:02d}-modular-diversity"
            )
    for pattern, population in sorted(by_pattern.items()):
        ordered = sorted(population, key=modular_sort_key)
        label = "_".join("p" if value > 0 else "m" for value in pattern)
        for row in ordered[:PER_RELATIVE_SIGN_TAIL_KEEP]:
            reasons[row["direction_id"]].add(f"sign-{label}-high-modular-tail")
        for row in ordered[-PER_RELATIVE_SIGN_TAIL_KEEP:]:
            reasons[row["direction_id"]].add(f"sign-{label}-low-modular-tail")

    selected = []
    for identifier, strata in reasons.items():
        row = dict(by_id[identifier])
        row["selection_strata"] = sorted(strata)
        selected.append(row)
    selected.sort(
        key=lambda row: (
            row["projective_height"], Q(row["quartic_x"]), row["direction_id"]
        )
    )
    return selected, {
        "selection_uses_weight4_specialized_rank": False,
        "selection_uses_weight4_conductor": False,
        "selection_uses_weight4_slice_outcomes": False,
        "modular_primes": list(MODULAR_PRIMES),
        "modular_lookup_sha256": lookup_digest(lookup),
        "modular_score": (
            "for each slope, exact average of (#C(F_p)-(p+1))/(p+1) over "
            "usable primes; max slope then sum slopes gives descending order"
        ),
        "lowest_height_fraction": (
            f"{LOW_HEIGHT_NUMERATOR}/{LOW_HEIGHT_DENOMINATOR} rounded up"
        ),
        "modular_high_tail_fraction": (
            f"{MODULAR_TAIL_NUMERATOR}/{MODULAR_TAIL_DENOMINATOR} rounded up"
        ),
        "modular_low_tail_fraction": (
            f"{MODULAR_TAIL_NUMERATOR}/{MODULAR_TAIL_DENOMINATOR} rounded up"
        ),
        "per_first_support_high_modular_keep": PER_FIRST_SUPPORT_HIGH_KEEP,
        "per_relative_sign_pattern_each_modular_tail_keep": (
            PER_RELATIVE_SIGN_TAIL_KEEP
        ),
        "full_modular_profile_sha256": modular_profile_digest(scored),
        "selected_direction_count": len(selected),
        "selected_direction_sha256": direction_digest(selected),
        "selection_and_strata_sha256": selection_and_strata_digest(selected),
        "unselected_direction_count": len(scored) - len(selected),
    }


def load_lower_weight_directions(
    artifact_directory: Path,
    primary: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], set[Fraction], dict[str, Any]]:
    weight2_path = artifact_directory / WEIGHT2_ARTIFACT
    weight3_path = artifact_directory / WEIGHT3_ARTIFACT
    remainder_path = artifact_directory / WEIGHT3_REMAINDER_ARTIFACT
    stream_path = artifact_directory / WEIGHT3_REMAINDER_STREAM
    observed = {
        WEIGHT2_ARTIFACT: sha256_file(weight2_path),
        WEIGHT3_ARTIFACT: sha256_file(weight3_path),
        WEIGHT3_REMAINDER_ARTIFACT: sha256_file(remainder_path),
        WEIGHT3_REMAINDER_STREAM: sha256_file(stream_path),
    }
    expected = {
        WEIGHT2_ARTIFACT: EXPECTED_WEIGHT2_ARTIFACT_SHA256,
        WEIGHT3_ARTIFACT: EXPECTED_WEIGHT3_ARTIFACT_SHA256,
        WEIGHT3_REMAINDER_ARTIFACT: EXPECTED_WEIGHT3_REMAINDER_ARTIFACT_SHA256,
        WEIGHT3_REMAINDER_STREAM: EXPECTED_WEIGHT3_REMAINDER_STREAM_SHA256,
    }
    if observed != expected:
        raise AssertionError("a frozen lower-weight source artifact changed")

    weight2_artifact = json.loads(weight2_path.read_text())
    weight3_artifact = json.loads(weight3_path.read_text())
    remainder_artifact = json.loads(remainder_path.read_text())
    weight2 = [
        {key: value for key, value in row.items() if key != "slice_searches"}
        for row in weight2_artifact["direction_searches"]
    ]
    selected_weight3 = [dict(row) for row in weight3_artifact["selected_directions"]]
    remainder_weight3 = []
    with stream_path.open() as handle:
        for line in handle:
            row = json.loads(line)
            remainder_weight3.append(
                {key: value for key, value in row.items() if key != "slice_searches"}
            )
    weight3 = selected_weight3 + remainder_weight3
    if (
        len(weight2) != EXPECTED_WEIGHT2_DIRECTION_COUNT
        or direction_digest(weight2) != EXPECTED_WEIGHT2_DIRECTION_SHA256
        or len(selected_weight3) != EXPECTED_WEIGHT3_SELECTED_COUNT
        or len(remainder_weight3) != EXPECTED_WEIGHT3_REMAINDER_COUNT
        or len(weight3) != EXPECTED_WEIGHT3_DIRECTION_COUNT
        or direction_digest(weight3) != EXPECTED_WEIGHT3_DIRECTION_SHA256
        or remainder_artifact["outcome"]["full_remainder_exhausted"] is not True
    ):
        raise AssertionError("the exact lower-weight direction population changed")

    known_x, known_record = known_record_abscissas(primary)
    lower_x = known_x | {Q(row["quartic_x"]) for row in weight2 + weight3}
    lower_digest = rational_digest(sorted(lower_x))
    if EXPECTED_LOWER_WEIGHT_X_COUNT is not None and (
        len(lower_x) != EXPECTED_LOWER_WEIGHT_X_COUNT
        or lower_digest != EXPECTED_LOWER_WEIGHT_X_SHA256
    ):
        raise AssertionError("the exact l1<=3 record-abscissa union changed")
    return weight2, weight3, lower_x, {
        "source_artifact_sha256": observed,
        "known_record_abscissas": known_record,
        "weight2_direction_count": len(weight2),
        "weight2_direction_sha256": direction_digest(weight2),
        "weight3_direction_count": len(weight3),
        "weight3_direction_sha256": direction_digest(weight3),
        "full_lower_weight_x_count": len(lower_x),
        "full_lower_weight_x_sha256": lower_digest,
    }


def selector_calibration(
    weight2: Sequence[dict[str, Any]],
    weight3: Sequence[dict[str, Any]],
    lookup: dict[int, dict[int, tuple[int | None, ...]]],
) -> dict[str, Any]:
    selected2, audit2 = select_directions(weight2, lookup)
    selected3, audit3 = select_directions(weight3, lookup)
    selected3_ids = {row["direction_id"] for row in selected3}
    height3 = sorted(
        weight3,
        key=lambda row: (
            row["projective_height"], Q(row["quartic_x"]), row["direction_id"]
        ),
    )
    scored3 = [dict(row) for row in weight3]
    for row in scored3:
        row["modular_square_profile"] = modular_profile(row, lookup)
    modular3 = sorted(scored3, key=modular_sort_key)
    hit_height_rank = next(
        index
        for index, row in enumerate(height3, start=1)
        if row["direction_id"] == KNOWN_WEIGHT3_SECOND_FIBRE_DIRECTION
    )
    hit_modular_rank = next(
        index
        for index, row in enumerate(modular3, start=1)
        if row["direction_id"] == KNOWN_WEIGHT3_SECOND_FIBRE_DIRECTION
    )
    if KNOWN_WEIGHT3_SECOND_FIBRE_DIRECTION not in selected3_ids:
        raise AssertionError("the declared lower-weight calibration was not retained")
    return {
        "fully_searched_weight2_population": len(weight2),
        "known_weight2_second_fibre_count": 0,
        "selector_weight2_retained_count": len(selected2),
        "selector_weight2_audit": audit2,
        "fully_searched_weight3_population": len(weight3),
        "known_weight3_second_fibre_count": 1,
        "known_weight3_second_fibre_direction": KNOWN_WEIGHT3_SECOND_FIBRE_DIRECTION,
        "known_weight3_second_fibre_parameter": "29771/78",
        "selector_weight3_retained_count": len(selected3),
        "known_hit_retained": True,
        "known_hit_height_rank_one_based": hit_height_rank,
        "known_hit_descending_modular_rank_one_based": hit_modular_rank,
        "selector_weight3_audit": audit3,
        "interpretation": (
            "the known hit is retained by the height stratum; the modular score "
            "alone ranks it poorly, so modular tails are used only as a diversified "
            "prefilter and are not claimed predictive"
        ),
        "selection_leakage": (
            "in-sample cross-weight calibration: selector quotas were fixed after "
            "the completed weight3 outcome was known; no weight4 slice, conductor, "
            "or rank outcome entered selection"
        ),
    }


def stream_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stream_result_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open() as handle:
        for line in handle:
            direction = json.loads(line)
            for row in direction["slice_searches"]:
                search = row["search"]
                digest.update(
                    (
                        f"{row['slice_id']}|{row['quartic_polynomial_sha256']}|"
                        f"{search['status']}|{search.get('signed_point_count')}|"
                        f"{row['record_T0_calibration_count']}|"
                        f"{len(row['incidences'])}\n"
                    ).encode()
                )
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--slice-height", type=int, default=SLICE_HEIGHT)
    parser.add_argument(
        "--slice-wall-cap", type=float, default=DEFAULT_SLICE_WALL_CAP_SECONDS
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--slice-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--h50-timeout", type=float, default=30.0)
    parser.add_argument("--h250-timeout", type=float, default=45.0)
    parser.add_argument("--h1m-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_fermigier_rank22_record_group_quads.json",
    )
    parser.add_argument(
        "--stream-output",
        type=Path,
        default=generated
        / "elliptic_fermigier_rank22_record_group_quads_stream.jsonl",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != SLICE_HEIGHT:
        raise SystemExit("this weight-four tranche is pinned at slice H=50000")
    if not 60 <= args.slice_wall_cap <= 3_600:
        raise SystemExit("the slice wall cap must lie in [60,3600] seconds")
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    caps = (
        args.slice_timeout,
        args.conductor_timeout,
        args.h50_timeout,
        args.h250_timeout,
        args.h1m_timeout,
        args.height_timeout,
        args.saturation_timeout,
    )
    if min(caps) <= 0 or max(caps) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")
    if args.output.exists():
        raise SystemExit("refusing to overwrite a weight-four artifact")
    if not args.selection_only and args.stream_output.exists():
        raise SystemExit("refusing to overwrite a weight-four stream")

    root = Path(__file__).resolve().parents[2]
    artifact_directory = root / "artifacts" / "generated-results"
    primary_path = artifact_directory / PRIMARY_ARTIFACT
    primary = json.loads(primary_path.read_text())
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published accidental preimages changed")
    auxiliary = RecordQuarticAuxiliary.construct()
    auxiliary_basis, short_basis, _ = load_transport_source(primary, auxiliary)
    if (
        transport_source_digest(primary["published_point_preimages"])
        != EXPECTED_TRANSPORT_SOURCE_SHA256
        or point_digest(short_basis) != EXPECTED_PUBLISHED_BASIS_POINT_SHA256
    ):
        raise AssertionError("the transported certified rank22 basis changed")

    started_preflight = time.monotonic()
    weight2, weight3, lower_x, lower_record = load_lower_weight_directions(
        artifact_directory, primary
    )
    directions, population = generate_quad_directions(
        auxiliary, auxiliary_basis, short_basis, lower_x
    )
    if (
        population["full_vector_count"] != EXPECTED_QUAD_VECTOR_COUNT
        or population["exceptional_inverse_count"] != 0
        or population["prior_quartic_abscissa_exclusions"] != 1
        or population["duplicate_new_abscissas"] != 0
        or len(directions) != EXPECTED_QUAD_VECTOR_COUNT - 1
    ):
        raise AssertionError(
            f"the exact weight-four direction population changed: {population}"
        )
    if (
        EXPECTED_FULL_QUAD_DIRECTION_SHA256 is not None
        and population["genuinely_new_direction_sha256"]
        != EXPECTED_FULL_QUAD_DIRECTION_SHA256
    ):
        raise AssertionError("the exact weight-four direction digest changed")

    lookup = modular_lookup()
    calibration = selector_calibration(weight2, weight3, lookup)
    selected, selection = select_directions(directions, lookup)
    if EXPECTED_SELECTED_QUAD_COUNT is not None and (
        len(selected) != EXPECTED_SELECTED_QUAD_COUNT
        or direction_digest(selected) != EXPECTED_SELECTED_QUAD_SHA256
        or selection_and_strata_digest(selected)
        != EXPECTED_SELECTION_AND_STRATA_SHA256
    ):
        raise AssertionError("the pinned weight-four selection changed")

    prior_parameters, prior_record = prior_fermigier_parameters(
        artifact_directory, args.output
    )
    if (
        len(prior_parameters) != EXPECTED_PRIOR_PARAMETER_COUNT
        or rational_digest(sorted(prior_parameters))
        != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the complete prior Fermigier parameter set changed")
    preflight_seconds = time.monotonic() - started_preflight

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "selection-only exact weight4 orbit prepared"
            if args.selection_only
            else "in-progress selected exact weight4 manufactured-slice search"
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
        },
        "lower_weight_direction_exclusion": lower_record,
        "full_weight4_population": population,
        "rank_and_conductor_blind_selection": selection,
        "lower_weight_selector_calibration": calibration,
        "selected_directions": selected,
        "prior_parameter_decontamination": prior_record,
        "parameters": {
            "exact_l1_norm": 4,
            "coefficient_alphabet": [-1, 0, 1],
            "global_sign_quotient": True,
            "selected_direction_count": len(selected),
            "unselected_direction_count": len(directions) - len(selected),
            "slopes_per_direction": [-1, 1],
            "declared_slice_call_count": 2 * len(selected),
            "slice_height": args.slice_height,
            "slice_wall_cap_seconds": args.slice_wall_cap,
            "wall_cap_checked_only_between_complete_parallel_batches": True,
            "parallel_workers": args.workers,
            "slice_timeout_seconds": args.slice_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "H50000_timeout_seconds": args.h50_timeout,
            "H250000_timeout_seconds": args.h250_timeout,
            "H1000000_timeout_seconds": args.h1m_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": [72, 120],
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "stack_bytes_per_worker": args.stack_bytes,
            "no_retries": True,
            "append_only_checkpoint_after_each_complete_parallel_batch": True,
        },
        "stream": None,
        "candidates": [],
        "execution": {
            "phase": (
                "selection-only-complete" if args.selection_only else "slice-search-in-progress"
            ),
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

    artifact["stream"] = {
        "path": str(args.stream_output),
        "format": "one canonical JSON completed direction record per line",
        "created_exclusive": True,
    }
    write_artifact(args.output, artifact)
    qualifying_all: list[dict[str, Any]] = []
    classifications = Counter()
    slice_completed = 0
    slice_failed = 0
    calibrated = 0
    completed = 0
    started_search = time.monotonic()
    args.stream_output.parent.mkdir(parents=True, exist_ok=True)
    with args.stream_output.open("x") as stream, ThreadPoolExecutor(
        max_workers=args.workers
    ) as executor:
        for batch_start in range(0, len(selected), args.workers):
            if time.monotonic() - started_search >= args.slice_wall_cap:
                artifact["execution"]["wall_cap_triggered"] = True
                break
            batch = selected[batch_start : batch_start + args.workers]
            futures = [
                executor.submit(
                    search_direction,
                    direction,
                    height_bound=args.slice_height,
                    timeout=args.slice_timeout,
                    stack_bytes=args.stack_bytes,
                    prior_parameters=prior_parameters,
                )
                for direction in batch
            ]
            results = [future.result() for future in futures]
            for direction_row, qualifying in results:
                stream.write(json.dumps(direction_row, sort_keys=True) + "\n")
                qualifying_all.extend(qualifying)
                completed += 1
                for row in direction_row["slice_searches"]:
                    if row["search"]["status"] == "completed":
                        slice_completed += 1
                    else:
                        slice_failed += 1
                    calibrated += row["record_T0_calibration_count"] == 1
                    classifications.update(
                        incidence["classification"] for incidence in row["incidences"]
                    )
            stream.flush()
            os.fsync(stream.fileno())
            artifact["execution"].update(
                {
                    "directions_completed": completed,
                    "slice_calls_completed_or_attempted": completed * 2,
                    "last_direction_id": batch[-1]["direction_id"],
                    "wall_seconds_so_far": time.monotonic() - started_search,
                    "wall_cap_triggered": False,
                }
            )
            if completed % 100 < args.workers or completed == len(selected):
                write_artifact(args.output, artifact)
                print(
                    f"weight4 selected {completed}/{len(selected)} "
                    f"qualifying_incidences={len(qualifying_all)}",
                    flush=True,
                )

    slice_phase_seconds = time.monotonic() - started_search
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
                candidate["staged_rank_triage"] = staged_rank_triage(
                    candidate,
                    h50_timeout=args.h50_timeout,
                    h250_timeout=args.h250_timeout,
                    h1m_timeout=args.h1m_timeout,
                    height_timeout=args.height_timeout,
                    stack_bytes=args.stack_bytes,
                    precisions=(72, 120),
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except subprocess.TimeoutExpired as error:
                candidate["staged_rank_triage"] = {
                    "status": "timeout-no-retry",
                    "error": str(error)[:1000],
                }
            except (RuntimeError, AssertionError, ValueError) as error:
                candidate["staged_rank_triage"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            write_artifact(args.output, artifact)

    open_count = len(selected) - completed
    rank_records = [
        candidate["staged_rank_triage"]
        for candidate in candidates
        if "maximum_stable_numerical_rank"
        in candidate.get("staged_rank_triage", {})
    ]
    artifact["stream"].update(
        {
            "completed_direction_line_count": completed,
            "sha256": stream_sha256(args.stream_output),
            "exact_slice_result_sha256": stream_result_digest(args.stream_output),
        }
    )
    artifact["outcome"] = {
        "full_exact_weight4_direction_count": len(directions),
        "declared_selected_direction_count": len(selected),
        "completed_selected_direction_count": completed,
        "open_selected_direction_count": open_count,
        "selected_tranche_exhausted": open_count == 0,
        "unselected_direction_count": len(directions) - len(selected),
        "unselected_directions_claimed_negative": False,
        "slice_calls_attempted": completed * 2,
        "slice_calls_completed": slice_completed,
        "slice_calls_timed_out_or_errored": slice_failed,
        "record_T0_calibrated_slices": calibrated,
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
                record["maximum_stable_numerical_rank"]
                for record in rank_records
            ),
            default=None,
        ),
        "slice_phase_wall_seconds": slice_phase_seconds,
    }
    artifact["target"]["hit"] = any(
        stage.get("finite_reduction_attempt", {}).get(
            "certified_algebraic_rank_lower_bound", 0
        )
        >= 21
        and candidate.get("conductor_probe", {}).get(
            "below_strict_log_conductor_target"
        )
        for candidate in candidates
        for stage in candidate.get("staged_rank_triage", {}).get("stages", [])
    )
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no completed selected weight4 direction produced a certified "
            "subtarget rank21 fibre"
        )
    artifact["status"] = (
        "completed selected exact weight4 manufactured-slice tranche"
        if open_count == 0
        else "safe-wall-cap checkpoint of selected exact weight4 tranche"
    )
    artifact["execution"].update(
        {
            "phase": "complete",
            "wall_seconds": time.monotonic() - started_search,
            "owned_processes_remaining": 0,
        }
    )
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()

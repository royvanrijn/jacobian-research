#!/usr/bin/env python3
"""One-pass simultaneous-square search on Fermigier exceptional fiber products.

The pinned structural artifact contains 80 affine genus-two covers

    y_i^2 = f_i(T),

and all 3,160 unordered genus-nine fiber products.  This program searches
each ``f_i`` once in the exact projective box ``H(T) <= 200000`` and then
intersects the returned parameter sets.  Consequently a retained incidence
has two separately checked rational square roots; a square root of
``f_i(T) f_j(T)`` is never used as a surrogate.

Before starting PARI, the program constructs explicit local-solubility
bitsets on P^1(F_p) for the fixed primes 13, 23, 37, 41, and 43.  The bitsets
are replayed against every returned point.  PARI's bounded search remains an
experiment, not a proof that no rational point exists outside the box.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from ecsearch.fermigier_rank import (
    certify_fermigier_rank_sections,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (
    IndependenceCertificate,
    select_independent_subset,
    verify_independence_certificate,
)
from ek_k3 import rational_square_root, rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_rank22_accidental_slices import (
    conductor_probe,
    search_polynomial,
)
from analyze_fermigier_exceptional_transport import family_expression


Q = Fraction
HEIGHT_BOUND = 200_000
SIEVE_PRIMES = (13, 23, 37, 41, 43)
WORKERS = 4
SEARCH_TIMEOUT_SECONDS = 30.0
PARI_STACK_BYTES = 512_000_000
SPECIALIZATION_HEIGHT_BOUND = 50_000
SPECIALIZATION_TIMEOUT_SECONDS = 30.0
CERTIFICATE_PRIME_BOUND = 2_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXPECTED_DIRECTION_COUNT = 80
EXPECTED_PAIR_COUNT = 3_160
E22_T = Q(39_508, 39)
RANK20_T = Q(28_917, 10)
TRANSPORT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_fermigier_exceptional_transport.json"
)
TRANSPORT_SHA256 = "a767e849119d4eb974eb8e85536031413c6d52a59151933239fa141235de5777"
TRANSPORT_RESULT_SHA256 = (
    "db07f28e39c73c4e66fdf29bf11d652bb8714525e93f05f819ec1477107a4d0c"
)

T_SCALAR_KEYS = {
    "T",
    "t",
    "candidate_t",
    "canonical_parameter_t",
    "escalated_parameter_t",
    "internal_normalized_parameter",
    "known_weight3_second_fibre_parameter",
    "literal_shift",
    "literal_shift_T",
    "normalized_parameter",
    "normalized_parameter_t",
    "normalized_record_parameter",
    "parameter",
    "parameter_t",
    "published_parameter",
    "raw_parameter_t",
    "record_parameter_normalized_T",
    "record_parameter_t",
    "signed_parameter_t",
    "specialized_quartic_frontier_parameter_t",
}
T_LIST_KEYS = {
    "H50000_seen_parameters",
    "baseline_Fermigier_only_parameters",
    "batch_rank_triage_parameters",
    "candidate_parameters",
    "canonical_candidate_parameters",
    "current_slice_parameters_including_record",
    "deep_tranche_parameters",
    "dense_excluded_parameters",
    "excluded_parameters",
    "fermigier_accidental_slice_parameters_including_benchmark",
    "genuinely_new_parameters",
    "legacy_accidental_parameters",
    "legacy_fermigier_accidental_slice_parameters",
    "legacy_fermigier_parameters",
    "legacy_superseded_slice_parameters",
    "multiple_source_parameters",
    "new_parameters",
    "parameters",
    "parameters_t",
    "prior_parameters",
    "prior_parameters_intersecting_population",
    "raw_parameters",
    "signed_parameters",
    "signed_slice_parameters",
    "unique_parameters",
}
U_SCALAR_KEYS = {
    "adapter_parameter",
    "adapter_u",
    "canonical_adapter_u",
    "canonical_parameter_u",
    "imported_seed_adapter_u",
}


@dataclass(frozen=True)
class Direction:
    identifier: str
    e22_label: str
    rank20_label: str
    e22_x: Fraction
    rank20_x: Fraction
    slope: Fraction
    intercept: Fraction
    polynomial: tuple[Fraction, ...]

    def x_value(self, parameter: Fraction) -> Fraction:
        return self.slope * Q(parameter) + self.intercept


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def rational_digest(values: Iterable[Fraction]) -> str:
    return sha256_lines(rational_to_string(value) for value in values)


def polynomial_digest(coefficients: Sequence[Fraction]) -> str:
    return rational_digest(Q(value) for value in coefficients)


def point_text(point: tuple[Fraction, Fraction]) -> dict[str, str | bool]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
        "exact_membership_checked": True,
    }


def parse_scalar(value: Any) -> Fraction | None:
    if not isinstance(value, (str, int)) or isinstance(value, bool):
        return None
    try:
        return Q(value)
    except (ValueError, ZeroDivisionError):
        return None


def extract_literal_parameters(value: Any, *, key: str | None = None) -> set[Fraction]:
    """Conservatively extract legacy literal-shift T values from prior artifacts."""

    answer: set[Fraction] = set()
    if key in T_SCALAR_KEYS:
        parsed = parse_scalar(value)
        if parsed is not None:
            answer.add(abs(parsed))
    elif key in U_SCALAR_KEYS:
        parsed = parse_scalar(value)
        if parsed is not None:
            answer.add(2 * abs(parsed))
    elif key in T_LIST_KEYS and isinstance(value, list):
        for item in value:
            parsed = parse_scalar(item)
            if parsed is not None:
                answer.add(abs(parsed))

    if isinstance(value, dict):
        coordinate = value.get("coordinate")
        alias_value = parse_scalar(value.get("value"))
        if alias_value is not None:
            if coordinate == "adapter_u":
                answer.add(2 * abs(alias_value))
            elif coordinate in {"literal_symmetric_shift_s", "literal_shift_T"}:
                answer.add(abs(alias_value))
        if value.get("name") == "u":
            canonical = parse_scalar(value.get("value"))
            if canonical is not None:
                answer.add(2 * abs(canonical))
        for child_key, child in value.items():
            answer.update(extract_literal_parameters(child, key=child_key))
    elif isinstance(value, list):
        for child in value:
            answer.update(extract_literal_parameters(child))
    return answer


def prior_parameter_snapshot(
    root: Path, output_path: Path
) -> tuple[set[Fraction], dict[str, Any]]:
    artifact_dir = root / "artifacts/generated-results"
    paths = {
        *artifact_dir.glob("elliptic_fermigier*.json"),
        *artifact_dir.glob("elliptic_curve_candidate_fermigier*.json"),
        *(artifact_dir / "elliptic-curves").glob("fermigier*.json"),
    }
    output_resolved = output_path.resolve()
    parameters = {abs(E22_T), abs(RANK20_T)}
    sources: dict[str, Any] = {}
    for path in sorted(paths):
        if path.resolve() == output_resolved:
            continue
        raw = path.read_bytes()
        data = json.loads(raw)
        extracted = extract_literal_parameters(data)
        parameters.update(extracted)
        sources[str(path.relative_to(root))] = {
            "sha256": sha256_bytes(raw),
            "extracted_literal_T_count": len(extracted),
        }
    ordered = tuple(sorted(parameters))
    return parameters, {
        "coordinate": "absolute legacy literal shift T=2*u",
        "source_count": len(sources),
        "sources": sources,
        "unique_prior_parameter_count": len(ordered),
        "prior_parameter_sha256": rational_digest(ordered),
    }


def load_transport(root: Path) -> tuple[dict[str, Any], bytes]:
    path = root / TRANSPORT_RELATIVE
    raw = path.read_bytes()
    if sha256_bytes(raw) != TRANSPORT_SHA256:
        raise AssertionError("the pinned exceptional-transport artifact changed")
    data = json.loads(raw)
    if data["result_sha256"] != TRANSPORT_RESULT_SHA256:
        raise AssertionError("the pinned exceptional-transport result digest changed")
    products = data["fiber_products"]
    if products["pair_count"] != EXPECTED_PAIR_COUNT:
        raise AssertionError("the pinned fiber-product population changed")
    if products["histogram"] != [
        {
            "common_branch_gcd_degree": 0,
            "count": EXPECTED_PAIR_COUNT,
            "fiber_product_genus": 9,
            "third_quotient_genus": 5,
            "third_quotient_squarefree_degree": 12,
        }
    ]:
        raise AssertionError("the pinned genus-nine classification changed")
    return data, raw


def _label_number(label: str) -> int:
    return int("".join(character for character in label if character.isdigit()))


def build_directions(transport: dict[str, Any]) -> tuple[Direction, ...]:
    e22_map = transport["exceptional_quotients"]["E22"]["accidental_source_x"]
    independent = transport["exceptional_quotients"]["E22"][
        "independent_exceptional_labels_modulo_generic"
    ]
    rank20_map = transport["exceptional_quotients"]["rank20"][
        "independent_exceptional_preimages"
    ]
    if independent != [f"P{index}" for index in range(13, 23)]:
        raise AssertionError("the independent E22 endpoint population changed")
    e22 = sorted(((label, Q(e22_map[label])) for label in independent), key=lambda row: _label_number(row[0]))
    rank20 = sorted(((label, Q(value)) for label, value in rank20_map.items()), key=lambda row: _label_number(row[0]))
    if len(e22) != 10 or len(rank20) != 8:
        raise AssertionError("the exceptional endpoint dimensions changed")

    T, X = sp.symbols("T X")
    family = family_expression(T, X)
    first = sp.Rational(E22_T.numerator, E22_T.denominator)
    second = sp.Rational(RANK20_T.numerator, RANK20_T.denominator)
    directions: list[Direction] = []
    for left_label, left_x_q in e22:
        for right_label, right_x_q in rank20:
            left_x = sp.Rational(left_x_q.numerator, left_x_q.denominator)
            right_x = sp.Rational(right_x_q.numerator, right_x_q.denominator)
            line = sp.cancel(left_x + (right_x - left_x) * (T - first) / (second - first))
            polynomial = sp.Poly(sp.expand(family.subs(X, line)), T, domain=sp.QQ)
            if polynomial.degree() != 6:
                raise AssertionError("an affine transport lost degree six")
            coefficients = tuple(
                Q(int(value.p), int(value.q))
                for value in reversed(polynomial.all_coeffs())
            )
            slope = Q(int(sp.diff(line, T).p), int(sp.diff(line, T).q))
            intercept_sp = line.subs(T, 0)
            intercept = Q(int(intercept_sp.p), int(intercept_sp.q))
            direction = Direction(
                identifier=f"{left_label}__{right_label}",
                e22_label=left_label,
                rank20_label=right_label,
                e22_x=left_x_q,
                rank20_x=right_x_q,
                slope=slope,
                intercept=intercept,
                polynomial=coefficients,
            )
            for parameter in (E22_T, RANK20_T):
                if rational_square_root(evaluate_polynomial(coefficients, parameter)) is None:
                    raise AssertionError("an affine equation lost an anchor square")
            directions.append(direction)
    directions.sort(key=lambda item: (_label_number(item.e22_label), _label_number(item.rank20_label)))
    if len(directions) != EXPECTED_DIRECTION_COUNT:
        raise AssertionError("the 80-direction population changed")
    if len(tuple(combinations(directions, 2))) != EXPECTED_PAIR_COUNT:
        raise AssertionError("the 3,160-pair population changed")
    return tuple(directions)


def evaluate_polynomial(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def fraction_mod(value: Fraction, prime: int) -> int:
    if value.denominator % prime == 0:
        raise ValueError(f"coefficient denominator is not a unit modulo {prime}")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def local_square_bitset(coefficients: Sequence[Fraction], prime: int) -> int:
    """Return allowed P^1(F_p) residues; index p denotes infinity."""

    degree = len(coefficients) - 1
    modular = tuple(fraction_mod(value, prime) for value in coefficients)
    mask = 0
    for index in range(prime + 1):
        if index == prime:
            a_value, b_value = 1, 0
        else:
            a_value, b_value = index, 1
        value = sum(
            coefficient
            * pow(a_value, exponent, prime)
            * pow(b_value, degree - exponent, prime)
            for exponent, coefficient in enumerate(modular)
        ) % prime
        if value == 0 or pow(value, (prime - 1) // 2, prime) == 1:
            mask |= 1 << index
    return mask


def projective_residue(value: Fraction, prime: int) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    if denominator == 0:
        if numerator == 0:
            raise AssertionError("a reduced rational vanished projectively")
        return prime
    return numerator * pow(denominator, -1, prime) % prime


def build_local_sieve(directions: Sequence[Direction]) -> tuple[dict[str, tuple[int, ...]], dict[str, Any]]:
    masks = {
        direction.identifier: tuple(
            local_square_bitset(direction.polynomial, prime)
            for prime in SIEVE_PRIMES
        )
        for direction in directions
    }
    for direction in directions:
        for anchor in (E22_T, RANK20_T):
            for prime, mask in zip(SIEVE_PRIMES, masks[direction.identifier], strict=True):
                if not (mask >> projective_residue(anchor, prime)) & 1:
                    raise AssertionError("an anchor failed the pre-search local sieve")

    pair_histogram: Counter[tuple[int, ...]] = Counter()
    manifest_lines = []
    for left, right in combinations(directions, 2):
        intersections = tuple(
            left_mask & right_mask
            for left_mask, right_mask in zip(
                masks[left.identifier], masks[right.identifier], strict=True
            )
        )
        counts = tuple(mask.bit_count() for mask in intersections)
        if any(count == 0 for count in counts):
            raise AssertionError("an anchor-calibrated pair became locally insoluble")
        pair_histogram[counts] += 1
        manifest_lines.append(f"{left.identifier}|{right.identifier}|{counts}")

    return masks, {
        "built_before_any_bounded_search": True,
        "definition": (
            "for every f_i and p, bit r is set iff the degree-six homogenization "
            "of f_i at r in P^1(F_p) is zero or a quadratic residue"
        ),
        "primes": list(SIEVE_PRIMES),
        "direction_bitsets_hex": {
            identifier: {
                str(prime): hex(mask)
                for prime, mask in zip(SIEVE_PRIMES, direction_masks, strict=True)
            }
            for identifier, direction_masks in sorted(masks.items())
        },
        "pair_intersection_count_histogram": [
            {"allowed_residue_counts": list(counts), "pair_count": count}
            for counts, count in sorted(pair_histogram.items())
        ],
        "pair_count": sum(pair_histogram.values()),
        "pair_manifest_sha256": sha256_lines(manifest_lines),
        "interpretation": (
            "a necessary local square test on both factors separately; it is not "
            "a product-square test and it does not prove global solubility"
        ),
    }


def unique_points_by_parameter(
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    by_parameter: dict[Fraction, Fraction] = {}
    for parameter, ordinate in points:
        parameter, ordinate = Q(parameter), abs(Q(ordinate))
        existing = by_parameter.get(parameter)
        if existing is not None and existing * existing != ordinate * ordinate:
            raise AssertionError("one parameter had inconsistent ordinates")
        by_parameter[parameter] = ordinate
    return tuple(sorted(by_parameter.items()))


def search_direction(
    direction: Direction,
    masks: tuple[int, ...],
) -> tuple[str, tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    points, record = search_polynomial(
        direction.polynomial,
        height_bound=HEIGHT_BOUND,
        timeout=SEARCH_TIMEOUT_SECONDS,
        stack_bytes=PARI_STACK_BYTES,
    )
    unique = unique_points_by_parameter(points)
    for parameter, ordinate in unique:
        if max(abs(parameter.numerator), parameter.denominator) > HEIGHT_BOUND:
            raise AssertionError("PARI returned a parameter outside the pinned box")
        if evaluate_polynomial(direction.polynomial, parameter) != ordinate * ordinate:
            raise AssertionError("PARI returned a point off an individual cover")
        for prime, mask in zip(SIEVE_PRIMES, masks, strict=True):
            if not (mask >> projective_residue(parameter, prime)) & 1:
                raise AssertionError("an exact square failed the pinned local sieve")
    return direction.identifier, unique, record


def select_and_certify(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[int, ...], IndependenceCertificate]:
    best: tuple[tuple[int, ...], IndependenceCertificate] | None = None
    failures = []
    for relation_prime in (5, 3, 7, 2):
        try:
            trial = select_independent_subset(
                coefficients,
                points,
                relation_prime=relation_prime,
                maximum_reduction_prime=CERTIFICATE_PRIME_BOUND,
            )
        except ArithmeticError as error:
            failures.append(f"ell={relation_prime}: {error}")
            continue
        if best is None or len(trial[0]) > len(best[0]):
            best = trial
        if len(trial[0]) == len(points):
            break
    if best is None:
        raise ArithmeticError("; ".join(failures))
    indices, certificate = best
    verify_independence_certificate(
        coefficients, tuple(points[index] for index in indices), certificate
    )
    return indices, certificate


def certify_third_parameter(
    signed_parameter: Fraction,
    directions: Sequence[Direction],
) -> dict[str, Any]:
    canonical_t = abs(signed_parameter)
    if canonical_t == 0 or FermigierMestreFamily.discriminant_factor(canonical_t) == 0:
        return {
            "status": "singular-or-zero-fiber",
            "signed_parameter_T": rational_to_string(signed_parameter),
            "canonical_parameter_T": rational_to_string(canonical_t),
        }
    specialization = specialize_fermigier_rank_sections(canonical_t / 2)
    generic_certificate = certify_fermigier_rank_sections(
        specialization, maximum_reduction_prime=CERTIFICATE_PRIME_BOUND
    )
    generic_x = {point[0] for point in specialization.quartic_points}
    forced_by_x: dict[Fraction, dict[str, Any]] = {}
    for direction in directions:
        value = evaluate_polynomial(direction.polynomial, signed_parameter)
        ordinate = rational_square_root(value)
        if ordinate is None:
            raise AssertionError("a simultaneous incidence lost an individual square")
        x_value = direction.x_value(signed_parameter)
        if evaluate_polynomial(specialization.quartic_model.quartic, x_value) != ordinate * ordinate:
            raise AssertionError("a forced point missed the canonical positive-T quartic")
        row = forced_by_x.setdefault(
            x_value,
            {
                "quartic_x": rational_to_string(x_value),
                "quartic_y": rational_to_string(ordinate),
                "direction_ids": [],
                "generic_abscissa_collision": x_value in generic_x,
                "exact_membership_checked": True,
            },
        )
        row["direction_ids"].append(direction.identifier)

    forced_quartic_points = tuple(
        point
        for x_value, row in forced_by_x.items()
        for point in (
            (x_value, Q(row["quartic_y"])),
            (x_value, -Q(row["quartic_y"])),
        )
    )
    searched, search_record = search_polynomial(
        specialization.quartic_model.quartic,
        height_bound=SPECIALIZATION_HEIGHT_BOUND,
        timeout=SPECIALIZATION_TIMEOUT_SECONDS,
        stack_bytes=PARI_STACK_BYTES,
    )
    pool = tuple(searched) + forced_quartic_points
    cloud = section_and_point_cloud_differences(specialization, pool)
    certificate_error = None
    try:
        selected_indices, certificate = select_and_certify(
            specialization.canonical_model, cloud
        )
        selected = tuple(cloud[index] for index in selected_indices)
        certificate_record: dict[str, Any] | None = certificate.to_json_object()
    except ArithmeticError as error:
        selected_indices, selected, certificate_record = (), (), None
        certificate_error = str(error)

    conductor = conductor_probe(
        canonical_t,
        timeout=30.0,
        stack_bytes=PARI_STACK_BYTES,
    )
    lower_bound = len(selected)
    return {
        "status": "finite-certified" if certificate_record is not None else "certificate-failed",
        "signed_parameter_T": rational_to_string(signed_parameter),
        "canonical_parameter_T": rational_to_string(canonical_t),
        "canonical_adapter_u": rational_to_string(canonical_t / 2),
        "forced_points": sorted(forced_by_x.values(), key=lambda row: Q(row["quartic_x"])),
        "distinct_forced_abscissa_count": len(forced_by_x),
        "generic_rank_12_certificate": generic_certificate.to_json_object(),
        "specialized_quartic_search": search_record,
        "specialized_quartic_signed_point_count": len(searched),
        "complete_discovered_group_point_count": len(cloud),
        "selected_indices_zero_based": list(selected_indices),
        "selected_points": [point_text(point) for point in selected],
        "independence_certificate": certificate_record,
        "certificate_error": certificate_error,
        "certified_rank_lower_bound": lower_bound if certificate_record is not None else None,
        "certified_gain_over_generic_12": (
            lower_bound - 12 if certificate_record is not None else None
        ),
        "conductor_probe": conductor,
        "rank_conductor_target_met": bool(
            certificate_record is not None
            and lower_bound >= 21
            and conductor.get("status") == "completed"
            and Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
        ),
    }


def stable_projection(value: Any) -> Any:
    """Remove timestamps and measured runtimes from a result-digest payload."""

    if isinstance(value, dict):
        return {
            key: stable_projection(child)
            for key, child in value.items()
            if key not in {"generated_at_utc", "wall_seconds", "search_wall_seconds", "pari_milliseconds"}
        }
    if isinstance(value, list):
        return [stable_projection(child) for child in value]
    return value


def stable_result_digest(artifact: dict[str, Any]) -> str:
    stable = {
        "schema_version": artifact["schema_version"],
        "source": artifact["source"],
        "search_box": artifact["search_box"],
        "local_sieve": artifact["local_sieve"],
        "prior_parameter_snapshot": artifact["prior_parameter_snapshot"],
        "directions": artifact["directions"],
        "individual_beyond_anchor_incidences": artifact[
            "individual_beyond_anchor_incidences"
        ],
        "pair_results": artifact["pair_results"],
        "third_parameter_certifications": artifact["third_parameter_certifications"],
        "outcome": artifact["outcome"],
    }
    return hashlib.sha256(
        json.dumps(stable_projection(stable), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def run(root: Path, output_path: Path, *, workers: int = WORKERS) -> dict[str, Any]:
    if workers != WORKERS:
        raise ValueError(f"this one-pass artifact pins workers={WORKERS}")
    transport, transport_raw = load_transport(root)
    directions = build_directions(transport)
    prior_parameters, prior_snapshot = prior_parameter_snapshot(root, output_path)

    # This ordering is intentional: all explicit local data are fixed before
    # the first bounded search subprocess is launched.
    masks, sieve_record = build_local_sieve(directions)
    if sieve_record["pair_count"] != EXPECTED_PAIR_COUNT:
        raise AssertionError("the local sieve did not cover all 3,160 pairs")

    started = time.monotonic()
    search_results: dict[str, tuple[tuple[Fraction, Fraction], ...]] = {}
    search_records: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(search_direction, direction, masks[direction.identifier]): direction
            for direction in directions
        }
        completed_count = 0
        for future in as_completed(futures):
            identifier, points, record = future.result()
            search_results[identifier] = points
            search_records[identifier] = record
            completed_count += 1
            if completed_count % 8 == 0:
                print(
                    f"completed {completed_count}/{EXPECTED_DIRECTION_COUNT} individual square searches",
                    flush=True,
                )
    search_wall_seconds = time.monotonic() - started

    by_identifier = {direction.identifier: direction for direction in directions}
    pair_rows = []
    simultaneous_by_parameter: dict[Fraction, set[str]] = defaultdict(set)
    for left, right in combinations(directions, 2):
        left_points = dict(search_results[left.identifier])
        right_points = dict(search_results[right.identifier])
        common = tuple(sorted(set(left_points) & set(right_points)))
        for parameter in common:
            if evaluate_polynomial(left.polynomial, parameter) != left_points[parameter] ** 2:
                raise AssertionError("left simultaneous factor ceased to be a square")
            if evaluate_polynomial(right.polynomial, parameter) != right_points[parameter] ** 2:
                raise AssertionError("right simultaneous factor ceased to be a square")
            simultaneous_by_parameter[parameter].update((left.identifier, right.identifier))
        beyond = tuple(
            parameter
            for parameter in common
            if abs(parameter) not in {abs(E22_T), abs(RANK20_T)}
        )
        pair_rows.append(
            {
                "first": left.identifier,
                "second": right.identifier,
                "simultaneous_square_parameters": [
                    rational_to_string(value) for value in common
                ],
                "beyond_anchor_parameters": [
                    rational_to_string(value) for value in beyond
                ],
                "both_factors_checked_separately": True,
            }
        )
    if len(pair_rows) != EXPECTED_PAIR_COUNT:
        raise AssertionError("the exact pair intersection did not cover all pairs")

    third_parameters = sorted(
        parameter
        for parameter in simultaneous_by_parameter
        if abs(parameter) not in {abs(E22_T), abs(RANK20_T)}
    )
    certifications = []
    for parameter in third_parameters:
        identifiers = sorted(simultaneous_by_parameter[parameter])
        record = certify_third_parameter(
            parameter, tuple(by_identifier[identifier] for identifier in identifiers)
        )
        record["direction_ids"] = identifiers
        record["direction_count"] = len(identifiers)
        record["prior_parameter"] = abs(parameter) in prior_parameters
        record["classification"] = (
            "prior-Fermigier-parameter"
            if record["prior_parameter"]
            else "new-third-parameter"
        )
        certifications.append(record)

    individual_beyond_anchor_incidences = []
    for direction in directions:
        for parameter, ordinate in search_results[direction.identifier]:
            if abs(parameter) in {abs(E22_T), abs(RANK20_T)}:
                continue
            individual_beyond_anchor_incidences.append(
                {
                    "direction_id": direction.identifier,
                    "signed_parameter_T": rational_to_string(parameter),
                    "canonical_parameter_T": rational_to_string(abs(parameter)),
                    "ordinate": rational_to_string(ordinate),
                    "prior_parameter": abs(parameter) in prior_parameters,
                    "individual_square_checked_exactly": True,
                    "rejection_reason": "no second affine cover is square at this parameter",
                }
            )
    individual_beyond_anchor_incidences.sort(
        key=lambda row: (Q(row["canonical_parameter_T"]), row["direction_id"])
    )

    direction_rows = []
    for direction in directions:
        record = search_records[direction.identifier]
        direction_rows.append(
            {
                "direction_id": direction.identifier,
                "E22_endpoint": {
                    "label": direction.e22_label,
                    "x": rational_to_string(direction.e22_x),
                },
                "rank20_endpoint": {
                    "label": direction.rank20_label,
                    "x": rational_to_string(direction.rank20_x),
                },
                "x_of_T": {
                    "slope": rational_to_string(direction.slope),
                    "intercept": rational_to_string(direction.intercept),
                },
                "polynomial_low_to_high": [
                    rational_to_string(value) for value in direction.polynomial
                ],
                "polynomial_sha256": polynomial_digest(direction.polynomial),
                "search": record,
                "points": [
                    {
                        "T": rational_to_string(parameter),
                        "y": rational_to_string(ordinate),
                        "individual_square_checked_exactly": True,
                    }
                    for parameter, ordinate in search_results[direction.identifier]
                ],
            }
        )

    all_completed = all(
        row["search"]["status"] == "completed" for row in direction_rows
    )
    new_certifications = [
        row for row in certifications if row["classification"] == "new-third-parameter"
    ]
    target_rows = [row for row in certifications if row["rank_conductor_target_met"]]
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-exceptional-pair-simultaneous-h200000.v1",
        "status": "complete bounded search" if all_completed else "incomplete bounded search",
        "claim_level": (
            "exact simultaneous-square computation inside one finite projective-height box; "
            "absence outside the box is not proved"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "transport_artifact": str(TRANSPORT_RELATIVE),
            "transport_artifact_sha256": sha256_bytes(transport_raw),
            "transport_result_sha256": transport["result_sha256"],
            "transport_fiber_product_manifest_sha256": transport["fiber_products"]["manifest_sha256"],
            "script_sha256": sha256_file(Path(__file__)),
        },
        "search_box": {
            "coordinate": "legacy literal shift T=2*u",
            "definition": "T=a/b in lowest terms, b>0, max(abs(a),b)<=200000",
            "projective_height_bound": HEIGHT_BOUND,
            "engine": "PARI/GP hyperellratpoints",
            "one_pass": True,
            "retries": 0,
            "direction_timeout_seconds": SEARCH_TIMEOUT_SECONDS,
            "workers": workers,
            "search_wall_seconds": search_wall_seconds,
        },
        "local_sieve": sieve_record,
        "prior_parameter_snapshot": prior_snapshot,
        "directions": direction_rows,
        "individual_beyond_anchor_incidences": individual_beyond_anchor_incidences,
        "pair_results": pair_rows,
        "third_parameter_certifications": certifications,
        "outcome": {
            "all_80_direction_searches_completed": all_completed,
            "direction_count": len(direction_rows),
            "fiber_product_pair_count": len(pair_rows),
            "anchor_calibration_parameters": [
                rational_to_string(E22_T),
                rational_to_string(RANK20_T),
            ],
            "signed_simultaneous_parameter_count_including_anchors": len(simultaneous_by_parameter),
            "individual_beyond_anchor_incidence_count": len(
                individual_beyond_anchor_incidences
            ),
            "individual_beyond_anchor_unique_parameter_count": len(
                {
                    row["canonical_parameter_T"]
                    for row in individual_beyond_anchor_incidences
                }
            ),
            "individual_beyond_anchor_prior_parameter_count": sum(
                row["prior_parameter"]
                for row in individual_beyond_anchor_incidences
            ),
            "signed_third_parameter_count": len(third_parameters),
            "prior_third_parameter_count": sum(row["prior_parameter"] for row in certifications),
            "new_third_parameter_count": len(new_certifications),
            "new_third_parameters": [row["canonical_parameter_T"] for row in new_certifications],
            "highest_certified_rank_lower_bound_among_third_parameters": max(
                (row["certified_rank_lower_bound"] or 0 for row in certifications),
                default=0,
            ),
            "rank_conductor_target_hits": len(target_rows),
            "target_met": bool(target_rows),
            "product_square_surrogate_used": False,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/search_fermigier_exceptional_pair_simultaneous_h200000.py"
        ),
    }
    artifact["result_sha256"] = stable_result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=WORKERS)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts/generated-results/elliptic_fermigier_exceptional_pair_simultaneous_h200000.json",
    )
    args = parser.parse_args()
    artifact = run(root, args.output, workers=args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()

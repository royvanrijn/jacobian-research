#!/usr/bin/env python3
"""Bounded pairwise fiber-product screen for Fermigier's published directions.

At Fermigier's rank-22 record fiber, published points P6 and P13--P22 have
exact non-generic quartic preimages.  For each preimage and each slope
``m in {-1,+1}``, the line ``x=m*T+n`` gives a genus-one quartic slice
``z^2=f_i(T)``.  Simultaneous persistence of two distinct published source
directions requires a rational point on

    w^2 = f_i(T) f_j(T).

This script searches all 220 cross-direction products once at H=5000.  A
pair advances once to H=50000 only if its pilot produces a genuinely new
parameter at which *both* factors are individually rational squares.  Every
survivor is checked on both specialized quartics, decontaminated against the
record fiber, both signs of the generic section catalog, and all stored prior
Fermigier parameter populations.  Exact conductor and H=50000 specialized
point/rank triage are reserved for new fibers forcing at least two distinct
source labels and lying below the strict conductor target.

The product search and conductor/rank calls are bounded computations.  A
numerical height rank is not a rank certificate; an exact finite-reduction
attempt is triggered only at stable numerical rank at least 21.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from alternate_quartic_covers import point_on_short_curve
from ek_k3 import rational_square_root, rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_rank22_accidental_slices import (
    T0,
    Slice,
    build_slices,
    canonical_signless_points,
    conductor_probe,
    finite_reduction_attempt,
    generic_group_seed_points,
    point_record,
    poly_evaluate,
    poly_multiply,
    quartic_group_pullback,
    search_polynomial,
    search_specialized_quartic,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    stable_height_rank,
)


TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXPECTED_LABELS = ("P6", *tuple(f"P{index}" for index in range(13, 23)))
EXPECTED_SLICE_COUNT = 22
EXPECTED_PAIR_COUNT = 220
PILOT_HEIGHT = 5_000
ESCALATION_HEIGHT = 50_000
SPECIALIZATION_HEIGHT = 50_000
PRIMARY_ARTIFACT = "elliptic_fermigier_rank22_accidental_slices.json"
EXPECTED_PUBLISHED_PREIMAGE_SHA256 = (
    "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15"
)
EXPECTED_PRIOR_PARAMETER_COUNT = 590
EXPECTED_PRIOR_PARAMETER_SHA256 = (
    "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291"
)
EXPECTED_H5000_RESULT_SHA256 = (
    "80413701447b6468a826fa2185528da74057faa02619bdf02f237e7efb8b1b8b"
)
LEGACY_ACCIDENTAL_PARAMETERS = tuple(
    Fraction(value)
    for value in ("19033/135", "22253/114", "31331/104", "38633/138")
)
PARAMETER_KEYS = {
    "T",
    "t",
    "candidate_t",
    "normalized_parameter_t",
    "normalized_record_parameter",
    "parameter",
    "parameter_t",
    "published_parameter",
    "record_parameter_normalized_T",
}
PARAMETER_LIST_KEYS = {
    "candidate_parameters",
    "canonical_candidate_parameters",
    "parameters_t",
}


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def rational_digest(values: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update((rational_to_string(value) + "\n").encode())
    return digest.hexdigest()


def polynomial_digest(coefficients: Sequence[Fraction]) -> str:
    return rational_digest(Fraction(value) for value in coefficients)


def load_primary(path: Path) -> tuple[dict[str, Any], bytes]:
    raw = path.read_bytes()
    data = json.loads(raw)
    if published_preimage_digest(data) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact 11-point published-preimage population changed")
    return data, raw


def published_preimage_digest(data: dict[str, Any]) -> str:
    digest = hashlib.sha256()
    for record in data["published_point_preimages"]:
        if record["classification"] != "accidental rational preimage":
            continue
        point = record["quartic_preimage"]
        digest.update(
            (
                f"{record['label']}|{Fraction(point['x'])}|{Fraction(point['z'])}\n"
            ).encode()
        )
    return digest.hexdigest()


def pair_result_digest(pair_rows: Sequence[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in pair_rows:
        pilot = row["pilot"]
        digest.update(
            (
                f"{row['pair_id']}|{pilot['product_polynomial_sha256']}|"
                f"{pilot['search'].get('signed_point_count')}|"
                f"{pilot['qualifying_new_parameter_count']}\n"
            ).encode()
        )
    return digest.hexdigest()


def published_accidentals(data: dict[str, Any]) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    rows = []
    for record in data["published_point_preimages"]:
        if record["classification"] != "accidental rational preimage":
            continue
        point = record["quartic_preimage"]
        if not point["exact_quartic_membership_checked"]:
            raise AssertionError("a published preimage lost exact membership")
        rows.append(
            (
                record["label"],
                (Fraction(point["x"]), Fraction(point["z"])),
            )
        )
    if tuple(label for label, _ in rows) != EXPECTED_LABELS:
        raise AssertionError("the 11 published accidental labels changed")
    return tuple(rows)


def exact_slices(accidentals: Sequence[tuple[str, tuple[Fraction, Fraction]]]) -> tuple[Slice, ...]:
    slices = build_slices(accidentals)
    if len(slices) != EXPECTED_SLICE_COUNT:
        raise AssertionError("expected 22 signed published-point slices")
    if any(poly_evaluate(slice_data.coefficients, T0) != slice_data.source_point[1] ** 2 for slice_data in slices):
        raise AssertionError("a published slice missed its record source point")
    return slices


def pair_population(slices: Sequence[Slice]) -> tuple[tuple[Slice, Slice], ...]:
    pairs = tuple(
        (left, right)
        for left, right in combinations(slices, 2)
        if left.accidental_label != right.accidental_label
    )
    if len(pairs) != EXPECTED_PAIR_COUNT:
        raise AssertionError("the cross-direction pair population changed")
    return pairs


def pair_identifier(left: Slice, right: Slice) -> str:
    return f"{left.identifier}__{right.identifier}"


def polynomial_coprime(left: Sequence[Fraction], right: Sequence[Fraction]) -> bool:
    symbol = sp.symbols("T")
    left_poly = sp.Poly(
        sum(sp.Rational(value.numerator, value.denominator) * symbol**index for index, value in enumerate(left)),
        symbol,
        domain=sp.QQ,
    )
    right_poly = sp.Poly(
        sum(sp.Rational(value.numerator, value.denominator) * symbol**index for index, value in enumerate(right)),
        symbol,
        domain=sp.QQ,
    )
    return sp.gcd(left_poly, right_poly).degree() == 0


def extract_parameter_values(value: Any, *, key: str | None = None) -> set[Fraction]:
    answer: set[Fraction] = set()
    if key in PARAMETER_KEYS and isinstance(value, (str, int)):
        try:
            answer.add(abs(Fraction(value)))
        except (ValueError, ZeroDivisionError):
            pass
    if key in PARAMETER_LIST_KEYS and isinstance(value, list):
        for item in value:
            if isinstance(item, (str, int)):
                try:
                    answer.add(abs(Fraction(item)))
                except (ValueError, ZeroDivisionError):
                    pass
    if isinstance(value, dict):
        for child_key, child in value.items():
            answer.update(extract_parameter_values(child, key=child_key))
    elif isinstance(value, list):
        for child in value:
            answer.update(extract_parameter_values(child))
    return answer


def prior_fermigier_parameters(
    artifact_dir: Path, output_path: Path
) -> tuple[set[Fraction], dict[str, Any]]:
    parameters = {abs(T0), *[abs(value) for value in LEGACY_ACCIDENTAL_PARAMETERS]}
    sources: dict[str, dict[str, Any]] = {}
    output_resolved = output_path.resolve()
    for path in sorted(artifact_dir.glob("elliptic_fermigier*.json")):
        if path.resolve() == output_resolved:
            continue
        raw = path.read_bytes()
        data = json.loads(raw)
        extracted = extract_parameter_values(data)
        parameters.update(extracted)
        sources[path.name] = {
            "sha256": sha256_bytes(raw),
            "extracted_parameter_count": len(extracted),
        }
    ordered = tuple(sorted(parameters))
    return parameters, {
        "canonicalization": "T -> abs(T), since the Fermigier family is even",
        "artifact_sources": sources,
        "legacy_accidental_parameters": [
            rational_to_string(abs(value)) for value in LEGACY_ACCIDENTAL_PARAMETERS
        ],
        "unique_prior_parameter_count": len(parameters),
        "prior_parameter_sha256": rational_digest(ordered),
    }


def generic_abscissas(parameter: Fraction) -> set[Fraction]:
    parameter = abs(Fraction(parameter))
    return {
        point[0]
        for signed_parameter in (parameter, -parameter)
        for point in FermigierMestreFamily.known_quartic_points(signed_parameter)
    }


def forced_point_record(
    slice_data: Slice,
    signed_parameter: Fraction,
    ordinate: Fraction,
    canonical_parameter: Fraction,
) -> tuple[dict[str, Any], tuple[Fraction, Fraction]]:
    x_value = slice_data.x_value(signed_parameter)
    point = (x_value, abs(ordinate))
    if point[1] ** 2 != FermigierMestreFamily.quartic_value(
        canonical_parameter, point[0]
    ):
        raise AssertionError("a factor square did not give an exact quartic point")
    direct = FermigierMestreFamily.quartic_point_to_jacobian(canonical_parameter, point)
    pullback = quartic_group_pullback(canonical_parameter, point)
    if pullback is None:
        raise AssertionError("a non-generic forced point pulled back to the origin")
    return (
        {
            "source_label": slice_data.accidental_label,
            "slice_id": slice_data.identifier,
            "quartic_x": rational_to_string(point[0]),
            "quartic_z": rational_to_string(point[1]),
            "direct_covariant_image": point_record(direct),
            "basepoint_group_pullback": point_record(pullback),
            "exact_membership_checked": True,
        },
        pullback,
    )


def exact_forced_quartic_record(
    slice_data: Slice,
    signed_parameter: Fraction,
    ordinate: Fraction,
    canonical_parameter: Fraction,
) -> dict[str, Any]:
    x_value = slice_data.x_value(signed_parameter)
    point = (x_value, abs(ordinate))
    if point[1] ** 2 != FermigierMestreFamily.quartic_value(
        canonical_parameter, point[0]
    ):
        raise AssertionError("an individual factor square missed its quartic fiber")
    return {
        "source_label": slice_data.accidental_label,
        "slice_id": slice_data.identifier,
        "quartic_x": rational_to_string(point[0]),
        "quartic_z": rational_to_string(point[1]),
        "exact_membership_checked": True,
    }


def classify_product_point(
    left: Slice,
    right: Slice,
    product_point: tuple[Fraction, Fraction],
    prior_parameters: set[Fraction],
) -> dict[str, Any]:
    signed_parameter, product_ordinate = map(Fraction, product_point)
    left_value = poly_evaluate(left.coefficients, signed_parameter)
    right_value = poly_evaluate(right.coefficients, signed_parameter)
    if product_ordinate**2 != left_value * right_value:
        raise AssertionError("PARI returned a point off the exact fiber product")
    left_root = rational_square_root(left_value)
    right_root = rational_square_root(right_value)
    record: dict[str, Any] = {
        "signed_parameter_t": rational_to_string(signed_parameter),
        "canonical_parameter_t": rational_to_string(abs(signed_parameter)),
        "product_ordinate": rational_to_string(abs(product_ordinate)),
        "left_factor_is_square": left_root is not None,
        "right_factor_is_square": right_root is not None,
        "classification": "product-square-only",
    }
    if left_root is None or right_root is None:
        return record
    canonical_parameter = abs(signed_parameter)
    exact_left = exact_forced_quartic_record(
        left, signed_parameter, left_root, canonical_parameter
    )
    exact_right = exact_forced_quartic_record(
        right, signed_parameter, right_root, canonical_parameter
    )
    record["exact_forced_quartic_points"] = [exact_left, exact_right]
    if canonical_parameter == 0:
        record["classification"] = "zero-parameter-excluded"
        return record
    if canonical_parameter == abs(T0):
        record["classification"] = "record-fiber-excluded"
        return record
    if FermigierMestreFamily.discriminant_factor(canonical_parameter) == 0:
        record["classification"] = "singular-fiber-excluded"
        return record
    left_x = Fraction(exact_left["quartic_x"])
    right_x = Fraction(exact_right["quartic_x"])
    generic_x = generic_abscissas(canonical_parameter)
    if left_x in generic_x or right_x in generic_x:
        record["classification"] = "generic-section-collision-excluded"
        return record
    if canonical_parameter in prior_parameters:
        record["classification"] = "prior-Fermigier-population-excluded"
        return record
    if left_x == right_x:
        record["classification"] = "coincident-forced-abscissa-excluded"
        return record
    left_record, left_pullback = forced_point_record(
        left, signed_parameter, left_root, canonical_parameter
    )
    right_record, right_pullback = forced_point_record(
        right, signed_parameter, right_root, canonical_parameter
    )
    if left_pullback[0] == right_pullback[0]:
        record["classification"] = "coincident-specialized-direction-excluded"
        return record
    record.update(
        {
            "classification": "genuinely-new-double-forced-fiber",
            "distinct_published_source_direction_count": 2,
            "forced_points": [left_record, right_record],
        }
    )
    return record


def unique_product_parameters(
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    by_parameter: dict[Fraction, Fraction] = {}
    for parameter, ordinate in points:
        if parameter in by_parameter and by_parameter[parameter] ** 2 != ordinate**2:
            raise AssertionError("one product parameter had inconsistent ordinates")
        by_parameter[parameter] = abs(ordinate)
    return tuple(sorted(by_parameter.items()))


def search_pair(
    left: Slice,
    right: Slice,
    *,
    height: int,
    timeout: float,
    stack_bytes: int,
    prior_parameters: set[Fraction],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    product = poly_multiply(left.coefficients, right.coefficients)
    if len(product) != 9:
        raise AssertionError("a pair product did not have degree eight")
    points, search = search_polynomial(
        product, height_bound=height, timeout=timeout, stack_bytes=stack_bytes
    )
    incidences = [
        classify_product_point(left, right, point, prior_parameters)
        for point in unique_product_parameters(points)
    ]
    qualifying = [
        record
        for record in incidences
        if record["classification"] == "genuinely-new-double-forced-fiber"
    ]
    return (
        {
            "height_bound": height,
            "search": search,
            "product_polynomial_sha256": polynomial_digest(product),
            "product_degree": 8,
            "incidences": incidences,
            "qualifying_new_parameter_count": len(
                {record["canonical_parameter_t"] for record in qualifying}
            ),
        },
        qualifying,
    )


def aggregate_candidates(
    qualifying: Sequence[tuple[str, dict[str, Any]]],
) -> dict[Fraction, dict[str, Any]]:
    candidates: dict[Fraction, dict[str, Any]] = {}
    for pair_id, incidence in qualifying:
        parameter = Fraction(incidence["canonical_parameter_t"])
        candidate = candidates.setdefault(
            parameter,
            {
                "parameter_t": rational_to_string(parameter),
                "pair_ids": set(),
                "signed_parameters": set(),
                "source_labels": set(),
                "forced_points_by_x": {},
            },
        )
        candidate["pair_ids"].add(pair_id)
        candidate["signed_parameters"].add(incidence["signed_parameter_t"])
        for point in incidence["forced_points"]:
            candidate["source_labels"].add(point["source_label"])
            existing = candidate["forced_points_by_x"].get(point["quartic_x"])
            if existing is None:
                candidate["forced_points_by_x"][point["quartic_x"]] = point
            elif existing["quartic_z"] != point["quartic_z"]:
                raise AssertionError("a forced abscissa had inconsistent ordinates")
    return candidates


def finalized_candidate_record(candidate: dict[str, Any]) -> dict[str, Any]:
    points = list(candidate["forced_points_by_x"].values())
    pullback_x = {
        point["basepoint_group_pullback"]["jacobian_x"] for point in points
    }
    return {
        "parameter_t": candidate["parameter_t"],
        "pair_ids": sorted(candidate["pair_ids"]),
        "signed_parameters": sorted(candidate["signed_parameters"]),
        "published_source_labels": sorted(candidate["source_labels"]),
        "distinct_published_source_direction_count": len(candidate["source_labels"]),
        "distinct_forced_quartic_abscissa_count": len(points),
        "distinct_group_pullback_classes_modulo_inverse": len(pullback_x),
        "forced_points": sorted(points, key=lambda point: (point["source_label"], Fraction(point["quartic_x"]))),
    }


def triage_specialization(
    candidate: dict[str, Any],
    *,
    search_timeout: float,
    height_timeout: float,
    precisions: tuple[int, ...],
    stack_bytes: int,
    saturation_timeout: float,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    parameter = Fraction(candidate["parameter_t"])
    coefficients = FermigierMestreFamily.coefficients(parameter)
    seeds = generic_group_seed_points(parameter)
    seed_runs = height_matrix_replay(
        coefficients,
        seeds,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    forced_points = tuple(
        (
            Fraction(record["basepoint_group_pullback"]["jacobian_x"]),
            Fraction(record["basepoint_group_pullback"]["jacobian_y"]),
        )
        for record in candidate["forced_points"]
    )
    pool = list(seeds)
    seen_x = {point[0] for point in pool}
    for point in forced_points:
        if point[0] not in seen_x:
            pool.append(point)
            seen_x.add(point[0])
    forced_pool = tuple(pool)
    forced_runs = height_matrix_replay(
        coefficients,
        forced_pool,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )

    raw_points, search = search_specialized_quartic(
        parameter,
        height_bound=SPECIALIZATION_HEIGHT,
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    new_search_points = []
    if search["status"] == "completed":
        generic_x = generic_abscissas(parameter)
        forced_x = {Fraction(record["quartic_x"]) for record in candidate["forced_points"]}
        for quartic_point in canonical_signless_points(raw_points):
            if quartic_point[0] in generic_x or quartic_point[0] in forced_x:
                continue
            pullback = quartic_group_pullback(parameter, quartic_point)
            if pullback is None or pullback[0] in seen_x:
                continue
            seen_x.add(pullback[0])
            pool.append(pullback)
            new_search_points.append(
                {
                    "quartic_x": rational_to_string(quartic_point[0]),
                    "quartic_z": rational_to_string(quartic_point[1]),
                    "basepoint_group_pullback": point_record(pullback),
                }
            )
    full_pool = tuple(pool)
    full_runs = height_matrix_replay(
        coefficients,
        full_pool,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    seed_rank = stable_height_rank(seed_runs)
    forced_rank = stable_height_rank(forced_runs)
    full_rank = stable_height_rank(full_runs)
    indices = tuple(full_runs[-1]["subset_indices_one_based"])
    selected = tuple(full_pool[index - 1] for index in indices)
    record: dict[str, Any] = {
        "specialized_quartic_search": search,
        "signed_points_found": len(raw_points),
        "new_search_group_pullbacks": new_search_points,
        "generic_seed_point_count": len(seeds),
        "forced_augmented_point_count": len(forced_pool),
        "full_pool_point_count": len(full_pool),
        "generic_seed_stable_numerical_rank": seed_rank,
        "forced_augmented_stable_numerical_rank": forced_rank,
        "full_pool_stable_numerical_rank": full_rank,
        "numerical_rank_gain_from_forced_directions": forced_rank - seed_rank,
        "numerical_rank_gain_after_H50000": full_rank - seed_rank,
        "seed_height_runs": list(seed_runs),
        "forced_height_runs": list(forced_runs),
        "full_height_runs": list(full_runs),
        "selected_point_sha256": point_digest(selected),
        "selected_points": [point_record(point) for point in selected],
        "scope_warning": "stable height ranks are numerical triage evidence only",
    }
    if full_rank >= 21:
        record["finite_reduction_attempt"] = finite_reduction_attempt(
            coefficients,
            selected,
            saturation_timeout=saturation_timeout,
            stack_bytes=stack_bytes,
            certificate_prime_bound=certificate_prime_bound,
        )
    return record


def parse_precisions(value: str) -> tuple[int, ...]:
    values = tuple(int(part) for part in value.split(",") if part)
    if len(values) < 2 or tuple(sorted(set(values))) != values:
        raise argparse.ArgumentTypeError("provide increasing distinct precisions")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-height", type=int, default=PILOT_HEIGHT)
    parser.add_argument("--escalation-height", type=int, default=ESCALATION_HEIGHT)
    parser.add_argument("--pair-timeout", type=float, default=5.0)
    parser.add_argument("--escalation-timeout", type=float, default=15.0)
    parser.add_argument("--pilot-wall-cap", type=float, default=180.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_published_pair_fiber_products.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.pilot_height != PILOT_HEIGHT or args.escalation_height != ESCALATION_HEIGHT:
        raise SystemExit("the declared pair search pins pilot H=5000 and escalation H=50000")
    if min(
        args.pair_timeout,
        args.escalation_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if max(args.pair_timeout, args.escalation_timeout) > 60:
        raise SystemExit("pair subprocess timeouts may not exceed 60 seconds")

    root = Path(__file__).resolve().parents[2]
    artifact_dir = root / "artifacts" / "generated-results"
    primary_path = artifact_dir / PRIMARY_ARTIFACT
    primary, primary_raw = load_primary(primary_path)
    primary_script = root / "elliptic-curves" / "cas" / "search_fermigier_rank22_accidental_slices.py"
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
    all_qualifying: list[tuple[str, dict[str, Any]]] = []
    productive_pairs = []
    pilot_started = time.monotonic()
    stopped_disproportionate = False
    for left, right in pairs:
        if time.monotonic() - pilot_started > args.pilot_wall_cap:
            stopped_disproportionate = True
            break
        if not polynomial_coprime(left.coefficients, right.coefficients):
            raise AssertionError("two distinct source slices shared a polynomial factor")
        pair_id = pair_identifier(left, right)
        pilot, qualifying = search_pair(
            left,
            right,
            height=args.pilot_height,
            timeout=args.pair_timeout,
            stack_bytes=args.stack_bytes,
            prior_parameters=prior_parameters,
        )
        row = {
            "pair_id": pair_id,
            "left_source_label": left.accidental_label,
            "left_slice_id": left.identifier,
            "right_source_label": right.accidental_label,
            "right_slice_id": right.identifier,
            "distinct_source_labels": True,
            "factor_polynomials_coprime": True,
            "pilot": pilot,
            "escalation": None,
        }
        pair_rows.append(row)
        all_qualifying.extend((pair_id, incidence) for incidence in qualifying)
        if qualifying:
            productive_pairs.append((left, right, row))

    if not stopped_disproportionate and len(pair_rows) != EXPECTED_PAIR_COUNT:
        raise AssertionError("the pair pilot did not cover all 220 declared pairs")
    pilot_wall_seconds = time.monotonic() - pilot_started
    for left, right, row in productive_pairs:
        escalation, qualifying = search_pair(
            left,
            right,
            height=args.escalation_height,
            timeout=args.escalation_timeout,
            stack_bytes=args.stack_bytes,
            prior_parameters=prior_parameters,
        )
        row["escalation"] = escalation
        all_qualifying.extend((row["pair_id"], incidence) for incidence in qualifying)

    aggregated = aggregate_candidates(all_qualifying)
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
    exact_pair_result_sha256 = pair_result_digest(pair_rows)
    if exact_pair_result_sha256 != EXPECTED_H5000_RESULT_SHA256:
        raise AssertionError("the exact H=5000 pair-result population changed")
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded pairwise necessary-condition search; factor-square and quartic "
            "membership checks are exact, conductors are PARI computations, and "
            "height ranks remain numerical unless a finite-reduction block says certified"
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
            "artifact": str(primary_path.relative_to(root)),
            "artifact_sha256_observed": sha256_bytes(primary_raw),
            "script": str(primary_script.relative_to(root)),
            "script_sha256_observed": sha256_file(primary_script),
            "published_accidental_preimage_sha256": published_preimage_digest(
                primary
            ),
            "record_parameter_t": rational_to_string(T0),
            "published_accidental_labels": list(EXPECTED_LABELS),
        },
        "slice_population": {
            "source_direction_count": len(accidentals),
            "signed_slice_count": len(slices),
            "slices": [
                {
                    "slice_id": slice_data.identifier,
                    "source_label": slice_data.accidental_label,
                    "slope": slice_data.slope,
                    "intercept": rational_to_string(slice_data.intercept),
                    "source_quartic_x": rational_to_string(slice_data.source_point[0]),
                    "source_quartic_z": rational_to_string(slice_data.source_point[1]),
                    "quartic_coefficients_ascending": [
                        rational_to_string(value) for value in slice_data.coefficients
                    ],
                }
                for slice_data in slices
            ],
        },
        "prior_decontamination": prior_record,
        "parameters": {
            "pilot_height": args.pilot_height,
            "escalation_height": args.escalation_height,
            "pair_timeout_seconds": args.pair_timeout,
            "escalation_timeout_seconds": args.escalation_timeout,
            "pilot_global_wall_cap_seconds": args.pilot_wall_cap,
            "conductor_timeout_seconds": args.conductor_timeout,
            "specialization_height": SPECIALIZATION_HEIGHT,
            "specialization_timeout_seconds": args.specialization_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": list(args.precisions),
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
        },
        "pair_searches": pair_rows,
        "candidates": candidates,
        "outcome": {
            "declared_pair_count": EXPECTED_PAIR_COUNT,
            "pilot_pairs_attempted": len(pair_rows),
            "pilot_pairs_completed": sum(
                row["pilot"]["search"]["status"] == "completed" for row in pair_rows
            ),
            "pilot_pairs_timed_out_or_errored": sum(
                row["pilot"]["search"]["status"] != "completed" for row in pair_rows
            ),
            "stopped_as_computationally_disproportionate": stopped_disproportionate,
            "pilot_wall_seconds": pilot_wall_seconds,
            "productive_pilot_pairs": len(productive_pairs),
            "escalation_pairs_attempted": len(productive_pairs),
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
            "exact_pair_result_sha256": exact_pair_result_sha256,
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
            "no new double-forced subtarget fiber received an exact rank-21 certificate"
        )
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()

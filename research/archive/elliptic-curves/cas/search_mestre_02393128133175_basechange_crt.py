#!/usr/bin/env python3
"""Exact multiple-root CRT search on a rank-at-least-13 Mestre base change.

The base family has centers ``(0,23,93,128,133,175)`` and the conic base
change ``T=(14406-u^2)/(2u)``.  The base change supplies thirteen exactly
independent generic Jacobian directions (eleven visible directions, one
linear-extra direction, and split infinity).

This search does not use an additive rank score.  It factors the arithmetic
search into exact local constraints on the primitive degree-40 homogeneous
base discriminant, CRT, and exact two-dimensional Gauss reduction.  The
chosen residue balls force clean multiplicative discriminant powers at
11,19,23,37,47.  Candidate ordering uses only lattice height and a radical
upper proxy.  PARI computes exact conductors/root numbers before any bounded
point search.  Every promoted point cloud is mapped exactly and receives a
finite-reduction mod-3 subgroup-rank certificate.

All searches are bounded; failure to find more points is not a rank upper
bound.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd, isqrt, lcm, log
from pathlib import Path
import platform
import sys
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, short_rational_representatives
from mestre_rank13_02393128133175 import (
    BASE_CHANGE_CONSTANT,
    CONSTRUCTION,
    ROOTS,
    base_changed_short_jacobian_coefficients,
    base_parameter,
    known_jacobian_points,
)
from multiple_root_lifting import (
    affine_variable_coefficients,
    all_roots_mod_prime_power,
    fixed_divisor_valuation,
)
from pari_bridge import minimal_curve_data
from search_mestre_rank13_multifamily_rational import derive_even_coefficients
from search_mestre_root_tuple_scale import (
    bounded_quartic_points,
    canonical_signless_points,
    point_digest,
    quartic_point_to_jacobian,
    sha256_file,
)
from search_mestre_root_tuple_scale_max200 import (
    gf_l_rank_and_pivots,
    mod3_independence_certificate,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
TARGET_LOG_CONDUCTOR = Decimal("182.72")
STACK_BYTES = 512_000_000
DEFAULT_OUTPUT_DIRECTORY = Path(
    "artifacts/local/elliptic-curves/mestre-02393128133175-basechange-crt-v1"
)


@dataclass(frozen=True)
class ConstraintGroup:
    name: str
    prime: int
    modulus: int
    residues: tuple[int, ...]
    forced_discriminant_valuation: int


GROUPS = {
    "p11v3": ConstraintGroup("p11v3", 11, 11, (3, 5, 6, 8), 3),
    "p19v4": ConstraintGroup("p19v4", 19, 19, (2, 3, 5, 14, 16, 17), 4),
    "p23v3": ConstraintGroup("p23v3", 23, 23, (10, 13), 3),
    "p37v3": ConstraintGroup("p37v3", 37, 37, (16, 17, 20, 21), 3),
    "p47v2": ConstraintGroup("p47v2", 47, 47, (11, 15, 20, 27, 32, 36), 2),
}

TIERS = (
    ("p19v4", "p23v3", "p37v3"),
    ("p11v3", "p19v4", "p23v3"),
    ("p19v4", "p23v3", "p47v2"),
    ("p11v3", "p19v4", "p23v3", "p37v3"),
    ("p19v4", "p23v3", "p37v3", "p47v2"),
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_text(point[0]), "y": rational_text(point[1])}


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def poly_multiply(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return tuple(answer)


def primitive_base_discriminant() -> tuple[tuple[int, ...], int]:
    """Return primitive F(u)=(2u)^20 D((14406-u^2)/(2u))/content."""

    direct = CONSTRUCTION.primitive_discriminant_polynomial
    if any(value.denominator != 1 for value in direct):
        raise AssertionError("the direct primitive discriminant ceased to be integral")
    constant = int(BASE_CHANGE_CONSTANT)
    numerator = (constant, 0, -1)
    denominator = (0, 2)
    numerator_powers = [(1,)]
    denominator_powers = [(1,)]
    for _ in range(20):
        numerator_powers.append(poly_multiply(numerator_powers[-1], numerator))
        denominator_powers.append(poly_multiply(denominator_powers[-1], denominator))
    answer = [0] * 41
    for power, coefficient in enumerate(direct):
        term = poly_multiply(
            numerator_powers[power], denominator_powers[20 - power]
        )
        for index, value in enumerate(term):
            answer[index] += int(coefficient) * value
    content = gcd(*(abs(value) for value in answer))
    primitive = tuple(value // content for value in answer)
    if content != 47_775_744 or gcd(*(abs(value) for value in primitive)) != 1:
        raise AssertionError("the base discriminant content changed")
    if len(primitive) != 41 or any(primitive[index] for index in range(1, 41, 2)):
        raise AssertionError("the primitive base discriminant ceased to be even degree 40")
    return primitive, content


BASE_DISCRIMINANT, BASE_DISCRIMINANT_CONTENT = primitive_base_discriminant()


def polynomial_value(coefficients: Sequence[int], value: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def homogeneous_value(coefficients: Sequence[int], numerator: int, denominator: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def valuation_integer(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    value = abs(value)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def primes_up_to(bound: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(bound) + 1):
        if sieve[prime]:
            sieve[prime * prime : bound + 1 : prime] = b"\x00" * (
                (bound - prime * prime) // prime + 1
            )
    return tuple(index for index, value in enumerate(sieve) if value)


TRIAL_PRIMES = primes_up_to(997)


def radical_upper(value: int) -> int:
    """Return an upper bound for rad(value), exact through trial prime 997."""

    value = abs(int(value))
    if value <= 1:
        return 1
    answer = 1
    for prime in TRIAL_PRIMES:
        if value % prime:
            continue
        answer *= prime
        while value % prime == 0:
            value //= prime
    if value > 1:
        answer *= value
    return answer


def log_integer(value: int) -> float:
    value = abs(int(value))
    if value <= 1:
        return 0.0
    bits = value.bit_length()
    if bits <= 1022:
        return log(value)
    shift = bits - 53
    return log(value >> shift) + shift * log(2.0)


def group_certificate(group: ConstraintGroup) -> dict[str, Any]:
    lift = all_roots_mod_prime_power(
        BASE_DISCRIMINANT,
        group.prime,
        group.forced_discriminant_valuation,
        max_roots=200_000,
    )
    maximal = lift.maximal_balls()
    residue_records = []
    for residue in group.residues:
        fixed = fixed_divisor_valuation(
            affine_variable_coefficients(
                BASE_DISCRIMINANT, residue, group.modulus
            ),
            group.prime,
        )
        if fixed < group.forced_discriminant_valuation:
            raise AssertionError(f"{group.name} lost its fixed valuation")
        if not any(
            ball.exponent <= 1
            and residue % ball.modulus == ball.residue
            for ball in maximal
        ):
            raise AssertionError(f"{group.name} residue is not a compressed root ball")
        residue_records.append(
            {"residue": residue, "fixed_discriminant_valuation": fixed}
        )
    return {
        **asdict(group),
        "complete_lift_level_counts": list(lift.level_counts),
        "complete_root_count_at_target_power": len(lift.roots),
        "complete_maximal_balls": [asdict(ball) for ball in maximal],
        "selected_efficient_unit_residues": residue_records,
        "selected_union_density": rational_text(Q(len(group.residues), group.modulus)),
    }


def enumerate_tier(
    group_names: Sequence[str], *, coefficient_radius: int, representatives_per_class: int
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    groups = tuple(GROUPS[name] for name in group_names)
    records = []
    class_count = 0
    for residues in product(*(group.residues for group in groups)):
        residue = 0
        modulus = 1
        for group, local_residue in zip(groups, residues, strict=True):
            residue, modulus = crt_pair(
                residue, modulus, local_residue, group.modulus
            )
        class_count += 1
        representatives = short_rational_representatives(
            residue,
            modulus,
            coefficient_radius=coefficient_radius,
            limit=representatives_per_class,
        )
        for representative in representatives:
            records.append(
                {
                    "numerator": representative.numerator,
                    "denominator": representative.denominator,
                    "u": rational_text(Q(representative.numerator, representative.denominator)),
                    "height": representative.height,
                    "tier": "_".join(group_names),
                    "crt_residue": residue,
                    "crt_modulus": modulus,
                    "groups": list(group_names),
                    "local_residues": {
                        group.name: local_residue
                        for group, local_residue in zip(groups, residues, strict=True)
                    },
                }
            )
    expected = 1
    for group in groups:
        expected *= len(group.residues)
    if class_count != expected:
        raise AssertionError("CRT tier class count changed")
    return records, {
        "tier": "_".join(group_names),
        "groups": list(group_names),
        "crt_modulus": records[0]["crt_modulus"],
        "class_count": class_count,
        "representatives_before_global_deduplication": len(records),
    }


def exact_candidate(record: dict[str, Any]) -> dict[str, Any] | None:
    record = {
        **record,
        "local_residues": dict(record["local_residues"]),
    }
    numerator = int(record["numerator"])
    denominator = int(record["denominator"])
    if numerator == 0:
        return None
    sign_normalized = numerator < 0
    if sign_normalized:
        numerator = -numerator
        record["numerator"] = numerator
        record["u"] = rational_text(Q(numerator, denominator))
        record["crt_residue"] = (-int(record["crt_residue"])) % int(
            record["crt_modulus"]
        )
        record["local_residues"] = {
            group_name: (-int(residue)) % GROUPS[group_name].modulus
            for group_name, residue in record["local_residues"].items()
        }
    value = homogeneous_value(BASE_DISCRIMINANT, numerator, denominator)
    if value == 0:
        return None
    valuations = {}
    for group_name in record["groups"]:
        group = GROUPS[group_name]
        if gcd(denominator, group.modulus) != 1:
            raise AssertionError("a lattice denominator is not a constrained unit")
        residue = numerator * pow(denominator, -1, group.modulus) % group.modulus
        if residue != record["local_residues"][group_name]:
            raise AssertionError("a lattice representative lost its local residue")
        actual = valuation_integer(value, group.prime)
        if actual < group.forced_discriminant_valuation:
            raise AssertionError("a lattice representative lost its forced valuation")
        valuations[str(group.prime)] = actual
    radical = radical_upper(value * numerator * denominator)
    result = {
        **record,
        "base_T": rational_text(base_parameter(Q(numerator, denominator))),
        "homogeneous_primitive_discriminant": str(value),
        "forced_prime_actual_valuations": valuations,
        "radical_upper_trial_bound": 997,
        "radical_upper": str(radical),
        "log_radical_upper": f"{log_integer(radical):.15f}",
        "sign_normalized_using_u_to_minus_u_curve_symmetry": sign_normalized,
    }
    return result


def merge_candidates(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[tuple[int, int], dict[str, Any]] = {}
    for raw in records:
        candidate = exact_candidate(raw)
        if candidate is None:
            continue
        key = (candidate["numerator"], candidate["denominator"])
        previous = merged.get(key)
        path = {
            key: candidate[key]
            for key in ("tier", "crt_residue", "crt_modulus", "groups", "local_residues")
        }
        if previous is None:
            candidate["forcing_paths"] = [path]
            for key_to_remove in path:
                candidate.pop(key_to_remove, None)
            merged[key] = candidate
        else:
            if path not in previous["forcing_paths"]:
                previous["forcing_paths"].append(path)
    return sorted(
        merged.values(),
        key=lambda item: (
            float(item["log_radical_upper"]),
            item["height"],
            abs(item["numerator"]),
            item["denominator"],
        ),
    )


def select_conductor_candidates(
    candidates: Sequence[dict[str, Any]], *, proxy_keep: int, height_keep: int, per_tier_keep: int
) -> list[dict[str, Any]]:
    selected: dict[tuple[int, int], dict[str, Any]] = {}

    def add(items: Iterable[dict[str, Any]]) -> None:
        for item in items:
            selected[(item["numerator"], item["denominator"])] = item

    add(candidates[:proxy_keep])
    add(sorted(candidates, key=lambda item: (item["height"], float(item["log_radical_upper"])))[:height_keep])
    for tier_names in TIERS:
        tier = "_".join(tier_names)
        members = [
            item
            for item in candidates
            if any(path["tier"] == tier for path in item["forcing_paths"])
        ]
        add(members[:per_tier_keep])
        add(sorted(members, key=lambda item: (item["height"], float(item["log_radical_upper"])))[:per_tier_keep])
    return sorted(
        selected.values(),
        key=lambda item: (
            float(item["log_radical_upper"]), item["height"], item["numerator"], item["denominator"]
        ),
    )


def conductor_worker(candidate: dict[str, Any]) -> dict[str, Any]:
    parameter_u = Q(candidate["numerator"], candidate["denominator"])
    primes = sorted(
        {
            GROUPS[group_name].prime
            for path in candidate["forcing_paths"]
            for group_name in path["groups"]
        }
    )
    global_data = minimal_curve_data(
        base_changed_short_jacobian_coefficients(parameter_u),
        timeout=90,
        local_primes=primes,
        stack_bytes=STACK_BYTES,
    )
    local_checks = {}
    for prime in primes:
        forced = max(
            GROUPS[group_name].forced_discriminant_valuation
            for path in candidate["forcing_paths"]
            for group_name in path["groups"]
            if GROUPS[group_name].prime == prime
        )
        local = global_data["local_reduction"][str(prime)]
        valid = (
            local["conductor_exponent"] == 1
            and local["minimal_c4_valuation"] == 0
            and local["minimal_discriminant_valuation"] >= forced
        )
        if not valid:
            raise AssertionError(f"engineered p={prime} class was not clean multiplicative")
        local_checks[str(prime)] = {
            "forced_minimal_discriminant_valuation": forced,
            "verified": True,
        }
    return {
        **candidate,
        "status": "completed exact conductor and local-reduction computation",
        "global_curve": global_data,
        "below_strict_log_conductor_182_72": Decimal(str(global_data["log_conductor"]))
        < TARGET_LOG_CONDUCTOR,
        "engineered_local_reduction_checks": local_checks,
    }


def rank_restricted_columns(certificate: dict[str, Any], column_count: int) -> int:
    rows = []
    for signature in certificate["signatures"]:
        rows.extend(tuple(row[:column_count]) for row in signature["rows"])
    rank, _ = gf_l_rank_and_pivots(rows, column_count, 3)
    return rank


def point_worker(candidate: dict[str, Any], height_bound: int) -> dict[str, Any]:
    parameter_u = Q(candidate["numerator"], candidate["denominator"])
    parameter_t = base_parameter(parameter_u)
    coefficients = base_changed_short_jacobian_coefficients(parameter_u)
    known = known_jacobian_points(parameter_u)
    by_x: dict[Fraction, tuple[tuple[Fraction, Fraction], str]] = {}
    for index, point in enumerate(known):
        by_x.setdefault(point[0], (point, f"known-{index:02d}"))
    raw = bounded_quartic_points(
        CONSTRUCTION.primitive_quartic_coefficients(parameter_t),
        height_bound=height_bound,
        timeout=120,
        stack_bytes=STACK_BYTES,
    )
    signless = canonical_signless_points(raw)
    quartic_records = []
    for index, quartic_point in enumerate(signless):
        jacobian_point = quartic_point_to_jacobian(
            CONSTRUCTION, parameter_t, quartic_point
        )
        by_x.setdefault(jacobian_point[0], (jacobian_point, f"searched-{index:04d}"))
        quartic_records.append(point_record(quartic_point))
    pool = tuple(item[0] for item in by_x.values())
    sources = tuple(item[1] for item in by_x.values())
    known_columns = sum(source.startswith("known-") for source in sources)
    certificate = mod3_independence_certificate(
        coefficients, pool, prime_bound=499
    )
    known_rank = rank_restricted_columns(certificate, known_columns)
    if known_rank != 13:
        raise AssertionError("the base-changed generic rank-13 subgroup changed")
    return {
        **candidate,
        "status": "completed bounded point search and exact finite-reduction certification",
        "point_search_height_bound": height_bound,
        "known_jacobian_point_count_before_inverse_deduplication": len(known),
        "known_jacobian_columns_modulo_inverse": known_columns,
        "exact_known_subgroup_dimension": known_rank,
        "signed_quartic_points_returned": len(raw),
        "signless_quartic_points_returned": len(signless),
        "searched_quartic_points": quartic_records,
        "pool_point_count_modulo_inverse": len(pool),
        "pool_point_sha256": point_digest(pool),
        "pool_sources": list(sources),
        "exact_specialization_rank_lower_bound": certificate["combined_exact_rank_over_F3"],
        "finite_reduction_certificate": certificate,
        "bounded_search_is_not_a_rank_upper_bound": True,
    }


def candidate_identifier(candidate: dict[str, Any]) -> str:
    return f"u{candidate['numerator']}_{candidate['denominator']}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    parser.add_argument("--coefficient-radius", type=int, default=12)
    parser.add_argument("--representatives-per-class", type=int, default=12)
    parser.add_argument("--proxy-keep", type=int, default=120)
    parser.add_argument("--height-keep", type=int, default=80)
    parser.add_argument("--per-tier-keep", type=int, default=24)
    parser.add_argument("--point-keep", type=int, default=160)
    parser.add_argument("--point-height", type=int, default=200_000)
    parser.add_argument("--workers", type=int, default=8)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 12:
        raise SystemExit("--workers must lie in [1,12]")
    output = ROOT / args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    conductor_directory = output / "conductor-records"
    point_directory = output / "point-certificates-h200000"
    conductor_directory.mkdir(exist_ok=True)
    point_directory.mkdir(exist_ok=True)

    group_certificates = [group_certificate(group) for group in GROUPS.values()]
    all_raw = []
    tier_summaries = []
    for tier in TIERS:
        records, summary = enumerate_tier(
            tier,
            coefficient_radius=args.coefficient_radius,
            representatives_per_class=args.representatives_per_class,
        )
        all_raw.extend(records)
        tier_summaries.append(summary)
        print(
            f"tier {summary['tier']}: classes={summary['class_count']} representatives={len(records)}",
            flush=True,
        )
    candidates = merge_candidates(all_raw)
    selected = select_conductor_candidates(
        candidates,
        proxy_keep=args.proxy_keep,
        height_keep=args.height_keep,
        per_tier_keep=args.per_tier_keep,
    )
    input_payload = {
        "scope": "frozen exact multiple-root CRT/Gauss promotion population",
        "base_discriminant_coefficients_ascending": list(BASE_DISCRIMINANT),
        "base_discriminant_content_removed": BASE_DISCRIMINANT_CONTENT,
        "groups": group_certificates,
        "tiers": tier_summaries,
        "raw_representatives_before_global_deduplication": len(all_raw),
        "unique_nonsingular_representatives": len(candidates),
        "selected_conductor_count": len(selected),
        "selection": {
            "proxy_keep": args.proxy_keep,
            "height_keep": args.height_keep,
            "per_tier_proxy_and_height_keep": args.per_tier_keep,
            "rank_score_used": False,
        },
        "candidates": selected,
    }
    input_payload["result_sha256"] = canonical_digest(input_payload)
    atomic_json(output / "candidate-input.json", input_payload)

    conductor_results = []
    pending = []
    for candidate in selected:
        path = conductor_directory / f"{candidate_identifier(candidate)}.json"
        if path.exists():
            conductor_results.append(json.loads(path.read_text()))
        else:
            pending.append(candidate)
    print(f"conductors cached={len(conductor_results)} pending={len(pending)}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(conductor_worker, candidate): candidate for candidate in pending}
        for index, future in enumerate(as_completed(futures), start=1):
            candidate = futures[future]
            try:
                record = future.result()
            except Exception as error:
                record = {**candidate, "status": "error", "error": repr(error)}
            record["result_sha256"] = canonical_digest(record)
            atomic_json(conductor_directory / f"{candidate_identifier(candidate)}.json", record)
            conductor_results.append(record)
            print(
                f"conductor {index}/{len(pending)} {candidate_identifier(candidate)} "
                f"lnN={record.get('global_curve', {}).get('log_conductor')} "
                f"W={record.get('global_curve', {}).get('root_number')} status={record['status']}",
                flush=True,
            )

    completed = [
        record for record in conductor_results if record["status"].startswith("completed")
    ]
    completed.sort(
        key=lambda item: (
            not item["below_strict_log_conductor_182_72"],
            item["global_curve"]["root_number"] != -1,
            Decimal(str(item["global_curve"]["log_conductor"])),
            float(item["log_radical_upper"]),
            item["height"],
        )
    )
    below = [item for item in completed if item["below_strict_log_conductor_182_72"]]
    point_selected = below[: args.point_keep]
    if len(point_selected) < args.point_keep:
        used = {(item["numerator"], item["denominator"]) for item in point_selected}
        point_selected.extend(
            item
            for item in completed
            if (item["numerator"], item["denominator"]) not in used
        )
        point_selected = point_selected[: args.point_keep]
    point_input = {
        "scope": "conductor-first bounded point population",
        "point_height_bound": args.point_height,
        "point_keep": args.point_keep,
        "target_qualified_available": len(below),
        "ordering": "target qualified, W=-1 first, exact lnN, radical upper, height",
        "candidates": point_selected,
    }
    point_input["result_sha256"] = canonical_digest(point_input)
    atomic_json(output / "point-search-input.json", point_input)

    point_results = []
    pending_points = []
    for candidate in point_selected:
        path = point_directory / f"{candidate_identifier(candidate)}.json"
        if path.exists():
            point_results.append(json.loads(path.read_text()))
        else:
            pending_points.append(candidate)
    print(f"points cached={len(point_results)} pending={len(pending_points)}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(point_worker, candidate, args.point_height): candidate
            for candidate in pending_points
        }
        for index, future in enumerate(as_completed(futures), start=1):
            candidate = futures[future]
            try:
                record = future.result()
            except Exception as error:
                record = {**candidate, "status": "error", "error": repr(error)}
            record["result_sha256"] = canonical_digest(record)
            atomic_json(point_directory / f"{candidate_identifier(candidate)}.json", record)
            point_results.append(record)
            rank = record.get("exact_specialization_rank_lower_bound")
            print(
                f"point {index}/{len(pending_points)} {candidate_identifier(candidate)} "
                f"rank={rank} lnN={candidate.get('global_curve', {}).get('log_conductor')} "
                f"W={candidate.get('global_curve', {}).get('root_number')} status={record['status']}",
                flush=True,
            )
            if rank is not None and rank > 16:
                print(f"ALERT exact LB>16 {candidate_identifier(candidate)} rank={rank}", flush=True)

    point_results.sort(
        key=lambda item: (
            -item.get("exact_specialization_rank_lower_bound", -1),
            Decimal(str(item.get("global_curve", {}).get("log_conductor", "1e999"))),
            item["height"],
        )
    )
    summary = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed bounded exact multiple-root CRT specialization search",
        "family": {
            "roots": list(ROOTS),
            "base_change": "T=(14406-u^2)/(2u)",
            "generic_exact_rank_lower_bound": 13,
        },
        "exact_base_changed_form": {
            "direct_jacobian_A_B_coefficients_ascending": [
                list(values) for values in derive_even_coefficients(CONSTRUCTION)
            ],
            "homogenization": {
                "h": "14406*b^2-a^2",
                "d": "2*a*b",
                "A_h": "d^8*A(h/d)",
                "B_h": "d^12*B(h/d)",
                "primitive_F_h": "d^20*D(h/d)/47775744",
                "scaled_model_discriminant": "16*3^12*47775744*d^4*F_h",
            },
            "primitive_base_discriminant_coefficients_ascending": list(BASE_DISCRIMINANT),
            "removed_content": BASE_DISCRIMINANT_CONTENT,
        },
        "candidate_input": {
            "path": str((output / "candidate-input.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "candidate-input.json"),
            "result_sha256": input_payload["result_sha256"],
        },
        "point_search_input": {
            "path": str((output / "point-search-input.json").relative_to(ROOT)),
            "sha256": sha256_file(output / "point-search-input.json"),
            "result_sha256": point_input["result_sha256"],
        },
        "conductor_summary": {
            "selected": len(selected),
            "completed": len(completed),
            "errors": len(conductor_results) - len(completed),
            "target_qualified": len(below),
            "target_qualified_W_minus_1": sum(
                item["global_curve"]["root_number"] == -1 for item in below
            ),
        },
        "point_summary": {
            "selected": len(point_selected),
            "completed": sum(item["status"].startswith("completed") for item in point_results),
            "errors": sum(item["status"] == "error" for item in point_results),
            "rank_distribution": dict(
                sorted(
                    Counter(
                        item.get("exact_specialization_rank_lower_bound")
                        for item in point_results
                        if item.get("exact_specialization_rank_lower_bound") is not None
                    ).items()
                )
            ),
            "rank_above_16": [
                {"u": item["u"], "rank": item["exact_specialization_rank_lower_bound"]}
                for item in point_results
                if item.get("exact_specialization_rank_lower_bound", -1) > 16
            ],
        },
        "point_results": point_results,
        "scope_warning": "bounded H=200000 quartic search; nonpromotion is not a rank upper bound",
        "provenance": {
            "script_path": str(Path(__file__).resolve().relative_to(ROOT)),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "formula_module_path": "elliptic-curves/cas/mestre_rank13_02393128133175.py",
            "formula_module_sha256": sha256_file(
                ROOT / "elliptic-curves/cas/mestre_rank13_02393128133175.py"
            ),
            "reproducing_command": "PYTHONPATH=elliptic-curves/cas python3 " + " ".join(sys.argv),
        },
        "software": {"python": platform.python_version(), "platform": platform.platform()},
    }
    summary["result_sha256"] = canonical_digest(
        {key: value for key, value in summary.items() if key != "generated_at_utc"}
    )
    atomic_json(output / "summary.json", summary)
    print(
        f"complete conductor={len(completed)} target={len(below)} "
        f"ranks={summary['point_summary']['rank_distribution']} "
        f"rank>16={summary['point_summary']['rank_above_16']} "
        f"sha={summary['result_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()

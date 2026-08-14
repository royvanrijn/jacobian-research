#!/usr/bin/env python3
"""Freeze a disjoint conductor population near Mestre T=2731/36.

Every emitted parameter is a primitive positive rational with reduced
denominator 257..2048, so this population is exactly disjoint from the prior
thirteen-family box whose reduced denominators were at most 256.  The union
combines a narrow numerator annulus, Gauss/Farey perturbation rays, simple
Hensel lifts of discriminant roots modulo p^2, and fresh favorable local-trace
residue rays.  Exact canonical values are deduplicated and manifest-hashed.

Fresh primes 919..971 rank the complete population.  Only the top 2048 are
then evaluated on the disjoint held band 977..1021 and by exact homogeneous
discriminant features.  Fixed rank-blind local-score/radical/power/height and
source/denominator-diversity quotas select the conductor population.  This
script makes no conductor, point, height, or rank call.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
import os
from pathlib import Path
import platform
import sys
from typing import Any, Sequence

from search_mestre_root_tuple_scale import sha256_file
from search_mestre_root_tuple_scale_max100 import stable_json_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 7, 93, 154, 161, 191)
CENTER = Q(2731, 36)
DENOMINATOR_MIN = 257
DENOMINATOR_MAX = 2048
ORDINARY_OFFSET_MAX = 24
POWER_AND_TRACE_OFFSET_MAX = 96
GAUSS_M_MIN = 8
GAUSS_M_MAX = 56
GAUSS_COEFFICIENT_BOUND = 4
POWER_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
DISCOVERY_PRIMES = (919, 929, 937, 941, 947, 953, 967, 971)
HELD_PRIMES = (977, 983, 991, 997, 1009, 1013, 1019, 1021)
DISCOVERY_KEEP = 2048
TRIAL_DIVISION_LIMIT = 997
FRONTIER_FILENAME = "elliptic_mestre_rank13_multifamily_rational.json"
EXPECTED_FRONTIER_SHA256 = (
    "0f664e937b9983bd7fa1cfb80269b5c734faddbf6d02dbe4dfca0e3b573ac41f"
)
CERTIFICATE_FILENAME = "elliptic_mestre_rank15_2731_36.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "5f91987e9fd21887afbe0cd376e7b56844a37e0e70ade6fc713aaa3121e87c1a"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_rank15_2731_36_disjoint_neighborhood_selection.json"
)


@dataclass(frozen=True)
class LocalSymbol:
    good: bool
    trace: int
    score_units: int


def nearest_integer(value: Fraction) -> int:
    """Round a nonnegative rational half-up, deterministically."""

    value = Q(value)
    if value < 0:
        return -nearest_integer(-value)
    quotient, remainder = divmod(value.numerator, value.denominator)
    return quotient + (2 * remainder >= value.denominator)


def evaluate_polynomial(coefficients: Sequence[int], value: int, modulus: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % modulus
    return answer


def derivative(coefficients: Sequence[int]) -> tuple[int, ...]:
    return tuple(index * coefficient for index, coefficient in enumerate(coefficients[1:], 1))


def simple_hensel_roots_mod_p2(
    coefficients: Sequence[int], prime: int
) -> tuple[int, ...]:
    derivative_coefficients = derivative(coefficients)
    roots = []
    for residue in range(prime):
        if evaluate_polynomial(coefficients, residue, prime) != 0:
            continue
        if evaluate_polynomial(derivative_coefficients, residue, prime) == 0:
            continue
        modulus = prime * prime
        lifts = [
            residue + prime * digit
            for digit in range(prime)
            if evaluate_polynomial(coefficients, residue + prime * digit, modulus) == 0
        ]
        if len(lifts) != 1:
            raise AssertionError("a simple root did not have one Hensel lift")
        roots.append(lifts[0])
    return tuple(roots)


def quadratic_characters(prime: int) -> tuple[int, ...]:
    answer = [-1] * prime
    answer[0] = 0
    for value in range(1, prime):
        answer[value * value % prime] = 1
    return tuple(answer)


def local_symbol(
    a_coefficients: Sequence[int], b_coefficients: Sequence[int],
    residue: int | None, prime: int, characters: Sequence[int],
) -> LocalSymbol:
    if residue is None:
        coefficient_a = a_coefficients[-1] % prime
        coefficient_b = b_coefficients[-1] % prime
    else:
        coefficient_a = evaluate_polynomial(a_coefficients, residue, prime)
        coefficient_b = evaluate_polynomial(b_coefficients, residue, prime)
    if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
        return LocalSymbol(False, 0, 0)
    character_sum = 0
    for x_value in range(prime):
        rhs = (x_value**3 + coefficient_a * x_value + coefficient_b) % prime
        character_sum += characters[rhs]
    trace = -character_sum
    score = ((2 - trace) / (prime + 1 - trace)) * log(float(prime))
    units = int(score * 1.0e12 + (0.5 if score >= 0 else -0.5))
    return LocalSymbol(True, trace, units)


def build_tables(
    a_coefficients: Sequence[int], b_coefficients: Sequence[int], primes: Sequence[int]
) -> dict[int, tuple[LocalSymbol, ...]]:
    tables = {}
    for prime in primes:
        characters = quadratic_characters(prime)
        tables[prime] = tuple(
            [
                local_symbol(a_coefficients, b_coefficients, residue, prime, characters)
                for residue in range(prime)
            ]
            + [local_symbol(a_coefficients, b_coefficients, None, prime, characters)]
        )
    return tables


def table_digest(tables: dict[int, tuple[LocalSymbol, ...]]) -> str:
    return stable_json_digest(
        [
            [prime, [[symbol.good, symbol.trace] for symbol in tables[prime]]]
            for prime in sorted(tables)
        ]
    )


def score_parameter(
    parameter: Fraction, primes: Sequence[int],
    tables: dict[int, tuple[LocalSymbol, ...]],
) -> tuple[int, int]:
    units = 0
    good = 0
    for prime in primes:
        if parameter.denominator % prime == 0:
            symbol = tables[prime][prime]
        else:
            residue = (
                parameter.numerator
                * pow(parameter.denominator % prime, prime - 2, prime)
            ) % prime
            symbol = tables[prime][residue]
        if symbol.good:
            units += symbol.score_units
            good += 1
    return units, good


def score_text(units: int) -> str:
    sign = "-" if units < 0 else ""
    absolute = abs(units)
    return (
        f"{sign}{absolute // 1_000_000_000_000}."
        f"{absolute % 1_000_000_000_000:012d}"
    )


def primes_up_to(bound: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(bound**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : bound + 1 : prime] = b"\x00" * (
                (bound - prime * prime) // prime + 1
            )
    return tuple(index for index, value in enumerate(sieve) if value)


TRIAL_PRIMES = primes_up_to(TRIAL_DIVISION_LIMIT)


def complete_radical(value: int) -> int:
    remaining = abs(value)
    radical = 1
    for prime in TRIAL_PRIMES:
        if prime * prime > remaining:
            break
        if remaining % prime == 0:
            radical *= prime
            while remaining % prime == 0:
                remaining //= prime
    if remaining > 1:
        radical *= remaining
    return radical


def homogeneous_value(coefficients: Sequence[int], numerator: int, denominator: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def discriminant_feature(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> dict[str, Any]:
    absolute = abs(homogeneous_value(coefficients, numerator, denominator))
    if absolute == 0:
        return {"singular": True, "absolute_homogeneous_discriminant": "0"}
    remaining = absolute
    valuations = []
    known_radical = 1
    known_powerful = 1
    for prime in TRIAL_PRIMES:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            valuations.append([prime, exponent])
            known_radical *= prime
            if exponent > 1:
                known_powerful *= prime ** (exponent - 1)
    denominator_radical = complete_radical(denominator)
    radical_upper = known_radical * remaining * denominator_radical
    return {
        "singular": False,
        "absolute_homogeneous_discriminant": str(absolute),
        "small_prime_valuations": valuations,
        "residual_cofactor": str(remaining),
        "residual_cofactor_bit_length": remaining.bit_length(),
        "known_discriminant_radical": str(known_radical),
        "known_powerful_part": str(known_powerful),
        "denominator_radical": denominator_radical,
        "combined_radical_upper_bound": str(radical_upper),
        "upper_bound_semantics": "rad(residual)<=residual exactly",
    }


def source_category(source: str) -> str:
    if source.startswith("ordinary"):
        return "ordinary-annulus"
    if source.startswith("gauss"):
        return "gauss-farey"
    if source.startswith("power"):
        return "discriminant-power"
    if source.startswith("trace"):
        return "local-trace"
    raise AssertionError("unknown generator source")


def generate_population(
    discriminant: Sequence[int], discovery_tables: dict[int, tuple[LocalSymbol, ...]]
) -> tuple[dict[Fraction, set[str]], dict[str, Any]]:
    population: dict[Fraction, set[str]] = {}
    attempted: dict[str, int] = {}
    rejected_nonprimitive_or_outside = 0

    def add(numerator: int, denominator: int, source: str) -> None:
        nonlocal rejected_nonprimitive_or_outside
        attempted[source] = attempted.get(source, 0) + 1
        if numerator <= 0 or denominator <= 0:
            rejected_nonprimitive_or_outside += 1
            return
        divisor = gcd(numerator, denominator)
        parameter = Q(numerator // divisor, denominator // divisor)
        if not DENOMINATOR_MIN <= parameter.denominator <= DENOMINATOR_MAX:
            rejected_nonprimitive_or_outside += 1
            return
        if parameter == CENTER:
            raise AssertionError("the fixed center leaked into the disjoint annulus")
        population.setdefault(parameter, set()).add(source)

    for denominator in range(DENOMINATOR_MIN, DENOMINATOR_MAX + 1):
        center_numerator = nearest_integer(CENTER * denominator)
        for offset in range(-ORDINARY_OFFSET_MAX, ORDINARY_OFFSET_MAX + 1):
            if offset:
                add(center_numerator + offset, denominator, "ordinary")

    for multiplier in range(GAUSS_M_MIN, GAUSS_M_MAX + 1):
        for c_value in range(-GAUSS_COEFFICIENT_BOUND, GAUSS_COEFFICIENT_BOUND + 1):
            for d_value in range(-GAUSS_COEFFICIENT_BOUND, GAUSS_COEFFICIENT_BOUND + 1):
                if c_value == d_value == 0:
                    continue
                add(
                    CENTER.numerator * multiplier + c_value,
                    CENTER.denominator * multiplier + d_value,
                    f"gauss_c{c_value}_d{d_value}",
                )

    power_lifts = {}
    for prime in POWER_PRIMES:
        roots = simple_hensel_roots_mod_p2(discriminant, prime)
        power_lifts[str(prime)] = list(roots)
        modulus = prime * prime
        for root in roots:
            for denominator in range(DENOMINATOR_MIN, DENOMINATOR_MAX + 1):
                residue = root * denominator % modulus
                quotient = nearest_integer((CENTER * denominator - residue) / modulus)
                numerator = residue + modulus * quotient
                if abs(CENTER.denominator * numerator - CENTER.numerator * denominator) > (
                    CENTER.denominator * POWER_AND_TRACE_OFFSET_MAX
                ):
                    continue
                add(numerator, denominator, f"power_p{prime}_r{root}")

    favorable_residues = {}
    for prime in DISCOVERY_PRIMES:
        ranked = sorted(
            (
                (-symbol.score_units, residue)
                for residue, symbol in enumerate(discovery_tables[prime][:-1])
                if symbol.good
            )
        )
        residues = tuple(residue for _, residue in ranked[:2])
        favorable_residues[str(prime)] = list(residues)
        for residue in residues:
            for denominator in range(DENOMINATOR_MIN, DENOMINATOR_MAX + 1):
                congruent = residue * denominator % prime
                quotient = nearest_integer((CENTER * denominator - congruent) / prime)
                numerator = congruent + prime * quotient
                if abs(CENTER.denominator * numerator - CENTER.numerator * denominator) > (
                    CENTER.denominator * POWER_AND_TRACE_OFFSET_MAX
                ):
                    continue
                add(numerator, denominator, f"trace_p{prime}_r{residue}")

    if not population or min(parameter.denominator for parameter in population) < DENOMINATOR_MIN:
        raise AssertionError("the disjoint denominator gate failed")
    digest = hashlib.sha256()
    for parameter in sorted(population, key=lambda value: (value.denominator, value.numerator)):
        digest.update(
            (
                f"{parameter.numerator}/{parameter.denominator}|"
                f"{','.join(sorted(population[parameter]))}\n"
            ).encode()
        )
    category_counts = {
        category: sum(
            any(source_category(source) == category for source in sources)
            for sources in population.values()
        )
        for category in (
            "ordinary-annulus", "gauss-farey", "discriminant-power", "local-trace"
        )
    }
    return population, {
        "raw_attempts_by_exact_source": attempted,
        "raw_attempt_count": sum(attempted.values()),
        "nonprimitive_or_reduced_outside_annulus_rejections": rejected_nonprimitive_or_outside,
        "unique_reduced_parameter_count": len(population),
        "parameters_by_source_category_after_dedup": category_counts,
        "simple_hensel_roots_mod_p2": power_lifts,
        "favorable_trace_residues_top_two_per_prime": favorable_residues,
        "canonical_parameter_source_manifest_sha256": digest.hexdigest(),
    }


def select_conductor_population(pool: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reasons: dict[str, set[str]] = {}

    def take(label: str, order: Sequence[dict[str, Any]], count: int) -> None:
        for record in order[:count]:
            reasons.setdefault(record["parameter"], set()).add(label)

    held_order = sorted(
        pool,
        key=lambda row: (
            -Decimal(row["held_score"]), -Decimal(row["discovery_score"]),
            row["denominator"], row["numerator"],
        ),
    )
    radical_order = sorted(
        pool,
        key=lambda row: (
            int(row["discriminant_feature"]["combined_radical_upper_bound"]),
            -Decimal(row["held_score"]), row["denominator"], row["numerator"],
        ),
    )
    powerful_order = sorted(
        pool,
        key=lambda row: (
            -int(row["discriminant_feature"]["known_powerful_part"]),
            int(row["discriminant_feature"]["combined_radical_upper_bound"]),
            row["denominator"], row["numerator"],
        ),
    )
    height_order = sorted(
        pool,
        key=lambda row: (
            max(row["numerator"], row["denominator"]),
            row["denominator"], row["numerator"],
        ),
    )
    take("highest-held-score", held_order, 16)
    take("smallest-radical-upper-bound", radical_order, 16)
    take("largest-known-powerful-part", powerful_order, 8)
    take("lowest-projective-height", height_order, 8)
    for category in (
        "ordinary-annulus", "gauss-farey", "discriminant-power", "local-trace"
    ):
        order = [row for row in held_order if category in row["source_categories"]]
        take(f"held-score-source-{category}", order, 4)
    for lower, upper in ((257, 704), (705, 1152), (1153, 1600), (1601, 2048)):
        order = [row for row in held_order if lower <= row["denominator"] <= upper]
        take(f"held-score-denominator-{lower}-{upper}", order, 2)

    by_parameter = {row["parameter"]: row for row in pool}
    selected = []
    digest = hashlib.sha256()
    for parameter in sorted(
        reasons,
        key=lambda text: (Q(text).denominator, Q(text).numerator),
    ):
        record = dict(by_parameter[parameter])
        record["conductor_selection_strata"] = sorted(reasons[parameter])
        selected.append(record)
        digest.update(
            f"{parameter}|{','.join(record['conductor_selection_strata'])}\n".encode()
        )
    return selected, {
        "selection_uses_conductor_or_point_or_rank_data": False,
        "discovery_population_closed_before_held_scores": True,
        "quota_rule": {
            "highest_held": 16,
            "smallest_radical_upper": 16,
            "largest_known_powerful": 8,
            "lowest_projective_height": 8,
            "held_per_source_category": 4,
            "held_per_denominator_quartile": 2,
        },
        "selected_population": len(selected),
        "selected_population_sha256": digest.hexdigest(),
    }


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=root / "artifacts/generated-results" / DEFAULT_OUTPUT.name,
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite the neighborhood selection artifact")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    generated = root / "artifacts/generated-results"
    frontier_path = generated / FRONTIER_FILENAME
    certificate_path = generated / CERTIFICATE_FILENAME
    if sha256_file(frontier_path) != EXPECTED_FRONTIER_SHA256:
        raise AssertionError("the pinned multifamily frontier changed")
    if sha256_file(certificate_path) != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("the pinned rank-15 certificate changed")
    frontier = json.loads(frontier_path.read_text())
    family = next(row for row in frontier["families"] if tuple(row["roots"]) == ROOTS)
    a_coefficients = tuple(int(value) for value in family["A_coefficients_ascending"])
    b_coefficients = tuple(int(value) for value in family["B_coefficients_ascending"])
    discriminant = tuple(
        int(value) for value in family["content_free_discriminant_coefficients_ascending"]
    )
    if len(a_coefficients) != 9 or len(b_coefficients) != 13 or len(discriminant) != 21:
        raise AssertionError("the exact family formulas changed")

    all_tables = build_tables(
        a_coefficients, b_coefficients, DISCOVERY_PRIMES + HELD_PRIMES
    )
    discovery_tables = {prime: all_tables[prime] for prime in DISCOVERY_PRIMES}
    held_tables = {prime: all_tables[prime] for prime in HELD_PRIMES}
    population, population_audit = generate_population(discriminant, discovery_tables)
    print(
        f"population closed count={len(population)} "
        f"manifest={population_audit['canonical_parameter_source_manifest_sha256']}",
        flush=True,
    )

    ranked = []
    for parameter, sources in population.items():
        discovery_units, discovery_good = score_parameter(
            parameter, DISCOVERY_PRIMES, discovery_tables
        )
        ranked.append((
            -discovery_units, -discovery_good,
            parameter.denominator, parameter.numerator, parameter, sources,
            discovery_units, discovery_good,
        ))
    ranked.sort()
    retained = ranked[:DISCOVERY_KEEP]
    if len(retained) != DISCOVERY_KEEP:
        raise AssertionError("the discovery survivor gate shrank")

    pool = []
    feature_digest = hashlib.sha256()
    singular = 0
    for _, _, _, _, parameter, sources, discovery_units, discovery_good in retained:
        held_units, held_good = score_parameter(parameter, HELD_PRIMES, held_tables)
        feature = discriminant_feature(
            discriminant, parameter.numerator, parameter.denominator
        )
        if feature["singular"]:
            singular += 1
            continue
        categories = sorted({source_category(source) for source in sources})
        record = {
            "numerator": parameter.numerator,
            "denominator": parameter.denominator,
            "parameter": str(parameter),
            "sources": sorted(sources),
            "source_categories": categories,
            "discovery_score": score_text(discovery_units),
            "held_score": score_text(held_units),
            "discovery_good": discovery_good,
            "held_good": held_good,
            "discriminant_feature": feature,
        }
        pool.append(record)
        feature_digest.update(
            (
                f"{parameter}|{record['discovery_score']}|{record['held_score']}|"
                f"{','.join(record['sources'])}|"
                f"{feature['absolute_homogeneous_discriminant']}|"
                f"{feature['combined_radical_upper_bound']}|"
                f"{feature['known_powerful_part']}\n"
            ).encode()
        )
    selected, selection = select_conductor_population(pool)
    print(
        f"feature pool={len(pool)} singular={singular} selected={len(selected)}",
        flush=True,
    )

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "closed rank-blind selection for disjoint T=2731/36 neighborhood",
        "scope": {
            "roots": list(ROOTS),
            "center": str(CENTER),
            "center_excluded": True,
            "reduced_denominator_interval": [DENOMINATOR_MIN, DENOMINATOR_MAX],
            "prior_common_box_reduced_denominator_maximum": 256,
            "exactly_disjoint_from_prior_common_box": True,
            "positive_T_only_under_exact_sign_quotient": True,
        },
        "family": {
            "A_coefficients_ascending": list(a_coefficients),
            "B_coefficients_ascending": list(b_coefficients),
            "content_free_discriminant_coefficients_ascending": list(discriminant),
        },
        "generator": {
            "ordinary_offset_range": [-ORDINARY_OFFSET_MAX, ORDINARY_OFFSET_MAX],
            "ordinary_zero_offset_excluded": True,
            "gauss_multiplier_interval": [GAUSS_M_MIN, GAUSS_M_MAX],
            "gauss_coefficient_interval": [-GAUSS_COEFFICIENT_BOUND, GAUSS_COEFFICIENT_BOUND],
            "power_primes": list(POWER_PRIMES),
            "power_and_trace_integer_offset_cap": POWER_AND_TRACE_OFFSET_MAX,
            **population_audit,
        },
        "local_trace_screen": {
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_relative_to_prior_multifamily_bands_through_911": min(
                DISCOVERY_PRIMES + HELD_PRIMES
            ) > 911,
            "discovery_table_digest": table_digest(discovery_tables),
            "held_table_digest": table_digest(held_tables),
            "complete_population_scored_on_discovery_band": len(population),
            "discovery_survivors_retained": len(retained),
            "held_scores_computed_only_after_discovery_survivors_closed": True,
        },
        "exact_discriminant_feature_screen": {
            "trial_division_prime_bound": TRIAL_DIVISION_LIMIT,
            "retained_before_singular_rejection": len(retained),
            "exact_singular_rejections": singular,
            "admissible_feature_pool_count": len(pool),
            "exact_feature_population_sha256": feature_digest.hexdigest(),
        },
        "conductor_selection": selection,
        "selected_records": selected,
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "frontier_path": str(frontier_path.relative_to(root)),
            "frontier_sha256": EXPECTED_FRONTIER_SHA256,
            "certificate_path": str(certificate_path.relative_to(root)),
            "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "external_conductor_point_height_rank_calls": 0,
        },
        "software": {"python": platform.python_version()},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "scope": artifact["scope"],
            "family": artifact["family"],
            "generator": artifact["generator"],
            "trace": artifact["local_trace_screen"],
            "features": artifact["exact_discriminant_feature_screen"],
            "selection": artifact["conductor_selection"],
            "records": artifact["selected_records"],
        }
    )
    exclusive_write(args.output, artifact)
    print(f"wrote {args.output} result={artifact['result_sha256']}", flush=True)


if __name__ == "__main__":
    main()

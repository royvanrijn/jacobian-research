#!/usr/bin/env python3
"""Closed rational scan of the thirteen max-root-200 H5000-rank-13 families.

The frozen max-root-200 artifact is the sole family input.  This script
reconstructs exactly the thirteen families whose panel follow-up had stable
numerical rank 13 at H=5000, and rejects the three rank-14 families and every
family that already received a rational scan.  It scans the common primitive
positive box 1<=a<=4096, 1<=b<=256, excluding the already-screened integer
panel T=1,...,8.  T and -T define the same curve in these even families.

An exact C++ scanner uses fresh discovery primes 811..857 and a disjoint
held-out band 859..911.  Per-family discovery survivors close before held
scores are read.  Equal rank-blind quotas then select local-score,
discriminant-radical, powerful-part, low-height, and denominator-quartile
leaders.  Every selected curve receives a conductor call before any point
search.  Fixed nested point tiers are H=5000, 50000, 250000, and 1000000.
Stable numerical rank at least 16 triggers an immediate exact mod-3
independence attempt.  The protocol stops after the declared H1m tier.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
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
import shlex
import shutil
import sys
import tempfile
import time
from typing import Any, Sequence

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_0430313946_frontier import exact_point_stage
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    capped_minimal_curve_data,
    point_digest,
    point_on_short_curve,
    run_capped_process,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

TARGET_LOG_CONDUCTOR = Decimal("182.72")
PLAUSIBLE_LOG_CONDUCTOR = Decimal("190")
FROZEN_MAX200_ARTIFACT = "elliptic_mestre_root_tuple_scale_max200.json"
EXPECTED_MAX200_ARTIFACT_SHA256 = (
    "5e1b53e187520735efba46fc8fd9cbdd4dfd4284545a815f6416baf3be84f342"
)
EXPECTED_MAX200_RESULT_SHA256 = (
    "546907f506326b4878c070ff33f0bea46716d2443cd758797a821aef71b09139"
)
EXPECTED_MAX200_SCRIPT_SHA256 = (
    "405a2b9f7653c89af0e3e6caf2e77765cb4bfc88fccf88edffa67d3435aebf24"
)

EXPECTED_RANK13_PANEL = (
    ((0, 12, 33, 142, 150, 169), 2),
    ((0, 12, 50, 93, 114, 131), 5),
    ((0, 21, 95, 100, 121, 155), 6),
    ((0, 23, 89, 124, 147, 181), 5),
    ((0, 23, 93, 128, 133, 175), 2),
    ((0, 25, 83, 124, 149, 183), 7),
    ((0, 26, 53, 70, 88, 117), 7),
    ((0, 32, 65, 97, 108, 148), 6),
    ((0, 40, 55, 100, 108, 151), 2),
    ((0, 5, 110, 111, 115, 133), 4),
    ((0, 7, 54, 127, 148, 166), 6),
    ((0, 7, 93, 154, 161, 191), 2),
    ((0, 8, 60, 93, 108, 125), 7),
)
EXPECTED_RANK14_PANEL = (
    ((0, 17, 142, 145, 162, 200), 7),
    ((0, 25, 57, 104, 116, 148), 1),
    ((0, 7, 121, 128, 183, 194), 1),
)

# These are all families with a closed rational-box frontier in this research
# directory before this run.  Two entries share the rank14-pair artifact.
PRIOR_RATIONAL_SCANS = (
    (
        (0, 4, 30, 31, 39, 46),
        "elliptic_mestre_0430313946_frontier.json",
        "546cfc676b28f6956808b2698260d3bab4f9490dab5f2efc195f487ab6a2e514",
    ),
    (
        (0, 6, 49, 73, 82, 96),
        "elliptic_mestre_0649738296_rational.json",
        "1c0e001c6c03e557722e16897f66ad2c90c93aa7d88f7cbbdd286700c66eaa78",
    ),
    (
        (0, 17, 142, 145, 162, 200),
        "elliptic_mestre_rank14_pair_rational_frontier.json",
        "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7",
    ),
    (
        (0, 25, 57, 104, 116, 148),
        "elliptic_mestre_02557104116148_direct_rational.json",
        "4874478c553c81ed69fffb49738b5975900a26a17d96f4dca9203a8244e75db6",
    ),
    (
        (0, 7, 121, 128, 183, 194),
        "elliptic_mestre_rank14_pair_rational_frontier.json",
        "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7",
    ),
)

DISCOVERY_PRIMES = (811, 821, 823, 827, 829, 839, 853, 857)
HELD_PRIMES = (859, 863, 877, 881, 883, 887, 907, 911)
PRIOR_PANEL_PARAMETERS = tuple(Q(value) for value in range(1, 9))
TRIAL_DIVISION_LIMIT = 997
FINITE_REDUCTION_TRIGGER = 16
STACK_BYTES = 512_000_000
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_rank13_multifamily_rational.json"
)


@dataclass(frozen=True)
class FamilySpec:
    index: int
    label: str
    roots: tuple[int, ...]
    calibration: Fraction
    a_coefficients: tuple[int, ...]
    b_coefficients: tuple[int, ...]
    source_record: dict[str, Any]


@dataclass(frozen=True)
class ScannerCandidate:
    numerator: int
    denominator: int
    discovery_score: str
    held_score: str
    discovery_good: int
    held_good: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)


@dataclass(frozen=True)
class ScannerResult:
    family_index: int
    numerator_bound: int
    denominator_bound: int
    keep: int
    primitive_population: int
    panel_excluded: int
    evaluated_population: int
    discovery_table_digest: str
    held_table_digest: str
    calibration: ScannerCandidate
    candidates: tuple[ScannerCandidate, ...]
    stdout_sha256: str


def poly_add(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    size = max(len(left), len(right))
    return tuple(
        (left[index] if index < len(left) else Q(0))
        + (right[index] if index < len(right) else Q(0))
        for index in range(size)
    )


def poly_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            answer[i + j] += a_value * b_value
    return tuple(answer)


def interpolate_consecutive_from_one(values: Sequence[Fraction]) -> tuple[int, ...]:
    differences = [Q(value) for value in values]
    newton = []
    while differences:
        newton.append(differences[0])
        differences = [right - left for left, right in zip(differences, differences[1:])]
    answer: tuple[Fraction, ...] = (Q(0),)
    basis: tuple[Fraction, ...] = (Q(1),)
    for index, coefficient in enumerate(newton):
        answer = poly_add(answer, tuple(coefficient * value for value in basis))
        basis = tuple(
            value / (index + 1)
            for value in poly_multiply(basis, (Q(-(index + 1)), Q(1)))
        )
    while len(answer) > 1 and answer[-1] == 0:
        answer = answer[:-1]
    if any(value.denominator != 1 for value in answer):
        raise AssertionError("an interpolated primitive coefficient was nonintegral")
    return tuple(value.numerator for value in answer)


def derive_even_coefficients(
    construction: SixRootMestreConstruction,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    values = [
        construction.primitive_jacobian_coefficients(Q(parameter))
        for parameter in range(1, 15)
    ]
    a_coefficients = interpolate_consecutive_from_one(
        [value[3] for value in values]
    )
    b_coefficients = interpolate_consecutive_from_one(
        [value[4] for value in values]
    )
    if len(a_coefficients) > 9 or len(b_coefficients) > 13:
        raise AssertionError("a Mestre Jacobian exceeded its degree bound")
    a_coefficients += (0,) * (9 - len(a_coefficients))
    b_coefficients += (0,) * (13 - len(b_coefficients))
    if any(a_coefficients[index] for index in range(1, 9, 2)) or any(
        b_coefficients[index] for index in range(1, 13, 2)
    ):
        raise AssertionError("the exact T-to-minus-T quotient changed")
    for parameter in (Q(17), Q(-17, 3), Q(11, 5)):
        expected = construction.primitive_jacobian_coefficients(parameter)
        actual = (
            Q(0), Q(0), Q(0),
            evaluate_polynomial(a_coefficients, parameter),
            evaluate_polynomial(b_coefficients, parameter),
        )
        if actual != expected:
            raise AssertionError("Jacobian coefficient interpolation failed")
    return a_coefficients, b_coefficients


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def load_families(max200: dict[str, Any]) -> tuple[FamilySpec, ...]:
    records = max200["leader_followup"]["records"]
    rank13 = [
        record for record in records
        if int(record["point_triage"].get("stable_numerical_rank", -1)) == 13
    ]
    rank14 = [
        record for record in records
        if int(record["point_triage"].get("stable_numerical_rank", -1)) == 14
    ]
    observed13 = tuple((tuple(record["roots"]), int(record["parameter"])) for record in rank13)
    observed14 = tuple((tuple(record["roots"]), int(record["parameter"])) for record in rank14)
    if observed13 != EXPECTED_RANK13_PANEL or observed14 != EXPECTED_RANK14_PANEL:
        raise AssertionError("the frozen rank-13/rank-14 frontier changed")
    prior_roots = {roots for roots, _, _ in PRIOR_RATIONAL_SCANS}
    if prior_roots & {roots for roots, _ in EXPECTED_RANK13_PANEL}:
        raise AssertionError("an included rank-13 family already had a rational scan")
    families = []
    for index, record in enumerate(rank13):
        roots = tuple(int(value) for value in record["roots"])
        construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
        if construction.quartic_condition or construction.is_reflection_symmetric:
            raise AssertionError("an included family became degenerate or reflection symmetric")
        a_coefficients, b_coefficients = derive_even_coefficients(construction)
        families.append(
            FamilySpec(
                index=index,
                label="r" + "_".join(map(str, roots)),
                roots=roots,
                calibration=Q(record["parameter"]),
                a_coefficients=a_coefficients,
                b_coefficients=b_coefficients,
                source_record=record,
            )
        )
    return tuple(families)


def family_coefficients(spec: FamilySpec, parameter: Fraction) -> tuple[Fraction, ...]:
    return (
        Q(0), Q(0), Q(0),
        evaluate_polynomial(spec.a_coefficients, parameter),
        evaluate_polynomial(spec.b_coefficients, parameter),
    )


def exact_local_trace_projective(
    spec: FamilySpec, numerator: int, denominator: int, prime: int
) -> int | None:
    def homogeneous(coefficients: Sequence[int]) -> int:
        degree = len(coefficients) - 1
        return sum(
            coefficient * pow(numerator, power, prime)
            * pow(denominator, degree - power, prime)
            for power, coefficient in enumerate(coefficients)
        ) % prime

    coefficient_a = homogeneous(spec.a_coefficients)
    coefficient_b = homogeneous(spec.b_coefficients)
    if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
        return None
    character_sum = 0
    for x_value in range(prime):
        rhs = (x_value**3 + coefficient_a * x_value + coefficient_b) % prime
        if rhs:
            symbol = pow(rhs, (prime - 1) // 2, prime)
            character_sum += 1 if symbol == 1 else -1
    return -character_sum


def llround(value: float) -> int:
    return int(value + 0.5) if value >= 0 else int(value - 0.5)


def score_text(
    spec: FamilySpec, parameter: Fraction, primes: Sequence[int]
) -> tuple[str, int]:
    units = 0
    good = 0
    for prime in primes:
        trace = exact_local_trace_projective(
            spec, parameter.numerator, parameter.denominator, prime
        )
        if trace is None:
            continue
        units += llround(
            ((2 - trace) / (prime + 1 - trace)) * log(float(prime)) * 1.0e12
        )
        good += 1
    sign = "-" if units < 0 else ""
    absolute = abs(units)
    return (
        f"{sign}{absolute // 1_000_000_000_000}."
        f"{absolute % 1_000_000_000_000:012d}",
        good,
    )


def parse_scanner_output(stdout: str) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_RANK13_MULTIFAMILY_SCAN_V1":
        raise AssertionError("the multifamily scanner omitted its header")
    family_index = int(lines[1].split()[1])
    if tuple(map(int, lines[2].split()[1:])) != DISCOVERY_PRIMES:
        raise AssertionError("the discovery band changed")
    if tuple(map(int, lines[3].split()[1:])) != HELD_PRIMES:
        raise AssertionError("the held-out band changed")
    digests = lines[4].split()
    calibration_line = lines[5].split()
    calibration = ScannerCandidate(
        int(calibration_line[1]), int(calibration_line[2]),
        calibration_line[3], calibration_line[4],
        int(calibration_line[5]), int(calibration_line[6]),
    )
    candidates = []
    for line in lines[6:-1]:
        fields = line.split()
        if fields[0] != "C" or len(fields) != 7:
            raise AssertionError("malformed multifamily scanner candidate")
        candidates.append(
            ScannerCandidate(
                int(fields[1]), int(fields[2]), fields[3], fields[4],
                int(fields[5]), int(fields[6]),
            )
        )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 8:
        raise AssertionError("malformed multifamily scanner summary")
    result = ScannerResult(
        family_index=family_index,
        numerator_bound=int(summary[1]),
        denominator_bound=int(summary[2]),
        keep=int(summary[3]),
        primitive_population=int(summary[4]),
        panel_excluded=int(summary[5]),
        evaluated_population=int(summary[6]),
        discovery_table_digest=digests[1],
        held_table_digest=digests[2],
        calibration=calibration,
        candidates=tuple(candidates),
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )
    if len(candidates) != int(summary[7]) or len(candidates) != result.keep:
        raise AssertionError("the retained-population gate changed")
    if len({candidate.parameter for candidate in candidates}) != len(candidates):
        raise AssertionError("the scanner emitted duplicate reduced parameters")
    return result


def coefficient_manifest(spec: FamilySpec) -> str:
    return "\n".join(
        ["MESTRE_EVEN_COEFFICIENTS_V1"]
        + [str(value) for value in spec.a_coefficients]
        + [str(value) for value in spec.b_coefficients]
    ) + "\n"


def run_scanners(
    families: Sequence[FamilySpec], source: Path, *, compiler: str,
    compile_timeout: float, scan_timeout: float, numerator_bound: int,
    denominator_bound: int, keep: int,
) -> tuple[ScannerResult, ...]:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError("a C++17 compiler is required")
    results = []
    with tempfile.TemporaryDirectory(prefix="mestre-rank13-multifamily-") as directory:
        temporary = Path(directory)
        binary = temporary / "scan"
        run_capped_process(
            (executable, "-std=c++17", "-O3", "-DNDEBUG", str(source), "-o", str(binary)),
            timeout=compile_timeout,
        )
        for spec in families:
            manifest = temporary / f"family-{spec.index}.txt"
            manifest.write_text(coefficient_manifest(spec))
            stdout, _ = run_capped_process(
                (
                    str(binary), str(manifest), str(spec.index),
                    str(spec.calibration.numerator), str(spec.calibration.denominator),
                    str(numerator_bound), str(denominator_bound), str(keep),
                ),
                timeout=scan_timeout,
            )
            result = parse_scanner_output(stdout)
            if result.family_index != spec.index:
                raise AssertionError("scanner family index changed")
            expected_discovery = score_text(spec, spec.calibration, DISCOVERY_PRIMES)
            expected_held = score_text(spec, spec.calibration, HELD_PRIMES)
            if (
                (result.calibration.discovery_score, result.calibration.discovery_good)
                != expected_discovery
                or (result.calibration.held_score, result.calibration.held_good)
                != expected_held
            ):
                raise AssertionError("exact Python/C++ calibration replay failed")
            results.append(result)
    return tuple(results)


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


def polynomial_content(coefficients: Sequence[Fraction]) -> int:
    nonzero = [abs(value.numerator) for value in coefficients if value]
    if any(value.denominator != 1 for value in coefficients) or not nonzero:
        raise AssertionError("the primitive discriminant polynomial changed")
    return gcd(*nonzero)


def homogeneous_value(coefficients: Sequence[int], numerator: int, denominator: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


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
        "trial_division_prime_bound": TRIAL_DIVISION_LIMIT,
        "small_prime_valuations": valuations,
        "residual_cofactor": str(remaining),
        "residual_cofactor_bit_length": remaining.bit_length(),
        "known_discriminant_radical": str(known_radical),
        "known_powerful_part": str(known_powerful),
        "denominator_radical": denominator_radical,
        "combined_radical_upper_bound": str(radical_upper),
        "upper_bound_semantics": (
            "rad(residual)<=residual exactly; valuations p<=997 and rad(b) are exact"
        ),
    }


def common_parameter_manifest(
    numerator_bound: int, denominator_bound: int
) -> tuple[int, int, int, str]:
    primitive = excluded = evaluated = 0
    digest = hashlib.sha256()
    for denominator in range(1, denominator_bound + 1):
        for numerator in range(1, numerator_bound + 1):
            if gcd(numerator, denominator) != 1:
                continue
            primitive += 1
            if denominator == 1 and numerator <= 8:
                excluded += 1
                continue
            evaluated += 1
            digest.update(f"{numerator}/{denominator}\n".encode())
    return primitive, excluded, evaluated, digest.hexdigest()


def build_feature_pools(
    families: Sequence[FamilySpec], scans: Sequence[ScannerResult],
    discriminants: dict[int, tuple[int, ...]],
) -> tuple[dict[int, list[dict[str, Any]]], dict[str, Any]]:
    pools: dict[int, list[dict[str, Any]]] = {}
    audits = {}
    for spec, scan in zip(families, scans):
        pool = []
        digest = hashlib.sha256()
        singular = 0
        for candidate in scan.candidates:
            if candidate.parameter in PRIOR_PANEL_PARAMETERS:
                raise AssertionError("an excluded fixed panel fiber leaked")
            feature = discriminant_feature(
                discriminants[spec.index], candidate.numerator, candidate.denominator
            )
            digest.update(
                (
                    f"{spec.label}|{candidate.numerator}/{candidate.denominator}|"
                    f"{candidate.discovery_score}|{candidate.held_score}|"
                    f"{feature['absolute_homogeneous_discriminant']}|"
                    f"{feature.get('combined_radical_upper_bound')}|"
                    f"{feature.get('known_powerful_part')}\n"
                ).encode()
            )
            if feature["singular"]:
                singular += 1
                continue
            pool.append(
                {
                    "family_index": spec.index,
                    "family_label": spec.label,
                    "roots": list(spec.roots),
                    "numerator": candidate.numerator,
                    "denominator": candidate.denominator,
                    "parameter": str(candidate.parameter),
                    "discovery_score": candidate.discovery_score,
                    "held_score": candidate.held_score,
                    "discovery_good": candidate.discovery_good,
                    "held_good": candidate.held_good,
                    "discriminant_feature": feature,
                }
            )
        pools[spec.index] = pool
        audits[spec.label] = {
            "retained_before_singular_rejection": len(scan.candidates),
            "exact_singular_rejections": singular,
            "admissible_feature_pool_count": len(pool),
            "exact_feature_population_sha256": digest.hexdigest(),
        }
    return pools, audits


def select_conductor_population(
    families: Sequence[FamilySpec], pools: dict[int, list[dict[str, Any]]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    selected = []
    digest = hashlib.sha256()
    selected_per_family = {}
    for spec in families:
        pool = pools[spec.index]
        reasons: dict[str, set[str]] = {}

        def take(label: str, order: Sequence[dict[str, Any]], count: int = 1) -> None:
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
        take("highest-held-score", held_order)
        take("smallest-radical-upper-bound", radical_order)
        take("largest-known-powerful-part", powerful_order)
        take("lowest-projective-height", height_order)
        for lower, upper in ((1, 64), (65, 128), (129, 192), (193, 256)):
            bucket = [row for row in held_order if lower <= row["denominator"] <= upper]
            if bucket:
                take(f"held-score-denominator-{lower}-{upper}", bucket)

        by_parameter = {record["parameter"]: record for record in pool}
        family_selected = []
        for parameter, labels in reasons.items():
            record = dict(by_parameter[parameter])
            record["conductor_selection_strata"] = sorted(labels)
            family_selected.append(record)
        family_selected.sort(key=lambda row: (row["denominator"], row["numerator"]))
        selected_per_family[spec.label] = len(family_selected)
        selected.extend(family_selected)
        for record in family_selected:
            digest.update(
                (
                    f"{record['family_label']}|{record['parameter']}|"
                    f"{','.join(record['conductor_selection_strata'])}\n"
                ).encode()
            )
    return selected, {
        "selection_uses_point_or_numerical_rank_data": False,
        "selection_uses_conductor_data": False,
        "equal_per_family_rule": True,
        "discovery_population_closed_before_held_scores": True,
        "quota_rule": (
            "one each held, radical, powerful, height, plus held leader in each "
            "of four denominator quartiles; overlaps retained once"
        ),
        "selected_per_family": selected_per_family,
        "selected_population": len(selected),
        "selected_population_sha256": digest.hexdigest(),
    }


def conductor_worker(
    spec: FamilySpec, numerator: int, denominator: int, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    coefficients = family_coefficients(spec, Q(numerator, denominator))
    try:
        result = capped_minimal_curve_data(
            coefficients, timeout=timeout, stack_bytes=stack_bytes
        )
        logarithm = Decimal(result["log_conductor"])
        return {
            "status": "completed exact PARI minimal-model/conductor computation",
            **result,
            "below_strict_log_conductor_target_numerically": logarithm < TARGET_LOG_CONDUCTOR,
            "below_plausible_log_conductor_190_numerically": logarithm < PLAUSIBLE_LOG_CONDUCTOR,
        }
    except CappedProcessTimeout:
        return {"status": "timeout", "timeout_seconds": timeout, "retried": False}
    except Exception as error:
        return {"status": "error", "error": str(error)[:1000], "retried": False}


def point_stage_worker(
    spec: FamilySpec, numerator: int, denominator: int, height_bound: int,
    point_timeout: float, height_timeout: float, ellrank_timeout: float,
    stack_bytes: int, mapping_cap: int, certificate_prime_bound: int,
) -> dict[str, Any]:
    construction = SixRootMestreConstruction(tuple(Q(root) for root in spec.roots))
    parameter = Q(numerator, denominator)
    coefficients = family_coefficients(spec, parameter)
    try:
        stage, subset = exact_point_stage(
            construction,
            parameter,
            coefficients,
            height_bound=height_bound,
            point_timeout=point_timeout,
            height_timeout=height_timeout,
            ellrank_timeout=ellrank_timeout,
            stack_bytes=stack_bytes,
            mapping_cap=mapping_cap,
        )
        rank = int(stage["stable_numerical_rank"])
        if rank >= FINITE_REDUCTION_TRIGGER:
            stage["finite_reduction_attempt"] = mod3_independence_certificate(
                coefficients, subset, prime_bound=certificate_prime_bound
            )
        else:
            stage["finite_reduction_attempt"] = {
                "status": "not triggered",
                "trigger_stable_numerical_rank": FINITE_REDUCTION_TRIGGER,
            }
        return stage
    except CappedProcessTimeout:
        return {
            "status": "timeout", "timeout_seconds": point_timeout,
            "same_height_retry": False,
        }
    except Exception as error:
        return {
            "status": "error", "error": str(error)[:1000],
            "same_height_retry": False,
        }


def stage_rank_key(record: dict[str, Any], stage_name: str) -> tuple[Any, ...]:
    stage = record["point_stages"][stage_name]
    return (
        -int(stage["stable_numerical_rank"]),
        -Decimal(record["held_score"]),
        int(record["discriminant_feature"]["combined_radical_upper_bound"]),
        record["denominator"], record["numerator"],
    )


def calibration_audit(spec: FamilySpec) -> dict[str, Any]:
    source = spec.source_record
    coefficients = family_coefficients(spec, spec.calibration)
    subset = tuple(
        (Q(point["x"]), Q(point["y"]))
        for point in source["point_triage"]["numerical_subset"]
    )
    if len(subset) != 13 or any(
        not point_on_short_curve(coefficients, point) for point in subset
    ):
        raise AssertionError("a frozen rank-13 calibration point set changed")
    return {
        "family_index": spec.index,
        "family_label": spec.label,
        "roots": list(spec.roots),
        "parameter": str(spec.calibration),
        "source_stable_numerical_rank": 13,
        "source_height_bound": source["point_triage"]["height_bound"],
        "point_count": len(subset),
        "point_sha256": point_digest(subset),
        "numerical_rank_is_not_an_independence_certificate": True,
        "excluded_from_common_rational_box": True,
    }


def scanner_record(spec: FamilySpec, scan: ScannerResult) -> dict[str, Any]:
    return {
        "family_index": spec.index,
        "family_label": spec.label,
        "numerator_bound": scan.numerator_bound,
        "denominator_bound": scan.denominator_bound,
        "keep": scan.keep,
        "primitive_population": scan.primitive_population,
        "fixed_panel_excluded": scan.panel_excluded,
        "evaluated_population": scan.evaluated_population,
        "discovery_table_digest_fnv64": scan.discovery_table_digest,
        "held_table_digest_fnv64": scan.held_table_digest,
        "stdout_sha256": scan.stdout_sha256,
        "retained_candidate_sha256": stable_json_digest(
            [
                [
                    candidate.numerator, candidate.denominator,
                    candidate.discovery_score, candidate.held_score,
                    candidate.discovery_good, candidate.held_good,
                ]
                for candidate in scan.candidates
            ]
        ),
        "calibration": {
            "parameter": str(scan.calibration.parameter),
            "discovery_score": scan.calibration.discovery_score,
            "held_score": scan.calibration.held_score,
            "discovery_good": scan.calibration.discovery_good,
            "held_good": scan.calibration.held_good,
            "exact_python_replay": True,
            "excluded": True,
        },
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
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--numerator-bound", type=int, default=4096)
    parser.add_argument("--denominator-bound", type=int, default=256)
    parser.add_argument("--keep-per-family", type=int, default=512)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=10.0)
    parser.add_argument("--h5000-timeout", type=float, default=15.0)
    parser.add_argument("--h50000-timeout", type=float, default=20.0)
    parser.add_argument("--h250000-timeout", type=float, default=30.0)
    parser.add_argument("--h1000000-timeout", type=float, default=45.0)
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--ellrank-timeout", type=float, default=12.0)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=499)
    parser.add_argument("--stack-bytes", type=int, default=STACK_BYTES)
    parser.add_argument(
        "--output", type=Path,
        default=root / "artifacts/generated-results" / DEFAULT_OUTPUT.name,
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if (
        args.numerator_bound != 4096 or args.denominator_bound != 256
        or args.keep_per_family != 512
    ):
        raise SystemExit("the common rational box and survivor count are pinned")
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.compile_timeout, args.scan_timeout, args.conductor_timeout,
        args.h5000_timeout, args.h50000_timeout, args.h250000_timeout,
        args.h1000000_timeout, args.height_timeout, args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")
    if not 32 <= args.mapping_cap <= 1024:
        raise SystemExit("mapping cap must lie in [32,1024]")
    if not 211 <= args.certificate_prime_bound <= 2000:
        raise SystemExit("certificate prime bound must lie in [211,2000]")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rank13-multifamily artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    scanner_path = script_path.with_name("scan_mestre_rank13_multifamily.cpp")
    artifact_dir = root / "artifacts/generated-results"
    max200_path = artifact_dir / FROZEN_MAX200_ARTIFACT
    max200_script = script_path.with_name("search_mestre_root_tuple_scale_max200.py")
    if (
        sha256_file(max200_path) != EXPECTED_MAX200_ARTIFACT_SHA256
        or sha256_file(max200_script) != EXPECTED_MAX200_SCRIPT_SHA256
    ):
        raise AssertionError("the frozen max-root-200 inputs changed")
    max200 = json.loads(max200_path.read_text())
    if max200["result_sha256"] != EXPECTED_MAX200_RESULT_SHA256:
        raise AssertionError("the frozen max-root-200 result digest changed")
    for _, artifact_name, expected_sha in PRIOR_RATIONAL_SCANS:
        if sha256_file(artifact_dir / artifact_name) != expected_sha:
            raise AssertionError(f"prior rational-scan artifact changed: {artifact_name}")

    families = load_families(max200)
    calibrations = [calibration_audit(spec) for spec in families]
    discriminants: dict[int, tuple[int, ...]] = {}
    family_records = []
    for spec in families:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in spec.roots))
        raw_discriminant = construction.primitive_discriminant_polynomial
        content = polynomial_content(raw_discriminant)
        normalized = tuple(value.numerator // content for value in raw_discriminant)
        if len(normalized) != 21:
            raise AssertionError("an included discriminant ceased to have degree 20")
        discriminants[spec.index] = normalized
        family_records.append(
            {
                "family_index": spec.index,
                "family_label": spec.label,
                "roots": list(spec.roots),
                "calibration_parameter": str(spec.calibration),
                "A_coefficients_ascending": list(spec.a_coefficients),
                "B_coefficients_ascending": list(spec.b_coefficients),
                "content_free_discriminant_coefficients_ascending": list(normalized),
                "removed_discriminant_polynomial_content": str(content),
                "exact_symmetry": "A(-T)=A(T), B(-T)=B(T), R_-T=R_T",
            }
        )

    primitive, excluded, evaluated, parameter_digest = common_parameter_manifest(
        args.numerator_bound, args.denominator_bound
    )
    scans = run_scanners(
        families,
        scanner_path,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
        keep=args.keep_per_family,
    )
    if any(
        (scan.primitive_population, scan.panel_excluded, scan.evaluated_population)
        != (primitive, excluded, evaluated)
        for scan in scans
    ):
        raise AssertionError("the exact common population replay failed")
    print(
        f"exact scans closed: families={len(families)} common_population={evaluated}",
        flush=True,
    )

    pools, pool_audits = build_feature_pools(families, scans, discriminants)
    selected, selection = select_conductor_population(families, pools)
    print(
        f"rank-blind diversity selection closed: conductors={len(selected)}",
        flush=True,
    )

    excluded_family_lines = "".join(
        ",".join(map(str, roots)) + "\n" for roots, _, _ in PRIOR_RATIONAL_SCANS
    )
    included_family_lines = "".join(
        ",".join(map(str, spec.roots)) + f"|{spec.calibration}\n" for spec in families
    )
    common: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "completed selection-only rank13 multifamily rational scan"
            if args.selection_only else "in-progress conductor-first rank13 multifamily rational scan"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "scope": {
            "included_families": [list(spec.roots) for spec in families],
            "included_family_anchor_manifest_sha256": hashlib.sha256(
                included_family_lines.encode()
            ).hexdigest(),
            "excluded_rank14_panel_families": [list(roots) for roots, _ in EXPECTED_RANK14_PANEL],
            "excluded_prior_rational_scan_families": [
                {
                    "roots": list(roots),
                    "artifact": artifact_name,
                    "artifact_sha256": artifact_sha,
                }
                for roots, artifact_name, artifact_sha in PRIOR_RATIONAL_SCANS
            ],
            "excluded_family_manifest_sha256": hashlib.sha256(
                excluded_family_lines.encode()
            ).hexdigest(),
            "included_and_prior_rational_families_disjoint": True,
            "T_sign_quotient": "primitive positive T=a/b only; T and -T are identical",
            "fixed_panel_parameters_excluded": [str(value) for value in PRIOR_PANEL_PARAMETERS],
            "fixed_fiber_policy": (
                "no T in the frozen integer panel 1..8 is searched, including any "
                "fiber conditionally closed by explicit-formula analysis"
            ),
        },
        "families": family_records,
        "frozen_calibrations": calibrations,
        "modular_scan": {
            "score": (
                "sum ((2-a_p)/(p+1-a_p))*log(p); exact traces, each term "
                "quantized to 1e-12"
            ),
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_relative_to_all_prior_rational_scanners_through_prime_809": (
                min(DISCOVERY_PRIMES + HELD_PRIMES) > 809
            ),
            "common_box": {
                "numerator": [1, args.numerator_bound],
                "denominator": [1, args.denominator_bound],
                "primitive_positive_rationals": primitive,
                "fixed_panel_excluded": excluded,
                "evaluated_per_family": evaluated,
                "family_count": len(families),
                "evaluated_family_parameter_pairs": evaluated * len(families),
                "ordered_parameter_manifest_sha256": parameter_digest,
            },
            "family_scans": [
                scanner_record(spec, scan) for spec, scan in zip(families, scans)
            ],
        },
        "exact_discriminant_feature_screen": {
            "trial_division_prime_bound": TRIAL_DIVISION_LIMIT,
            "content_free_homogeneous_degree": 20,
            "pool_audits": pool_audits,
        },
        "conductor_selection": selection,
        "selected_records": selected,
        "parameters": {
            key: value for key, value in vars(args).items()
            if key not in {"output", "selection_only"}
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "scanner_path": str(scanner_path.relative_to(root)),
            "scanner_sha256": sha256_file(scanner_path),
            "frozen_max200_artifact": str(max200_path.relative_to(root)),
            "frozen_max200_artifact_sha256": EXPECTED_MAX200_ARTIFACT_SHA256,
            "frozen_max200_result_sha256": EXPECTED_MAX200_RESULT_SHA256,
            "frozen_max200_script_sha256": EXPECTED_MAX200_SCRIPT_SHA256,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "temporary_scanner_binary_and_manifests_removed": True,
            "external_calls_use_foreground_process_groups": True,
            "same_stage_retries": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "compiler": shutil.which(args.compiler),
        },
        "timings": {"pre_conductor_wall_seconds": time.monotonic() - started},
    }
    if args.selection_only:
        common["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
        common["result_sha256"] = stable_json_digest(
            {
                "scope": common["scope"], "families": common["families"],
                "scan": common["modular_scan"],
                "features": common["exact_discriminant_feature_screen"],
                "selection": common["conductor_selection"],
            }
        )
        exclusive_write(args.output, common)
        return

    conductor_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                conductor_worker,
                families[record["family_index"]], record["numerator"],
                record["denominator"], args.conductor_timeout, args.stack_bytes,
            )
            for record in selected
        ]
        for position, (record, future) in enumerate(zip(selected, futures), start=1):
            record["conductor_phase"] = future.result()
            if position % 20 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    eligible = [
        record for record in selected
        if record["conductor_phase"]["status"].startswith("completed")
    ]
    print(f"conductor-first boundary closed: completed={len(eligible)}", flush=True)

    stages = (
        ("H5000", 5_000, None, args.h5000_timeout),
        ("H50000", 50_000, 2, args.h50000_timeout),
        ("H250000", 250_000, 1, args.h250000_timeout),
        ("H1000000", 1_000_000, 1, args.h1000000_timeout),
    )
    current = eligible
    for stage_index, (name, height_bound, keep_per_family, timeout) in enumerate(stages):
        if stage_index:
            prior_name = stages[stage_index - 1][0]
            next_current = []
            for spec in families:
                completed = [
                    record for record in current
                    if record["family_index"] == spec.index
                    and record.get("point_stages", {}).get(prior_name, {}).get("status")
                    == "completed"
                ]
                completed.sort(key=lambda record: stage_rank_key(record, prior_name))
                next_current.extend(completed[: int(keep_per_family)])
            current = next_current
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    point_stage_worker,
                    families[record["family_index"]], record["numerator"],
                    record["denominator"], height_bound, timeout,
                    args.height_timeout, args.ellrank_timeout, args.stack_bytes,
                    args.mapping_cap, args.certificate_prime_bound,
                )
                for record in current
            ]
            for record, future in zip(current, futures):
                stage = future.result()
                record.setdefault("point_stages", {})[name] = stage
                rank = int(stage.get("stable_numerical_rank", -1))
                if rank >= FINITE_REDUCTION_TRIGGER:
                    certificate = stage.get("finite_reduction_attempt", {})
                    print(
                        f"EARLY_SIGNAL family={record['family_label']} "
                        f"T={record['parameter']} stage={name} rank={rank} "
                        f"certified={certificate.get('certified_algebraic_rank_lower_bound')}",
                        flush=True,
                    )
        maximum = max(
            (
                int(record["point_stages"][name].get("stable_numerical_rank", -1))
                for record in current
            ),
            default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={maximum}", flush=True)

    completed_stages = [
        (record, stage_name, stage)
        for record in selected
        for stage_name, stage in record.get("point_stages", {}).items()
        if stage.get("status") == "completed"
    ]
    maximum_rank = max(
        (int(stage["stable_numerical_rank"]) for _, _, stage in completed_stages),
        default=None,
    )
    finite_attempts = []
    target_hits = []
    for record, stage_name, stage in completed_stages:
        certificate = stage.get("finite_reduction_attempt", {})
        certified = certificate.get("certified_algebraic_rank_lower_bound")
        if certified is None:
            continue
        finite_attempts.append(
            {
                "family_label": record["family_label"],
                "parameter": record["parameter"],
                "stage": stage_name,
                "certified_rank_lower_bound": certified,
            }
        )
        below = record["conductor_phase"].get(
            "below_strict_log_conductor_target_numerically", False
        )
        if certified >= 30 or (certified >= 21 and below):
            target_hits.append(
                {
                    "family_label": record["family_label"],
                    "parameter": record["parameter"],
                    "stage": stage_name,
                    "certified_rank_lower_bound": certified,
                    "conductor": record["conductor_phase"]["conductor"],
                    "log_conductor": record["conductor_phase"]["log_conductor"],
                }
            )

    common["status"] = (
        "completed fixed rank13 multifamily rational scan; stopped without broadening"
    )
    common["target"]["hits"] = target_hits
    common["conductor_first_screen"] = {
        "population_closed_before_any_point_or_rank_call": True,
        "selected_population": len(selected),
        "completed": len(eligible),
        "timeouts": sum(
            record["conductor_phase"]["status"] == "timeout" for record in selected
        ),
        "errors": sum(
            record["conductor_phase"]["status"] == "error" for record in selected
        ),
        "strict_subtarget": sum(
            record["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically", False
            )
            for record in selected
        ),
        "plausible_below_190": sum(
            record["conductor_phase"].get(
                "below_plausible_log_conductor_190_numerically", False
            )
            for record in selected
        ),
        "all_selected_received_exact_conductor_attempt": True,
    }
    common["point_search_protocol"] = {
        "stages": [
            {
                "name": name,
                "height_bound": height_bound,
                "keep_per_family_after_previous_stage": keep,
                "attempted": sum(
                    name in record.get("point_stages", {}) for record in selected
                ),
            }
            for name, height_bound, keep, _ in stages
        ],
        "increasing_heights_are_not_retries": True,
        "same_height_retries": 0,
        "finite_reduction_trigger_stable_rank": FINITE_REDUCTION_TRIGGER,
        "finite_reduction_attempts": finite_attempts,
        "maximum_stable_numerical_rank": maximum_rank,
        "completed_stage_calls": len(completed_stages),
        "stop_rule": "stop after the fixed H1m per-family diversity tier",
        "stop_rule_fired": True,
        "broadening_calls_after_fixed_protocol": 0,
    }
    common["timings"].update(
        {
            "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
            "total_wall_seconds": time.monotonic() - started,
        }
    )
    common["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    common["provenance"]["owned_processes_remaining"] = 0
    common["result_sha256"] = stable_json_digest(
        {
            "scope": common["scope"], "families": common["families"],
            "calibrations": common["frozen_calibrations"],
            "scan": common["modular_scan"],
            "features": common["exact_discriminant_feature_screen"],
            "selection": common["conductor_selection"],
            "records": selected,
            "conductor": common["conductor_first_screen"],
            "points": common["point_search_protocol"],
            "target": common["target"],
        }
    )
    exclusive_write(args.output, common)
    print(
        f"complete max_rank={maximum_rank} target_hits={len(target_hits)} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

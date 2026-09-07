#!/usr/bin/env python3
"""Leakage-controlled rational search in two max-root-200 rank-14 families.

This standalone continuation treats exactly the two new families

* roots ``(0,17,142,145,162,200)``, calibration ``T=7``;
* roots ``(0,7,121,128,183,194)``, calibration ``T=1``.

The distinct max-root-200 family ``(0,25,57,104,116,148)`` is deliberately
out of scope.  Both selected families are even in T, so the search exhausts
one representative ``T=a/b>0`` in a primitive 30,000 by 1,000 box.

A compiled scanner uses a fresh discovery-prime band (211..281) and a
disjoint held-out band (283..373).  Discovery scores close three nested
survivor populations before held-out scores are read.  For their union this
driver evaluates the exact content-free degree-20 homogeneous discriminant,
extracts all prime valuations through 997, and records an exact radical upper
bound and exact known powerful part.  A fixed union of held-score, radical,
powerful-part, and low-height quotas closes before any conductor call.

The complete retained population receives conductor computations first.
Only after that boundary do exact quartic searches and numerical height
triage run at H=5,000, 50,000, 250,000, and 1,000,000 in fixed nested tiers.
Increasing heights are not retries.  Stable numerical rank at least 18
immediately receives an exact finite-reduction independence attempt.  If the
fixed protocol never exceeds rank 15, the run stops without broadening.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import ceil, gcd, log
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

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

DISCOVERY_PRIMES = (
    211, 223, 227, 229, 233, 239, 241,
    251, 257, 263, 269, 271, 277, 281,
)
HELD_PRIMES = (
    283, 293, 307, 311, 313, 317, 331,
    337, 347, 349, 353, 359, 367, 373,
)
TRIAL_DIVISION_LIMIT = 997
PRIOR_PANEL_PARAMETERS = tuple(Q(value) for value in range(1, 9))
STACK_BYTES = 512_000_000
FINITE_REDUCTION_TRIGGER = 18
DEFAULT_OUTPUT = Path(
    "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_rank14_pair_rational_frontier.json"
)


@dataclass(frozen=True)
class FamilySpec:
    index: int
    label: str
    roots: tuple[int, ...]
    calibration: Fraction
    a_coefficients: tuple[int, ...]
    b_coefficients: tuple[int, ...]


FAMILIES = (
    FamilySpec(
        0,
        "r0_17_142_145_162_200",
        (0, 17, 142, 145, 162, 200),
        Q(7),
        (
            -314_966_562_629_647_516_875,
            0,
            -6_216_335_576_543_616_300,
            0,
            853_623_006_080_112,
            0,
            10_204_530_336,
            0,
            -7_641_648,
        ),
        (
            -1_574_645_096_287_482_927_205_008_656_250,
            0,
            218_635_455_077_380_331_614_991_902_500,
            0,
            -43_358_567_227_166_411_264_450_700,
            0,
            13_182_804_975_172_483_234_512,
            0,
            -1_378_693_790_037_559_296,
            0,
            -16_286_430_416_256,
            0,
            8_130_713_472,
        ),
    ),
    FamilySpec(
        1,
        "r0_7_121_128_183_194",
        (0, 7, 121, 128, 183, 194),
        Q(1),
        (
            -23_176_885_126_717_023_443_712,
            0,
            -30_716_437_929_772_416_480,
            0,
            5_349_668_724_315_261,
            0,
            -44_248_113_000,
            0,
            -29_428_272,
        ),
        (
            867_772_568_274_887_992_065_902_700_355_584,
            0,
            4_115_037_444_689_925_520_425_570_654_720,
            0,
            -945_618_543_766_059_469_223_433_696,
            0,
            142_515_300_626_499_263_347_710,
            0,
            -16_911_444_014_164_334_952,
            0,
            138_585_089_916_000,
            0,
            61_446_231_936,
        ),
    ),
)


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
    stratum: str
    numerator_bound: int
    denominator_bound: int
    keep: int
    primitive_population: int
    prior_excluded: int
    evaluated_population: int
    discovery_table_digest: str
    held_table_digest: str
    calibration: ScannerCandidate
    candidates: tuple[ScannerCandidate, ...]
    stdout_sha256: str


def family_spec(index: int) -> FamilySpec:
    return FAMILIES[index]


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + coefficient
    return answer


def family_coefficients(index: int, parameter: Fraction) -> tuple[Fraction, ...]:
    spec = family_spec(index)
    parameter = Q(parameter)
    return (
        Q(0), Q(0), Q(0),
        evaluate_polynomial(spec.a_coefficients, parameter),
        evaluate_polynomial(spec.b_coefficients, parameter),
    )


def exact_local_trace_projective(
    index: int, numerator: int, denominator: int, prime: int
) -> int | None:
    spec = family_spec(index)

    def homogeneous(coefficients: Sequence[int]) -> int:
        degree = len(coefficients) - 1
        return sum(
            coefficient
            * pow(numerator, power, prime)
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


def exact_table_digest(index: int, primes: Sequence[int]) -> str:
    digest = 1_469_598_103_934_665_603

    def mix(value: int) -> None:
        nonlocal digest
        value &= (1 << 64) - 1
        for offset in range(8):
            digest ^= (value >> (8 * offset)) & 255
            digest = (digest * 1_099_511_628_211) & ((1 << 64) - 1)

    for prime in primes:
        mix(prime)
        for numerator, denominator in (
            *((residue, 1) for residue in range(prime)),
            (1, 0),
        ):
            trace = exact_local_trace_projective(
                index, numerator, denominator, prime
            )
            mix(trace is not None)
            mix(0 if trace is None else trace)
    return str(digest)


def llround(value: float) -> int:
    return int(value + 0.5) if value >= 0 else int(value - 0.5)


def score_text(index: int, parameter: Fraction, primes: Sequence[int]) -> tuple[str, int]:
    units = 0
    good = 0
    for prime in primes:
        trace = exact_local_trace_projective(
            index, parameter.numerator, parameter.denominator, prime
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


def parse_scanner_output(stratum: str, stdout: str) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_RANK14_PAIR_SCAN_V1":
        raise AssertionError("the rank14-pair scanner omitted its header")
    family_index = int(lines[1].split()[1])
    discovery = tuple(int(value) for value in lines[2].split()[1:])
    held = tuple(int(value) for value in lines[3].split()[1:])
    if discovery != DISCOVERY_PRIMES or held != HELD_PRIMES:
        raise AssertionError("the fresh scanner prime bands changed")
    local = lines[4].split()
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
            raise AssertionError("malformed rank14-pair candidate")
        candidates.append(
            ScannerCandidate(
                int(fields[1]), int(fields[2]), fields[3], fields[4],
                int(fields[5]), int(fields[6]),
            )
        )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 8:
        raise AssertionError("malformed rank14-pair summary")
    result = ScannerResult(
        family_index=family_index,
        stratum=stratum,
        numerator_bound=int(summary[1]),
        denominator_bound=int(summary[2]),
        keep=int(summary[3]),
        primitive_population=int(summary[4]),
        prior_excluded=int(summary[5]),
        evaluated_population=int(summary[6]),
        discovery_table_digest=local[1],
        held_table_digest=local[2],
        calibration=calibration,
        candidates=tuple(candidates),
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )
    if len(candidates) != int(summary[7]) or len(candidates) != result.keep:
        raise AssertionError("the scanner retained-count gate changed")
    if len({candidate.parameter for candidate in candidates}) != len(candidates):
        raise AssertionError("the scanner emitted a duplicate rational parameter")
    return result


def run_scanners(
    source: Path,
    *,
    compiler: str,
    compile_timeout: float,
    scan_timeout: float,
    strata: Sequence[tuple[str, int, int]],
    denominator_bound: int,
) -> tuple[ScannerResult, ...]:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError("a C++17 compiler is required")
    results = []
    with tempfile.TemporaryDirectory(prefix="mestre-rank14-pair-") as directory:
        binary = Path(directory) / "scan"
        run_capped_process(
            (
                executable, "-std=c++17", "-O3", "-DNDEBUG",
                str(source), "-o", str(binary),
            ),
            timeout=compile_timeout,
        )
        for spec in FAMILIES:
            for name, numerator_bound, keep in strata:
                stdout, _ = run_capped_process(
                    (
                        str(binary), str(spec.index), str(numerator_bound),
                        str(denominator_bound), str(keep),
                    ),
                    timeout=scan_timeout,
                )
                results.append(parse_scanner_output(name, stdout))
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
    if any(value.denominator != 1 for value in coefficients):
        raise AssertionError("the primitive discriminant polynomial became rational")
    return gcd(*(abs(value.numerator) for value in coefficients if value))


def homogeneous_value(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def complete_radical(value: int) -> int:
    value = abs(value)
    radical = 1
    remaining = value
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
            "rad(residual)<=residual exactly; all valuations p<=997 and rad(b) "
            "are exact, so this is a rigorous radical upper bound"
        ),
    }


def pool_and_features(
    spec: FamilySpec,
    scans: Sequence[ScannerResult],
    discriminant_coefficients: Sequence[int],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_parameter: dict[Fraction, dict[str, Any]] = {}
    for scan in scans:
        if scan.family_index != spec.index:
            continue
        for candidate in scan.candidates:
            if candidate.parameter in PRIOR_PANEL_PARAMETERS:
                raise AssertionError("a prior panel parameter escaped scanner exclusion")
            record = by_parameter.setdefault(
                candidate.parameter,
                {
                    "family_index": spec.index,
                    "family_label": spec.label,
                    "numerator": candidate.numerator,
                    "denominator": candidate.denominator,
                    "parameter": str(candidate.parameter),
                    "discovery_score": candidate.discovery_score,
                    "held_score": candidate.held_score,
                    "discovery_good": candidate.discovery_good,
                    "held_good": candidate.held_good,
                    "discovery_survivor_strata": set(),
                },
            )
            if (
                record["discovery_score"] != candidate.discovery_score
                or record["held_score"] != candidate.held_score
            ):
                raise AssertionError("one parameter acquired inconsistent local scores")
            record["discovery_survivor_strata"].add(scan.stratum)
    pool = []
    singular = 0
    digest = hashlib.sha256()
    for parameter, record in sorted(by_parameter.items()):
        feature = discriminant_feature(
            discriminant_coefficients, parameter.numerator, parameter.denominator
        )
        record["discovery_survivor_strata"] = sorted(
            record["discovery_survivor_strata"]
        )
        record["discriminant_feature"] = feature
        digest.update(
            (
                f"{spec.label}|{parameter}|{record['discovery_score']}|"
                f"{record['held_score']}|"
                f"{','.join(record['discovery_survivor_strata'])}|"
                f"{feature['absolute_homogeneous_discriminant']}|"
                f"{feature.get('combined_radical_upper_bound')}|"
                f"{feature.get('known_powerful_part')}\n"
            ).encode()
        )
        if feature["singular"]:
            singular += 1
        else:
            pool.append(record)
    return pool, {
        "discovery_survivor_union_count_before_singular_rejection": len(by_parameter),
        "exact_singular_rejections": singular,
        "admissible_feature_pool_count": len(pool),
        "exact_feature_population_sha256": digest.hexdigest(),
    }


def select_conductor_population(
    pools: dict[int, list[dict[str, Any]]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reasons: dict[tuple[int, str], set[str]] = defaultdict(set)
    by_key = {
        (record["family_index"], record["parameter"]): record
        for pool in pools.values()
        for record in pool
    }
    quotas = {
        "highest-held-score": 32,
        "smallest-exact-radical-upper-bound": 32,
        "largest-exact-known-powerful-part": 16,
        "lowest-projective-parameter-height": 8,
    }
    for spec in FAMILIES:
        pool = pools[spec.index]
        orders = {
            "highest-held-score": sorted(
                pool,
                key=lambda row: (
                    -Decimal(row["held_score"]),
                    -Decimal(row["discovery_score"]),
                    row["denominator"], row["numerator"],
                ),
            ),
            "smallest-exact-radical-upper-bound": sorted(
                pool,
                key=lambda row: (
                    int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                    -Decimal(row["held_score"]), row["denominator"], row["numerator"],
                ),
            ),
            "largest-exact-known-powerful-part": sorted(
                pool,
                key=lambda row: (
                    -int(row["discriminant_feature"]["known_powerful_part"]),
                    int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                    row["denominator"], row["numerator"],
                ),
            ),
            "lowest-projective-parameter-height": sorted(
                pool,
                key=lambda row: (
                    max(row["numerator"], row["denominator"]),
                    row["denominator"], row["numerator"],
                ),
            ),
        }
        for label, order in orders.items():
            for record in order[: quotas[label]]:
                reasons[spec.index, record["parameter"]].add(label)
    selected = []
    for key, strata in reasons.items():
        record = dict(by_key[key])
        record["conductor_selection_strata"] = sorted(strata)
        selected.append(record)
    selected.sort(
        key=lambda row: (row["family_index"], row["denominator"], row["numerator"])
    )
    digest = hashlib.sha256()
    for record in selected:
        digest.update(
            (
                f"{record['family_label']}|{record['parameter']}|"
                f"{','.join(record['conductor_selection_strata'])}\n"
            ).encode()
        )
    return selected, {
        "selection_uses_conductor": False,
        "selection_uses_point_or_rank_data": False,
        "discovery_survivors_closed_before_held_scores": True,
        "held_scores_rank_only_discovery_survivors": True,
        "exact_discriminant_features_use_no_conductor_or_rank_data": True,
        "quota_per_family": quotas,
        "selected_per_family": {
            spec.label: sum(row["family_index"] == spec.index for row in selected)
            for spec in FAMILIES
        },
        "selected_population": len(selected),
        "selected_population_sha256": digest.hexdigest(),
    }


def calibration_records(artifact: dict[str, Any]) -> list[dict[str, Any]]:
    records = []
    for spec in FAMILIES:
        source = next(
            record
            for record in artifact["leader_followup"]["records"]
            if tuple(record["roots"]) == spec.roots
            and Q(record["parameter"]) == spec.calibration
        )
        construction = SixRootMestreConstruction(tuple(Q(root) for root in spec.roots))
        coefficients = construction.primitive_jacobian_coefficients(spec.calibration)
        subset = tuple(
            (Q(point["x"]), Q(point["y"]))
            for point in source["point_triage"]["numerical_subset"]
        )
        if (
            len(subset) != 14
            or any(not point_on_short_curve(coefficients, point) for point in subset)
        ):
            raise AssertionError("a frozen rank14 calibration point set changed")
        replay = mod3_independence_certificate(coefficients, subset, prime_bound=499)
        source_certificate = source["immediate_exact_gain_attempt"]["mod3"]
        if json.loads(json.dumps(replay)) != source_certificate:
            raise AssertionError("a frozen rank14 finite-reduction certificate changed")
        records.append(
            {
                "family_index": spec.index,
                "family_label": spec.label,
                "roots": list(spec.roots),
                "parameter": str(spec.calibration),
                "point_sha256": point_digest(subset),
                "finite_reduction_certificate": replay,
                "certified_algebraic_rank_lower_bound": 14,
                "conductor": source["conductor_phase"],
                "excluded_from_every_scanner_stratum": True,
            }
        )
    return records


def conductor_worker(
    index: int,
    numerator: int,
    denominator: int,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    coefficients = family_coefficients(index, Q(numerator, denominator))
    try:
        result = capped_minimal_curve_data(
            coefficients, timeout=timeout, stack_bytes=stack_bytes
        )
        return {
            "status": "completed exact PARI minimal-model/conductor computation",
            **result,
            "below_strict_log_conductor_target_numerically": (
                Decimal(result["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
        }
    except CappedProcessTimeout:
        return {"status": "timeout", "timeout_seconds": timeout, "retried": False}
    except Exception as error:
        return {"status": "error", "error": str(error)[:1000], "retried": False}


def point_stage_worker(
    index: int,
    numerator: int,
    denominator: int,
    height_bound: int,
    point_timeout: float,
    height_timeout: float,
    ellrank_timeout: float,
    stack_bytes: int,
    mapping_cap: int,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    spec = family_spec(index)
    construction = SixRootMestreConstruction(tuple(Q(root) for root in spec.roots))
    parameter = Q(numerator, denominator)
    coefficients = family_coefficients(index, parameter)
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
        if int(stage["stable_numerical_rank"]) >= FINITE_REDUCTION_TRIGGER:
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
            "status": "timeout",
            "timeout_seconds": point_timeout,
            "same_height_retry": False,
        }
    except Exception as error:
        return {
            "status": "error",
            "error": str(error)[:1000],
            "same_height_retry": False,
        }


def stage_rank_key(record: dict[str, Any], prior_stage: str) -> tuple[Any, ...]:
    stage = record["point_stages"][prior_stage]
    return (
        -int(stage["stable_numerical_rank"]),
        -Decimal(record["held_score"]),
        int(record["discriminant_feature"]["combined_radical_upper_bound"]),
        record["denominator"],
        record["numerator"],
    )


def scanner_record(scan: ScannerResult) -> dict[str, Any]:
    return {
        "family_index": scan.family_index,
        "family_label": family_spec(scan.family_index).label,
        "stratum": scan.stratum,
        "numerator_bound": scan.numerator_bound,
        "denominator_bound": scan.denominator_bound,
        "keep": scan.keep,
        "primitive_population": scan.primitive_population,
        "prior_panel_excluded": scan.prior_excluded,
        "evaluated_population": scan.evaluated_population,
        "discovery_table_digest": scan.discovery_table_digest,
        "held_table_digest": scan.held_table_digest,
        "stdout_sha256": scan.stdout_sha256,
        "calibration": {
            "parameter": str(scan.calibration.parameter),
            "discovery_score": scan.calibration.discovery_score,
            "held_score": scan.calibration.held_score,
            "discovery_good": scan.calibration.discovery_good,
            "held_good": scan.calibration.held_good,
            "excluded": True,
        },
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
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--numerator-bound", type=int, default=30_000)
    parser.add_argument("--denominator-bound", type=int, default=1_000)
    parser.add_argument("--global-keep", type=int, default=4_096)
    parser.add_argument("--medium-numerator-bound", type=int, default=5_000)
    parser.add_argument("--medium-keep", type=int, default=2_048)
    parser.add_argument("--low-numerator-bound", type=int, default=1_000)
    parser.add_argument("--low-keep", type=int, default=1_024)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=12.0)
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
    if args.numerator_bound != 30_000 or args.denominator_bound != 1_000:
        raise SystemExit("the declared broad box is pinned at 30000 by 1000")
    if (
        args.global_keep != 4_096
        or args.medium_numerator_bound != 5_000
        or args.medium_keep != 2_048
        or args.low_numerator_bound != 1_000
        or args.low_keep != 1_024
    ):
        raise SystemExit("the three discovery survivor strata are pinned")
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.compile_timeout, args.scan_timeout, args.conductor_timeout,
        args.h5000_timeout, args.h50000_timeout, args.h250000_timeout,
        args.h1000000_timeout, args.height_timeout, args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")
    if args.mapping_cap < 32 or args.mapping_cap > 1024:
        raise SystemExit("mapping cap must lie in [32,1024]")
    if args.certificate_prime_bound < 211 or args.certificate_prime_bound > 2000:
        raise SystemExit("certificate prime bound must lie in [211,2000]")


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rank14-pair frontier artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    scanner_path = script_path.with_name("scan_mestre_rank14_pair.cpp")
    max200_path = root / "artifacts/generated-results" / FROZEN_MAX200_ARTIFACT
    max200_script = script_path.with_name("search_mestre_root_tuple_scale_max200.py")
    if (
        sha256_file(max200_path) != EXPECTED_MAX200_ARTIFACT_SHA256
        or sha256_file(max200_script) != EXPECTED_MAX200_SCRIPT_SHA256
    ):
        raise AssertionError("the frozen max-root-200 inputs changed")
    max200 = json.loads(max200_path.read_text())
    if max200["result_sha256"] != EXPECTED_MAX200_RESULT_SHA256:
        raise AssertionError("the frozen max-root-200 result digest changed")

    calibrations = calibration_records(max200)
    constructions: dict[int, SixRootMestreConstruction] = {}
    discriminants: dict[int, tuple[int, ...]] = {}
    family_records = []
    for spec in FAMILIES:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in spec.roots))
        constructions[spec.index] = construction
        if construction.quartic_condition or construction.is_reflection_symmetric:
            raise AssertionError("a selected rank14 family geometry changed")
        for parameter in (Q(1), Q(7), Q(11, 5), Q(-17, 3)):
            if construction.primitive_jacobian_coefficients(parameter) != family_coefficients(
                spec.index, parameter
            ):
                raise AssertionError("a pinned rank14 A(T),B(T) formula changed")
            if family_coefficients(spec.index, parameter) != family_coefficients(
                spec.index, -parameter
            ):
                raise AssertionError("the exact T <-> -T quotient changed")
        raw_discriminant = construction.primitive_discriminant_polynomial
        content = polynomial_content(raw_discriminant)
        normalized = tuple(value.numerator // content for value in raw_discriminant)
        discriminants[spec.index] = normalized
        family_records.append(
            {
                "family_index": spec.index,
                "family_label": spec.label,
                "roots": list(spec.roots),
                "calibration_parameter": str(spec.calibration),
                "quartic_square_scale": str(construction.quartic_square_scale),
                "A_coefficients_ascending": list(spec.a_coefficients),
                "B_coefficients_ascending": list(spec.b_coefficients),
                "content_free_discriminant_coefficients_ascending": list(normalized),
                "removed_discriminant_polynomial_content": str(content),
                "exact_symmetry": "A(-T)=A(T), B(-T)=B(T), R_-T=R_T",
            }
        )

    discovery_digest = {
        spec.index: exact_table_digest(spec.index, DISCOVERY_PRIMES)
        for spec in FAMILIES
    }
    held_digest = {
        spec.index: exact_table_digest(spec.index, HELD_PRIMES)
        for spec in FAMILIES
    }
    scans = run_scanners(
        scanner_path,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
        denominator_bound=args.denominator_bound,
        strata=(
            ("global", args.numerator_bound, args.global_keep),
            ("medium", args.medium_numerator_bound, args.medium_keep),
            ("low", args.low_numerator_bound, args.low_keep),
        ),
    )
    for scan in scans:
        spec = family_spec(scan.family_index)
        if (
            scan.discovery_table_digest != discovery_digest[spec.index]
            or scan.held_table_digest != held_digest[spec.index]
            or scan.prior_excluded != 8
            or scan.calibration.parameter != spec.calibration
        ):
            raise AssertionError("a scanner exact-local or exclusion gate changed")
        expected_discovery = score_text(
            spec.index, spec.calibration, DISCOVERY_PRIMES
        )
        expected_held = score_text(spec.index, spec.calibration, HELD_PRIMES)
        if (
            (scan.calibration.discovery_score, scan.calibration.discovery_good)
            != expected_discovery
            or (scan.calibration.held_score, scan.calibration.held_good)
            != expected_held
        ):
            raise AssertionError("the scanner calibration score failed exact replay")
    global_scans = [scan for scan in scans if scan.stratum == "global"]
    if len(global_scans) != 2 or any(
        scan.primitive_population < 18_000_000 for scan in global_scans
    ):
        raise AssertionError("the declared broad rational population shrank")
    print(
        "fresh-prime scans closed: "
        + ", ".join(
            f"{family_spec(scan.family_index).label}={scan.evaluated_population}"
            for scan in global_scans
        ),
        flush=True,
    )

    pools: dict[int, list[dict[str, Any]]] = {}
    pool_audits = {}
    for spec in FAMILIES:
        pool, audit = pool_and_features(
            spec, scans, discriminants[spec.index]
        )
        pools[spec.index] = pool
        pool_audits[spec.label] = audit
    selected, selection = select_conductor_population(pools)
    print(
        f"exact discriminant features closed; conductor population={len(selected)}",
        flush=True,
    )

    common_artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "selection-only leakage-controlled rank14-pair rational frontier"
            if args.selection_only
            else "in-progress conductor-first rank14-pair rational frontier"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "scope": {
            "included_families": [list(spec.roots) for spec in FAMILIES],
            "explicitly_excluded_family": [0, 25, 57, 104, 116, 148],
            "T_sign_quotient": "primitive positive T=a/b only; T and -T are identical",
            "prior_panel_parameters_excluded": [str(value) for value in PRIOR_PANEL_PARAMETERS],
        },
        "families": family_records,
        "frozen_calibrations": calibrations,
        "modular_scan": {
            "score": (
                "sum ((2-a_p)/(p+1-a_p))*log(p), each exact trace term "
                "quantized to 1e-12"
            ),
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_relative_to_prior_broad_scanners_through_prime_199": True,
            "strata": [scanner_record(scan) for scan in scans],
            "global_boxes": [
                {
                    "family_label": family_spec(scan.family_index).label,
                    "numerator": [1, scan.numerator_bound],
                    "denominator": [1, scan.denominator_bound],
                    "primitive_positive_rationals": scan.primitive_population,
                    "prior_panel_excluded": scan.prior_excluded,
                    "evaluated": scan.evaluated_population,
                }
                for scan in global_scans
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
            key: value
            for key, value in vars(args).items()
            if key not in {"output", "selection_only"}
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "scanner_path": str(scanner_path.relative_to(root)),
            "scanner_sha256": sha256_file(scanner_path),
            "frozen_max200_artifact": str(max200_path.relative_to(root)),
            "frozen_max200_artifact_sha256": EXPECTED_MAX200_ARTIFACT_SHA256,
            "frozen_max200_script_sha256": EXPECTED_MAX200_SCRIPT_SHA256,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "temporary_scanner_binary_removed": True,
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
        common_artifact["status"] = "completed selection-only rank14-pair frontier"
        common_artifact["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
        common_artifact["result_sha256"] = stable_json_digest(
            {
                "scope": common_artifact["scope"],
                "scan": common_artifact["modular_scan"],
                "features": common_artifact["exact_discriminant_feature_screen"],
                "selection": selection,
            }
        )
        exclusive_write(args.output, common_artifact)
        return

    conductor_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                conductor_worker,
                record["family_index"], record["numerator"], record["denominator"],
                args.conductor_timeout, args.stack_bytes,
            )
            for record in selected
        ]
        for position, (record, future) in enumerate(zip(selected, futures), start=1):
            record["conductor_phase"] = future.result()
            if position % 32 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    conductor_population_closed = True
    eligible = [
        record
        for record in selected
        if record["conductor_phase"]["status"].startswith("completed")
    ]

    stages = (
        ("H5000", 5_000, None, args.h5000_timeout),
        ("H50000", 50_000, 16, args.h50000_timeout),
        ("H250000", 250_000, 4, args.h250000_timeout),
        ("H1000000", 1_000_000, 1, args.h1000000_timeout),
    )
    current = eligible
    for stage_index, (name, height, keep_per_family, timeout) in enumerate(stages):
        if stage_index:
            prior_name = stages[stage_index - 1][0]
            next_current = []
            for spec in FAMILIES:
                completed = [
                    record
                    for record in current
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
                    record["family_index"], record["numerator"], record["denominator"],
                    height, timeout, args.height_timeout, args.ellrank_timeout,
                    args.stack_bytes, args.mapping_cap, args.certificate_prime_bound,
                )
                for record in current
            ]
            for record, future in zip(current, futures):
                record.setdefault("point_stages", {})[name] = future.result()
        max_rank = max(
            (
                record["point_stages"][name].get("stable_numerical_rank", -1)
                for record in current
            ),
            default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={max_rank}", flush=True)
        if max_rank >= 15:
            print(f"EARLY_SIGNAL {name} stable_rank={max_rank}", flush=True)

    completed_stages = [
        stage
        for record in selected
        for stage in record.get("point_stages", {}).values()
        if stage.get("status") == "completed"
    ]
    maximum_rank = max(
        (int(stage["stable_numerical_rank"]) for stage in completed_stages),
        default=None,
    )
    target_hits = []
    finite_attempts = []
    for record in selected:
        for stage_name, stage in record.get("point_stages", {}).items():
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

    common_artifact["status"] = (
        "completed fixed rank14-pair rational frontier; stopped without broadening"
    )
    common_artifact["target"]["hits"] = target_hits
    common_artifact["conductor_first_screen"] = {
        "population_closed_before_any_point_or_rank_call": conductor_population_closed,
        "selected_population": len(selected),
        "completed": len(eligible),
        "timeouts": sum(
            record["conductor_phase"]["status"] == "timeout" for record in selected
        ),
        "errors": sum(
            record["conductor_phase"]["status"] == "error" for record in selected
        ),
        "subtarget": sum(
            record["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically"
            )
            is True
            for record in selected
        ),
    }
    common_artifact["point_search_protocol"] = {
        "stages": [
            {
                "name": name,
                "height_bound": height,
                "keep_per_family_after_previous_stage": keep,
                "attempted": sum(
                    name in record.get("point_stages", {}) for record in selected
                ),
            }
            for name, height, keep, _ in stages
        ],
        "increasing_heights_are_not_retries": True,
        "same_height_retries": 0,
        "finite_reduction_trigger_stable_rank": FINITE_REDUCTION_TRIGGER,
        "finite_reduction_attempts": finite_attempts,
        "maximum_stable_numerical_rank": maximum_rank,
        "completed_stage_calls": len(completed_stages),
        "stop_rule": "stop without broadening after fixed H1m tier if maximum rank <=15",
        "stop_rule_fired": maximum_rank is not None and maximum_rank <= 15,
        "broadening_calls_after_fixed_protocol": 0,
    }
    common_artifact["timings"].update(
        {
            "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
            "total_wall_seconds": time.monotonic() - started,
        }
    )
    common_artifact["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    common_artifact["provenance"]["owned_processes_remaining"] = 0
    digest_payload = {
        "scope": common_artifact["scope"],
        "families": common_artifact["families"],
        "calibrations": common_artifact["frozen_calibrations"],
        "scan": common_artifact["modular_scan"],
        "features": common_artifact["exact_discriminant_feature_screen"],
        "selection": common_artifact["conductor_selection"],
        "records": selected,
        "conductor": common_artifact["conductor_first_screen"],
        "points": common_artifact["point_search_protocol"],
        "target": common_artifact["target"],
    }
    common_artifact["result_sha256"] = stable_json_digest(digest_payload)
    exclusive_write(args.output, common_artifact)
    print(
        f"complete max_rank={maximum_rank} target_hits={len(target_hits)} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

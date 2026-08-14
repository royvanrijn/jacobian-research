#!/usr/bin/env python3
"""Leakage-controlled rational neighborhood of the new exact rank-15 fiber.

For the six-root Mestre family ``(0,2,136,217,261,290)``, the specialization
``T=2`` has an unconditional finite-reduction rank lower bound 15 and is
conditionally closed at rank 15 by a Delta=11/5 explicit-formula diagnostic.
This continuation therefore searches nearby specializations, not hidden
points on the fixed fiber.  The exact symmetry ``T <-> -T`` lets us retain
only positive T.

The two raw populations are disjoint:

* near strip: ``2 <= b <= 20000``, ``a=2b+d``, ``1 <= |d| <= 32``;
* ordinary window: ``2 <= b <= 2000``, ``3b/2 <= a <= 5b/2``,
  ``|a-2b| >= 33``.

Only primitive pairs are evaluated, so the calibration ``T=2`` and every
previous integer panel fiber are excluded.  A compiled exact-local scanner
uses fresh-for-this-family discovery primes 587..647.  A disjoint 653..733
band is read only after each discovery survivor set closes.  The survivor
union receives exact degree-20 discriminant/radical/powerful-part features.
A fixed rank-aware/diversity conductor population closes before any conductor
or point call.  Every completed conductor receives H=5000; fixed nested tiers
then receive H=50000, H=250000, and H=1000000.  Stable numerical rank at
least 15 triggers immediate exact finite-reduction replay.  No stage retries.
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

ROOTS = (0, 2, 136, 217, 261, 290)
ANCHOR = Q(2)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DISCOVERY_PRIMES = (587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647)
HELD_PRIMES = (653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733)
TRIAL_DIVISION_LIMIT = 997
NEAR_DENOMINATOR_BOUND = 20_000
NEAR_OFFSET_BOUND = 32
NEAR_KEEP = 8_192
NEAR_LOW_DENOMINATOR_BOUND = 512
NEAR_LOW_KEEP = 4_096
ORDINARY_DENOMINATOR_BOUND = 2_000
ORDINARY_KEEP = 8_192
FINITE_REDUCTION_TRIGGER = 15
CERTIFICATE_PRIME_BOUND = 499
STACK_BYTES = 512_000_000
EXPECTED_CERTIFICATE_SHA256 = (
    "35abefefab42b19f49fad074f0c2cd65b039e8f36c398fbe7b46f68a0c2f09ea"
)
EXPECTED_CERTIFICATE_RESULT_SHA256 = (
    "c1de5071cf9ac8bb993345804bb0ab6f96656c72912c294f8e5fe097d002a77b"
)
EXPECTED_EXPLICIT_FORMULA_SHA256 = (
    "6092f3b547e53275cd16842a595488ba97230f299fafae2305d75f92239b070a"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_mestre_02136217261290_t2_rational_neighborhood.json"
)

A_COEFFICIENTS = (
    -542_564_766_112_960_201_029_552,
    0,
    169_665_434_230_038_056_352,
    0,
    -17_668_537_074_117_936,
    0,
    742_147_660_800,
    0,
    -21_676_032,
)
B_COEFFICIENTS = (
    153_797_356_987_326_659_323_597_489_852_634_496,
    0,
    -71_778_413_690_052_019_619_903_289_482_112,
    0,
    12_434_077_250_691_567_747_253_858_944,
    0,
    -867_706_293_306_445_539_767_424,
    0,
    3_733_083_805_979_308_032,
    0,
    1_994_892_912_230_400,
    0,
    -38_843_449_344,
)
DISCRIMINANT_COEFFICIENTS = (
    33_017_927_500_670_460_546_272_292_860_628_938_363_953_258_054_740_480_000,
    0,
    -472_052_647_835_131_008_114_470_174_820_135_715_149_770_726_056_566_928,
    0,
    1_093_181_429_009_327_291_065_198_487_556_339_542_453_831_845_498_465,
    0,
    -848_446_288_998_889_233_947_509_911_294_041_690_038_257_415_110,
    0,
    334_502_116_786_881_401_365_043_998_448_770_646_025_906_761,
    0,
    -77_650_782_783_130_359_146_604_331_928_149_714_065_296,
    0,
    11_241_321_655_351_862_394_306_739_260_262_298_208,
    0,
    -1_023_327_376_153_748_695_480_017_822_390_016,
    0,
    56_504_459_266_007_679_447_127_542_016,
    0,
    -1_704_471_769_394_756_528_947_200,
    0,
    20_987_012_285_890_560_000,
)
REMOVED_DISCRIMINANT_CONTENT = 12_845_056


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
    label: str
    mode: str
    denominator_bound: int
    keep: int
    proposed_population: int
    primitive_population: int
    discovery_table_digest: str
    held_table_digest: str
    calibration_discovery_score: str
    calibration_held_score: str
    calibration_discovery_good: int
    calibration_held_good: int
    candidates: tuple[ScannerCandidate, ...]
    stdout_sha256: str


CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + coefficient
    return answer


def family_coefficients(parameter: Fraction) -> tuple[Fraction, ...]:
    parameter = Q(parameter)
    return (
        Q(0),
        Q(0),
        Q(0),
        evaluate_polynomial(A_COEFFICIENTS, parameter),
        evaluate_polynomial(B_COEFFICIENTS, parameter),
    )


def exact_local_trace_projective(
    numerator: int, denominator: int, prime: int
) -> int | None:
    def homogeneous(coefficients: Sequence[int]) -> int:
        degree = len(coefficients) - 1
        return sum(
            coefficient
            * pow(numerator, power, prime)
            * pow(denominator, degree - power, prime)
            for power, coefficient in enumerate(coefficients)
        ) % prime

    coefficient_a = homogeneous(A_COEFFICIENTS)
    coefficient_b = homogeneous(B_COEFFICIENTS)
    if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
        return None
    character_sum = 0
    for x_value in range(prime):
        rhs = (x_value**3 + coefficient_a * x_value + coefficient_b) % prime
        if rhs:
            symbol = pow(rhs, (prime - 1) // 2, prime)
            character_sum += 1 if symbol == 1 else -1
    return -character_sum


def exact_table_digest(primes: Sequence[int]) -> str:
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
            trace = exact_local_trace_projective(numerator, denominator, prime)
            mix(trace is not None)
            mix(0 if trace is None else trace)
    return str(digest)


def llround(value: float) -> int:
    return int(value + 0.5) if value >= 0 else int(value - 0.5)


def score_text(parameter: Fraction, primes: Sequence[int]) -> tuple[str, int]:
    units = 0
    good = 0
    for prime in primes:
        trace = exact_local_trace_projective(
            parameter.numerator, parameter.denominator, prime
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


def parse_scanner_output(label: str, stdout: str) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_02136217261290_T2_NEIGHBORHOOD_SCAN_V1":
        raise AssertionError("the neighborhood scanner omitted its V1 header")
    discovery = tuple(int(value) for value in lines[1].split()[1:])
    held = tuple(int(value) for value in lines[2].split()[1:])
    if discovery != DISCOVERY_PRIMES or held != HELD_PRIMES:
        raise AssertionError("the neighborhood scanner prime bands changed")
    local = lines[3].split()
    calibration = lines[4].split()
    if local[0] != "L" or calibration[0] != "A" or len(calibration) != 5:
        raise AssertionError("the neighborhood scanner audit header changed")
    candidates = []
    for line in lines[5:-1]:
        fields = line.split()
        if fields[0] != "C" or len(fields) != 7:
            raise AssertionError("a malformed neighborhood candidate escaped")
        candidates.append(
            ScannerCandidate(
                int(fields[1]),
                int(fields[2]),
                fields[3],
                fields[4],
                int(fields[5]),
                int(fields[6]),
            )
        )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 7:
        raise AssertionError("the neighborhood scanner summary changed")
    result = ScannerResult(
        label=label,
        mode=summary[1],
        denominator_bound=int(summary[2]),
        keep=int(summary[3]),
        proposed_population=int(summary[4]),
        primitive_population=int(summary[5]),
        discovery_table_digest=local[1],
        held_table_digest=local[2],
        calibration_discovery_score=calibration[1],
        calibration_held_score=calibration[2],
        calibration_discovery_good=int(calibration[3]),
        calibration_held_good=int(calibration[4]),
        candidates=tuple(candidates),
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )
    if len(candidates) != int(summary[6]) or len(candidates) != result.keep:
        raise AssertionError("the neighborhood scanner retained count changed")
    if len({candidate.parameter for candidate in candidates}) != len(candidates):
        raise AssertionError("the neighborhood scanner emitted a duplicate")
    return result


def run_scanners(
    source: Path,
    *,
    compiler: str,
    compile_timeout: float,
    scan_timeout: float,
) -> tuple[ScannerResult, ...]:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError("a C++17 compiler is required")
    specifications = (
        ("near-global", "near", NEAR_DENOMINATOR_BOUND, NEAR_KEEP),
        ("near-low", "near", NEAR_LOW_DENOMINATOR_BOUND, NEAR_LOW_KEEP),
        ("ordinary-global", "ordinary", ORDINARY_DENOMINATOR_BOUND, ORDINARY_KEEP),
    )
    with tempfile.TemporaryDirectory(prefix="mestre-02136217261290-neighborhood-") as directory:
        binary = Path(directory) / "scan"
        run_capped_process(
            (
                executable,
                "-std=c++17",
                "-O3",
                "-DNDEBUG",
                str(source),
                "-o",
                str(binary),
            ),
            timeout=compile_timeout,
        )
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = [
                executor.submit(
                    run_capped_process,
                    (str(binary), mode, str(bound), str(keep)),
                    timeout=scan_timeout,
                )
                for _, mode, bound, keep in specifications
            ]
            outputs = [future.result()[0] for future in futures]
    return tuple(
        parse_scanner_output(label, stdout)
        for (label, _, _, _), stdout in zip(specifications, outputs)
    )


def iter_scope(mode: str, denominator_bound: int) -> Iterable[tuple[int, int]]:
    for denominator in range(2, denominator_bound + 1):
        if mode == "near":
            for offset in range(-NEAR_OFFSET_BOUND, NEAR_OFFSET_BOUND + 1):
                if offset == 0:
                    continue
                numerator = 2 * denominator + offset
                if numerator > 0 and gcd(numerator, denominator) == 1:
                    yield numerator, denominator
        elif mode == "ordinary":
            first = (3 * denominator + 1) // 2
            last = 5 * denominator // 2
            for numerator in range(first, last + 1):
                if (
                    abs(numerator - 2 * denominator) >= NEAR_OFFSET_BOUND + 1
                    and gcd(numerator, denominator) == 1
                ):
                    yield numerator, denominator
        else:
            raise ValueError("unknown neighborhood mode")


def scope_audit(mode: str, denominator_bound: int) -> dict[str, Any]:
    digest = hashlib.sha256()
    count = 0
    minimum = None
    maximum = None
    for numerator, denominator in iter_scope(mode, denominator_bound):
        parameter = Q(numerator, denominator)
        if parameter == ANCHOR or parameter.denominator != denominator:
            raise AssertionError("the rational scope canonicalization failed")
        digest.update(f"{numerator}/{denominator}\n".encode())
        count += 1
        minimum = parameter if minimum is None or parameter < minimum else minimum
        maximum = parameter if maximum is None or parameter > maximum else maximum
    return {
        "mode": mode,
        "denominator_bound": denominator_bound,
        "primitive_parameter_count": count,
        "parameter_manifest_sha256": digest.hexdigest(),
        "minimum_parameter": str(minimum),
        "maximum_parameter": str(maximum),
    }


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


def homogeneous_value(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def discriminant_feature(numerator: int, denominator: int) -> dict[str, Any]:
    absolute = abs(
        homogeneous_value(DISCRIMINANT_COEFFICIENTS, numerator, denominator)
    )
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
    scans: Sequence[ScannerResult],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_parameter: dict[Fraction, dict[str, Any]] = {}
    for scan in scans:
        for candidate in scan.candidates:
            if candidate.parameter == ANCHOR or candidate.parameter.denominator == 1:
                raise AssertionError("a prior integer panel parameter escaped exclusion")
            expected_mode = (
                "near"
                if abs(candidate.numerator - 2 * candidate.denominator)
                <= NEAR_OFFSET_BOUND
                else "ordinary"
            )
            if scan.mode != expected_mode:
                raise AssertionError("a scanner survivor escaped its raw stratum")
            record = by_parameter.setdefault(
                candidate.parameter,
                {
                    "numerator": candidate.numerator,
                    "denominator": candidate.denominator,
                    "parameter": str(candidate.parameter),
                    "raw_stratum": expected_mode,
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
                or record["raw_stratum"] != expected_mode
            ):
                raise AssertionError("one parameter acquired inconsistent scanner data")
            record["discovery_survivor_strata"].add(scan.label)
    pool = []
    singular = 0
    digest = hashlib.sha256()
    for parameter, record in sorted(by_parameter.items()):
        feature = discriminant_feature(parameter.numerator, parameter.denominator)
        record["discovery_survivor_strata"] = sorted(
            record["discovery_survivor_strata"]
        )
        record["discriminant_feature"] = feature
        digest.update(
            (
                f"{parameter}|{record['raw_stratum']}|"
                f"{record['discovery_score']}|{record['held_score']}|"
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
        "discovery_survivor_union_before_singular_rejection": len(by_parameter),
        "exact_singular_rejections": singular,
        "admissible_feature_pool_count": len(pool),
        "exact_feature_population_sha256": digest.hexdigest(),
    }


def select_conductor_population(
    pool: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reasons: dict[str, set[str]] = defaultdict(set)
    by_parameter = {record["parameter"]: record for record in pool}
    quotas_per_raw_stratum = {
        "highest-held-score": 24,
        "smallest-exact-radical-upper-bound": 20,
        "largest-exact-known-powerful-part": 8,
    }
    for raw_stratum in ("near", "ordinary"):
        stratum_pool = [row for row in pool if row["raw_stratum"] == raw_stratum]
        orders = {
            "highest-held-score": sorted(
                stratum_pool,
                key=lambda row: (
                    -Decimal(row["held_score"]),
                    -Decimal(row["discovery_score"]),
                    row["denominator"],
                    row["numerator"],
                ),
            ),
            "smallest-exact-radical-upper-bound": sorted(
                stratum_pool,
                key=lambda row: (
                    int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                    -Decimal(row["held_score"]),
                    row["denominator"],
                    row["numerator"],
                ),
            ),
            "largest-exact-known-powerful-part": sorted(
                stratum_pool,
                key=lambda row: (
                    -int(row["discriminant_feature"]["known_powerful_part"]),
                    int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                    row["denominator"],
                    row["numerator"],
                ),
            ),
        }
        for label, order in orders.items():
            for record in order[: quotas_per_raw_stratum[label]]:
                reasons[record["parameter"]].add(f"{raw_stratum}-{label}")

    global_orders = {
        "lowest-projective-height": sorted(
            pool,
            key=lambda row: (
                max(row["numerator"], row["denominator"]),
                row["denominator"],
                row["numerator"],
            ),
        ),
        "closest-to-anchor": sorted(
            pool,
            key=lambda row: (
                abs(Q(row["numerator"], row["denominator"]) - ANCHOR),
                row["denominator"],
                row["numerator"],
            ),
        ),
    }
    for label, order in global_orders.items():
        for record in order[:8]:
            reasons[record["parameter"]].add(label)

    selected = []
    for parameter, strata in reasons.items():
        record = dict(by_parameter[parameter])
        record["conductor_selection_strata"] = sorted(strata)
        selected.append(record)
    selected.sort(key=lambda row: (row["denominator"], row["numerator"]))
    digest = hashlib.sha256()
    for record in selected:
        digest.update(
            (
                f"{record['parameter']}|{record['raw_stratum']}|"
                f"{','.join(record['conductor_selection_strata'])}\n"
            ).encode()
        )
    return selected, {
        "selection_uses_conductor": False,
        "selection_uses_point_or_rank_search_data": False,
        "discovery_survivors_closed_before_held_scores": True,
        "held_band_is_holdout_from_discovery_but_is_used_for_conductor_selection": True,
        "exact_discriminant_features_use_no_conductor_or_point_data": True,
        "quota_per_raw_stratum": quotas_per_raw_stratum,
        "global_diversity_quotas": {"lowest-projective-height": 8, "closest-to-anchor": 8},
        "selected_per_raw_stratum": dict(Counter(row["raw_stratum"] for row in selected)),
        "selected_population": len(selected),
        "selected_population_sha256": digest.hexdigest(),
    }


def conductor_worker(
    numerator: int, denominator: int, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    coefficients = family_coefficients(Q(numerator, denominator))
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
    parameter = Q(numerator, denominator)
    coefficients = family_coefficients(parameter)
    try:
        stage, subset = exact_point_stage(
            CONSTRUCTION,
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
        "label": scan.label,
        "mode": scan.mode,
        "denominator_bound": scan.denominator_bound,
        "keep": scan.keep,
        "proposed_population": scan.proposed_population,
        "primitive_population": scan.primitive_population,
        "discovery_table_digest": scan.discovery_table_digest,
        "held_table_digest": scan.held_table_digest,
        "stdout_sha256": scan.stdout_sha256,
        "calibration": {
            "parameter": "2",
            "excluded": True,
            "discovery_score": scan.calibration_discovery_score,
            "held_score": scan.calibration_held_score,
            "discovery_good": scan.calibration_discovery_good,
            "held_good": scan.calibration_held_good,
        },
        "retained_candidate_sha256": stable_json_digest(
            [
                [
                    candidate.numerator,
                    candidate.denominator,
                    candidate.discovery_score,
                    candidate.held_score,
                    candidate.discovery_good,
                    candidate.held_good,
                ]
                for candidate in scan.candidates
            ]
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=45.0)
    parser.add_argument("--conductor-timeout", type=float, default=12.0)
    parser.add_argument("--h5000-timeout", type=float, default=15.0)
    parser.add_argument("--h50000-timeout", type=float, default=20.0)
    parser.add_argument("--h250000-timeout", type=float, default=30.0)
    parser.add_argument("--h1000000-timeout", type=float, default=45.0)
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--ellrank-timeout", type=float, default=12.0)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=CERTIFICATE_PRIME_BOUND)
    parser.add_argument("--stack-bytes", type=int, default=STACK_BYTES)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.compile_timeout,
        args.scan_timeout,
        args.conductor_timeout,
        args.h5000_timeout,
        args.h50000_timeout,
        args.h250000_timeout,
        args.h1000000_timeout,
        args.height_timeout,
        args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 45:
        raise SystemExit("all subprocess caps must lie in (0,45]")
    if args.mapping_cap != 512:
        raise SystemExit("the exact-point mapping cap is pinned at 512")
    if args.certificate_prime_bound != CERTIFICATE_PRIME_BOUND:
        raise SystemExit("the finite-reduction prime bound is pinned at 499")
    if args.stack_bytes != STACK_BYTES:
        raise SystemExit("the PARI stack is pinned at 512000000 bytes")


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
        raise SystemExit("refusing to overwrite the rational-neighborhood artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    scanner_path = script_path.with_name(
        "scan_mestre_02136217261290_t2_neighborhood.cpp"
    )
    certificate_path = (
        root
        / "artifacts/generated-results/"
        "elliptic_mestre_02136217261290_t2_rank15_certificate.json"
    )
    formula_path = (
        root
        / "artifacts/generated-results/"
        "elliptic_mestre_02136217261290_t2_explicit_formula_delta22.json"
    )
    if (
        sha256_file(certificate_path) != EXPECTED_CERTIFICATE_SHA256
        or sha256_file(formula_path) != EXPECTED_EXPLICIT_FORMULA_SHA256
    ):
        raise AssertionError("the exact rank-15 anchor inputs changed")
    certificate = json.loads(certificate_path.read_text())
    formula = json.loads(formula_path.read_text())
    if (
        certificate["result_sha256"] != EXPECTED_CERTIFICATE_RESULT_SHA256
        or certificate["theorem"]["certified_algebraic_rank_lower_bound"] != 15
        or not formula["explicit_formula"][
            "conservative_upper_strictly_below_17"
        ]
    ):
        raise AssertionError("the exact/conditional anchor conclusion changed")

    if (
        CONSTRUCTION.quartic_condition
        or CONSTRUCTION.is_reflection_symmetric
        or CONSTRUCTION.quartic_square_scale != 1260
    ):
        raise AssertionError("the selected six-root family geometry changed")
    for parameter in (Q(1), Q(2), Q(11, 5), Q(-17, 3)):
        if (
            CONSTRUCTION.primitive_jacobian_coefficients(parameter)
            != family_coefficients(parameter)
            or family_coefficients(parameter) != family_coefficients(-parameter)
        ):
            raise AssertionError("the exact A(T),B(T) formula or sign quotient changed")
    raw_discriminant = CONSTRUCTION.primitive_discriminant_polynomial
    if any(value.denominator != 1 for value in raw_discriminant):
        raise AssertionError("the primitive discriminant polynomial became rational")
    content = gcd(*(abs(value.numerator) for value in raw_discriminant if value))
    normalized_discriminant = tuple(
        value.numerator // content for value in raw_discriminant
    )
    if (
        content != REMOVED_DISCRIMINANT_CONTENT
        or normalized_discriminant != DISCRIMINANT_COEFFICIENTS
    ):
        raise AssertionError("the pinned degree-20 discriminant formula changed")

    audits = {
        "near": scope_audit("near", NEAR_DENOMINATOR_BOUND),
        "near_low": scope_audit("near", NEAR_LOW_DENOMINATOR_BOUND),
        "ordinary": scope_audit("ordinary", ORDINARY_DENOMINATOR_BOUND),
    }
    # The definitions are disjoint by the exact offset partition
    # |a-2b|<=32 versus |a-2b|>=33.  Stream the union digest without retaining
    # roughly two million Python tuples in memory.
    raw_union_digest = hashlib.sha256()
    raw_union_count = 0
    for mode, bound in (
        ("near", NEAR_DENOMINATOR_BOUND),
        ("ordinary", ORDINARY_DENOMINATOR_BOUND),
    ):
        for numerator, denominator in iter_scope(mode, bound):
            raw_union_digest.update(f"{mode}|{numerator}/{denominator}\n".encode())
            raw_union_count += 1

    discovery_digest = exact_table_digest(DISCOVERY_PRIMES)
    held_digest = exact_table_digest(HELD_PRIMES)
    calibration_discovery = score_text(ANCHOR, DISCOVERY_PRIMES)
    calibration_held = score_text(ANCHOR, HELD_PRIMES)
    scans = run_scanners(
        scanner_path,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    expected_scans = {
        "near-global": ("near", NEAR_DENOMINATOR_BOUND, NEAR_KEEP, audits["near"]),
        "near-low": ("near", NEAR_LOW_DENOMINATOR_BOUND, NEAR_LOW_KEEP, audits["near_low"]),
        "ordinary-global": (
            "ordinary",
            ORDINARY_DENOMINATOR_BOUND,
            ORDINARY_KEEP,
            audits["ordinary"],
        ),
    }
    for scan in scans:
        mode, bound, keep, audit = expected_scans[scan.label]
        if (
            scan.mode != mode
            or scan.denominator_bound != bound
            or scan.keep != keep
            or scan.primitive_population != audit["primitive_parameter_count"]
            or scan.discovery_table_digest != discovery_digest
            or scan.held_table_digest != held_digest
            or (
                scan.calibration_discovery_score,
                scan.calibration_discovery_good,
            )
            != calibration_discovery
            or (scan.calibration_held_score, scan.calibration_held_good)
            != calibration_held
        ):
            raise AssertionError("a scanner population/local replay gate changed")
    print(
        "fresh-prime scopes closed: "
        f"near={audits['near']['primitive_parameter_count']} "
        f"ordinary={audits['ordinary']['primitive_parameter_count']}",
        flush=True,
    )

    pool, pool_audit = pool_and_features(scans)
    selected, selection = select_conductor_population(pool)
    print(
        f"exact discriminant features closed: pool={len(pool)} "
        f"conductors={len(selected)}",
        flush=True,
    )

    conductor_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_records = {
            executor.submit(
                conductor_worker,
                record["numerator"],
                record["denominator"],
                args.conductor_timeout,
                args.stack_bytes,
            ): record
            for record in selected
        }
        completed = 0
        for future in as_completed(future_records):
            record = future_records[future]
            record["conductor_phase"] = future.result()
            completed += 1
            if completed % 24 == 0 or completed == len(selected):
                print(f"conductors {completed}/{len(selected)}", flush=True)
    conductor_population_closed = True
    eligible = [
        record
        for record in selected
        if record["conductor_phase"]["status"].startswith("completed")
    ]

    stages = (
        ("H5000", 5_000, None, args.h5000_timeout),
        ("H50000", 50_000, 32, args.h50000_timeout),
        ("H250000", 250_000, 8, args.h250000_timeout),
        ("H1000000", 1_000_000, 2, args.h1000000_timeout),
    )
    current = eligible
    for stage_index, (name, height, keep, timeout) in enumerate(stages):
        if stage_index:
            prior_name = stages[stage_index - 1][0]
            completed_prior = [
                record
                for record in current
                if record.get("point_stages", {}).get(prior_name, {}).get("status")
                == "completed"
            ]
            completed_prior.sort(key=lambda record: stage_rank_key(record, prior_name))
            current = completed_prior[: int(keep)]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_records = {
                executor.submit(
                    point_stage_worker,
                    record["numerator"],
                    record["denominator"],
                    height,
                    timeout,
                    args.height_timeout,
                    args.ellrank_timeout,
                    args.stack_bytes,
                    args.mapping_cap,
                    args.certificate_prime_bound,
                ): record
                for record in current
            }
            completed = 0
            for future in as_completed(future_records):
                record = future_records[future]
                stage = future.result()
                record.setdefault("point_stages", {})[name] = stage
                completed += 1
                stable_rank = stage.get("stable_numerical_rank")
                if stable_rank is not None and int(stable_rank) >= FINITE_REDUCTION_TRIGGER:
                    print(
                        f"EARLY_SIGNAL {record['parameter']} {name} "
                        f"stable_rank={stable_rank}",
                        flush=True,
                    )
                    certified = stage.get("finite_reduction_attempt", {}).get(
                        "certified_algebraic_rank_lower_bound"
                    )
                    if certified is not None:
                        print(
                            f"EXACT_SIGNAL {record['parameter']} {name} "
                            f"certified_rank={certified}",
                            flush=True,
                        )
                if completed % 16 == 0 or completed == len(current):
                    print(f"{name} {completed}/{len(current)}", flush=True)
        maximum = max(
            (
                record["point_stages"][name].get("stable_numerical_rank", -1)
                for record in current
            ),
            default=-1,
        )
        print(f"{name} complete max_rank={maximum}", flush=True)

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
        exact = stage.get("finite_reduction_attempt", {})
        certified = exact.get("certified_algebraic_rank_lower_bound")
        if certified is None:
            continue
        finite_attempts.append(
            {
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
                    "parameter": record["parameter"],
                    "stage": stage_name,
                    "certified_rank_lower_bound": certified,
                    "conductor": record["conductor_phase"]["conductor"],
                    "log_conductor": record["conductor_phase"]["log_conductor"],
                }
            )

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete fixed rational-T neighborhood; stopped without broadening",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "family": {
            "roots": list(ROOTS),
            "anchor_parameter": "2",
            "quartic_square_scale": str(CONSTRUCTION.quartic_square_scale),
            "A_coefficients_ascending": list(A_COEFFICIENTS),
            "B_coefficients_ascending": list(B_COEFFICIENTS),
            "content_free_discriminant_coefficients_ascending": list(
                DISCRIMINANT_COEFFICIENTS
            ),
            "removed_discriminant_polynomial_content": str(
                REMOVED_DISCRIMINANT_CONTENT
            ),
            "exact_symmetry": "A(-T)=A(T), B(-T)=B(T), R_-T=R_T",
        },
        "anchor": {
            "excluded_from_every_searchable_population": True,
            "rank15_certificate": str(certificate_path.relative_to(root)),
            "rank15_certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "certified_algebraic_rank_lower_bound": 15,
            "conductor": certificate["curve"]["conductor"],
            "log_conductor": certificate["curve"]["log_conductor"],
            "root_number": certificate["curve"]["root_number"],
            "explicit_formula_diagnostic": str(formula_path.relative_to(root)),
            "explicit_formula_diagnostic_sha256": EXPECTED_EXPLICIT_FORMULA_SHA256,
            "conditional_upper_under_grh": formula["explicit_formula"][
                "conservative_explicit_formula_upper"
            ],
            "fixed_fiber_search_priority": "conditionally closed; search neighbors",
        },
        "scope": {
            "T_sign_quotient": "primitive positive T=a/b only; T and -T are identical",
            "near_strip": {
                **audits["near"],
                "definition": "b=2..20000, a=2b+delta, 1<=|delta|<=32, gcd(a,b)=1",
            },
            "ordinary_window": {
                **audits["ordinary"],
                "definition": (
                    "b=2..2000, ceil(3b/2)<=a<=floor(5b/2), "
                    "|a-2b|>=33, gcd(a,b)=1"
                ),
            },
            "raw_populations_are_disjoint": True,
            "raw_union_primitive_parameter_count": raw_union_count,
            "raw_union_manifest_sha256": raw_union_digest.hexdigest(),
            "prior_integer_panel_parameters_excluded": [str(value) for value in range(1, 9)],
            "unsearched_population": (
                "all other positive rational parameters outside the two declared raw "
                "populations, plus discovery non-survivors within them"
            ),
        },
        "modular_scan": {
            "score": (
                "sum ((2-a_p)/(p+1-a_p))*log(p), each exact trace term "
                "quantized to 1e-12"
            ),
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_for_this_family_relative_to_the_max300_panel": True,
            "discovery_table_digest": discovery_digest,
            "held_table_digest": held_digest,
            "strata": [scanner_record(scan) for scan in scans],
        },
        "exact_discriminant_feature_screen": {
            "trial_division_prime_bound": TRIAL_DIVISION_LIMIT,
            "content_free_homogeneous_degree": 20,
            **pool_audit,
        },
        "conductor_selection": selection,
        "selected_records": selected,
        "conductor_first_screen": {
            "population_closed_before_any_conductor_or_point_call": True,
            "all_conductor_calls_completed_before_any_point_call": conductor_population_closed,
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
        },
        "point_search_protocol": {
            "stages": [
                {
                    "name": name,
                    "height_bound": height,
                    "keep_after_previous_stage": keep,
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
            "broadening_calls_after_fixed_protocol": 0,
        },
        "parameters": {
            key: value for key, value in vars(args).items() if key != "output"
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "scanner_path": str(scanner_path.relative_to(root)),
            "scanner_sha256": sha256_file(scanner_path),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "temporary_scanner_binary_removed": True,
            "external_calls_use_foreground_process_groups": True,
            "same_stage_retries": 0,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "compiler": shutil.which(args.compiler),
        },
        "timings": {
            "pre_conductor_wall_seconds": conductor_started - started,
            "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
            "total_wall_seconds": time.monotonic() - started,
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "target": artifact["target"],
            "family": artifact["family"],
            "anchor": artifact["anchor"],
            "scope": artifact["scope"],
            "scan": artifact["modular_scan"],
            "features": artifact["exact_discriminant_feature_screen"],
            "selection": artifact["conductor_selection"],
            "records": selected,
            "conductor": artifact["conductor_first_screen"],
            "points": artifact["point_search_protocol"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        f"complete max_rank={maximum_rank} exact_attempts={len(finite_attempts)} "
        f"target_hits={len(target_hits)} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

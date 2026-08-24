#!/usr/bin/env python3
"""Leakage-controlled rational search in the rank-14 Mestre family.

The family attached to roots ``(0,6,49,73,82,96)`` has an exact even model

``y^2 = x^3 + A(T)x + B(T)``.

Thus ``T`` and ``-T`` define the same curve, and this search exhausts the
positive reduced representatives ``T=a/b`` in a declared rational box.  A
compiled scanner ranks all parameters using exact local traces at a discovery
prime band.  A disjoint held band ranks only the fixed discovery survivors.
Exact degree-20 discriminants and fixed-prime radical valuations provide two
independent conductor-blind selection strata.  The selected population is
fixed before any conductor call, and every conductor call closes before any
point or rank computation.

The complete conducted population receives H=5000 point searches.  Fixed
rank-driven tranches then receive H=50000, H=250000, and H=1000000.  These are
increasing stages, never retries.  Stable numerical rank at least 18 triggers
small-prime saturation followed immediately by an exact finite-reduction
attempt.  Numerical height ranks are triage evidence only.
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
from math import gcd, log
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

from extend_nagao_u42_frontier import saturate_exact_basis
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    TARGET_LOG_CONDUCTOR,
    bounded_quartic_points,
    canonical_signless_points,
    capped_minimal_curve_data,
    finite_reduction_attempt,
    height_matrix_replay,
    numerical_subset,
    pari_version_capped,
    point_digest,
    point_record,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    run_capped_process,
    sha256_file,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 6, 49, 73, 82, 96)
CALIBRATION_PARAMETER = Q(2)
PRIOR_PARAMETERS = tuple(Q(value) for value in range(1, 9))
DISCOVERY_PRIMES = (
    401,
    409,
    419,
    421,
    431,
    433,
    439,
    443,
    449,
    457,
    461,
    463,
    467,
    479,
    487,
    491,
    499,
)
HELD_PRIMES = (
    503,
    509,
    521,
    523,
    541,
    547,
    557,
    563,
    569,
    571,
    577,
    587,
    593,
    599,
)
A_COEFFICIENTS = (
    -5_657_385_312_367_322_049_797_427,
    0,
    9_743_173_073_723_224_641_600,
    0,
    -5_372_891_384_935_240_224,
    0,
    1_062_430_937_568_000,
    0,
    -130_657_052_592,
)
B_COEFFICIENTS = (
    5_176_687_041_098_134_891_263_619_737_397_215_246,
    0,
    -13_244_429_677_973_820_639_898_089_667_582_800,
    0,
    12_416_239_107_455_669_125_977_672_074_196,
    0,
    -4_698_896_181_839_363_925_068_659_200,
    0,
    230_904_028_542_265_215_172_992,
    0,
    221_720_837_222_941_056_000,
    0,
    -18_178_054_413_019_776,
)
DISCRIMINANT_CONSTANT = 14_517_450_288
DISCRIMINANT_Q16_COEFFICIENTS = (
    23_165_770_029_325_322_571_826_973_172_320_138_800_411_677_520_896,
    0,
    -1_163_925_863_273_107_648_241_445_795_210_860_755_572_313_596_900,
    0,
    4_830_411_827_954_034_723_871_093_515_078_324_644_658_541_485,
    0,
    -8_675_920_665_336_586_393_032_512_192_859_889_372_078_800,
    0,
    8_464_284_895_868_434_215_223_214_736_698_986_851_387,
    0,
    -4_766_812_593_706_309_656_183_919_466_422_516_200,
    0,
    1_519_365_706_898_882_251_653_231_974_688_432,
    0,
    -245_280_338_708_354_476_792_442_003_200,
    0,
    14_336_730_898_511_653_632_000_000,
)
STACK_BYTES = 512_000_000
SIGNAL_CERTIFICATE_TRIGGER = 18
SIGNAL_CERTIFICATE_PRIME_BOUND = 2_000
FROZEN_COMPLETE_SCRIPT_SHA256 = (
    "5fefa372dca3563bedef08a55b1adad8c7db26c5161db38faa89a2e12d59ef6c"
)
FROZEN_COMPLETE_ARTIFACT_SHA256 = (
    "c2cb2e68a54cc1625224cf1aed1cce0196582c16d21373bdbc86c67a1e91d24c"
)
SELECTION_QUOTAS = {
    "held_local_score": 80,
    "exact_fixed_prime_powerful_part": 64,
    "exact_discriminant_height": 32,
    "balanced_held_and_power_rank": 32,
}
POINT_STAGES = (
    ("H5000", 5_000, None, 128),
    ("H50000", 50_000, 32, 256),
    ("H250000", 250_000, 8, 384),
    ("H1000000", 1_000_000, 2, 512),
)


@dataclass(frozen=True)
class Candidate:
    numerator: int
    denominator: int
    discovery_score: str
    held_score: str
    discovery_good: int
    held_good: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)

    @property
    def identifier(self) -> str:
        return f"t{self.numerator}_{self.denominator}"


@dataclass(frozen=True)
class ScannerResult:
    numerator_bound: int
    denominator_bound: int
    keep: int
    primitive_population: int
    prior_excluded: int
    evaluated_population: int
    discovery_table_digest: str
    held_table_digest: str
    calibration: Candidate
    candidates: tuple[Candidate, ...]
    stdout_sha256: str


def stable_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def evaluate_polynomial(
    coefficients: Sequence[int], value: Fraction
) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def homogeneous_value(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def family_coefficients(parameter: Fraction) -> tuple[Fraction, ...]:
    parameter = Q(parameter)
    return (
        Q(0),
        Q(0),
        Q(0),
        evaluate_polynomial(A_COEFFICIENTS, parameter),
        evaluate_polynomial(B_COEFFICIENTS, parameter),
    )


def exact_discriminant(parameter: Fraction) -> Fraction:
    parameter = Q(parameter)
    numerator = parameter.numerator
    denominator = parameter.denominator
    raw = (
        DISCRIMINANT_CONSTANT
        * (2 * numerator - 45 * denominator) ** 2
        * (2 * numerator + 45 * denominator) ** 2
        * homogeneous_value(
            DISCRIMINANT_Q16_COEFFICIENTS, numerator, denominator
        )
    )
    return Q(raw, denominator**20)


def primes_up_to(bound: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(bound**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : bound + 1 : prime] = b"\x00" * (
                (bound - prime * prime) // prime + 1
            )
    return tuple(index for index, flag in enumerate(sieve) if flag)


RADICAL_PRIMES = primes_up_to(397)


def discriminant_features(candidate: Candidate) -> dict[str, Any]:
    discriminant = exact_discriminant(candidate.parameter)
    numerator = abs(discriminant.numerator)
    if numerator == 0:
        return {
            "exact_nonsingular": False,
            "discriminant_numerator": "0",
            "discriminant_denominator": str(discriminant.denominator),
        }
    residual = numerator
    valuations: dict[str, int] = {}
    radical = 1
    powerful_part = 1
    for prime in RADICAL_PRIMES:
        valuation = 0
        while residual % prime == 0:
            residual //= prime
            valuation += 1
        if valuation:
            valuations[str(prime)] = valuation
            radical *= prime
            powerful_part *= prime ** (valuation - 1)
    return {
        "exact_nonsingular": True,
        "discriminant_numerator": str(discriminant.numerator),
        "discriminant_denominator": str(discriminant.denominator),
        "discriminant_numerator_bits": numerator.bit_length(),
        "discriminant_denominator_bits": discriminant.denominator.bit_length(),
        "fixed_prime_panel": [2, 397],
        "fixed_prime_valuations": valuations,
        "fixed_prime_distinct_divisor_count": len(valuations),
        "fixed_prime_radical": str(radical),
        "fixed_prime_powerful_part": str(powerful_part),
        "fixed_prime_powerful_savings_log": sum(
            (valuation - 1) * log(int(prime))
            for prime, valuation in valuations.items()
        ),
        "residual_after_fixed_prime_removal_bits": residual.bit_length(),
        "forced_linear_square_factors": [
            str(abs(2 * candidate.numerator - 45 * candidate.denominator)),
            str(abs(2 * candidate.numerator + 45 * candidate.denominator)),
        ],
        "features_are_exact_but_not_a_complete_factorization": True,
    }


def homogeneous_coefficient(
    coefficients: Sequence[int], numerator: int, denominator: int, prime: int
) -> int:
    return homogeneous_value(coefficients, numerator, denominator) % prime


def exact_local_trace_projective(
    numerator: int, denominator: int, prime: int
) -> int | None:
    coefficient_a = homogeneous_coefficient(
        A_COEFFICIENTS, numerator, denominator, prime
    )
    coefficient_b = homogeneous_coefficient(
        B_COEFFICIENTS, numerator, denominator, prime
    )
    if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
        return None
    character_sum = 0
    for x_value in range(prime):
        rhs = (x_value**3 + coefficient_a * x_value + coefficient_b) % prime
        if rhs:
            symbol = pow(rhs, (prime - 1) // 2, prime)
            character_sum += 1 if symbol == 1 else -1
    return -character_sum


def exact_local_table_digest(primes: Sequence[int]) -> str:
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


def quantized_score(parameter: Fraction, primes: Sequence[int]) -> tuple[str, int]:
    parameter = Q(parameter)
    units = 0
    good = 0
    for prime in primes:
        trace = exact_local_trace_projective(
            parameter.numerator, parameter.denominator, prime
        )
        if trace is None:
            continue
        units += int(
            (2 - trace) / (prime + 1 - trace) * log(prime) * 1_000_000_000_000
            + 0.5
        )
        good += 1
    sign = "-" if units < 0 else ""
    absolute = abs(units)
    return (
        f"{sign}{absolute // 1_000_000_000_000}."
        f"{absolute % 1_000_000_000_000:012d}",
        good,
    )


def parse_candidate(fields: Sequence[str]) -> Candidate:
    if len(fields) != 7:
        raise AssertionError("malformed scanner candidate")
    return Candidate(
        numerator=int(fields[1]),
        denominator=int(fields[2]),
        discovery_score=fields[3],
        held_score=fields[4],
        discovery_good=int(fields[5]),
        held_good=int(fields[6]),
    )


def parse_scanner_output(stdout: str) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_0649738296_RATIONAL_SCAN_V1":
        raise AssertionError("the scanner omitted its exact format header")
    if tuple(map(int, lines[1].split()[1:])) != DISCOVERY_PRIMES:
        raise AssertionError("the discovery prime band changed")
    if tuple(map(int, lines[2].split()[1:])) != HELD_PRIMES:
        raise AssertionError("the held prime band changed")
    local = lines[3].split()
    if local[0] != "L" or len(local) != 3:
        raise AssertionError("the scanner omitted table digests")
    calibration_fields = lines[4].split()
    if calibration_fields[0] != "K":
        raise AssertionError("the scanner omitted T=2 calibration")
    calibration = parse_candidate(("C", *calibration_fields[1:]))
    candidates = tuple(
        parse_candidate(line.split()) for line in lines[5:-1] if line.startswith("C ")
    )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 8:
        raise AssertionError("the scanner omitted its population summary")
    numerator_bound, denominator_bound, keep = map(int, summary[1:4])
    primitive, prior, evaluated, retained = map(int, summary[4:])
    if retained != len(candidates) or retained != min(keep, evaluated):
        raise AssertionError("the scanner retained count changed")
    if calibration.parameter != CALIBRATION_PARAMETER:
        raise AssertionError("the scanner calibration parameter changed")
    if any(
        gcd(candidate.numerator, candidate.denominator) != 1
        or not (1 <= candidate.numerator <= numerator_bound)
        or not (1 <= candidate.denominator <= denominator_bound)
        or candidate.parameter in PRIOR_PARAMETERS
        for candidate in candidates
    ):
        raise AssertionError("a scanner survivor escaped the reduced/exclusion box")
    expected_order = sorted(
        candidates,
        key=lambda candidate: (
            -Decimal(candidate.discovery_score),
            -candidate.discovery_good,
            candidate.denominator,
            candidate.numerator,
        ),
    )
    if list(candidates) != expected_order:
        raise AssertionError("the discovery survivor ordering changed")
    return ScannerResult(
        numerator_bound=numerator_bound,
        denominator_bound=denominator_bound,
        keep=keep,
        primitive_population=primitive,
        prior_excluded=prior,
        evaluated_population=evaluated,
        discovery_table_digest=local[1],
        held_table_digest=local[2],
        calibration=calibration,
        candidates=candidates,
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )


def run_scanner(
    source: Path,
    *,
    numerator_bound: int,
    denominator_bound: int,
    keep: int,
    compile_timeout: float,
    scan_timeout: float,
) -> ScannerResult:
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="mestre-0649738296-") as directory:
        binary = Path(directory) / "scanner"
        run_capped_process(
            (
                compiler,
                "-std=c++17",
                "-O3",
                "-DNDEBUG",
                str(source),
                "-o",
                str(binary),
            ),
            timeout=compile_timeout,
        )
        stdout, _ = run_capped_process(
            (
                str(binary),
                str(numerator_bound),
                str(denominator_bound),
                str(keep),
            ),
            timeout=scan_timeout,
        )
    result = parse_scanner_output(stdout)
    if result.discovery_table_digest != exact_local_table_digest(DISCOVERY_PRIMES):
        raise AssertionError("the discovery table missed exact Python replay")
    if result.held_table_digest != exact_local_table_digest(HELD_PRIMES):
        raise AssertionError("the held table missed exact Python replay")
    expected_discovery = quantized_score(CALIBRATION_PARAMETER, DISCOVERY_PRIMES)
    expected_held = quantized_score(CALIBRATION_PARAMETER, HELD_PRIMES)
    if (
        abs(
            Decimal(result.calibration.discovery_score)
            - Decimal(expected_discovery[0])
        )
        > Decimal("0.000000000100")
        or result.calibration.discovery_good != expected_discovery[1]
        or abs(Decimal(result.calibration.held_score) - Decimal(expected_held[0]))
        > Decimal("0.000000000100")
        or result.calibration.held_good != expected_held[1]
    ):
        raise AssertionError("the T=2 local calibration missed exact replay")
    return result


def candidate_record(candidate: Candidate, features: dict[str, Any]) -> dict[str, Any]:
    return {
        "identifier": candidate.identifier,
        "numerator": candidate.numerator,
        "denominator": candidate.denominator,
        "parameter": str(candidate.parameter),
        "discovery_score": candidate.discovery_score,
        "held_score": candidate.held_score,
        "discovery_good_primes": candidate.discovery_good,
        "held_good_primes": candidate.held_good,
        "exact_discriminant_features": features,
    }


def select_population(
    candidates: Sequence[Candidate], features: dict[str, dict[str, Any]]
) -> tuple[tuple[Candidate, tuple[str, ...]], ...]:
    nonsingular = [
        candidate
        for candidate in candidates
        if features[candidate.identifier]["exact_nonsingular"]
    ]
    held = sorted(
        nonsingular,
        key=lambda candidate: (
            -Decimal(candidate.held_score),
            -Decimal(candidate.discovery_score),
            candidate.denominator,
            candidate.numerator,
        ),
    )
    power = sorted(
        nonsingular,
        key=lambda candidate: (
            -int(
                features[candidate.identifier]["fixed_prime_powerful_part"]
            ),
            features[candidate.identifier]["discriminant_numerator_bits"],
            -Decimal(candidate.held_score),
            candidate.denominator,
            candidate.numerator,
        ),
    )
    height = sorted(
        nonsingular,
        key=lambda candidate: (
            features[candidate.identifier]["discriminant_numerator_bits"],
            features[candidate.identifier]["discriminant_denominator_bits"],
            -Decimal(candidate.held_score),
            candidate.denominator,
            candidate.numerator,
        ),
    )
    held_rank = {candidate.identifier: index for index, candidate in enumerate(held)}
    power_rank = {candidate.identifier: index for index, candidate in enumerate(power)}
    balanced = sorted(
        nonsingular,
        key=lambda candidate: (
            held_rank[candidate.identifier] + power_rank[candidate.identifier],
            max(held_rank[candidate.identifier], power_rank[candidate.identifier]),
            candidate.denominator,
            candidate.numerator,
        ),
    )
    rankings = {
        "held_local_score": held,
        "exact_fixed_prime_powerful_part": power,
        "exact_discriminant_height": height,
        "balanced_held_and_power_rank": balanced,
    }
    selected: dict[str, tuple[Candidate, set[str]]] = {}
    for name, quota in SELECTION_QUOTAS.items():
        added = 0
        for candidate in rankings[name]:
            if candidate.identifier in selected:
                selected[candidate.identifier][1].add(name)
                continue
            selected[candidate.identifier] = (candidate, {name})
            added += 1
            if added == quota:
                break
        if added != quota:
            raise AssertionError(f"selection stratum {name} failed to fill")
    return tuple(
        (candidate, tuple(sorted(sources)))
        for candidate, sources in selected.values()
    )


def exact_point_stage(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    coefficients: Sequence[Fraction],
    *,
    height_bound: int,
    mapping_cap: int,
    point_timeout: float,
    height_timeout: float,
    saturation_timeout: float,
) -> dict[str, Any]:
    parameter = Q(parameter)
    visible_quartic = primitive_visible_points(construction, parameter)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point)
        for point in visible_quartic
    )
    quartic_coefficients = construction.primitive_quartic_coefficients(parameter)
    raw = bounded_quartic_points(
        quartic_coefficients,
        height_bound=height_bound,
        timeout=point_timeout,
        stack_bytes=STACK_BYTES,
    )
    signless = canonical_signless_points(raw)
    retained = signless[:mapping_cap]
    if any(
        point[1] ** 2 != quartic_value(quartic_coefficients, point[0])
        for point in retained
    ):
        raise AssertionError("a bounded search returned a point off the quartic")
    searched_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point)
        for point in retained
    )
    pool_by_x = {point[0]: point for point in visible_jacobian}
    for point in searched_jacobian:
        pool_by_x.setdefault(point[0], point)
    pool = tuple(pool_by_x.values())
    height = height_matrix_replay(
        coefficients,
        pool,
        precisions=(72, 120),
        timeout=height_timeout,
        stack_bytes=STACK_BYTES,
    )
    stable_rank = int(height[-1]["numerical_rank"])
    subset = numerical_subset(pool, height)
    result: dict[str, Any] = {
        "status": "completed exact point membership and numerical height triage",
        "height_bound": height_bound,
        "signed_quartic_points_returned": len(raw),
        "distinct_nonzero_ordinate_abscissas": len(signless),
        "retained_abscissas": len(retained),
        "mapping_cap": mapping_cap,
        "mapping_truncated": len(signless) > len(retained),
        "visible_quartic_point_count": len(visible_quartic),
        "visible_quartic_point_sha256": point_digest(visible_quartic),
        "visible_jacobian_point_sha256": point_digest(visible_jacobian),
        "pool_point_count_modulo_inverse": len(pool),
        "pool_point_sha256": point_digest(pool),
        "exact_quartic_and_jacobian_membership_checked": True,
        "height_matrix_runs": list(height),
        "stable_numerical_rank": stable_rank,
        "numerical_subset": [point_record(point) for point in subset],
        "numerical_rank_is_not_an_independence_certificate": True,
    }
    if stable_rank >= SIGNAL_CERTIFICATE_TRIGGER:
        certificate_basis = subset
        try:
            saturated, saturation = saturate_exact_basis(
                coefficients,
                subset,
                prime_bound=50,
                timeout=saturation_timeout,
                stack_bytes=STACK_BYTES,
            )
            if len(saturated) != len(subset):
                raise AssertionError("saturation changed the proposed basis length")
            certificate_basis = saturated
            result["small_prime_saturation"] = saturation
        except CappedProcessTimeout:
            result["small_prime_saturation"] = {
                "status": "timeout-no-retry",
                "timeout_seconds": saturation_timeout,
            }
        except Exception as error:
            result["small_prime_saturation"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
        result["finite_reduction_attempt"] = finite_reduction_attempt(
            coefficients,
            certificate_basis,
            prime_bound=SIGNAL_CERTIFICATE_PRIME_BOUND,
        )
    else:
        result["finite_reduction_attempt"] = {
            "status": "not triggered",
            "trigger_stable_numerical_rank": SIGNAL_CERTIFICATE_TRIGGER,
        }
    return result


def phase_rank(phase: dict[str, Any]) -> int:
    return int(phase.get("stable_numerical_rank", -1))


def stage_key(item: dict[str, Any], prior_stage: str) -> tuple[Any, ...]:
    phase = item["record"]["point_stages"][prior_stage]
    rank = phase_rank(phase)
    conductor = item["record"]["conductor_phase"]
    expected_parity = 1 if conductor["root_number"] == -1 else 0
    parity_mismatch = rank >= 0 and rank % 2 != expected_parity
    return (
        -rank,
        -int(parity_mismatch),
        not conductor["below_strict_log_conductor_target_numerically"],
        -Decimal(item["candidate"].held_score),
        -int(item["features"]["fixed_prime_powerful_part"]),
        item["candidate"].denominator,
        item["candidate"].numerator,
    )


def result_digest(artifact: dict[str, Any]) -> str:
    compact = []
    for record in artifact["conductor_first_screen"]["records"]:
        stages = []
        for name in sorted(record.get("point_stages", {})):
            phase = record["point_stages"][name]
            certificate = phase.get("finite_reduction_attempt", {})
            stages.append(
                [
                    name,
                    phase["status"],
                    phase.get("stable_numerical_rank"),
                    phase.get("pool_point_sha256"),
                    certificate.get("certified_algebraic_rank_lower_bound"),
                    certificate.get("point_sha256"),
                ]
            )
        compact.append(
            [
                record["identifier"],
                record["selection_strata"],
                record["conductor_phase"]["status"],
                record["conductor_phase"].get("conductor"),
                stages,
            ]
        )
    return stable_digest(
        {
            "family": artifact["family"],
            "calibration": artifact["calibration_T2"],
            "scan": artifact["modular_scan"],
            "selection": artifact["selection"],
            "population": artifact["conductor_first_screen"]["population"],
            "records": compact,
            "protocol": artifact["point_search_protocol"],
            "target": artifact["target"],
        }
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numerator-bound", type=int, default=30_000)
    parser.add_argument("--denominator-bound", type=int, default=1_000)
    parser.add_argument("--discovery-keep", type=int, default=8_192)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=10.0)
    parser.add_argument("--point-timeout", type=float, default=20.0)
    parser.add_argument("--deep-point-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=15.0)
    parser.add_argument("--saturation-timeout", type=float, default=25.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_mestre_0649738296_rational.json",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.numerator_bound < 30_000 or args.numerator_bound > 100_000:
        raise SystemExit("--numerator-bound must lie in [30000,100000]")
    if args.denominator_bound < 1_000 or args.denominator_bound > 5_000:
        raise SystemExit("--denominator-bound must lie in [1000,5000]")
    if args.discovery_keep < 1_024 or args.discovery_keep > 20_000:
        raise SystemExit("--discovery-keep must lie in [1024,20000]")
    caps = (
        args.compile_timeout,
        args.scan_timeout,
        args.conductor_timeout,
        args.point_timeout,
        args.deep_point_timeout,
        args.height_timeout,
        args.saturation_timeout,
    )
    if min(caps) <= 0 or max(caps) > 30:
        raise SystemExit("all process caps must lie in (0,30]")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rational-family artifact")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    root = Path(__file__).resolve().parents[2]
    script = Path(__file__).resolve()
    scanner_source = script.with_name("scan_mestre_0649738296_rational.cpp")
    complete_script = script.with_name(
        "search_mestre_root_tuple_scale_max100_complete.py"
    )
    complete_artifact_path = (
        root
        / "artifacts"
        / "generated-results"
        / "elliptic_mestre_root_tuple_scale_max100_complete.json"
    )
    if sha256_file(complete_script) != FROZEN_COMPLETE_SCRIPT_SHA256:
        raise AssertionError("the rank-14 calibration producer changed")
    if sha256_file(complete_artifact_path) != FROZEN_COMPLETE_ARTIFACT_SHA256:
        raise AssertionError("the rank-14 calibration artifact changed")
    complete_artifact = json.loads(complete_artifact_path.read_text())
    calibration_source = next(
        record
        for record in complete_artifact["fiber_records"]
        if record["identifier"] == "r0_6_49_73_82_96_t2"
    )

    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    if construction.quartic_condition != 0 or construction.is_reflection_symmetric:
        raise AssertionError("the pinned family geometry changed")
    if construction.quartic_square_scale != 8:
        raise AssertionError("the primitive quartic square scale changed")
    for parameter in (Q(1), Q(2), Q(7, 3), Q(29, 11), Q(-13, 7)):
        if construction.primitive_jacobian_coefficients(parameter) != family_coefficients(
            parameter
        ):
            raise AssertionError("the exact family polynomial changed")
        if family_coefficients(parameter) != family_coefficients(-parameter):
            raise AssertionError("the T to -T Jacobian quotient changed")
        if construction.primitive_quartic_coefficients(
            parameter
        ) != construction.primitive_quartic_coefficients(-parameter):
            raise AssertionError("the T to -T quartic quotient changed")
        if construction.primitive_quartic_discriminant(
            parameter
        ) != exact_discriminant(parameter):
            raise AssertionError("the factored discriminant formula changed")

    calibration_coefficients = family_coefficients(CALIBRATION_PARAMETER)
    calibration_phase = calibration_source["H5000_phase"]
    calibration_basis = tuple(
        (Q(point["jacobian_x"]), Q(point["jacobian_y"]))
        for point in calibration_phase["small_prime_saturation"]["saturated_basis"]
    )
    calibration_certificate = finite_reduction_attempt(
        calibration_coefficients,
        calibration_basis,
        prime_bound=SIGNAL_CERTIFICATE_PRIME_BOUND,
    )
    if calibration_certificate["certified_algebraic_rank_lower_bound"] != 14:
        raise AssertionError("the exact T=2 rank-14 calibration failed")
    calibration_features = discriminant_features(
        Candidate(2, 1, "0", "0", 0, 0)
    )

    started = time.monotonic()
    scan_started = time.monotonic()
    scan = run_scanner(
        scanner_source,
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
        keep=args.discovery_keep,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    scan_wall_seconds = time.monotonic() - scan_started
    if scan.prior_excluded != 8 or scan.evaluated_population < 18_000_000:
        raise AssertionError("the declared broad reduced population shrank")

    feature_started = time.monotonic()
    features = {
        candidate.identifier: discriminant_features(candidate)
        for candidate in scan.candidates
    }
    feature_records = [
        candidate_record(candidate, features[candidate.identifier])
        for candidate in scan.candidates
    ]
    selection = select_population(scan.candidates, features)
    feature_wall_seconds = time.monotonic() - feature_started
    print(
        f"scan evaluated={scan.evaluated_population} retained={len(scan.candidates)} "
        f"conductor_population={len(selection)}",
        flush=True,
    )

    conductor_started = time.monotonic()
    records: list[dict[str, Any]] = []
    runtime: dict[str, dict[str, Any]] = {}
    for position, (candidate, strata) in enumerate(selection, 1):
        coefficients = family_coefficients(candidate.parameter)
        record: dict[str, Any] = {
            **candidate_record(candidate, features[candidate.identifier]),
            "selection_strata": list(strata),
        }
        try:
            conductor = capped_minimal_curve_data(
                coefficients,
                timeout=args.conductor_timeout,
                stack_bytes=STACK_BYTES,
            )
            record["conductor_phase"] = {
                "status": "completed exact PARI minimal-model/conductor computation",
                **conductor,
                "below_strict_log_conductor_target_numerically": Decimal(
                    conductor["log_conductor"]
                )
                < TARGET_LOG_CONDUCTOR,
            }
            runtime[candidate.identifier] = {
                "candidate": candidate,
                "features": features[candidate.identifier],
                "coefficients": coefficients,
                "record": record,
            }
        except CappedProcessTimeout:
            record["conductor_phase"] = {
                "status": "timeout-no-retry",
                "timeout_seconds": args.conductor_timeout,
            }
        except Exception as error:
            record["conductor_phase"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
        records.append(record)
        if position % 32 == 0:
            print(f"conductor {position}/{len(selection)}", flush=True)
    conductor_wall_seconds = time.monotonic() - conductor_started

    point_started = time.monotonic()
    current = list(runtime.values())
    target_hits: list[dict[str, Any]] = []
    completed_stage_calls: Counter[str] = Counter()
    stage_rank_histograms: dict[str, Counter[str]] = {}
    prior_stage = ""
    for stage_name, height_bound, keep, mapping_cap in POINT_STAGES:
        if keep is not None:
            current = [
                item
                for item in current
                if item["record"]
                .get("point_stages", {})
                .get(prior_stage, {})
                .get("status", "")
                .startswith("completed")
            ]
            current.sort(key=lambda item: stage_key(item, prior_stage))
            current = current[:keep]
        stage_ranks: Counter[str] = Counter()
        for position, item in enumerate(current, 1):
            record = item["record"]
            record.setdefault("point_stages", {})
            timeout = (
                args.point_timeout if height_bound <= 50_000 else args.deep_point_timeout
            )
            try:
                phase = exact_point_stage(
                    construction,
                    item["candidate"].parameter,
                    item["coefficients"],
                    height_bound=height_bound,
                    mapping_cap=mapping_cap,
                    point_timeout=timeout,
                    height_timeout=args.height_timeout,
                    saturation_timeout=args.saturation_timeout,
                )
                record["point_stages"][stage_name] = phase
                completed_stage_calls[stage_name] += 1
                stage_ranks[str(phase_rank(phase))] += 1
                certificate = phase.get("finite_reduction_attempt", {})
                certified_rank = certificate.get(
                    "certified_algebraic_rank_lower_bound"
                )
                below_target = record["conductor_phase"][
                    "below_strict_log_conductor_target_numerically"
                ]
                if certified_rank is not None and (
                    certified_rank >= 30
                    or (certified_rank >= 21 and below_target)
                ):
                    target_hits.append(
                        {
                            "identifier": record["identifier"],
                            "parameter": record["parameter"],
                            "stage": stage_name,
                            "certified_algebraic_rank_lower_bound": certified_rank,
                            "conductor": record["conductor_phase"]["conductor"],
                            "log_conductor": record["conductor_phase"][
                                "log_conductor"
                            ],
                        }
                    )
            except CappedProcessTimeout:
                record["point_stages"][stage_name] = {
                    "status": "timeout-no-retry",
                    "height_bound": height_bound,
                    "timeout_seconds": timeout,
                }
            except Exception as error:
                record["point_stages"][stage_name] = {
                    "status": "error-no-retry",
                    "height_bound": height_bound,
                    "error": str(error)[:1000],
                }
            if stage_name == "H5000" and position % 32 == 0:
                print(f"H5000 {position}/{len(current)}", flush=True)
        stage_rank_histograms[stage_name] = stage_ranks
        print(
            f"{stage_name} attempted={len(current)} "
            f"max_rank={max((int(rank) for rank in stage_ranks), default=-1)}",
            flush=True,
        )
        prior_stage = stage_name
    point_wall_seconds = time.monotonic() - point_started

    all_completed = [
        (record, stage_name, phase)
        for record in records
        for stage_name, phase in record.get("point_stages", {}).items()
        if phase.get("status", "").startswith("completed")
    ]
    numerical_leaders = sorted(
        (
            {
                "identifier": record["identifier"],
                "parameter": record["parameter"],
                "stage": stage_name,
                "stable_numerical_rank": phase_rank(phase),
                "conductor": record["conductor_phase"],
                "finite_reduction_attempt": phase.get("finite_reduction_attempt"),
            }
            for record, stage_name, phase in all_completed
        ),
        key=lambda item: (
            -item["stable_numerical_rank"],
            Decimal(item["conductor"]["log_conductor"]),
            item["identifier"],
            item["stage"],
        ),
    )[:32]
    maximum_rank = max(
        (item["stable_numerical_rank"] for item in numerical_leaders), default=-1
    )
    stage_status_histograms = {
        stage_name: dict(
            sorted(
                Counter(
                    record.get("point_stages", {})[stage_name]["status"]
                    for record in records
                    if stage_name in record.get("point_stages", {})
                ).items()
            )
        )
        for stage_name, _, _, _ in POINT_STAGES
    }
    conductor_statuses = Counter(
        record["conductor_phase"]["status"] for record in records
    )
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "complete leakage-controlled rational-box screen in the exact rank-14 "
            "Mestre family; numerical ranks are triage only"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "family": {
            "roots": list(ROOTS),
            "primitive_quartic_square_scale": 8,
            "A_coefficients_ascending": list(A_COEFFICIENTS),
            "B_coefficients_ascending": list(B_COEFFICIENTS),
            "primitive_discriminant_factorization": {
                "constant": DISCRIMINANT_CONSTANT,
                "linear_square_factors": ["2*T-45", "2*T+45"],
                "degree_16_factor_coefficients_ascending": list(
                    DISCRIMINANT_Q16_COEFFICIENTS
                ),
            },
            "exact_symmetry": "A(-T)=A(T), B(-T)=B(T), and R_-T=R_T",
            "search_quotient": "one positive reduced representative T=a/b; T=0 omitted",
        },
        "calibration_T2": {
            "parameter": "2",
            "excluded_from_selection": True,
            "frozen_complete_artifact_sha256": FROZEN_COMPLETE_ARTIFACT_SHA256,
            "short_weierstrass_coefficients": [
                str(value) for value in calibration_coefficients
            ],
            "conductor_phase": calibration_source["conductor_phase"],
            "basis_point_sha256": point_digest(calibration_basis),
            "finite_reduction_certificate": calibration_certificate,
            "certified_algebraic_rank_lower_bound": 14,
            "exact_discriminant_features": calibration_features,
            "discovery_score": scan.calibration.discovery_score,
            "held_score": scan.calibration.held_score,
            "discovery_good_primes": scan.calibration.discovery_good,
            "held_good_primes": scan.calibration.held_good,
        },
        "modular_scan": {
            "numerator_bound": scan.numerator_bound,
            "denominator_bound": scan.denominator_bound,
            "primitive_positive_reduced_population": scan.primitive_population,
            "exact_prior_integer_parameters_excluded": scan.prior_excluded,
            "evaluated_population": scan.evaluated_population,
            "discovery_keep": scan.keep,
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "discovery_table_digest": scan.discovery_table_digest,
            "held_table_digest": scan.held_table_digest,
            "scanner_stdout_sha256": scan.stdout_sha256,
            "score_definition": (
                "sum ((2-a_p)/(p+1-a_p))*log(p), exact traces and good-reduction "
                "gates, each numerical term quantized to 10^-12"
            ),
            "retained_candidate_feature_sha256": stable_digest(feature_records),
            "retained_candidates_with_exact_features": feature_records,
        },
        "selection": {
            "fixed_before_any_conductor_call": True,
            "uses_no_conductor_point_or_rank_data": True,
            "quotas": SELECTION_QUOTAS,
            "selected_count": len(selection),
            "selected_identifier_sha256": hashlib.sha256(
                "\n".join(sorted(candidate.identifier for candidate, _ in selection)).encode()
            ).hexdigest(),
            "exact_discriminant_features_are_fixed_prime_radical_data_not_complete_factorizations": True,
        },
        "conductor_first_screen": {
            "population_closed_before_any_point_or_rank_call": True,
            "population": {
                "selected": len(selection),
                "status_histogram": dict(sorted(conductor_statuses.items())),
                "completed": len(runtime),
                "subtarget": sum(
                    record["conductor_phase"].get(
                        "below_strict_log_conductor_target_numerically", False
                    )
                    for record in records
                ),
            },
            "records": records,
        },
        "point_search_protocol": {
            "stages": [
                {
                    "name": name,
                    "height_bound": height,
                    "keep": "all conducted" if keep is None else keep,
                    "mapping_cap": cap,
                }
                for name, height, keep, cap in POINT_STAGES
            ],
            "increasing_stages_are_not_retries": True,
            "same_height_retries": 0,
            "completed_stage_calls": dict(completed_stage_calls),
            "stage_status_histograms": stage_status_histograms,
            "stable_numerical_rank_histograms": {
                name: dict(
                    sorted(histogram.items(), key=lambda item: int(item[0]))
                )
                for name, histogram in stage_rank_histograms.items()
            },
            "maximum_stable_numerical_rank": maximum_rank,
            "finite_reduction_trigger_stable_numerical_rank": SIGNAL_CERTIFICATE_TRIGGER,
            "adaptive_broadening_performed": False,
            "no_broadening_rule": (
                "the declared population is not broadened when the maximum stable rank "
                "remains at most 15"
            ),
        },
        "numerical_leaders": numerical_leaders,
        "parameters": {
            **{
                key: value
                for key, value in vars(args).items()
                if key != "output"
            },
            "stack_bytes": STACK_BYTES,
            "radical_prime_bound": 397,
            "signal_certificate_prime_bound": SIGNAL_CERTIFICATE_PRIME_BOUND,
        },
        "timings": {
            "scan_and_exact_table_replay_wall_seconds": scan_wall_seconds,
            "exact_discriminant_feature_and_selection_wall_seconds": feature_wall_seconds,
            "conductor_phase_wall_seconds": conductor_wall_seconds,
            "point_phase_wall_seconds": point_wall_seconds,
            "total_post_calibration_wall_seconds": time.monotonic() - started,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pari_gp": pari_version_capped(),
            "compiler": shutil.which("c++"),
        },
        "provenance": {
            "script": str(script.relative_to(root)),
            "script_sha256": sha256_file(script),
            "scanner": str(scanner_source.relative_to(root)),
            "scanner_sha256": sha256_file(scanner_source),
            "frozen_complete_script_sha256": FROZEN_COMPLETE_SCRIPT_SHA256,
            "frozen_complete_artifact_sha256": FROZEN_COMPLETE_ARTIFACT_SHA256,
            "reproducing_command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/search_mestre_0649738296_rational.py"
            ),
            "all_external_processes_foreground_and_capped": True,
            "whole_process_groups_killed_and_reaped_on_timeout": True,
            "no_retries": True,
            "temporary_scanner_binary_removed": True,
            "owned_processes_remaining": 0,
        },
    }
    artifact["result_sha256"] = result_digest(artifact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(args.output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "evaluated": scan.evaluated_population,
                "selected": len(selection),
                "conducted": len(runtime),
                "maximum_stable_numerical_rank": maximum_rank,
                "target_hits": target_hits,
                "result_sha256": artifact["result_sha256"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Certified T=5 baseline and rational frontier for a new Mestre family.

The family comes from the primitive normalized roots ``(0,4,30,31,39,46)``.
Its Jacobian has the exact even model

``y^2 = x^3 + A(T)x + B(T)``

with ``deg(A)=8`` and ``deg(B)=12``.  Thus ``T`` and ``-T`` define the same
curve, and the rational scan uses the exact quotient ``T=a/b>0``.

The frozen max-root-50 artifact supplies only the discovery provenance for the
``T=5`` point set.  This script independently saturates those twelve points,
checks them exactly, proves their independence by finite reductions, and
replays the exact minimal model and conductor.  It then compiles a standalone
modular scanner and exhausts a primitive ``30000 x 1000`` rational box.  The
scanner uses disjoint discovery and held prime bands; its local traces are
exact exhaustive finite-field point counts.

Candidate selection is leakage-controlled: the discovery band fixes bounded
survivor lists, the held band ranks only those survivors, and the entire final
population receives exact conductor computations before point/rank triage.
Successful fibers enter one ``H=50000`` search, with fixed smaller survivor
sets entering ``H=250000`` and ``H=1000000``.  These are increasing stages,
never retries.  Numerical height rank is triage only; stable rank at least 21
immediately triggers saturation followed by an exact finite-reduction test.
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
import re
import shutil
import sys
import tempfile
from typing import Any, Iterable, Sequence

from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from mestre_root_tuples import SixRootMestreConstruction
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    bounded_quartic_points,
    canonical_signless_points,
    capped_gp,
    capped_minimal_curve_data,
    ellrank_probe,
    height_matrix_replay,
    numerical_subset,
    pari_version_capped,
    point_digest,
    point_on_short_curve,
    point_record,
    primitive_visible_points,
    quartic_point_to_jacobian,
    run_capped_process,
    sha256_file,
)


Q = Fraction
ROOTS = (0, 4, 30, 31, 39, 46)
CALIBRATION_PARAMETER = Q(5)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
FROZEN_SCALE_ARTIFACT_SHA256 = (
    "fd2dccb1fd08aad70857df7ca19df77bd521e2be017b98f5579a748fd26cfc14"
)
FROZEN_SCALE_RESULT_SHA256 = (
    "7c3f451a92f208d241955d2500cdcf416d772e919bb54f7181f5c40fd8f53def"
)
FROZEN_SCALE_SCRIPT_SHA256 = (
    "5e7228b95ae995019fbc50b9f7667de41e06a86b4490f0feacff5702bb5cc174"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_0430313946_frontier.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_mestre_0430313946_frontier.py"
)

DISCOVERY_PRIMES = (
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)
HELD_PRIMES = (
    101,
    103,
    107,
    109,
    113,
    127,
    131,
    137,
    139,
    149,
    151,
    157,
    163,
    167,
    173,
    179,
    181,
    191,
    193,
    197,
    199,
)

# Ascending exact coefficients of A(T) and B(T).  Odd terms vanish.
A_COEFFICIENTS = (
    -4_840_084_292_061_123,
    0,
    -611_317_954_173_024,
    0,
    1_593_386_668_512,
    0,
    178_536_960,
    0,
    -4_762_800,
)
B_COEFFICIENTS = (
    -84_220_770_768_445_491_421_122,
    0,
    56_408_643_410_668_151_386_896,
    0,
    -193_116_296_662_521_546_636,
    0,
    930_870_844_049_215_872,
    0,
    -2_013_991_695_596_160,
    0,
    -224_956_569_600,
    0,
    4_000_752_000,
)


@dataclass(frozen=True)
class ScannerCandidate:
    numerator: int
    denominator: int
    discovery_score: str
    held_score: str
    discovery_good_primes: int
    held_good_primes: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)


@dataclass(frozen=True)
class ScannerResult:
    name: str
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


def evaluate_polynomial(
    coefficients: Sequence[int], value: Fraction
) -> Fraction:
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


def homogeneous_coefficient(
    coefficients: Sequence[int], numerator: int, denominator: int, prime: int
) -> int:
    degree = len(coefficients) - 1
    answer = 0
    for power, coefficient in enumerate(coefficients):
        answer += (
            coefficient
            * pow(numerator, power, prime)
            * pow(denominator, degree - power, prime)
        )
    return answer % prime


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


def exact_local_trace(parameter: Fraction, prime: int) -> int | None:
    parameter = Q(parameter)
    return exact_local_trace_projective(
        parameter.numerator, parameter.denominator, prime
    )


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
    units = 0
    good = 0
    for prime in primes:
        trace = exact_local_trace(parameter, prime)
        if trace is None:
            continue
        group_order = prime + 1 - trace
        term = (2 - trace) / group_order * log(prime)
        units += int(term * 1_000_000_000_000 + 0.5)
        good += 1
    sign = "-" if units < 0 else ""
    absolute = abs(units)
    return (
        f"{sign}{absolute // 1_000_000_000_000}."
        f"{absolute % 1_000_000_000_000:012d}",
        good,
    )


def parse_scanner_output(name: str, stdout: str) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_0430313946_SCAN_V1":
        raise AssertionError("the modular scanner omitted its format header")
    discovery = tuple(int(value) for value in lines[1].split()[1:])
    held = tuple(int(value) for value in lines[2].split()[1:])
    if lines[1].split()[0] != "D" or discovery != DISCOVERY_PRIMES:
        raise AssertionError("the modular scanner changed the discovery band")
    if lines[2].split()[0] != "H" or held != HELD_PRIMES:
        raise AssertionError("the modular scanner changed the held band")
    local_fields = lines[3].split()
    if local_fields[0] != "L" or len(local_fields) != 3:
        raise AssertionError("the modular scanner omitted its table digests")

    def candidate_from_fields(fields: Sequence[str]) -> ScannerCandidate:
        if len(fields) != 7:
            raise AssertionError("malformed modular-scanner candidate")
        return ScannerCandidate(
            numerator=int(fields[1]),
            denominator=int(fields[2]),
            discovery_score=fields[3],
            held_score=fields[4],
            discovery_good_primes=int(fields[5]),
            held_good_primes=int(fields[6]),
        )

    calibration_fields = lines[4].split()
    if calibration_fields[0] != "K":
        raise AssertionError("the modular scanner omitted T=5 calibration")
    calibration = candidate_from_fields(
        ("C", *calibration_fields[1:])
    )
    candidates = tuple(
        candidate_from_fields(line.split())
        for line in lines[5:-1]
        if line.startswith("C ")
    )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 8:
        raise AssertionError("the modular scanner omitted its summary")
    numerator_bound, denominator_bound, keep = map(int, summary[1:4])
    primitive, prior, evaluated, retained = map(int, summary[4:])
    if retained != len(candidates) or retained != min(keep, evaluated):
        raise AssertionError("the scanner retained-count summary changed")
    if any(
        gcd(candidate.numerator, candidate.denominator) != 1
        or not (1 <= candidate.numerator <= numerator_bound)
        or not (1 <= candidate.denominator <= denominator_bound)
        for candidate in candidates
    ):
        raise AssertionError("a modular candidate escaped the primitive box")
    if any(candidate.parameter in {Q(value) for value in range(1, 9)} for candidate in candidates):
        raise AssertionError("the modular scanner leaked a prior parameter")
    return ScannerResult(
        name=name,
        numerator_bound=numerator_bound,
        denominator_bound=denominator_bound,
        keep=keep,
        primitive_population=primitive,
        prior_excluded=prior,
        evaluated_population=evaluated,
        discovery_table_digest=local_fields[1],
        held_table_digest=local_fields[2],
        calibration=calibration,
        candidates=candidates,
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )


def run_scanner_strata(
    source: Path,
    *,
    compiler: str,
    compile_timeout: float,
    scan_timeout: float,
    denominator_bound: int,
    strata: Sequence[tuple[str, int, int]],
) -> tuple[ScannerResult, ...]:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError(f"C++ compiler {compiler!r} was not found")
    results: list[ScannerResult] = []
    with tempfile.TemporaryDirectory(prefix="mestre-0430313946-") as directory:
        binary = Path(directory) / "scanner"
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
        for name, numerator_bound, keep in strata:
            stdout, _ = run_capped_process(
                (
                    str(binary),
                    str(numerator_bound),
                    str(denominator_bound),
                    str(keep),
                ),
                timeout=scan_timeout,
            )
            results.append(parse_scanner_output(name, stdout))
    if len({result.discovery_table_digest for result in results}) != 1:
        raise AssertionError("discovery lookup tables changed between strata")
    if len({result.held_table_digest for result in results}) != 1:
        raise AssertionError("held lookup tables changed between strata")
    if results[0].discovery_table_digest != exact_local_table_digest(
        DISCOVERY_PRIMES
    ):
        raise AssertionError("the C++ discovery trace table missed exact Python replay")
    if results[0].held_table_digest != exact_local_table_digest(HELD_PRIMES):
        raise AssertionError("the C++ held trace table missed exact Python replay")
    expected_discovery = quantized_score(CALIBRATION_PARAMETER, DISCOVERY_PRIMES)
    expected_held = quantized_score(CALIBRATION_PARAMETER, HELD_PRIMES)
    for result in results:
        if result.calibration.discovery_good_primes != expected_discovery[1] or abs(
            Decimal(result.calibration.discovery_score)
            - Decimal(expected_discovery[0])
        ) > Decimal("0.000000000001"):
            raise AssertionError("the C++ discovery calibration missed exact Python traces")
        if result.calibration.held_good_primes != expected_held[1] or abs(
            Decimal(result.calibration.held_score) - Decimal(expected_held[0])
        ) > Decimal("0.000000000001"):
            raise AssertionError("the C++ held calibration missed exact Python traces")
    return tuple(results)


def parse_gp_points(text: str) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(
        (Q(x_value), Q(y_value))
        for x_value, y_value in re.findall(
            r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", text
        )
    )


def saturate_basis(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    prime_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a saturation input point missed the exact curve")
    curve = ",".join(f"({value})" for value in coefficients)
    point_vector = ",".join(f"[({x_value}),({y_value})]" for x_value, y_value in points)
    program = "\n".join(
        (
            "default(realprecision,100);",
            f"E=ellinit([{curve}]);P=[{point_vector}];",
            f"S=ellsaturation(E,P,{prime_bound});",
            'print("COUNT_BEGIN");print(#S);print("COUNT_END");',
            'print("POINTS_BEGIN");print(S);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    stdout, _ = capped_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    count_match = re.search(r"COUNT_BEGIN\n(\d+)\nCOUNT_END", stdout)
    point_match = re.search(r"POINTS_BEGIN\n(.*?)\nPOINTS_END", stdout, re.DOTALL)
    if count_match is None or point_match is None:
        raise AssertionError("PARI omitted saturation output")
    saturated = parse_gp_points(point_match.group(1))
    if len(saturated) != int(count_match.group(1)):
        raise AssertionError("the saturated point vector parsed incompletely")
    if any(not point_on_short_curve(coefficients, point) for point in saturated):
        raise AssertionError("a saturated point missed the exact curve")
    return saturated, {
        "status": "completed PARI saturation proposal; exact independence is separate",
        "input_point_count": len(points),
        "returned_point_count": len(saturated),
        "prime_bound_strict_upper_limit": prime_bound,
        "saturated_point_sha256": point_digest(saturated),
        "saturated_basis": [point_record(point) for point in saturated],
    }


def finite_certificate(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    saturation_bound: int,
    saturation_timeout: float,
    certificate_prime_bound: int,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    saturated, saturation = saturate_basis(
        coefficients,
        points,
        prime_bound=saturation_bound,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, saturated, prime_bound=certificate_prime_bound
    )
    exact_rank = combined_mod2_rank(signatures, len(saturated))
    certified = exact_rank == len(saturated)
    two_torsion_prime = (
        find_two_torsion_certificate_prime(coefficients, prime_bound=200)
        if certified
        else None
    )
    return saturated, {
        "status": "certified" if certified else "bounded-search-rank-deficient",
        "saturation": saturation,
        "point_count": len(saturated),
        "point_sha256": point_digest(saturated),
        "certificate_prime_bound": certificate_prime_bound,
        "certificate_primes": [signature.prime for signature in signatures],
        "combined_exact_rank_over_F2": exact_rank,
        "two_torsion_certificate_prime": two_torsion_prime,
        "signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "certified_algebraic_rank_lower_bound": len(saturated) if certified else None,
        "proof_semantics": (
            "full F2 rank in a product of E(F_p)/2E(F_p), together with the "
            "displayed irreducible 2-division-cubic prime, proves Z-independence"
        ),
    }


def frozen_scale_calibration(path: Path) -> dict[str, Any]:
    if sha256_file(path) != FROZEN_SCALE_ARTIFACT_SHA256:
        raise AssertionError("the frozen max-root-50 artifact hash changed")
    artifact = json.loads(path.read_text())
    if artifact["result_sha256"] != FROZEN_SCALE_RESULT_SHA256:
        raise AssertionError("the frozen max-root-50 result digest changed")
    if artifact["provenance"]["script_sha256"] != FROZEN_SCALE_SCRIPT_SHA256:
        raise AssertionError("the frozen max-root-50 producer hash changed")
    prior_parameters = {
        Q(record["parameter"])
        for record in artifact["specialization_screen"][
            "conductor_first_fiber_records"
        ]
        if tuple(record["roots"]) == ROOTS
    }
    prior_parameters.update(
        Q(record["parameter"])
        for record in artifact["specialization_screen"]["inadmissible_fibers"]
        if tuple(record["roots"]) == ROOTS
    )
    if prior_parameters != {Q(value) for value in range(1, 9)}:
        raise AssertionError("the exact prior parameter manifest changed")
    record = next(
        item
        for item in artifact["specialization_screen"]["point_search_finalists"]
        if tuple(item["roots"]) == ROOTS
        and Q(item["parameter"]) == CALIBRATION_PARAMETER
    )
    escalation = record["single_strongest_signal_escalation"]
    points = tuple(
        (Q(item["x"]), Q(item["y"])) for item in escalation["numerical_subset"]
    )
    if len(points) != 12:
        raise AssertionError("the frozen T=5 numerical subset changed")
    return {
        "artifact": artifact,
        "record": record,
        "points": points,
        "prior_parameters": tuple(sorted(prior_parameters)),
    }


def select_conductor_population(
    scans: Sequence[ScannerResult], quotas: dict[str, int]
) -> tuple[tuple[ScannerCandidate, tuple[str, ...]], ...]:
    selected: dict[Fraction, tuple[ScannerCandidate, set[str]]] = {}
    for scan in scans:
        ranked = sorted(
            scan.candidates,
            key=lambda candidate: (
                -Decimal(candidate.held_score),
                -Decimal(candidate.discovery_score),
                candidate.denominator,
                candidate.numerator,
            ),
        )
        added = 0
        for candidate in ranked:
            if candidate.parameter in selected:
                selected[candidate.parameter][1].add(scan.name)
                continue
            selected[candidate.parameter] = (candidate, {scan.name})
            added += 1
            if added == quotas[scan.name]:
                break
        if added != quotas[scan.name]:
            raise AssertionError(f"stratum {scan.name} could not fill its quota")
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
    point_timeout: float,
    height_timeout: float,
    ellrank_timeout: float,
    stack_bytes: int,
    mapping_cap: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    visible_quartic = primitive_visible_points(construction, parameter)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point)
        for point in visible_quartic
    )
    raw = bounded_quartic_points(
        construction.primitive_quartic_coefficients(parameter),
        height_bound=height_bound,
        timeout=point_timeout,
        stack_bytes=stack_bytes,
    )
    signless = canonical_signless_points(raw)
    retained = signless[:mapping_cap]
    searched_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point) for point in retained
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
        stack_bytes=stack_bytes,
    )
    subset = numerical_subset(pool, height)
    try:
        ellrank = ellrank_probe(
            coefficients,
            subset,
            timeout=ellrank_timeout,
            stack_bytes=stack_bytes,
        )
    except CappedProcessTimeout:
        ellrank = {"status": "timeout", "timeout_seconds": ellrank_timeout}
    except Exception as error:
        ellrank = {"status": "error", "error": str(error)[:1000]}
    return {
        "status": "completed",
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
        "stable_numerical_rank": int(height[-1]["numerical_rank"]),
        "numerical_subset": [point_record(point) for point in subset],
        "effort_zero_ellrank": ellrank,
        "numerical_rank_is_not_an_independence_certificate": True,
    }, subset


def scanner_record(scan: ScannerResult) -> dict[str, Any]:
    return {
        "name": scan.name,
        "numerator_bound": scan.numerator_bound,
        "denominator_bound": scan.denominator_bound,
        "keep": scan.keep,
        "primitive_population": scan.primitive_population,
        "prior_excluded": scan.prior_excluded,
        "evaluated_population": scan.evaluated_population,
        "discovery_table_digest": scan.discovery_table_digest,
        "held_table_digest": scan.held_table_digest,
        "stdout_sha256": scan.stdout_sha256,
        "calibration": {
            "parameter": "5",
            "discovery_score": scan.calibration.discovery_score,
            "held_score": scan.calibration.held_score,
            "discovery_good_primes": scan.calibration.discovery_good_primes,
            "held_good_primes": scan.calibration.held_good_primes,
            "exactly_excluded_from_population": True,
        },
        "retained_candidates": [
            {
                "numerator": candidate.numerator,
                "denominator": candidate.denominator,
                "parameter": str(candidate.parameter),
                "discovery_score": candidate.discovery_score,
                "held_score": candidate.held_score,
                "discovery_good_primes": candidate.discovery_good_primes,
                "held_good_primes": candidate.held_good_primes,
            }
            for candidate in scan.candidates
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scale-artifact",
        type=Path,
        default=root
        / "artifacts/generated-results/elliptic_mestre_root_tuple_scale.json",
    )
    parser.add_argument("--numerator-bound", type=int, default=30_000)
    parser.add_argument("--denominator-bound", type=int, default=1_000)
    parser.add_argument("--global-keep", type=int, default=4_096)
    parser.add_argument("--medium-numerator-bound", type=int, default=5_000)
    parser.add_argument("--medium-keep", type=int, default=2_048)
    parser.add_argument("--low-numerator-bound", type=int, default=1_000)
    parser.add_argument("--low-keep", type=int, default=1_024)
    parser.add_argument("--conductor-global-quota", type=int, default=48)
    parser.add_argument("--conductor-medium-quota", type=int, default=48)
    parser.add_argument("--conductor-low-quota", type=int, default=32)
    parser.add_argument("--h50000-keep", type=int, default=24)
    parser.add_argument("--h250000-keep", type=int, default=6)
    parser.add_argument("--h1000000-keep", type=int, default=2)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--h50000-timeout", type=float, default=20.0)
    parser.add_argument("--h250000-timeout", type=float, default=30.0)
    parser.add_argument("--h1000000-timeout", type=float, default=45.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--ellrank-timeout", type=float, default=15.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-bound", type=int, default=20)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.numerator_bound < 30_000 or args.numerator_bound > 100_000:
        raise SystemExit("--numerator-bound must lie in [30000,100000]")
    if args.denominator_bound < 1_000 or args.denominator_bound > 5_000:
        raise SystemExit("--denominator-bound must lie in [1000,5000]")
    if not (1 <= args.low_numerator_bound <= args.medium_numerator_bound <= args.numerator_bound):
        raise SystemExit("low/medium numerator bounds must be nested")
    keeps = (args.global_keep, args.medium_keep, args.low_keep)
    if min(keeps) < 128 or max(keeps) > 20_000:
        raise SystemExit("scanner keep counts must lie in [128,20000]")
    quotas = (
        args.conductor_global_quota,
        args.conductor_medium_quota,
        args.conductor_low_quota,
    )
    if min(quotas) < 1 or sum(quotas) > 256:
        raise SystemExit("conductor quotas must be positive and total at most 256")
    if not (1 <= args.h1000000_keep <= args.h250000_keep <= args.h50000_keep <= sum(quotas)):
        raise SystemExit("point-stage keep counts must be nested")
    if args.mapping_cap < 12 or args.mapping_cap > 1024:
        raise SystemExit("--mapping-cap must lie in [12,1024]")
    timeouts = (
        args.compile_timeout,
        args.scan_timeout,
        args.conductor_timeout,
        args.h50000_timeout,
        args.h250000_timeout,
        args.h1000000_timeout,
        args.height_timeout,
        args.ellrank_timeout,
        args.saturation_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all subprocess timeouts must lie in (0,60]")
    if args.stack_bytes < 8_000_000 or args.stack_bytes > 1_000_000_000:
        raise SystemExit("--stack-bytes must lie in [8000000,1000000000]")
    if args.saturation_bound < 3 or args.saturation_bound > 100:
        raise SystemExit("--saturation-bound must lie in [3,100]")
    if args.certificate_prime_bound < 53 or args.certificate_prime_bound > 2_000:
        raise SystemExit("--certificate-prime-bound must lie in [53,2000]")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    scanner_path = script_path.with_name("scan_mestre_0430313946.cpp")
    scale_script = script_path.with_name("search_mestre_root_tuple_scale.py")
    output = args.output if args.output.is_absolute() else repo_root / args.output
    if sha256_file(scale_script) != FROZEN_SCALE_SCRIPT_SHA256:
        raise AssertionError("the frozen scale helper script changed")

    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    if construction.quartic_condition != 0 or construction.is_reflection_symmetric:
        raise AssertionError("the pinned family geometry changed")
    if construction.quartic_square_scale != 40:
        raise AssertionError("the pinned primitive quartic scale changed")
    for parameter in (Q(1), Q(5), Q(7, 3), Q(-11, 4)):
        if construction.primitive_jacobian_coefficients(parameter) != family_coefficients(
            parameter
        ):
            raise AssertionError("the pinned A(T),B(T) polynomial changed")
        if family_coefficients(parameter) != family_coefficients(-parameter):
            raise AssertionError("the exact T -> -T quotient failed")
        if construction.primitive_quartic_coefficients(
            parameter
        ) != construction.primitive_quartic_coefficients(-parameter):
            raise AssertionError("the primitive quartic is not even in T")

    frozen = frozen_scale_calibration(args.scale_artifact)
    calibration_coefficients = family_coefficients(CALIBRATION_PARAMETER)
    saturated, calibration_certificate = finite_certificate(
        calibration_coefficients,
        frozen["points"],
        saturation_bound=args.saturation_bound,
        saturation_timeout=args.saturation_timeout,
        certificate_prime_bound=args.certificate_prime_bound,
        stack_bytes=args.stack_bytes,
    )
    if calibration_certificate["certified_algebraic_rank_lower_bound"] != 12:
        raise AssertionError("finite reductions did not certify the rank-12 baseline")
    conductor_replay = capped_minimal_curve_data(
        calibration_coefficients,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    frozen_conductor = frozen["record"]["conductor"]
    for key in ("minimal_model", "conductor", "root_number"):
        if conductor_replay[key] != frozen_conductor[key]:
            raise AssertionError(f"the T=5 conductor replay changed at {key}")
    calibration_record = {
        "roots": list(ROOTS),
        "parameter": "5",
        "short_weierstrass_coefficients": [str(value) for value in calibration_coefficients],
        "input_numerical_subset_from_frozen_artifact": {
            "point_count": len(frozen["points"]),
            "point_sha256": point_digest(frozen["points"]),
            "selection_is_not_used_as_proof": True,
        },
        "finite_reduction_certificate": calibration_certificate,
        "certified_algebraic_rank_lower_bound": 12,
        "conductor_replay": conductor_replay,
        "exact_log_conductor_bound": exact_log_conductor_certificate(
            int(conductor_replay["conductor"])
        ),
        "strictly_below_log_conductor_target": True,
        "excluded_from_every_rational_scan_stratum": True,
    }
    print(
        "T=5 exact rank>=12 certificate and conductor replay complete",
        flush=True,
    )

    scans = run_scanner_strata(
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
    if scans[0].primitive_population < 18_000_000:
        raise AssertionError("the declared broad primitive population shrank")
    if any(scan.prior_excluded != 8 for scan in scans):
        raise AssertionError("the exact prior-exclusion count changed")
    quotas = {
        "global": args.conductor_global_quota,
        "medium": args.conductor_medium_quota,
        "low": args.conductor_low_quota,
    }
    selected = select_conductor_population(scans, quotas)
    print(
        f"scanner closed global primitive population={scans[0].primitive_population}; "
        f"conductor population={len(selected)}",
        flush=True,
    )

    conductor_records: list[dict[str, Any]] = []
    runtime: dict[str, dict[str, Any]] = {}
    for position, (candidate, sources) in enumerate(selected, 1):
        parameter = candidate.parameter
        identifier = f"t{candidate.numerator}_{candidate.denominator}"
        coefficients = family_coefficients(parameter)
        record: dict[str, Any] = {
            "identifier": identifier,
            "parameter": str(parameter),
            "numerator": candidate.numerator,
            "denominator": candidate.denominator,
            "selection_strata": list(sources),
            "discovery_score": candidate.discovery_score,
            "held_score": candidate.held_score,
            "discovery_good_primes": candidate.discovery_good_primes,
            "held_good_primes": candidate.held_good_primes,
        }
        if parameter in frozen["prior_parameters"]:
            raise AssertionError("a prior parameter escaped exact exclusion")
        if construction.quartic_discriminant(parameter) == 0:
            record["conductor_phase"] = {
                "status": "not attempted: exact singular specialization"
            }
        else:
            try:
                conductor = capped_minimal_curve_data(
                    coefficients,
                    timeout=args.conductor_timeout,
                    stack_bytes=args.stack_bytes,
                )
                record["conductor_phase"] = {
                    "status": "completed exact PARI minimal-model/conductor computation",
                    **conductor,
                    "below_strict_log_conductor_target_numerically": (
                        Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                    ),
                }
                runtime[identifier] = {
                    "candidate": candidate,
                    "parameter": parameter,
                    "coefficients": coefficients,
                    "record": record,
                }
            except CappedProcessTimeout:
                record["conductor_phase"] = {
                    "status": "timeout",
                    "timeout_seconds": args.conductor_timeout,
                }
            except Exception as error:
                record["conductor_phase"] = {
                    "status": "error",
                    "error": str(error)[:1000],
                }
        conductor_records.append(record)
        if position % 16 == 0:
            print(
                f"conductors {position}/{len(selected)} complete={len(runtime)}",
                flush=True,
            )

    # The full selected conductor population is now closed.  Point/rank data
    # has not been computed for any scan candidate before this boundary.
    conductor_population_closed = True
    conducted = [item for item in runtime.values()]
    conducted.sort(
        key=lambda item: (
            not item["record"]["conductor_phase"][
                "below_strict_log_conductor_target_numerically"
            ],
            -Decimal(item["candidate"].held_score),
            -Decimal(item["candidate"].discovery_score),
            item["candidate"].denominator,
            item["candidate"].numerator,
        )
    )
    tiers = (
        ("H50000", 50_000, args.h50000_keep, args.h50000_timeout),
        ("H250000", 250_000, args.h250000_keep, args.h250000_timeout),
        ("H1000000", 1_000_000, args.h1000000_keep, args.h1000000_timeout),
    )
    current = conducted[: args.h50000_keep]
    target_hits: list[dict[str, Any]] = []
    for tier_index, (tier_name, height_bound, keep, point_timeout) in enumerate(tiers):
        if tier_index:
            completed_previous = [
                item
                for item in current
                if item["record"].get("point_stages", {}).get(
                    tiers[tier_index - 1][0], {}
                ).get("status")
                == "completed"
            ]
            completed_previous.sort(
                key=lambda item: (
                    -int(
                        item["record"]["point_stages"][tiers[tier_index - 1][0]][
                            "stable_numerical_rank"
                        ]
                    ),
                    -Decimal(item["candidate"].held_score),
                    -Decimal(item["candidate"].discovery_score),
                    item["candidate"].denominator,
                    item["candidate"].numerator,
                )
            )
            current = completed_previous[:keep]
        for item in current:
            record = item["record"]
            record.setdefault("point_stages", {})
            try:
                stage, subset = exact_point_stage(
                    construction,
                    item["parameter"],
                    item["coefficients"],
                    height_bound=height_bound,
                    point_timeout=point_timeout,
                    height_timeout=args.height_timeout,
                    ellrank_timeout=args.ellrank_timeout,
                    stack_bytes=args.stack_bytes,
                    mapping_cap=args.mapping_cap,
                )
                record["point_stages"][tier_name] = stage
                if int(stage["stable_numerical_rank"]) >= 21:
                    _, certificate = finite_certificate(
                        item["coefficients"],
                        subset,
                        saturation_bound=args.saturation_bound,
                        saturation_timeout=args.saturation_timeout,
                        certificate_prime_bound=args.certificate_prime_bound,
                        stack_bytes=args.stack_bytes,
                    )
                    stage["finite_reduction_attempt"] = certificate
                    certified_rank = certificate[
                        "certified_algebraic_rank_lower_bound"
                    ]
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
                                "tier": tier_name,
                                "certified_rank_lower_bound": certified_rank,
                                "conductor": record["conductor_phase"]["conductor"],
                                "log_conductor": record["conductor_phase"][
                                    "log_conductor"
                                ],
                            }
                        )
                else:
                    stage["finite_reduction_attempt"] = {
                        "status": "not triggered",
                        "trigger_stable_numerical_rank": 21,
                    }
            except CappedProcessTimeout:
                record["point_stages"][tier_name] = {
                    "status": "timeout",
                    "timeout_seconds": point_timeout,
                    "no_retry_at_same_height": True,
                }
            except Exception as error:
                record["point_stages"][tier_name] = {
                    "status": "error",
                    "error": str(error)[:1000],
                    "no_retry_at_same_height": True,
                }
        print(
            f"{tier_name} stage attempted={len(current)} "
            f"max_rank={max((item['record']['point_stages'][tier_name].get('stable_numerical_rank', -1) for item in current), default=-1)}",
            flush=True,
        )

    completed_stage_records = [
        stage
        for record in conductor_records
        for stage in record.get("point_stages", {}).values()
        if stage.get("status") == "completed"
    ]
    source_paths = {
        "script": script_path,
        "scanner": scanner_path,
        "frozen_scale_script": scale_script,
        "frozen_scale_artifact": args.scale_artifact,
    }
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "exact rank-12 baseline certificate plus bounded leakage-controlled "
            "rational specialization frontier; numerical ranks are triage only"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
            "explanation": (
                "no exact finite-reduction target certificate was produced"
                if not target_hits
                else "at least one exact finite-reduction target certificate was produced"
            ),
        },
        "family": {
            "roots": list(ROOTS),
            "quartic_condition": str(construction.quartic_condition),
            "reflection_symmetric": construction.is_reflection_symmetric,
            "removed_quartic_square_scale": str(construction.quartic_square_scale),
            "A_coefficients_ascending": list(A_COEFFICIENTS),
            "B_coefficients_ascending": list(B_COEFFICIENTS),
            "exact_symmetry": "A(-T)=A(T), B(-T)=B(T), and primitive R_-T=R_T",
            "search_quotient": "one representative T=a/b>0; T=0 excluded",
        },
        "frozen_scale_boundary": {
            "artifact": str(args.scale_artifact),
            "artifact_sha256": FROZEN_SCALE_ARTIFACT_SHA256,
            "result_sha256": FROZEN_SCALE_RESULT_SHA256,
            "producer_script_sha256": FROZEN_SCALE_SCRIPT_SHA256,
            "exact_prior_parameters": [str(value) for value in frozen["prior_parameters"]],
            "exact_prior_parameter_sha256": hashlib.sha256(
                "\n".join(str(value) for value in frozen["prior_parameters"]).encode()
            ).hexdigest(),
        },
        "calibration_T5": calibration_record,
        "modular_scan": {
            "score": (
                "sum ((2-a_p)/(p+1-a_p))*log(p) over exact good local traces; "
                "each term quantized to 10^-12 before summation"
            ),
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "strata": [scanner_record(scan) for scan in scans],
            "global_box_exhausted": {
                "numerator": [1, scans[0].numerator_bound],
                "denominator": [1, scans[0].denominator_bound],
                "primitive_positive_rationals": scans[0].primitive_population,
                "exact_prior_excluded": scans[0].prior_excluded,
                "evaluated": scans[0].evaluated_population,
            },
            "conductor_selection_quotas": quotas,
        },
        "conductor_first_screen": {
            "population_closed_before_point_or_rank_triage": conductor_population_closed,
            "selected_population": len(selected),
            "completed": sum(
                record["conductor_phase"]["status"].startswith("completed")
                for record in conductor_records
            ),
            "timeouts": sum(
                record["conductor_phase"]["status"] == "timeout"
                for record in conductor_records
            ),
            "errors": sum(
                record["conductor_phase"]["status"] == "error"
                for record in conductor_records
            ),
            "singular_rejections": sum(
                record["conductor_phase"]["status"].startswith("not attempted")
                for record in conductor_records
            ),
            "records": conductor_records,
        },
        "point_search_protocol": {
            "stages": [
                {"name": "H50000", "height_bound": 50_000, "keep": args.h50000_keep},
                {"name": "H250000", "height_bound": 250_000, "keep": args.h250000_keep},
                {"name": "H1000000", "height_bound": 1_000_000, "keep": args.h1000000_keep},
            ],
            "increasing_stages_are_not_retries": True,
            "same_height_retries": 0,
            "mapping_cap": args.mapping_cap,
            "height_precisions": [72, 120],
            "finite_reduction_trigger_rank": 21,
            "maximum_stable_numerical_rank": max(
                (int(stage["stable_numerical_rank"]) for stage in completed_stage_records),
                default=None,
            ),
            "completed_stage_calls": len(completed_stage_records),
        },
        "parameters": {
            key: value
            for key, value in vars(args).items()
            if key not in {"scale_artifact", "output"}
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pari_gp": pari_version_capped(),
            "compiler": shutil.which(args.compiler),
        },
        "provenance": {
            **{
                f"{name}_path": str(path.relative_to(repo_root))
                for name, path in source_paths.items()
            },
            **{
                f"{name}_sha256": sha256_file(path)
                for name, path in source_paths.items()
            },
            "reproducing_command": REPRODUCING_COMMAND,
            "temporary_scanner_binary_removed": True,
            "all_external_calls_are_foreground_process_groups": True,
            "whole_process_group_killed_and_reaped_on_timeout": True,
            "same_stage_retries": 0,
        },
    }
    digest_payload = {
        "target": artifact["target"],
        "family": artifact["family"],
        "frozen_scale_boundary": artifact["frozen_scale_boundary"],
        "calibration_T5": artifact["calibration_T5"],
        "modular_scan": artifact["modular_scan"],
        "conductor_first_screen": artifact["conductor_first_screen"],
        "point_search_protocol": artifact["point_search_protocol"],
    }
    artifact["result_sha256"] = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"max_rank={artifact['point_search_protocol']['maximum_stable_numerical_rank']} "
        f"target_hits={len(target_hits)} wrote={output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

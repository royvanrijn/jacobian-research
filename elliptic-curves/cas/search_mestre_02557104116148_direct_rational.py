#!/usr/bin/env python3
"""Leakage-controlled direct rational-T scan in Nagao's rank-13 root family.

This lane searches the six-root Mestre family ``(0,25,57,104,116,148)``
directly in its even parameter ``T``.  It is deliberately distinct from the
older quadratic base-change searches in ``u``.  Every positive reduced
``T=a/b`` in the pinned 30,000 by 1,000 box is considered once, except for
``T=1,...,8`` and every exact ``parameter_t`` discoverable in the frozen
Nagao rank-13 artifacts.  The latter exclusion is reconstructed and hashed at
runtime, so the populations are auditable rather than described informally.

Fresh discovery primes 601..701 close three survivor tails before the
disjoint held band 709..809 is read.  Exact degree-20 discriminant/radical
features then close a fixed conductor population.  All conductor calls finish
before any point search.  Exact quartic searches run at H=5,000, 50,000,
250,000, and 1,000,000 on fixed nested tiers.  Stable numerical rank at least
18 immediately triggers an exact mod-3 finite-reduction attempt.  Numerical
height rank remains triage evidence, and the search does not broaden when the
fixed protocol remains at rank at most 15.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd
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
import search_mestre_rank14_pair_rational_frontier as engine
from search_mestre_root_tuple_scale import (
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

ROOTS = (0, 25, 57, 104, 116, 148)
FAMILY_LABEL = "r0_25_57_104_116_148"
CALIBRATION_PARAMETER = Q(1)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DISCOVERY_PRIMES = (
    601, 607, 613, 617, 619, 631, 641, 643, 647,
    653, 659, 661, 673, 677, 683, 691, 701,
)
HELD_PRIMES = (
    709, 719, 727, 733, 739, 743, 751,
    757, 761, 769, 773, 787, 797, 809,
)
A_COEFFICIENTS = (
    -27_546_462_334_108_146_267,
    0,
    -1_530_834_958_134_000,
    0,
    3_106_634_557_536,
    0,
    -359_251_200,
    0,
    -34_992,
)
B_COEFFICIENTS = (
    55_624_396_621_360_883_431_446_459_126,
    0,
    5_508_337_118_494_541_464_141_200,
    0,
    -13_510_101_631_102_695_979_884,
    0,
    4_296_246_530_145_998_400,
    0,
    -634_269_830_133_888,
    0,
    38_799_129_600,
    0,
    2_519_424,
)
FAMILY = engine.FamilySpec(
    0,
    FAMILY_LABEL,
    ROOTS,
    CALIBRATION_PARAMETER,
    A_COEFFICIENTS,
    B_COEFFICIENTS,
)
FINITE_REDUCTION_TRIGGER = 18
TRIAL_DIVISION_LIMIT = 997
STACK_BYTES = 512_000_000
SELECTION_QUOTAS = {
    "highest-held-score": 80,
    "smallest-exact-radical-upper-bound": 64,
    "largest-exact-known-powerful-part": 32,
    "balanced-held-and-power-rank": 32,
}
POINT_STAGES = (
    ("H5000", 5_000, None),
    ("H50000", 50_000, 32),
    ("H250000", 250_000, 8),
    ("H1000000", 1_000_000, 2),
)

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

# These are the complete frozen artifact-name pattern
# elliptic_nagao_rank13*.json union elliptic_nagao_u*.json at lane creation.
PRIOR_ARTIFACT_SHA256 = {
    "elliptic_nagao_rank13_finalist_triage.json": "ba936299554ab578e61b732996983a8106f9bf76f0a66399f01ba073170b259b",
    "elliptic_nagao_rank13_generation3.json": "924879ff9e9202473dab7813870889a18df86a54aae89653a0473168935caefe",
    "elliptic_nagao_rank13_integer_u.json": "c238b04983c637359838f3a874e4695bf6b470dbc7e7c71fff5d78b9a2a90535",
    "elliptic_nagao_rank13_integer_u2000.json": "4153ebe30b72c8f583be0df8324f0f9320c9ed2186729d2a0d88a3e63243da51",
    "elliptic_nagao_rank13_local_candidate_triage.json": "9dd967422ed5b0cb14e85d673389ec2e8256c7089b5562d7e1186a4b8e6f6860",
    "elliptic_nagao_rank13_local_crt.json": "272fbf6eb8b227dc05b5ea82ebf0a8a3057bd786937ca1b917060d9212e16894",
    "elliptic_nagao_rank13_rank_gain_mutations.json": "b039dba6b629182f56ab8ee8777f19f0fb674696a3411f390213aa824f6dd0a4",
    "elliptic_nagao_rank13_rank_gain_search.json": "5f55e0f35368760c65dc2ee66da8edee9ca1153f88e70123e1844d4e2bade559",
    "elliptic_nagao_u118_height_1000000.json": "f526328e4710b4c7a2742982e6746dd902ec3d19bdbd69d72ded5cdc3a137f06",
    "elliptic_nagao_u135_alternate_covers.json": "c47c91155588ccaf6f6606f7f76115aedb3b61e32a4c1907b44fbc5ae4e1c83b",
    "elliptic_nagao_u135_ell2cover.json": "61566039b04722164dda4b2b02d459c7426cde8ac452e6ae4ea9dcc068914630",
    "elliptic_nagao_u135_skew_height.json": "4182d87746ac1903e0e67babb21398ed35a90d769dbecf7159d129aa88a0f001",
    "elliptic_nagao_u42_descent_toolchain.json": "902f83f22b1685a1f19f3c550079c2d4c2eea48fd58c69ae3b2361af5f4b1fcf",
    "elliptic_nagao_u42_height_10000000.json": "4fea0207fd637988bcc1147143657cbec5c2404cb81b4c4a487e2dde20cc43b8",
    "elliptic_nagao_u42_magma_probe.json": "5fefc8e6f0a4b9353e58542087849d73620748f71e2ed3620520741ae03bd8e3",
    "elliptic_nagao_u42_rank17_certificate.json": "ef942e81ed0a2e7bd1bfa8d5bf3f549b80ad71a100cc55333390aad614056437",
    "elliptic_nagao_u42_skew_height.json": "e0abb85d618cde90d3681efb6cef2ea4e7f22a89d59fcc5eda882f17f652b140",
    "elliptic_nagao_u471_11_alternate_covers.json": "3edbc41aa1faff001af37c6979a78b8b9ca12a4f15f6d9154019146457805152",
    "elliptic_nagao_u75_alternate_covers.json": "0cc7cd65f8adb5a01f770b3f8b383ff1a993bff3dc9aeb410a3d1508b6d5dc5c",
}
EXPECTED_PRIOR_PARAMETER_T_ABSOLUTE_COUNT = 13_119
EXPECTED_PRIOR_PARAMETER_T_BOX_COUNT = 660
EXPECTED_EXCLUSION_COUNT = 668
EXPECTED_EXCLUSION_SHA256 = (
    "4d8356b8c9031093fd0cb2a3c3f8eef83ecb3b7ae770b855c13a53926bcf4dc8"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_02557104116148_direct_rational.json"
)


def configure_engine(exclusions: Iterable[Fraction]) -> None:
    """Bind the generic frozen rank-14 rational APIs to this one family."""

    engine.FAMILIES = (FAMILY,)
    engine.DISCOVERY_PRIMES = DISCOVERY_PRIMES
    engine.HELD_PRIMES = HELD_PRIMES
    engine.PRIOR_PANEL_PARAMETERS = tuple(exclusions)
    engine.FINITE_REDUCTION_TRIGGER = FINITE_REDUCTION_TRIGGER
    engine.TARGET_LOG_CONDUCTOR = TARGET_LOG_CONDUCTOR


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + coefficient
    return answer


def family_coefficients(parameter: Fraction) -> tuple[Fraction, ...]:
    parameter = Q(parameter)
    return (
        Q(0), Q(0), Q(0),
        evaluate_polynomial(A_COEFFICIENTS, parameter),
        evaluate_polynomial(B_COEFFICIENTS, parameter),
    )


def _collect_parameter_t(value: Any, answer: set[Fraction]) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "parameter_t" and isinstance(child, (str, int)):
                try:
                    parameter = abs(Q(child))
                except (ValueError, ZeroDivisionError):
                    pass
                else:
                    if parameter:
                        answer.add(parameter)
            _collect_parameter_t(child, answer)
    elif isinstance(value, list):
        for child in value:
            _collect_parameter_t(child, answer)


def exclusion_lines(exclusions: Iterable[Fraction]) -> str:
    ordered = sorted(set(map(Q, exclusions)), key=lambda value: (
        value.denominator, value.numerator
    ))
    return "".join(
        f"{value.numerator}\t{value.denominator}\n" for value in ordered
    )


def reconstruct_exclusions(
    generated: Path, numerator_bound: int, denominator_bound: int
) -> tuple[tuple[Fraction, ...], dict[str, Any]]:
    observed_files = {
        path.name
        for pattern in ("elliptic_nagao_rank13*.json", "elliptic_nagao_u*.json")
        for path in generated.glob(pattern)
    }
    if observed_files != set(PRIOR_ARTIFACT_SHA256):
        raise AssertionError("the discoverable prior Nagao artifact set changed")
    all_parameters: set[Fraction] = set()
    source_records = []
    for name, expected_sha in sorted(PRIOR_ARTIFACT_SHA256.items()):
        path = generated / name
        actual_sha = sha256_file(path)
        if actual_sha != expected_sha:
            raise AssertionError(f"frozen prior artifact changed: {name}")
        local: set[Fraction] = set()
        _collect_parameter_t(json.loads(path.read_text()), local)
        all_parameters.update(local)
        source_records.append(
            {
                "artifact": name,
                "sha256": actual_sha,
                "unique_absolute_parameter_t_count": len(local),
                "inside_declared_box": sum(
                    value.numerator <= numerator_bound
                    and value.denominator <= denominator_bound
                    for value in local
                ),
            }
        )
    inside = {
        value
        for value in all_parameters
        if value.numerator <= numerator_bound
        and value.denominator <= denominator_bound
    }
    panel = {Q(value) for value in range(1, 9)}
    exclusions = tuple(sorted(inside | panel))
    text = exclusion_lines(exclusions)
    digest = hashlib.sha256(text.encode()).hexdigest()
    if (
        len(all_parameters) != EXPECTED_PRIOR_PARAMETER_T_ABSOLUTE_COUNT
        or len(inside) != EXPECTED_PRIOR_PARAMETER_T_BOX_COUNT
        or len(exclusions) != EXPECTED_EXCLUSION_COUNT
        or inside & panel
        or digest != EXPECTED_EXCLUSION_SHA256
    ):
        raise AssertionError("the exact prior direct-T exclusion census changed")
    return exclusions, {
        "method": (
            "recursive exact Fraction extraction of every `parameter_t`, "
            "canonicalized by absolute value because A and B are even"
        ),
        "source_artifacts": source_records,
        "unique_absolute_parameter_t_all_heights": len(all_parameters),
        "unique_absolute_parameter_t_inside_declared_box": len(inside),
        "panel_T_1_through_8_added": 8,
        "panel_overlap_with_artifact_parameters": 0,
        "total_excluded_inside_declared_box": len(exclusions),
        "canonical_exclusion_lines_sha256": digest,
    }


ScannerCandidate = engine.ScannerCandidate
ScannerResult = engine.ScannerResult


def _parse_candidate(fields: Sequence[str]) -> ScannerCandidate:
    if len(fields) != 7:
        raise AssertionError("malformed direct-T scanner candidate")
    return ScannerCandidate(
        int(fields[1]), int(fields[2]), fields[3], fields[4],
        int(fields[5]), int(fields[6]),
    )


def parse_scanner_output(
    stratum: str, stdout: str, exclusions: set[Fraction]
) -> ScannerResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_02557104116148_DIRECT_SCAN_V1":
        raise AssertionError("the direct-T scanner omitted its header")
    if tuple(map(int, lines[1].split()[1:])) != DISCOVERY_PRIMES:
        raise AssertionError("the fresh discovery band changed")
    if tuple(map(int, lines[2].split()[1:])) != HELD_PRIMES:
        raise AssertionError("the held band changed")
    local = lines[3].split()
    calibration = _parse_candidate(("C", *lines[4].split()[1:]))
    candidates = tuple(
        _parse_candidate(line.split()) for line in lines[5:-1]
    )
    summary = lines[-1].split()
    if local[0] != "L" or len(local) != 3 or summary[0] != "S" or len(summary) != 9:
        raise AssertionError("the scanner population/digest footer changed")
    numerator_bound, denominator_bound, keep = map(int, summary[1:4])
    primitive, prior, evaluated, retained, exclusion_file_count = map(
        int, summary[4:]
    )
    if retained != len(candidates) or retained != min(keep, evaluated):
        raise AssertionError("the scanner retained-count gate changed")
    if exclusion_file_count != len(exclusions):
        raise AssertionError("the scanner read an incomplete exclusion file")
    if calibration.parameter != CALIBRATION_PARAMETER:
        raise AssertionError("the T=1 calibration changed")
    if any(
        gcd(candidate.numerator, candidate.denominator) != 1
        or not (1 <= candidate.numerator <= numerator_bound)
        or not (1 <= candidate.denominator <= denominator_bound)
        or candidate.parameter in exclusions
        for candidate in candidates
    ):
        raise AssertionError("a survivor escaped the reduced/excluded box")
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
        family_index=0,
        stratum=stratum,
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


def run_scanners(
    source: Path,
    exclusions: Sequence[Fraction],
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
    exclusion_set = set(exclusions)
    results = []
    with tempfile.TemporaryDirectory(prefix="mestre-02557104116148-") as directory:
        directory_path = Path(directory)
        binary = directory_path / "scanner"
        exclusion_path = directory_path / "exclusions.tsv"
        exclusion_path.write_text(exclusion_lines(exclusions))
        run_capped_process(
            (
                executable, "-std=c++17", "-O3", "-DNDEBUG",
                str(source), "-o", str(binary),
            ),
            timeout=compile_timeout,
        )
        for name, numerator_bound, keep in strata:
            stdout, _ = run_capped_process(
                (
                    str(binary), str(numerator_bound), str(denominator_bound),
                    str(keep), str(exclusion_path),
                ),
                timeout=scan_timeout,
            )
            results.append(parse_scanner_output(name, stdout, exclusion_set))
    return tuple(results)


def select_conductor_population(
    pool: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    orders = {
        "highest-held-score": sorted(
            pool,
            key=lambda row: (
                -Decimal(row["held_score"]), -Decimal(row["discovery_score"]),
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
    }
    held_rank = {
        row["parameter"]: index
        for index, row in enumerate(orders["highest-held-score"])
    }
    power_rank = {
        row["parameter"]: index
        for index, row in enumerate(orders["largest-exact-known-powerful-part"])
    }
    orders["balanced-held-and-power-rank"] = sorted(
        pool,
        key=lambda row: (
            held_rank[row["parameter"]] + power_rank[row["parameter"]],
            max(held_rank[row["parameter"]], power_rank[row["parameter"]]),
            row["denominator"], row["numerator"],
        ),
    )
    selected: dict[str, dict[str, Any]] = {}
    reasons: dict[str, set[str]] = {}
    for label, quota in SELECTION_QUOTAS.items():
        added = 0
        for row in orders[label]:
            if row["parameter"] in selected:
                reasons[row["parameter"]].add(label)
                continue
            selected[row["parameter"]] = row
            reasons[row["parameter"]] = {label}
            added += 1
            if added == quota:
                break
        if added != quota:
            raise AssertionError(f"selection quota failed to fill: {label}")
    answer = []
    for parameter, row in selected.items():
        record = dict(row)
        record["conductor_selection_strata"] = sorted(reasons[parameter])
        answer.append(record)
    answer.sort(key=lambda row: (row["denominator"], row["numerator"]))
    digest = hashlib.sha256()
    for row in answer:
        digest.update(
            f"{row['parameter']}|{','.join(row['conductor_selection_strata'])}\n".encode()
        )
    return answer, {
        "selection_uses_conductor": False,
        "selection_uses_point_or_rank_data": False,
        "discovery_survivors_closed_before_held_scores": True,
        "held_scores_rank_only_discovery_survivors": True,
        "exact_discriminant_features_use_no_conductor_or_rank_data": True,
        "novel_candidates_added_per_stratum": SELECTION_QUOTAS,
        "selected_population": len(answer),
        "selected_population_sha256": digest.hexdigest(),
    }


def scanner_record(scan: ScannerResult) -> dict[str, Any]:
    return {
        "stratum": scan.stratum,
        "numerator_bound": scan.numerator_bound,
        "denominator_bound": scan.denominator_bound,
        "keep": scan.keep,
        "primitive_population": scan.primitive_population,
        "prior_specializations_excluded": scan.prior_excluded,
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
        default=root / DEFAULT_OUTPUT,
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.numerator_bound != 30_000 or args.denominator_bound != 1_000:
        raise SystemExit("the direct rational box is pinned at 30000 by 1000")
    if (
        args.global_keep != 4_096
        or args.medium_numerator_bound != 5_000
        or args.medium_keep != 2_048
        or args.low_numerator_bound != 1_000
        or args.low_keep != 1_024
    ):
        raise SystemExit("the three discovery tails are pinned")
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
    if not 211 <= args.certificate_prime_bound <= 2000:
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
        raise SystemExit("refusing to overwrite the direct-T artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    generated = root / "artifacts/generated-results"
    scanner_path = script_path.with_name(
        "scan_mestre_02557104116148_direct_rational.cpp"
    )
    max200_path = generated / FROZEN_MAX200_ARTIFACT
    max200_script = script_path.with_name(
        "search_mestre_root_tuple_scale_max200.py"
    )
    if (
        sha256_file(max200_path) != EXPECTED_MAX200_ARTIFACT_SHA256
        or sha256_file(max200_script) != EXPECTED_MAX200_SCRIPT_SHA256
    ):
        raise AssertionError("the frozen max-root-200 inputs changed")
    max200 = json.loads(max200_path.read_text())
    if max200["result_sha256"] != EXPECTED_MAX200_RESULT_SHA256:
        raise AssertionError("the frozen max-root-200 result digest changed")

    exclusions, exclusion_audit = reconstruct_exclusions(
        generated, args.numerator_bound, args.denominator_bound
    )
    configure_engine(exclusions)
    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    if construction.quartic_condition or construction.is_reflection_symmetric:
        raise AssertionError("the selected family geometry changed")
    for parameter in (Q(1), Q(7), Q(11, 5), Q(-17, 3)):
        if construction.primitive_jacobian_coefficients(parameter) != family_coefficients(
            parameter
        ):
            raise AssertionError("the pinned A(T),B(T) formula changed")
        if family_coefficients(parameter) != family_coefficients(-parameter):
            raise AssertionError("the T <-> -T quotient changed")
    raw_discriminant = construction.primitive_discriminant_polynomial
    content = engine.polynomial_content(raw_discriminant)
    discriminant = tuple(value.numerator // content for value in raw_discriminant)

    calibration = engine.calibration_records(max200)
    if len(calibration) != 1 or calibration[0][
        "certified_algebraic_rank_lower_bound"
    ] != 14:
        raise AssertionError("the frozen rank-14 calibration changed")

    exact_discovery_digest = engine.exact_table_digest(0, DISCOVERY_PRIMES)
    exact_held_digest = engine.exact_table_digest(0, HELD_PRIMES)
    scans = run_scanners(
        scanner_path,
        exclusions,
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
    exclusion_set = set(exclusions)
    for scan in scans:
        expected_prior = sum(
            value.numerator <= scan.numerator_bound
            and value.denominator <= scan.denominator_bound
            for value in exclusion_set
        )
        if (
            scan.discovery_table_digest != exact_discovery_digest
            or scan.held_table_digest != exact_held_digest
            or scan.prior_excluded != expected_prior
        ):
            raise AssertionError("a scanner local-table/exclusion gate changed")
        expected_discovery = engine.score_text(0, CALIBRATION_PARAMETER, DISCOVERY_PRIMES)
        expected_held = engine.score_text(0, CALIBRATION_PARAMETER, HELD_PRIMES)
        if (
            abs(Decimal(scan.calibration.discovery_score) - Decimal(expected_discovery[0]))
            > Decimal("0.000000000100")
            or scan.calibration.discovery_good != expected_discovery[1]
            or abs(Decimal(scan.calibration.held_score) - Decimal(expected_held[0]))
            > Decimal("0.000000000100")
            or scan.calibration.held_good != expected_held[1]
        ):
            raise AssertionError("the T=1 exact local-score replay failed")
    global_scan = next(scan for scan in scans if scan.stratum == "global")
    if (
        global_scan.primitive_population != 18_244_819
        or global_scan.prior_excluded != 668
        or global_scan.evaluated_population != 18_244_151
    ):
        raise AssertionError("the declared global direct-T population changed")
    print(
        "fresh-prime direct-T scan closed: "
        f"primitive={global_scan.primitive_population} "
        f"excluded={global_scan.prior_excluded} "
        f"evaluated={global_scan.evaluated_population}",
        flush=True,
    )

    pool, pool_audit = engine.pool_and_features(
        FAMILY, scans, discriminant
    )
    selected, selection = select_conductor_population(pool)
    print(
        f"exact discriminant features closed; conductor population={len(selected)}",
        flush=True,
    )

    common: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "completed selection-only direct-T rational scan"
            if args.selection_only
            else "in-progress conductor-first direct-T rational scan"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "scope": {
            "included_roots": list(ROOTS),
            "direct_parameter": "primitive positive T=a/b",
            "exact_sign_quotient": "A(-T)=A(T), B(-T)=B(T)",
            "distinct_from_u_base_change": True,
            "declared_box": {
                "numerator": [1, args.numerator_bound],
                "denominator": [1, args.denominator_bound],
            },
        },
        "family": {
            "label": FAMILY_LABEL,
            "roots": list(ROOTS),
            "quartic_square_scale": str(construction.quartic_square_scale),
            "A_coefficients_ascending": list(A_COEFFICIENTS),
            "B_coefficients_ascending": list(B_COEFFICIENTS),
            "content_free_discriminant_coefficients_ascending": list(discriminant),
            "removed_discriminant_polynomial_content": str(content),
        },
        "frozen_calibration_T1": calibration[0],
        "prior_specialization_exclusion_audit": exclusion_audit,
        "modular_scan": {
            "score": (
                "sum ((2-a_p)/(p+1-a_p))*log(p), each exact trace term "
                "quantized to 1e-12"
            ),
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_relative_to_all_predecessor_bands_through_prime_599": True,
            "strata": [scanner_record(scan) for scan in scans],
        },
        "exact_discriminant_feature_screen": {
            "trial_division_prime_bound": TRIAL_DIVISION_LIMIT,
            "content_free_homogeneous_degree": 20,
            "pool_audit": pool_audit,
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
            "temporary_scanner_and_exclusion_file_removed": True,
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
                "scope": common["scope"],
                "family": common["family"],
                "calibration": common["frozen_calibration_T1"],
                "exclusions": common["prior_specialization_exclusion_audit"],
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
                engine.conductor_worker,
                0, row["numerator"], row["denominator"],
                args.conductor_timeout, args.stack_bytes,
            )
            for row in selected
        ]
        for position, (row, future) in enumerate(zip(selected, futures), start=1):
            row["conductor_phase"] = future.result()
            if position % 32 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    eligible = [
        row
        for row in selected
        if row["conductor_phase"]["status"].startswith("completed")
    ]

    stage_timeouts = {
        "H5000": args.h5000_timeout,
        "H50000": args.h50000_timeout,
        "H250000": args.h250000_timeout,
        "H1000000": args.h1000000_timeout,
    }
    current = eligible
    for stage_index, (name, height, keep) in enumerate(POINT_STAGES):
        if stage_index:
            prior_name = POINT_STAGES[stage_index - 1][0]
            current = [
                row
                for row in current
                if row.get("point_stages", {}).get(prior_name, {}).get("status")
                == "completed"
            ]
            current.sort(key=lambda row: engine.stage_rank_key(row, prior_name))
            current = current[: int(keep)]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    engine.point_stage_worker,
                    0, row["numerator"], row["denominator"],
                    height, stage_timeouts[name], args.height_timeout,
                    args.ellrank_timeout, args.stack_bytes, args.mapping_cap,
                    args.certificate_prime_bound,
                )
                for row in current
            ]
            for row, future in zip(current, futures):
                row.setdefault("point_stages", {})[name] = future.result()
        maximum = max(
            (
                row["point_stages"][name].get("stable_numerical_rank", -1)
                for row in current
            ),
            default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={maximum}", flush=True)
        if maximum >= 15:
            print(f"EARLY_SIGNAL {name} stable_rank={maximum}", flush=True)

    completed_stages = [
        stage
        for row in selected
        for stage in row.get("point_stages", {}).values()
        if stage.get("status") == "completed"
    ]
    maximum_rank = max(
        (int(stage["stable_numerical_rank"]) for stage in completed_stages),
        default=None,
    )
    finite_attempts = []
    target_hits = []
    for row in selected:
        for stage_name, stage in row.get("point_stages", {}).items():
            certificate = stage.get("finite_reduction_attempt", {})
            certified = certificate.get("certified_algebraic_rank_lower_bound")
            if certified is None:
                continue
            finite_attempts.append(
                {
                    "parameter": row["parameter"],
                    "stage": stage_name,
                    "certified_rank_lower_bound": certified,
                }
            )
            below = row["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically", False
            )
            if certified >= 30 or (certified >= 21 and below):
                target_hits.append(
                    {
                        "parameter": row["parameter"],
                        "stage": stage_name,
                        "certified_rank_lower_bound": certified,
                        "conductor": row["conductor_phase"]["conductor"],
                        "log_conductor": row["conductor_phase"]["log_conductor"],
                    }
                )

    common["status"] = (
        "completed fixed direct-T rational scan; stopped without broadening"
    )
    common["target"]["hits"] = target_hits
    common["conductor_first_screen"] = {
        "population_closed_before_any_point_or_rank_call": True,
        "selected_population": len(selected),
        "completed": len(eligible),
        "timeouts": sum(
            row["conductor_phase"]["status"] == "timeout" for row in selected
        ),
        "errors": sum(
            row["conductor_phase"]["status"] == "error" for row in selected
        ),
        "subtarget": sum(
            row["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically"
            ) is True
            for row in selected
        ),
    }
    common["point_search_protocol"] = {
        "stages": [
            {
                "name": name,
                "height_bound": height,
                "keep_after_previous_stage": keep,
                "attempted": sum(
                    name in row.get("point_stages", {}) for row in selected
                ),
            }
            for name, height, keep in POINT_STAGES
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
            "scope": common["scope"],
            "family": common["family"],
            "calibration": common["frozen_calibration_T1"],
            "exclusions": common["prior_specialization_exclusion_audit"],
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

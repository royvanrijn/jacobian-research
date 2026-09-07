#!/usr/bin/env python3
"""Leakage-controlled neighboring annulus around the exact T=490/9 lead.

The sole parameter family is the even Mestre family with roots
``(0,7,121,128,183,194)``.  The population is the primitive ordinary annulus

``4097 <= b <= 16000, a=nearest(490*b/9)+delta, 5 <= |delta| <= 16``.

It is disjoint from the companion neighborhood search by construction: its
two exact Farey rays are removed and every possible reduced denominator of
its two CRT grids is conservatively removed.  Fresh discovery primes
587--659 close an 8192-candidate survivor set before scores at the disjoint
661--757 band or any exact discriminant feature are used.  A fixed union of
held-score, radical, powerful-part, offset-diversity, and low-denominator
quotas closes before every conductor call.  Exact quartic/height stages then
run at H=5k, 50k, 250k, and 1m with strict one-pass caps.  Stable rank at
least 15 immediately receives an exact finite-reduction attempt.
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
from math import gcd, isqrt
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
from search_mestre_rank14_pair_rational_frontier import (
    discriminant_feature,
    exact_table_digest,
    family_coefficients,
    polynomial_content,
    score_text,
)
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

ROOTS = (0, 7, 121, 128, 183, 194)
ANCHOR = Q(490, 9)
DENOMINATOR_MIN = 4_097
DENOMINATOR_MAX = 16_000
OFFSETS = tuple(range(-16, -4)) + tuple(range(5, 17))
DISCOVERY_PRIMES = (
    587, 593, 599, 601, 607, 613, 617,
    619, 631, 641, 643, 647, 653, 659,
)
HELD_PRIMES = (
    661, 673, 677, 683, 691, 701, 709,
    719, 727, 733, 739, 743, 751, 757,
)
SURVIVOR_COUNT = 8_192
TRIGGER_RANK = 15
TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXPECTED_SCANNER_SHA256 = (
    "545bee442bedf6750f51e36d727e1202a7e3477ac1de9a515f551930e4ba4479"
)
RANK15_CERTIFICATE_FILENAME = "elliptic_mestre_rank15_490_9.json"
EXPECTED_RANK15_CERTIFICATE_SHA256 = (
    "50b2b9c8bd24bcb5533534446af6404f3a9a761b5f33e0e28e04dc572227f950"
)
EXPECTED_RANK15_CERTIFICATE_RESULT_SHA256 = (
    "31ff386e802c368efdcdd027168883e8edd459921eff431de71bfcfa5e517648"
)
EXPECTED_PARAMETER_MANIFEST_SHA256 = (
    "6c2fd81b95bb0d7c0531bbccaeff6bcb7fecd1353420febab32492dd04aef4f2"
)
EXPECTED_FAREY_MANIFEST_SHA256 = (
    "f4cf84c3882563743fb6069cd7699752e7b4cbb401850696d0fb1394d4c28f00"
)
EXPECTED_GRID_DIVISOR_MANIFEST_SHA256 = (
    "5fe65a2ce0d54db5ba7065afb7f53045a2882651e281a45c191bed7cce8bb0ea"
)


@dataclass(frozen=True)
class AnnulusCandidate:
    numerator: int
    denominator: int
    offset: int
    discovery_score: str
    held_score: str
    discovery_good: int
    held_good: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)


@dataclass(frozen=True)
class ScanResult:
    candidates: tuple[AnnulusCandidate, ...]
    discovery_table_digest: str
    held_table_digest: str
    raw_population: int
    nonprimitive: int
    grid_excluded: int
    farey_excluded: int
    evaluated: int
    stdout_sha256: str


def nearest_anchor_numerator(denominator: int) -> int:
    # The reduced fractional denominator is 9, so a half-integer tie cannot
    # occur.  Adding 4 implements nearest-integer rounding exactly.
    return (490 * denominator + 4) // 9


def peer_farey_manifest() -> set[Fraction]:
    return {
        value
        for m in range(61, 2_049)
        for value in (
            Q(490 * m + 381, 9 * m + 7),
            Q(490 * m + 109, 9 * m + 2),
        )
    }


def peer_grid_divisors() -> set[int]:
    answer: set[int] = set()
    for j in range(256):
        for raw in (552 + 37 * j, 553 + 41 * j):
            for divisor in range(1, isqrt(raw) + 1):
                if raw % divisor:
                    continue
                answer.add(divisor)
                answer.add(raw // divisor)
    return answer


def population_audit() -> dict[str, Any]:
    farey = peer_farey_manifest()
    divisors = peer_grid_divisors()
    digest = hashlib.sha256()
    counts = Counter()
    for denominator in range(DENOMINATOR_MIN, DENOMINATOR_MAX + 1):
        nearest = nearest_anchor_numerator(denominator)
        for offset in OFFSETS:
            counts["raw"] += 1
            numerator = nearest + offset
            if gcd(numerator, denominator) != 1:
                counts["nonprimitive"] += 1
                continue
            if denominator in divisors:
                counts["grid_divisor_excluded"] += 1
                continue
            if Q(numerator, denominator) in farey:
                counts["farey_excluded"] += 1
                continue
            counts["evaluated"] += 1
            digest.update(f"{numerator}/{denominator}|{offset}\n".encode())
    farey_digest = hashlib.sha256()
    for parameter in sorted(farey):
        farey_digest.update(f"{parameter}\n".encode())
    relevant_divisors = sorted(
        value
        for value in divisors
        if DENOMINATOR_MIN <= value <= DENOMINATOR_MAX
    )
    divisor_digest = hashlib.sha256()
    for value in relevant_divisors:
        divisor_digest.update(f"{value}\n".encode())
    record = {
        "denominator_interval": [DENOMINATOR_MIN, DENOMINATOR_MAX],
        "offsets": list(OFFSETS),
        "raw_population": counts["raw"],
        "nonprimitive_rejections": counts["nonprimitive"],
        "peer_grid_divisor_rejections": counts["grid_divisor_excluded"],
        "peer_farey_rejections": counts["farey_excluded"],
        "evaluated_unique_primitive_parameters": counts["evaluated"],
        "ordered_parameter_manifest_sha256": digest.hexdigest(),
        "peer_farey_manifest_count": len(farey),
        "peer_farey_manifest_sha256": farey_digest.hexdigest(),
        "peer_grid_divisor_count_within_annulus": len(relevant_divisors),
        "peer_grid_divisor_manifest_sha256": divisor_digest.hexdigest(),
        "prior_rectangle_exclusion": (
            "primitive b>=4097 is disjoint from every reduced prior parameter "
            "with denominator<=1000"
        ),
        "anchor_excluded": ANCHOR not in {
            Q(nearest_anchor_numerator(b) + offset, b)
            for b in range(DENOMINATOR_MIN, DENOMINATOR_MAX + 1)
            for offset in OFFSETS
            if gcd(nearest_anchor_numerator(b) + offset, b) == 1
        },
    }
    expected = {
        "raw_population": 285_696,
        "nonprimitive_rejections": 112_027,
        "peer_grid_divisor_rejections": 5_572,
        "peer_farey_rejections": 0,
        "evaluated_unique_primitive_parameters": 168_097,
    }
    for key, value in expected.items():
        if record[key] != value:
            raise AssertionError(f"the annulus population count {key} changed")
    if record["ordered_parameter_manifest_sha256"] != EXPECTED_PARAMETER_MANIFEST_SHA256:
        raise AssertionError("the ordered annulus manifest changed")
    if record["peer_farey_manifest_sha256"] != EXPECTED_FAREY_MANIFEST_SHA256:
        raise AssertionError("the peer Farey exclusion manifest changed")
    if record["peer_grid_divisor_manifest_sha256"] != EXPECTED_GRID_DIVISOR_MANIFEST_SHA256:
        raise AssertionError("the peer grid-divisor exclusion manifest changed")
    if not record["anchor_excluded"]:
        raise AssertionError("the exact T=490/9 anchor escaped the annulus exclusion")
    return record


def parse_scan(stdout: str) -> ScanResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_RANK15_ANNULUS_SCAN_V1":
        raise AssertionError("the annulus scanner omitted its header")
    if tuple(map(int, lines[1].split()[1:])) != DISCOVERY_PRIMES:
        raise AssertionError("the discovery prime band changed")
    if tuple(map(int, lines[2].split()[1:])) != HELD_PRIMES:
        raise AssertionError("the held prime band changed")
    local = lines[3].split()
    candidates = []
    for line in lines[4:-1]:
        fields = line.split()
        if fields[0] != "C" or len(fields) != 8:
            raise AssertionError("a malformed annulus candidate was emitted")
        candidates.append(
            AnnulusCandidate(
                int(fields[1]), int(fields[2]), int(fields[3]), fields[4],
                fields[5], int(fields[6]), int(fields[7]),
            )
        )
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 10:
        raise AssertionError("the annulus scanner summary changed")
    if (int(summary[1]), int(summary[2]), int(summary[3])) != (
        DENOMINATOR_MIN, DENOMINATOR_MAX, SURVIVOR_COUNT
    ):
        raise AssertionError("the scanner bounds changed")
    result = ScanResult(
        candidates=tuple(candidates),
        discovery_table_digest=local[1],
        held_table_digest=local[2],
        raw_population=int(summary[4]),
        nonprimitive=int(summary[5]),
        grid_excluded=int(summary[6]),
        farey_excluded=int(summary[7]),
        evaluated=int(summary[8]),
        stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
    )
    if len(candidates) != int(summary[9]) or len(candidates) != SURVIVOR_COUNT:
        raise AssertionError("the annulus retained count changed")
    if len({candidate.parameter for candidate in candidates}) != len(candidates):
        raise AssertionError("the annulus scanner emitted a duplicate parameter")
    return result


def run_scan(
    source: Path, *, compiler: str, compile_timeout: float, scan_timeout: float
) -> ScanResult:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="mestre-rank15-annulus-") as directory:
        binary = Path(directory) / "scan"
        run_capped_process(
            (
                executable, "-std=c++17", "-O3", "-DNDEBUG",
                str(source), "-o", str(binary),
            ),
            timeout=compile_timeout,
        )
        stdout, _ = run_capped_process(
            (str(binary), str(SURVIVOR_COUNT)), timeout=scan_timeout
        )
    return parse_scan(stdout)


def exact_feature_pool(
    scan: ScanResult, discriminant_coefficients: Sequence[int]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    pool = []
    digest = hashlib.sha256()
    singular = 0
    for candidate in scan.candidates:
        feature = discriminant_feature(
            discriminant_coefficients,
            candidate.numerator,
            candidate.denominator,
        )
        record = {
            "numerator": candidate.numerator,
            "denominator": candidate.denominator,
            "parameter": str(candidate.parameter),
            "anchor_offset": candidate.offset,
            "discovery_score": candidate.discovery_score,
            "held_score": candidate.held_score,
            "discovery_good": candidate.discovery_good,
            "held_good": candidate.held_good,
            "discriminant_feature": feature,
        }
        digest.update(
            (
                f"{candidate.parameter}|{candidate.offset}|"
                f"{candidate.discovery_score}|{candidate.held_score}|"
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
        "survivor_count": len(scan.candidates),
        "singular_rejections": singular,
        "admissible_count": len(pool),
        "feature_population_sha256": digest.hexdigest(),
    }


def select_conductor_population(
    pool: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reasons: dict[str, set[str]] = defaultdict(set)
    by_parameter = {record["parameter"]: record for record in pool}
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
        "lowest-denominator": sorted(
            pool,
            key=lambda row: (
                row["denominator"], abs(row["anchor_offset"]), row["numerator"],
            ),
        ),
    }
    quotas = {
        "highest-held-score": 32,
        "smallest-exact-radical-upper-bound": 32,
        "largest-exact-known-powerful-part": 16,
        "lowest-denominator": 8,
    }
    for label, order in orders.items():
        for record in order[: quotas[label]]:
            reasons[record["parameter"]].add(label)
    for absolute_offset in range(5, 17):
        order = sorted(
            (
                row for row in pool
                if abs(int(row["anchor_offset"])) == absolute_offset
            ),
            key=lambda row: (
                -Decimal(row["held_score"]), -Decimal(row["discovery_score"]),
                row["denominator"], row["numerator"],
            ),
        )
        for record in order[:2]:
            reasons[record["parameter"]].add("two-per-absolute-offset")
    selected = []
    for parameter, labels in reasons.items():
        record = dict(by_parameter[parameter])
        record["conductor_selection_strata"] = sorted(labels)
        selected.append(record)
    selected.sort(key=lambda row: (row["denominator"], row["numerator"]))
    digest = hashlib.sha256()
    for record in selected:
        digest.update(
            (
                f"{record['parameter']}|{record['anchor_offset']}|"
                f"{','.join(record['conductor_selection_strata'])}\n"
            ).encode()
        )
    return selected, {
        "selection_uses_conductor": False,
        "selection_uses_point_or_rank_data": False,
        "discovery_survivors_closed_before_held_scores": True,
        "exact_discriminant_features_use_no_conductor_or_rank_data": True,
        "fixed_quotas": quotas,
        "absolute_offset_diversity_quota_each": 2,
        "selected_population": len(selected),
        "selected_population_sha256": digest.hexdigest(),
    }


def conductor_worker(
    numerator: int, denominator: int, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    coefficients = family_coefficients(1, Q(numerator, denominator))
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


def point_worker(
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
    construction = SixRootMestreConstruction(tuple(map(Q, ROOTS)))
    parameter = Q(numerator, denominator)
    coefficients = family_coefficients(1, parameter)
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
        if int(stage["stable_numerical_rank"]) >= TRIGGER_RANK:
            stage["finite_reduction_attempt"] = mod3_independence_certificate(
                coefficients, subset, prime_bound=certificate_prime_bound
            )
        else:
            stage["finite_reduction_attempt"] = {
                "status": "not triggered",
                "trigger_stable_numerical_rank": TRIGGER_RANK,
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


def stage_key(record: dict[str, Any], prior: str) -> tuple[Any, ...]:
    stage = record["point_stages"][prior]
    conductor = record["conductor_phase"]
    return (
        -int(stage["stable_numerical_rank"]),
        not conductor["below_strict_log_conductor_target_numerically"],
        -Decimal(record["held_score"]),
        int(record["discriminant_feature"]["combined_radical_upper_bound"]),
        record["denominator"], record["numerator"],
    )


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=12.0)
    parser.add_argument("--h5000-timeout", type=float, default=15.0)
    parser.add_argument("--h50000-timeout", type=float, default=20.0)
    parser.add_argument("--h250000-timeout", type=float, default=30.0)
    parser.add_argument("--h1000000-timeout", type=float, default=45.0)
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--ellrank-timeout", type=float, default=8.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=499)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root / "artifacts/generated-results/elliptic_mestre_rank15_annulus.json"
        ),
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    timeouts = (
        args.compile_timeout, args.scan_timeout, args.conductor_timeout,
        args.h5000_timeout, args.h50000_timeout, args.h250000_timeout,
        args.h1000000_timeout, args.height_timeout, args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    if not 32 <= args.mapping_cap <= 1024:
        raise SystemExit("mapping cap must lie in [32,1024]")
    if not 211 <= args.certificate_prime_bound <= 2_000:
        raise SystemExit("certificate prime bound must lie in [211,2000]")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the annulus artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    scanner_path = script_path.with_name("scan_mestre_rank15_annulus.cpp")
    certificate_path = (
        root / "artifacts/generated-results" / RANK15_CERTIFICATE_FILENAME
    )
    if sha256_file(scanner_path) != EXPECTED_SCANNER_SHA256:
        raise SystemExit("the pinned annulus scanner changed")
    if sha256_file(certificate_path) != EXPECTED_RANK15_CERTIFICATE_SHA256:
        raise SystemExit("the pinned rank-15 anchor certificate changed")
    anchor_certificate = json.loads(certificate_path.read_text())
    if (
        anchor_certificate["result_sha256"]
        != EXPECTED_RANK15_CERTIFICATE_RESULT_SHA256
        or tuple(anchor_certificate["curve"]["roots"]) != ROOTS
        or Q(anchor_certificate["curve"]["parameter"]) != ANCHOR
        or anchor_certificate["claim"]["certified_algebraic_rank_lower_bound"] != 15
    ):
        raise AssertionError("the exact rank-15 anchor certificate changed")

    population = population_audit()
    print(
        f"annulus manifest closed: evaluated="
        f"{population['evaluated_unique_primitive_parameters']}", flush=True
    )
    scan = run_scan(
        scanner_path,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    if (
        scan.raw_population != population["raw_population"]
        or scan.nonprimitive != population["nonprimitive_rejections"]
        or scan.grid_excluded != population["peer_grid_divisor_rejections"]
        or scan.farey_excluded != population["peer_farey_rejections"]
        or scan.evaluated != population["evaluated_unique_primitive_parameters"]
    ):
        raise AssertionError("the independent scanner/population counts disagree")
    python_discovery_digest = exact_table_digest(1, DISCOVERY_PRIMES)
    python_held_digest = exact_table_digest(1, HELD_PRIMES)
    if (
        scan.discovery_table_digest != python_discovery_digest
        or scan.held_table_digest != python_held_digest
    ):
        raise AssertionError("the C++/Python exact local tables disagree")
    anchor_discovery = score_text(1, ANCHOR, DISCOVERY_PRIMES)
    anchor_held = score_text(1, ANCHOR, HELD_PRIMES)
    retained_digest = hashlib.sha256()
    for candidate in scan.candidates:
        retained_digest.update(
            (
                f"{candidate.parameter}|{candidate.offset}|"
                f"{candidate.discovery_score}|{candidate.held_score}|"
                f"{candidate.discovery_good}|{candidate.held_good}\n"
            ).encode()
        )
    print("fresh local tables and 8192 survivors replayed exactly", flush=True)

    construction = SixRootMestreConstruction(tuple(map(Q, ROOTS)))
    if construction.primitive_jacobian_coefficients(ANCHOR) != family_coefficients(
        1, ANCHOR
    ):
        raise AssertionError("the family model changed at the rank-15 anchor")
    raw_discriminant = construction.primitive_discriminant_polynomial
    content = polynomial_content(raw_discriminant)
    discriminant_coefficients = tuple(
        value.numerator // content for value in raw_discriminant
    )
    pool, feature_audit = exact_feature_pool(scan, discriminant_coefficients)
    selected, selection = select_conductor_population(pool)
    print(
        f"exact features closed; conductor population={len(selected)}", flush=True
    )

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "selection closed; downstream phases pending",
        "scope": {
            "family_roots": list(ROOTS),
            "T_sign_quotient": "primitive positive T=a/b; T and -T are identical",
            "anchor_parameter": str(ANCHOR),
            "anchor_excluded": True,
            "population": population,
            "companion_search_nonoverlap": {
                "ordinary_low_denominator_lane": "b>=4097 separates peer b<=4096 lane",
                "farey_rays": "all 3976 exact canonical values explicitly excluded",
                "CRT_grids": (
                    "every possible canonical denominator dividing the peer raw "
                    "C/D grid denominators is excluded"
                ),
            },
        },
        "anchor_certificate": {
            "path": str(certificate_path.relative_to(root)),
            "sha256": EXPECTED_RANK15_CERTIFICATE_SHA256,
            "result_sha256": EXPECTED_RANK15_CERTIFICATE_RESULT_SHA256,
            "certified_algebraic_rank_lower_bound": 15,
            "root_number": anchor_certificate["curve"]["root_number"],
            "log_conductor": anchor_certificate["curve"]["log_conductor"],
        },
        "fresh_local_scan": {
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": True,
            "bands_disjoint_from_companion_search_through_prime_577": True,
            "discovery_table_digest": scan.discovery_table_digest,
            "held_table_digest": scan.held_table_digest,
            "Python_table_replay_matches": True,
            "anchor_discovery_score": anchor_discovery[0],
            "anchor_discovery_good_primes": anchor_discovery[1],
            "anchor_held_score": anchor_held[0],
            "anchor_held_good_primes": anchor_held[1],
            "survivor_count": len(scan.candidates),
            "survivor_population_sha256": retained_digest.hexdigest(),
            "stdout_sha256": scan.stdout_sha256,
        },
        "exact_discriminant_features": {
            "content_free_homogeneous_degree": 20,
            "removed_polynomial_content": str(content),
            "trial_division_prime_bound": 997,
            **feature_audit,
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
            "scanner_sha256": EXPECTED_SCANNER_SHA256,
            "external_calls_use_foreground_process_groups": True,
            "temporary_scanner_binary_removed": True,
            "same_stage_retries": 0,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "compiler": shutil.which(args.compiler),
        },
        "timings": {"selection_wall_seconds": time.monotonic() - started},
    }
    if args.selection_only:
        artifact["status"] = "completed selection-only annulus"
        artifact["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
        artifact["result_sha256"] = stable_json_digest(
            {
                "scope": artifact["scope"],
                "scan": artifact["fresh_local_scan"],
                "features": artifact["exact_discriminant_features"],
                "selection": selection,
            }
        )
        exclusive_write(args.output, artifact)
        return

    conductor_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                conductor_worker,
                record["numerator"], record["denominator"],
                args.conductor_timeout, args.stack_bytes,
            )
            for record in selected
        ]
        for position, (record, future) in enumerate(zip(selected, futures), start=1):
            record["conductor_phase"] = future.result()
            if position % 32 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    eligible = [
        record for record in selected
        if record["conductor_phase"]["status"].startswith("completed")
    ]
    stages = (
        ("H5000", 5_000, None, args.h5000_timeout),
        ("H50000", 50_000, 24, args.h50000_timeout),
        ("H250000", 250_000, 6, args.h250000_timeout),
        ("H1000000", 1_000_000, 2, args.h1000000_timeout),
    )
    current = eligible
    for index, (name, height, keep, timeout) in enumerate(stages):
        if index:
            prior = stages[index - 1][0]
            completed = [
                record for record in current
                if record.get("point_stages", {}).get(prior, {}).get("status")
                == "completed"
            ]
            completed.sort(key=lambda record: stage_key(record, prior))
            current = completed[: int(keep)]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    point_worker,
                    record["numerator"], record["denominator"], height, timeout,
                    args.height_timeout, args.ellrank_timeout, args.stack_bytes,
                    args.mapping_cap, args.certificate_prime_bound,
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
        if max_rank >= TRIGGER_RANK:
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
    certificates = []
    target_hits = []
    for record in selected:
        for stage_name, stage in record.get("point_stages", {}).items():
            certificate = stage.get("finite_reduction_attempt", {})
            exact_rank = certificate.get("certified_algebraic_rank_lower_bound")
            if exact_rank is None:
                continue
            entry = {
                "parameter": record["parameter"],
                "stage": stage_name,
                "certified_algebraic_rank_lower_bound": exact_rank,
                "log_conductor": record["conductor_phase"]["log_conductor"],
            }
            certificates.append(entry)
            if (
                (exact_rank >= 21 and Decimal(entry["log_conductor"]) < TARGET_LOG_CONDUCTOR)
                or exact_rank >= 30
            ):
                target_hits.append(entry)

    artifact["status"] = (
        "completed fixed annulus; stopped without broadening"
        if maximum_rank is None or maximum_rank <= 15
        else "completed fixed annulus with rank signal"
    )
    artifact["conductor_first_screen"] = {
        "population_closed_before_any_conductor_or_point_call": True,
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
                "below_strict_log_conductor_target_numerically", False
            )
            for record in selected
        ),
    }
    artifact["point_search_protocol"] = {
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
        "finite_reduction_trigger_stable_rank": TRIGGER_RANK,
        "completed_stage_calls": len(completed_stages),
        "maximum_stable_numerical_rank": maximum_rank,
        "finite_reduction_certificates": certificates,
        "broadening_calls_after_fixed_protocol": 0,
    }
    artifact["target"] = {
        "rank_at_least": 21,
        "strict_log_conductor_upper_bound": "182.72",
        "alternative_rank_at_least": 30,
        "hits": target_hits,
    }
    artifact["timings"].update(
        {
            "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
            "total_wall_seconds": time.monotonic() - started,
        }
    )
    artifact["provenance"]["owned_processes_remaining"] = 0
    artifact["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    artifact["result_sha256"] = stable_json_digest(
        {
            "scope": artifact["scope"],
            "scan": artifact["fresh_local_scan"],
            "features": artifact["exact_discriminant_features"],
            "selection": selection,
            "records": selected,
            "conductors": artifact["conductor_first_screen"],
            "points": artifact["point_search_protocol"],
            "target": artifact["target"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        f"complete max_rank={maximum_rank} exact_certificates={len(certificates)} "
        f"target_hits={len(target_hits)} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

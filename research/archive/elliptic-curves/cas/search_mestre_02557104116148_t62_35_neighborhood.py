#!/usr/bin/env python3
"""Disjoint local/trace neighborhood search around exact-rank-16 T=62/35.

The preceding direct scan exhausts denominators at most 1,000.  This search
exhausts the disjoint Farey annulus ``1001 <= b <= 5000`` in the exact real
window ``11/7 <= T=a/b <= 69/35`` around 62/35.  Every prior Nagao
``parameter_t`` in the annulus is reconstructed from pinned artifacts and
excluded exactly.

Fresh primes 811..919 close independent local-rank-score and Frobenius-trace
neighborhood tails.  Their fixed union alone is read at held primes 929..1021.
Exact degree-20 discriminant features close a fixed conductor population
before any point search.  Staged exact quartic searches use H=5k,50k,250k,1m;
stable numerical rank at least 15 immediately receives an exact mod-3
finite-reduction attempt.  Numerical ranks remain triage evidence.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, sqrt
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
import search_mestre_rank14_pair_rational_frontier as engine
import search_mestre_02557104116148_direct_rational as direct
from search_mestre_root_tuple_scale import run_capped_process, sha256_file
from search_mestre_root_tuple_scale_max100 import stable_json_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ANCHOR = Q(62, 35)
DENOMINATOR_MIN = 1_001
DENOMINATOR_MAX = 5_000
WINDOW_MIN = Q(11, 7)
WINDOW_MAX = Q(69, 35)
DISCOVERY_PRIMES = (
    811, 821, 823, 827, 829, 839, 853, 857, 859,
    863, 877, 881, 883, 887, 907, 911, 919,
)
HELD_PRIMES = (
    929, 937, 941, 947, 953, 967, 971, 977,
    983, 991, 997, 1009, 1013, 1019, 1021,
)
SCORE_KEEP = 4_096
TRACE_KEEP = 4_096
SELECTION_QUOTAS = {
    "highest-held-local-score": 64,
    "closest-held-anchor-traces": 32,
    "smallest-exact-radical-upper-bound": 40,
    "largest-exact-known-powerful-part": 24,
}
POINT_STAGES = (
    ("H5000", 5_000, None),
    ("H50000", 50_000, 32),
    ("H250000", 250_000, 8),
    ("H1000000", 1_000_000, 2),
)
FINITE_REDUCTION_TRIGGER = 15
TARGET_LOG_CONDUCTOR = Decimal("182.72")
EXPECTED_PRIOR_ALL_COUNT = 13_119
EXPECTED_PRIOR_WINDOW = (
    Q(2065, 1296), Q(4351, 2478), Q(4223, 2530),
    Q(6065, 3674), Q(7093, 3790), Q(8351, 4856),
)
EXPECTED_GLOBAL_PRIMITIVE_POPULATION = 2_918_494
EXPECTED_GLOBAL_EVALUATED_POPULATION = 2_918_488
FROZEN_RANK16_CERTIFICATE = (
    "elliptic_mestre_02557104116148_t62_35_rank16_certificate.json"
)
EXPECTED_RANK16_CERTIFICATE_SHA256 = (
    "2c6d918546548227ac8f83287b3242e8d4261a98facd2665d506a8308f4c9fc7"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_mestre_02557104116148_t62_35_neighborhood.json"
)


class Candidate:
    __slots__ = (
        "numerator", "denominator", "discovery_score", "held_score",
        "discovery_good", "held_good", "discovery_similarity",
        "held_similarity", "discovery_matches", "held_matches",
    )

    def __init__(self, fields: Sequence[str]):
        if len(fields) != 11 or fields[0] not in {"C", "K"}:
            raise AssertionError("malformed neighborhood candidate")
        self.numerator = int(fields[1])
        self.denominator = int(fields[2])
        self.discovery_score = fields[3]
        self.held_score = fields[4]
        self.discovery_good = int(fields[5])
        self.held_good = int(fields[6])
        self.discovery_similarity = fields[7]
        self.held_similarity = fields[8]
        self.discovery_matches = int(fields[9])
        self.held_matches = int(fields[10])

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)

    def digest_row(self) -> list[Any]:
        return [
            self.numerator, self.denominator,
            self.discovery_score, self.held_score,
            self.discovery_good, self.held_good,
            self.discovery_similarity, self.held_similarity,
            self.discovery_matches, self.held_matches,
        ]


def configure_engine(prior: Sequence[Fraction]) -> None:
    engine.FAMILIES = (direct.FAMILY,)
    engine.DISCOVERY_PRIMES = DISCOVERY_PRIMES
    engine.HELD_PRIMES = HELD_PRIMES
    engine.PRIOR_PANEL_PARAMETERS = tuple(prior)
    engine.FINITE_REDUCTION_TRIGGER = FINITE_REDUCTION_TRIGGER
    engine.TARGET_LOG_CONDUCTOR = TARGET_LOG_CONDUCTOR


def prior_window_exclusions(generated: Path) -> tuple[Fraction, ...]:
    observed_files = {
        path.name
        for pattern in ("elliptic_nagao_rank13*.json", "elliptic_nagao_u*.json")
        for path in generated.glob(pattern)
    }
    if observed_files != set(direct.PRIOR_ARTIFACT_SHA256):
        raise AssertionError("the frozen prior artifact census changed")
    all_parameters: set[Fraction] = set()
    for name, expected in direct.PRIOR_ARTIFACT_SHA256.items():
        path = generated / name
        if sha256_file(path) != expected:
            raise AssertionError(f"frozen prior artifact changed: {name}")
        direct._collect_parameter_t(json.loads(path.read_text()), all_parameters)
    window = tuple(
        sorted(
            value for value in all_parameters
            if DENOMINATOR_MIN <= value.denominator <= DENOMINATOR_MAX
            and WINDOW_MIN <= value <= WINDOW_MAX
        )
    )
    if len(all_parameters) != EXPECTED_PRIOR_ALL_COUNT or window != tuple(
        sorted(EXPECTED_PRIOR_WINDOW)
    ):
        raise AssertionError("the prior annulus exclusion set changed")
    return window


def exclusion_text(values: Sequence[Fraction]) -> str:
    return "".join(
        f"{value.numerator}\t{value.denominator}\n"
        for value in sorted(values, key=lambda item: (
            item.denominator, item.numerator
        ))
    )


def parse_scanner_output(
    stdout: str, prior: set[Fraction]
) -> tuple[Candidate, tuple[Candidate, ...], dict[str, int | str]]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_02557104116148_T62_35_NEIGHBORHOOD_V1":
        raise AssertionError("the neighborhood scanner header changed")
    if tuple(map(int, lines[1].split()[1:])) != DISCOVERY_PRIMES:
        raise AssertionError("the neighborhood discovery band changed")
    if tuple(map(int, lines[2].split()[1:])) != HELD_PRIMES:
        raise AssertionError("the neighborhood held band changed")
    local = lines[3].split()
    anchor = Candidate(lines[4].split())
    candidates = tuple(Candidate(line.split()) for line in lines[5:-1])
    summary = lines[-1].split()
    if local[0] != "L" or len(local) != 3 or summary[0] != "S" or len(summary) != 10:
        raise AssertionError("the scanner digest/population footer changed")
    values = list(map(int, summary[1:]))
    keys = (
        "denominator_min", "denominator_max", "score_keep", "trace_keep",
        "primitive_population", "prior_excluded", "evaluated_population",
        "retained_union", "exclusion_file_count",
    )
    audit: dict[str, int | str] = dict(zip(keys, values))
    audit["discovery_table_digest"] = local[1]
    audit["held_table_digest"] = local[2]
    if (
        anchor.parameter != ANCHOR
        or Decimal(anchor.discovery_similarity) != 0
        or Decimal(anchor.held_similarity) != 0
        or anchor.discovery_matches != len(DISCOVERY_PRIMES)
        or anchor.held_matches != len(HELD_PRIMES)
    ):
        raise AssertionError("the exact anchor trace fingerprint changed")
    if (
        audit["denominator_min"] != DENOMINATOR_MIN
        or audit["denominator_max"] != DENOMINATOR_MAX
        or audit["score_keep"] != SCORE_KEEP
        or audit["trace_keep"] != TRACE_KEEP
        or audit["primitive_population"] != EXPECTED_GLOBAL_PRIMITIVE_POPULATION
        or audit["prior_excluded"] != len(prior)
        or audit["evaluated_population"] != EXPECTED_GLOBAL_EVALUATED_POPULATION
        or audit["retained_union"] != len(candidates)
        or audit["exclusion_file_count"] != len(prior)
    ):
        raise AssertionError("the exact neighborhood population changed")
    if len({candidate.parameter for candidate in candidates}) != len(candidates):
        raise AssertionError("the scanner emitted a duplicate")
    if any(
        candidate.parameter in prior
        or not DENOMINATOR_MIN <= candidate.denominator <= DENOMINATOR_MAX
        or not WINDOW_MIN <= candidate.parameter <= WINDOW_MAX
        or gcd(candidate.numerator, candidate.denominator) != 1
        for candidate in candidates
    ):
        raise AssertionError("a neighborhood survivor escaped scope")
    return anchor, candidates, audit


def run_scanner(
    source: Path,
    prior: Sequence[Fraction],
    *, compiler: str, compile_timeout: float, scan_timeout: float,
) -> tuple[Candidate, tuple[Candidate, ...], dict[str, int | str]]:
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="mestre-t62-35-neighbor-") as directory:
        directory_path = Path(directory)
        binary = directory_path / "scanner"
        exclusion_path = directory_path / "prior.tsv"
        exclusion_path.write_text(exclusion_text(prior))
        run_capped_process(
            (
                executable, "-std=c++17", "-O3", "-DNDEBUG",
                str(source), "-o", str(binary),
            ), timeout=compile_timeout,
        )
        stdout, _ = run_capped_process(
            (
                str(binary), str(DENOMINATOR_MIN), str(DENOMINATOR_MAX),
                str(SCORE_KEEP), str(TRACE_KEEP), str(exclusion_path),
            ), timeout=scan_timeout,
        )
    anchor, candidates, audit = parse_scanner_output(stdout, set(prior))
    audit["stdout_sha256"] = hashlib.sha256(stdout.encode()).hexdigest()
    return anchor, candidates, audit


def local_table_digest(primes: Sequence[int]) -> str:
    digest = 1_469_598_103_934_665_603

    def mix(value: int) -> None:
        nonlocal digest
        value &= (1 << 64) - 1
        for offset in range(8):
            digest ^= (value >> (8 * offset)) & 255
            digest = (digest * 1_099_511_628_211) & ((1 << 64) - 1)

    for prime in primes:
        mix(prime)
        traces = []
        for numerator, denominator in (
            *((residue, 1) for residue in range(prime)), (1, 0)
        ):
            trace = engine.exact_local_trace_projective(
                0, numerator, denominator, prime
            )
            traces.append(trace)
            mix(trace is not None)
            mix(0 if trace is None else trace)
        anchor_trace = engine.exact_local_trace_projective(
            0, ANCHOR.numerator, ANCHOR.denominator, prime
        )
        if anchor_trace is None:
            raise AssertionError("the rank-16 anchor became bad in a score band")
        mix(anchor_trace)
    return str(digest)


def pool_with_features(
    candidates: Sequence[Candidate], discriminant: Sequence[int]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    pool = []
    digest = hashlib.sha256()
    singular = 0
    for candidate in candidates:
        feature = engine.discriminant_feature(
            discriminant, candidate.numerator, candidate.denominator
        )
        record = {
            "numerator": candidate.numerator,
            "denominator": candidate.denominator,
            "parameter": str(candidate.parameter),
            "discovery_score": candidate.discovery_score,
            "held_score": candidate.held_score,
            "discovery_good": candidate.discovery_good,
            "held_good": candidate.held_good,
            "discovery_trace_similarity": candidate.discovery_similarity,
            "held_trace_similarity": candidate.held_similarity,
            "discovery_exact_trace_matches": candidate.discovery_matches,
            "held_exact_trace_matches": candidate.held_matches,
            "discriminant_feature": feature,
        }
        digest.update(
            (
                f"{candidate.parameter}|{candidate.discovery_score}|"
                f"{candidate.held_score}|{candidate.discovery_similarity}|"
                f"{candidate.held_similarity}|"
                f"{feature['absolute_homogeneous_discriminant']}|"
                f"{feature.get('combined_radical_upper_bound')}\n"
            ).encode()
        )
        if feature["singular"]:
            singular += 1
        else:
            pool.append(record)
    return pool, {
        "retained_union_before_singular_rejection": len(candidates),
        "exact_singular_rejections": singular,
        "admissible_feature_pool_count": len(pool),
        "exact_feature_population_sha256": digest.hexdigest(),
    }


def select_population(
    pool: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    orders = {
        "highest-held-local-score": sorted(
            pool, key=lambda row: (
                -Decimal(row["held_score"]),
                -Decimal(row["discovery_score"]),
                row["denominator"], row["numerator"],
            )
        ),
        "closest-held-anchor-traces": sorted(
            pool, key=lambda row: (
                -Decimal(row["held_trace_similarity"]),
                -row["held_exact_trace_matches"],
                -Decimal(row["held_score"]),
                row["denominator"], row["numerator"],
            )
        ),
        "smallest-exact-radical-upper-bound": sorted(
            pool, key=lambda row: (
                int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                -Decimal(row["held_score"]),
                row["denominator"], row["numerator"],
            )
        ),
        "largest-exact-known-powerful-part": sorted(
            pool, key=lambda row: (
                -int(row["discriminant_feature"]["known_powerful_part"]),
                int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                row["denominator"], row["numerator"],
            )
        ),
    }
    selected: dict[str, dict[str, Any]] = {}
    reasons: dict[str, set[str]] = {}
    for label, quota in SELECTION_QUOTAS.items():
        added = 0
        for row in orders[label]:
            parameter = row["parameter"]
            if parameter in selected:
                reasons[parameter].add(label)
                continue
            selected[parameter] = row
            reasons[parameter] = {label}
            added += 1
            if added == quota:
                break
        if added != quota:
            raise AssertionError(f"selection quota did not fill: {label}")
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
        "novel_candidates_added_per_stratum": SELECTION_QUOTAS,
        "selected_population": len(answer),
        "selected_population_sha256": digest.hexdigest(),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-only", action="store_true")
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
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.compile_timeout, args.scan_timeout, args.conductor_timeout,
        args.h5000_timeout, args.h50000_timeout, args.h250000_timeout,
        args.h1000000_timeout, args.height_timeout, args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all process caps must lie in (0,60]")
    if args.mapping_cap < 32 or args.mapping_cap > 1024:
        raise SystemExit("mapping cap out of range")
    if args.certificate_prime_bound != 1_000:
        raise SystemExit("the finite-reduction prime bound is pinned at 1000")


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
        raise SystemExit("refusing to overwrite the neighborhood artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    generated = root / "artifacts/generated-results"
    scanner_path = script_path.with_name(
        "scan_mestre_02557104116148_t62_35_neighborhood.cpp"
    )
    certificate_path = generated / FROZEN_RANK16_CERTIFICATE
    if sha256_file(certificate_path) != EXPECTED_RANK16_CERTIFICATE_SHA256:
        raise AssertionError("the frozen exact rank-16 anchor changed")
    certificate = json.loads(certificate_path.read_text())
    if (
        certificate["theorem"]["certified_algebraic_rank_lower_bound"] != 16
        or Q(certificate["curve"]["parameter_T"]) != ANCHOR
    ):
        raise AssertionError("the anchor exact-rank checkpoint changed")
    prior = prior_window_exclusions(generated)
    configure_engine(prior)
    construction = SixRootMestreConstruction(tuple(Q(root) for root in direct.ROOTS))
    raw_discriminant = construction.primitive_discriminant_polynomial
    content = engine.polynomial_content(raw_discriminant)
    discriminant = tuple(value.numerator // content for value in raw_discriminant)

    anchor, candidates, scan_audit = run_scanner(
        scanner_path,
        prior,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    if (
        scan_audit["discovery_table_digest"] != local_table_digest(DISCOVERY_PRIMES)
        or scan_audit["held_table_digest"] != local_table_digest(HELD_PRIMES)
    ):
        raise AssertionError("the C++ exact local tables missed Python replay")
    expected_discovery = engine.score_text(0, ANCHOR, DISCOVERY_PRIMES)
    expected_held = engine.score_text(0, ANCHOR, HELD_PRIMES)
    if (
        abs(Decimal(anchor.discovery_score) - Decimal(expected_discovery[0]))
        > Decimal("0.000000000100")
        or anchor.discovery_good != expected_discovery[1]
        or abs(Decimal(anchor.held_score) - Decimal(expected_held[0]))
        > Decimal("0.000000000100")
        or anchor.held_good != expected_held[1]
    ):
        raise AssertionError("the anchor exact local score missed replay")
    print(
        f"disjoint neighborhood closed: primitive={scan_audit['primitive_population']} "
        f"excluded={scan_audit['prior_excluded']} union={len(candidates)}",
        flush=True,
    )
    pool, pool_audit = pool_with_features(candidates, discriminant)
    selected, selection = select_population(pool)
    print(
        f"exact features closed; conductor population={len(selected)}",
        flush=True,
    )

    common: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "completed selection-only disjoint T=62/35 neighborhood"
            if args.selection_only
            else "in-progress conductor-first T=62/35 neighborhood"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "scope": {
            "family_roots": list(direct.ROOTS),
            "anchor_T": str(ANCHOR),
            "anchor_exact_rank_lower_bound": 16,
            "denominator_interval": [DENOMINATOR_MIN, DENOMINATOR_MAX],
            "exact_parameter_window": [str(WINDOW_MIN), str(WINDOW_MAX)],
            "disjoint_from_complete_denominator_at_most_1000_scan": True,
            "prior_parameter_t_exclusions": [str(value) for value in prior],
        },
        "modular_scan": {
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "fresh_relative_to_predecessor_score_bands_through_809": True,
            "score_tail_keep": SCORE_KEEP,
            "trace_similarity_tail_keep": TRACE_KEEP,
            "trace_similarity": "-sum((a_p(T)-a_p(anchor))^2/(4p))",
            "anchor": anchor.digest_row(),
            "audit": scan_audit,
            "retained_union_sha256": stable_json_digest(
                [candidate.digest_row() for candidate in candidates]
            ),
        },
        "exact_discriminant_feature_screen": {
            "content_free_homogeneous_degree": 20,
            "removed_content": str(content),
            "trial_division_prime_bound": 997,
            "pool_audit": pool_audit,
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
            "rank16_anchor_certificate": str(certificate_path.relative_to(root)),
            "rank16_anchor_certificate_sha256": EXPECTED_RANK16_CERTIFICATE_SHA256,
            "command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "temporary_scanner_and_exclusion_file_removed": True,
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
            ) for row in selected
        ]
        for position, (row, future) in enumerate(zip(selected, futures), start=1):
            row["conductor_phase"] = future.result()
            if position % 32 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    eligible = [
        row for row in selected
        if row["conductor_phase"]["status"].startswith("completed")
    ]
    timeouts = {
        "H5000": args.h5000_timeout,
        "H50000": args.h50000_timeout,
        "H250000": args.h250000_timeout,
        "H1000000": args.h1000000_timeout,
    }
    current = eligible
    for index, (name, height, keep) in enumerate(POINT_STAGES):
        if index:
            prior_name = POINT_STAGES[index - 1][0]
            current = [
                row for row in current
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
                    height, timeouts[name], args.height_timeout,
                    args.ellrank_timeout, args.stack_bytes, args.mapping_cap,
                    args.certificate_prime_bound,
                ) for row in current
            ]
            for row, future in zip(current, futures):
                row.setdefault("point_stages", {})[name] = future.result()
        maximum = max(
            (
                row["point_stages"][name].get("stable_numerical_rank", -1)
                for row in current
            ), default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={maximum}", flush=True)
        if maximum >= 15:
            print(f"EARLY_SIGNAL {name} stable_rank={maximum}", flush=True)

    completed = [
        stage for row in selected
        for stage in row.get("point_stages", {}).values()
        if stage.get("status") == "completed"
    ]
    maximum_rank = max(
        (int(stage["stable_numerical_rank"]) for stage in completed),
        default=None,
    )
    finite_attempts, target_hits = [], []
    for row in selected:
        for stage_name, stage in row.get("point_stages", {}).items():
            certificate = stage.get("finite_reduction_attempt", {})
            certified = certificate.get("certified_algebraic_rank_lower_bound")
            if certified is None:
                continue
            finite_attempts.append({
                "parameter": row["parameter"], "stage": stage_name,
                "certified_rank_lower_bound": certified,
            })
            below = row["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically", False
            )
            if certified >= 30 or (certified >= 21 and below):
                target_hits.append({
                    "parameter": row["parameter"], "stage": stage_name,
                    "certified_rank_lower_bound": certified,
                    "conductor": row["conductor_phase"]["conductor"],
                    "log_conductor": row["conductor_phase"]["log_conductor"],
                })
    common["status"] = "completed fixed disjoint T=62/35 neighborhood"
    common["target"]["hits"] = target_hits
    common["conductor_first_screen"] = {
        "population_closed_before_any_point_or_rank_call": True,
        "selected_population": len(selected),
        "completed": len(eligible),
        "timeouts": sum(row["conductor_phase"]["status"] == "timeout" for row in selected),
        "errors": sum(row["conductor_phase"]["status"] == "error" for row in selected),
        "subtarget": sum(
            row["conductor_phase"].get("below_strict_log_conductor_target_numerically")
            is True for row in selected
        ),
    }
    common["point_search_protocol"] = {
        "stages": [
            {"name": name, "height_bound": height, "keep": keep,
             "attempted": sum(name in row.get("point_stages", {}) for row in selected)}
            for name, height, keep in POINT_STAGES
        ],
        "finite_reduction_trigger_stable_rank": FINITE_REDUCTION_TRIGGER,
        "finite_reduction_attempts": finite_attempts,
        "maximum_stable_numerical_rank": maximum_rank,
        "completed_stage_calls": len(completed),
        "same_height_retries": 0,
    }
    common["timings"].update({
        "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
        "total_wall_seconds": time.monotonic() - started,
    })
    common["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    common["provenance"]["owned_processes_remaining"] = 0
    common["result_sha256"] = stable_json_digest({
        "scope": common["scope"], "scan": common["modular_scan"],
        "features": common["exact_discriminant_feature_screen"],
        "selection": common["conductor_selection"], "records": selected,
        "conductor": common["conductor_first_screen"],
        "points": common["point_search_protocol"], "target": common["target"],
    })
    exclusive_write(args.output, common)
    print(
        f"complete max_rank={maximum_rank} target_hits={len(target_hits)} "
        f"output={args.output}", flush=True,
    )


if __name__ == "__main__":
    main()

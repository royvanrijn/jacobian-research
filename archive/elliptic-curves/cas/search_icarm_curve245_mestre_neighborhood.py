#!/usr/bin/env python3
"""Conductor-first rational neighborhood scan around ICARM curve #245.

The exact parent is the canonical six-root Mestre family
``(0,106,344,475,594,731)`` at ``T=5801/10``.  This bounded experiment scans
reduced positive ``T=a/b`` in a declared interval and denominator range.  A
compiled exact-local sieve closes its survivor pool before held-out scores,
exact discriminant features, conductors, or point searches are read.

Only a fixed union of high held-score and low-conductor fibers reaches the
quartic point stages.  Numerical height rank is triage; rank claims are made
only by the finite-reduction certificate.  A bounded miss is not an upper
bound on this family or on elliptic-curve rank.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile
import time
from typing import Any, Sequence

import search_mestre_02557104116148_direct_rational as direct
import search_mestre_0430313946_frontier as frontier
import search_mestre_rank14_pair_rational_frontier as engine
from icarm_curve245_mestre import (
    A_COEFFICIENTS,
    B_COEFFICIENTS,
    CANONICAL_PARAMETER,
    CANONICAL_ROOTS,
    CONSTRUCTION,
    extra_quartic_point,
    primitive_short_model,
)
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    point_digest,
    quartic_point_to_jacobian,
    run_capped_process,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

FAMILY_LABEL = "icarm245_mestre_u3_2_v2_r0_106_344_475_594_731"
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DISCOVERY_PRIMES = (211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281)
HELD_PRIMES = (283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373)
STACK_BYTES = 512_000_000
FINITE_REDUCTION_TRIGGER = 18
DEFAULT_OUTPUT = Path(
    "artifacts/local/elliptic-curves/icarm245-mestre-neighborhood-v1.json"
)

FAMILY = engine.FamilySpec(
    0,
    FAMILY_LABEL,
    CANONICAL_ROOTS,
    CANONICAL_PARAMETER,
    A_COEFFICIENTS,
    B_COEFFICIENTS,
)


def configure_engines() -> None:
    engine.FAMILIES = (FAMILY,)
    engine.DISCOVERY_PRIMES = DISCOVERY_PRIMES
    engine.HELD_PRIMES = HELD_PRIMES
    engine.PRIOR_PANEL_PARAMETERS = (CANONICAL_PARAMETER,)
    engine.FINITE_REDUCTION_TRIGGER = FINITE_REDUCTION_TRIGGER
    engine.TARGET_LOG_CONDUCTOR = TARGET_LOG_CONDUCTOR

    direct.FAMILY_LABEL = FAMILY_LABEL
    direct.ROOTS = CANONICAL_ROOTS
    direct.CALIBRATION_PARAMETER = CANONICAL_PARAMETER
    direct.A_COEFFICIENTS = A_COEFFICIENTS
    direct.B_COEFFICIENTS = B_COEFFICIENTS
    direct.FAMILY = FAMILY
    direct.DISCOVERY_PRIMES = DISCOVERY_PRIMES
    direct.HELD_PRIMES = HELD_PRIMES
    direct.SELECTION_QUOTAS = {
        "highest-held-score": 48,
        "smallest-exact-radical-upper-bound": 48,
        "largest-exact-known-powerful-part": 24,
        "balanced-held-and-power-rank": 24,
    }
    direct.configure_engine((CANONICAL_PARAMETER,))


def _cpp_integer_array(values: Sequence[int]) -> str:
    return ", ".join(f'"{int(value)}"' for value in values)


def render_scanner_source(
    template: str,
    *,
    parameter_min: Fraction,
    parameter_max: Fraction,
) -> str:
    """Adapt the audited direct-T scanner to this family and interval."""

    source = template
    source = re.sub(
        r"static constexpr std::array<int, 17> DISCOVERY_PRIMES\{.*?\};",
        "static constexpr std::array<int, 14> DISCOVERY_PRIMES{" +
        ", ".join(map(str, DISCOVERY_PRIMES)) + "};",
        source,
        count=1,
        flags=re.S,
    )
    source = re.sub(
        r"static constexpr std::array<int, 14> HELD_PRIMES\{.*?\};",
        "static constexpr std::array<int, 14> HELD_PRIMES{" +
        ", ".join(map(str, HELD_PRIMES)) + "};",
        source,
        count=1,
        flags=re.S,
    )
    source = re.sub(
        r"static const std::array<std::string, 9> A_COEFFICIENTS\{\{.*?\}\};",
        "static const std::array<std::string, 9> A_COEFFICIENTS{{" +
        _cpp_integer_array(A_COEFFICIENTS) + "}};",
        source,
        count=1,
        flags=re.S,
    )
    source = re.sub(
        r"static const std::array<std::string, 13> B_COEFFICIENTS\{\{.*?\}\};",
        "static const std::array<std::string, 13> B_COEFFICIENTS{{" +
        _cpp_integer_array(B_COEFFICIENTS) + "}};",
        source,
        count=1,
        flags=re.S,
    )
    source = source.replace(
        "calibration.numerator = 1;\n  calibration.denominator = 1;",
        f"calibration.numerator = {CANONICAL_PARAMETER.numerator};\n"
        f"  calibration.denominator = {CANONICAL_PARAMETER.denominator};",
        1,
    )
    source = source.replace(
        '<< "\\nL " << table_digest(discovery)',
        '<< "\\nL " << table_digest(discovery)',
        1,
    )
    source = source.replace(
        '<< "\\nK 1 1 " << score_text(calibration.discovery_score)',
        f'<< "\\nK {CANONICAL_PARAMETER.numerator} {CANONICAL_PARAMETER.denominator} " '
        '<< score_text(calibration.discovery_score)',
        1,
    )
    lower = (
        f"static_cast<std::int64_t>(numerator) * {parameter_min.denominator} "
        f"< static_cast<std::int64_t>(denominator) * {parameter_min.numerator}"
    )
    upper = (
        f"static_cast<std::int64_t>(numerator) * {parameter_max.denominator} "
        f"> static_cast<std::int64_t>(denominator) * {parameter_max.numerator}"
    )
    source = source.replace(
        "      if (std::gcd(numerator, denominator) != 1) continue;",
        f"      if ({lower} || {upper}) continue;\n"
        "      if (std::gcd(numerator, denominator) != 1) continue;",
        1,
    )
    source = source.replace(
        "numerator_bound > 100000", "numerator_bound > 1000000", 1
    ).replace("require NUM<=100000", "require NUM<=1000000", 1)
    required = (
        str(A_COEFFICIENTS[0]),
        str(B_COEFFICIENTS[0]),
        f"K {CANONICAL_PARAMETER.numerator} {CANONICAL_PARAMETER.denominator}",
        lower,
        upper,
    )
    if any(value not in source for value in required):
        raise AssertionError("the temporary scanner adaptation was incomplete")
    return source


def augmented_point_stage(
    numerator: int,
    denominator: int,
    *,
    height_bound: int,
    point_timeout: float,
    height_timeout: float,
    ellrank_timeout: float,
    mapping_cap: int,
    certificate_prime_bound: int,
    forced_quartic_points: Sequence[tuple[Fraction, Fraction]] = (),
) -> dict[str, Any]:
    """Run the standard stage and inject Fermigier's exact extra section."""

    parameter = Q(numerator, denominator)
    coefficients = primitive_short_model(parameter)
    try:
        stage, subset = frontier.exact_point_stage(
            CONSTRUCTION,
            parameter,
            coefficients,
            height_bound=height_bound,
            point_timeout=point_timeout,
            height_timeout=height_timeout,
            ellrank_timeout=ellrank_timeout,
            stack_bytes=STACK_BYTES,
            mapping_cap=mapping_cap,
        )
        quartic_extra = extra_quartic_point(parameter)
        jacobian_extra = quartic_point_to_jacobian(
            CONSTRUCTION, parameter, quartic_extra
        )
        forced_jacobian = tuple(
            quartic_point_to_jacobian(CONSTRUCTION, parameter, point)
            for point in forced_quartic_points
        )
        pool_by_x = {point[0]: point for point in subset}
        pool_by_x.setdefault(jacobian_extra[0], jacobian_extra)
        for point in forced_jacobian:
            pool_by_x.setdefault(point[0], point)
        pool = tuple(pool_by_x.values())
        height = frontier.height_matrix_replay(
            coefficients,
            pool,
            precisions=(72, 120),
            timeout=height_timeout,
            stack_bytes=STACK_BYTES,
        )
        augmented_subset = frontier.numerical_subset(pool, height)
        stage["fermigier_extra_quartic_point"] = frontier.point_record(quartic_extra)
        stage["fermigier_extra_jacobian_point"] = frontier.point_record(jacobian_extra)
        stage["fermigier_extra_injected_exactly"] = True
        stage["forced_quartic_points"] = [
            frontier.point_record(point) for point in forced_quartic_points
        ]
        stage["forced_jacobian_points"] = [
            frontier.point_record(point) for point in forced_jacobian
        ]
        stage["augmented_pool_sha256"] = point_digest(pool)
        stage["height_matrix_runs"] = list(height)
        stage["stable_numerical_rank"] = int(height[-1]["numerical_rank"])
        stage["numerical_subset"] = [
            frontier.point_record(point) for point in augmented_subset
        ]
        if stage["stable_numerical_rank"] >= FINITE_REDUCTION_TRIGGER:
            stage["finite_reduction_attempt"] = mod3_independence_certificate(
                coefficients,
                augmented_subset,
                prime_bound=certificate_prime_bound,
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
            "height_bound": height_bound,
            "timeout_seconds": point_timeout,
            "same_height_retry": False,
        }
    except Exception as error:
        return {
            "status": "error",
            "height_bound": height_bound,
            "error": str(error)[:1000],
            "same_height_retry": False,
        }


def point_population(records: Sequence[dict[str, Any]], keep: int) -> list[dict[str, Any]]:
    eligible = [
        row for row in records
        if row["conductor_phase"].get("below_strict_log_conductor_target_numerically")
    ]
    orders = (
        sorted(eligible, key=lambda row: (-Decimal(row["held_score"]), row["denominator"], row["numerator"])),
        sorted(eligible, key=lambda row: (Decimal(row["conductor_phase"]["log_conductor"]), -Decimal(row["held_score"]))),
        sorted(eligible, key=lambda row: (-int(row["discriminant_feature"]["known_powerful_part"]), -Decimal(row["held_score"]))),
    )
    selected: dict[str, dict[str, Any]] = {}
    quota = max(1, keep // len(orders))
    for order in orders:
        added = 0
        for row in order:
            if row["parameter"] in selected:
                continue
            selected[row["parameter"]] = row
            added += 1
            if added == quota:
                break
    for row in orders[0]:
        selected.setdefault(row["parameter"], row)
        if len(selected) == min(keep, len(eligible)):
            break
    return sorted(selected.values(), key=lambda row: (row["denominator"], row["numerator"]))


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-min", type=Q, default=Q(500))
    parser.add_argument("--parameter-max", type=Q, default=Q(660))
    parser.add_argument("--numerator-bound", type=int)
    parser.add_argument("--denominator-bound", type=int, default=320)
    parser.add_argument("--keep", type=int, default=4096)
    parser.add_argument("--conductor-keep", type=int, default=144)
    parser.add_argument("--h5000-keep", type=int, default=36)
    parser.add_argument("--h50000-keep", type=int, default=12)
    parser.add_argument("--h250000-keep", type=int, default=3)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--compile-timeout", type=float, default=30)
    parser.add_argument("--scan-timeout", type=float, default=45)
    parser.add_argument("--conductor-timeout", type=float, default=12)
    parser.add_argument("--h5000-timeout", type=float, default=15)
    parser.add_argument("--h50000-timeout", type=float, default=25)
    parser.add_argument("--h250000-timeout", type=float, default=40)
    parser.add_argument("--height-timeout", type=float, default=30)
    parser.add_argument("--ellrank-timeout", type=float, default=8)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=499)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not (Q(0) < args.parameter_min < CANONICAL_PARAMETER < args.parameter_max):
        raise SystemExit("the interval must be positive and contain the anchor")
    if not (1 <= args.denominator_bound <= 2000 and 128 <= args.keep <= 20000):
        raise SystemExit("invalid denominator/keep bound")
    if not (1 <= args.workers <= 8):
        raise SystemExit("workers must lie in [1,8]")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the neighborhood artifact")
    configure_engines()
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    template_path = script_path.with_name(
        "scan_mestre_02557104116148_direct_rational.cpp"
    )
    numerator_bound = (
        args.numerator_bound
        if args.numerator_bound is not None
        else int(args.parameter_max * args.denominator_bound) + 1
    )
    if not 1 <= numerator_bound <= 1_000_000:
        raise SystemExit("the numerator bound must lie in [1,1000000]")
    exclusions = (CANONICAL_PARAMETER,)
    template = template_path.read_text()
    source = render_scanner_source(
        template,
        parameter_min=args.parameter_min,
        parameter_max=args.parameter_max,
    )
    with tempfile.TemporaryDirectory(prefix="icarm245-neighborhood-source-") as directory:
        source_path = Path(directory) / "scanner.cpp"
        source_path.write_text(source)
        scans = direct.run_scanners(
            source_path,
            exclusions,
            compiler="c++",
            compile_timeout=args.compile_timeout,
            scan_timeout=args.scan_timeout,
            denominator_bound=args.denominator_bound,
            strata=(("local-window", numerator_bound, args.keep),),
        )
    scan = scans[0]
    print(
        f"local scan closed: primitive={scan.primitive_population} "
        f"evaluated={scan.evaluated_population} retained={len(scan.candidates)}",
        flush=True,
    )

    raw_discriminant = CONSTRUCTION.primitive_discriminant_polynomial
    content = engine.polynomial_content(raw_discriminant)
    discriminant = tuple(value.numerator // content for value in raw_discriminant)
    pool, pool_audit = engine.pool_and_features(FAMILY, scans, discriminant)
    selected, selection = direct.select_conductor_population(pool)
    selected = selected[: args.conductor_keep]
    print(f"conductor population closed: {len(selected)}", flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                engine.conductor_worker,
                0,
                row["numerator"],
                row["denominator"],
                args.conductor_timeout,
                STACK_BYTES,
            )
            for row in selected
        ]
        for position, (row, future) in enumerate(zip(selected, futures), start=1):
            row["conductor_phase"] = future.result()
            if position % 24 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    completed = [
        row for row in selected
        if row["conductor_phase"]["status"].startswith("completed")
    ]
    subtarget = [
        row for row in completed
        if row["conductor_phase"]["below_strict_log_conductor_target_numerically"]
    ]
    print(
        f"conductors complete={len(completed)} subtarget={len(subtarget)}",
        flush=True,
    )

    stages = (
        ("H5000", 5_000, args.h5000_keep, args.h5000_timeout),
        ("H50000", 50_000, args.h50000_keep, args.h50000_timeout),
        ("H250000", 250_000, args.h250000_keep, args.h250000_timeout),
    )
    current = point_population(completed, args.h5000_keep)
    for stage_index, (name, height, keep, timeout) in enumerate(stages):
        if stage_index:
            prior = stages[stage_index - 1][0]
            current = [
                row for row in current
                if row.get("point_stages", {}).get(prior, {}).get("status") == "completed"
            ]
            current.sort(
                key=lambda row: (
                    -int(row["point_stages"][prior]["stable_numerical_rank"]),
                    -Decimal(row["held_score"]),
                    Decimal(row["conductor_phase"]["log_conductor"]),
                )
            )
            current = current[:keep]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    augmented_point_stage,
                    row["numerator"],
                    row["denominator"],
                    height_bound=height,
                    point_timeout=timeout,
                    height_timeout=args.height_timeout,
                    ellrank_timeout=args.ellrank_timeout,
                    mapping_cap=args.mapping_cap,
                    certificate_prime_bound=args.certificate_prime_bound,
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

    target_hits = []
    finite_attempts = []
    maximum_rank = -1
    for row in selected:
        for name, stage in row.get("point_stages", {}).items():
            maximum_rank = max(maximum_rank, int(stage.get("stable_numerical_rank", -1)))
            certificate = stage.get("finite_reduction_attempt", {})
            rank = certificate.get("certified_algebraic_rank_lower_bound")
            if rank is None:
                continue
            finite_attempts.append({"parameter": row["parameter"], "stage": name, "rank": rank})
            if rank >= 30 or (
                rank >= 21
                and row["conductor_phase"].get("below_strict_log_conductor_target_numerically")
            ):
                target_hits.append(
                    {
                        "parameter": row["parameter"],
                        "stage": name,
                        "certified_rank_lower_bound": rank,
                        "conductor": row["conductor_phase"]["conductor"],
                        "log_conductor": row["conductor_phase"]["log_conductor"],
                    }
                )

    artifact = {
        "schema_version": 1,
        "status": "completed bounded conductor-first neighborhood experiment",
        "mathematical_status": (
            "exact local scores, discriminant features, conductors, quartic points, "
            "and finite-reduction certificates; numerical height ranks are triage"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "family": {
            "label": FAMILY_LABEL,
            "roots": list(CANONICAL_ROOTS),
            "anchor_parameter": str(CANONICAL_PARAMETER),
            "A_coefficients_ascending": list(A_COEFFICIENTS),
            "B_coefficients_ascending": list(B_COEFFICIENTS),
            "content_free_discriminant_coefficients_ascending": list(discriminant),
            "removed_discriminant_content": content,
        },
        "scope": {
            "parameter_interval": [str(args.parameter_min), str(args.parameter_max)],
            "denominator_bound": args.denominator_bound,
            "numerator_bound_implied": numerator_bound,
            "anchor_excluded": True,
            "primitive_population": scan.primitive_population,
            "evaluated_population": scan.evaluated_population,
            "discovery_survivors": len(scan.candidates),
        },
        "modular_scan": {
            "discovery_primes": list(DISCOVERY_PRIMES),
            "held_primes": list(HELD_PRIMES),
            "bands_disjoint": not set(DISCOVERY_PRIMES) & set(HELD_PRIMES),
            "scan": direct.scanner_record(scan),
        },
        "exact_discriminant_feature_screen": pool_audit,
        "conductor_selection": selection,
        "conductor_screen": {
            "selected": len(selected),
            "completed": len(completed),
            "subtarget": len(subtarget),
        },
        "point_search": {
            "stages": [
                {"name": name, "height": height, "keep": keep, "timeout": timeout}
                for name, height, keep, timeout in stages
            ],
            "maximum_stable_numerical_rank": maximum_rank,
            "finite_reduction_attempts": finite_attempts,
            "same_height_retries": 0,
        },
        "records": selected,
        "parameters": {
            key: str(value) if isinstance(value, Q) else value
            for key, value in vars(args).items()
            if key != "output"
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "scanner_template": str(template_path.relative_to(root)),
            "scanner_template_sha256": sha256_file(template_path),
            "rendered_scanner_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "reproducing_command": " ".join(sys.argv),
        },
        "software": {
            "python": platform.python_version(),
            "compiler": shutil.which("c++"),
        },
        "timings": {"total_wall_seconds": time.monotonic() - started},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "family": artifact["family"],
            "scope": artifact["scope"],
            "scan": artifact["modular_scan"],
            "conductors": artifact["conductor_screen"],
            "points": artifact["point_search"],
            "records": selected,
            "target": artifact["target"],
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"complete max_rank={maximum_rank} target_hits={len(target_hits)} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

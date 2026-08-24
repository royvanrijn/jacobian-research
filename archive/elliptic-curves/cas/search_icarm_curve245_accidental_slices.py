#!/usr/bin/env python3
"""Promote the eight accidental rank directions on ICARM curve #245.

An exact mod-3 certificate first replays rank 20 on the recovered Mestre
quartic and selects eight non-generic pivot points.  Through each pivot the
two lines ``x=+T+n`` and ``x=-T+n`` cut the Mestre surface in a genus-one
quartic.  A bounded rational-point search on these 16 slices therefore finds
nearby fibers carrying a prescribed non-generic point, rather than sampling
unconditioned fibers.

The candidate set closes before local scores, discriminant features,
conductors, or additional point searches are used.  Numerical height ranks
remain triage; only finite-reduction certificates are rank lower bounds.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Sequence

import search_mestre_rank14_pair_rational_frontier as engine
from icarm_curve245_mestre import (
    CANONICAL_PARAMETER,
    CANONICAL_ROOTS,
    CONSTRUCTION,
    extra_quartic_point,
    primitive_short_model,
)
from search_icarm_curve245_mestre_neighborhood import (
    DISCOVERY_PRIMES,
    FAMILY,
    FINITE_REDUCTION_TRIGGER,
    HELD_PRIMES,
    STACK_BYTES,
    TARGET_LOG_CONDUCTOR,
    augmented_point_stage,
    configure_engines,
)
from search_mestre_root_tuple_scale import (
    bounded_quartic_points,
    canonical_signless_points,
    point_digest,
    primitive_visible_points,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from search_six_root_low_conductor_centers import (
    Seed,
    build_slice_records,
    pool_with_sources,
)


Q = Fraction
DEFAULT_OUTPUT = Path(
    "artifacts/local/elliptic-curves/icarm245-accidental-slices-v1.json"
)


def anchor_accidental_pivots() -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    raw = bounded_quartic_points(
        CONSTRUCTION.primitive_quartic_coefficients(CANONICAL_PARAMETER),
        height_bound=50_000,
        timeout=30,
        stack_bytes=STACK_BYTES,
    )
    searched = canonical_signless_points(raw)
    generic_extra = extra_quartic_point(CANONICAL_PARAMETER)
    points, sources, visible_count, prescribed_count = pool_with_sources(
        CONSTRUCTION,
        CANONICAL_PARAMETER,
        searched,
        (generic_extra,),
    )
    certificate = mod3_independence_certificate(
        primitive_short_model(CANONICAL_PARAMETER),
        points,
        prime_bound=499,
    )
    if certificate.get("certified_algebraic_rank_lower_bound") != 20:
        raise AssertionError("the exact anchor rank-20 replay changed")
    generic_x = {
        point[0] for point in primitive_visible_points(CONSTRUCTION, CANONICAL_PARAMETER)
    }
    generic_x.add(generic_extra[0])
    pivots = []
    accidentals = []
    for one_based in certificate["independent_subset_indices_one_based"]:
        source = sources[one_based - 1]
        quartic_point = (
            Q(source["quartic_point"]["x"]),
            Q(source["quartic_point"]["y"]),
        )
        generic = quartic_point[0] in generic_x
        pivots.append(
            {
                "pool_index_one_based": one_based,
                **source,
                "generic_abscissa": generic,
            }
        )
        if not generic:
            accidentals.append(quartic_point)
    if len(accidentals) != 8:
        raise AssertionError("the anchor no longer has eight accidental pivots")
    return (
        {
            "parameter": str(CANONICAL_PARAMETER),
            "search_height": 50_000,
            "signed_quartic_points": len(raw),
            "signless_quartic_points": len(searched),
            "pool_points_modulo_inverse": len(points),
            "visible_columns": visible_count,
            "prescribed_visible_plus_extra_columns": prescribed_count,
            "pool_point_sha256": point_digest(points),
            "finite_reduction_certificate": certificate,
            "pivot_records": pivots,
            "accidental_pivot_points": [
                {"x": str(point[0]), "y": str(point[1])} for point in accidentals
            ],
        },
        tuple(accidentals),
    )


def select_conductors(
    candidates: Sequence[dict[str, Any]], keep: int
) -> list[dict[str, Any]]:
    orders = (
        sorted(candidates, key=lambda row: (-Decimal(row["held_score"]), row["parameter_height"])),
        sorted(candidates, key=lambda row: (int(row["discriminant_feature"]["combined_radical_upper_bound"]), -Decimal(row["held_score"]))),
        sorted(candidates, key=lambda row: (row["parameter_height"], Q(row["distance_from_anchor"]))),
        sorted(candidates, key=lambda row: (Q(row["distance_from_anchor"]), row["parameter_height"])),
        sorted(candidates, key=lambda row: (-row["slice_source_count"], -len(row["forced_quartic_points"]), -Decimal(row["held_score"]))),
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
        if len(selected) == min(keep, len(candidates)):
            break
    return sorted(selected.values(), key=lambda row: (row["denominator"], row["numerator"]))


def select_points(records: Sequence[dict[str, Any]], keep: int) -> list[dict[str, Any]]:
    eligible = [
        row for row in records
        if row["conductor_phase"].get("below_strict_log_conductor_target_numerically")
    ]
    orders = (
        sorted(eligible, key=lambda row: (-row["slice_source_count"], -len(row["forced_quartic_points"]), -Decimal(row["held_score"]))),
        sorted(eligible, key=lambda row: (-Decimal(row["held_score"]), Decimal(row["conductor_phase"]["log_conductor"]))),
        sorted(eligible, key=lambda row: (Decimal(row["conductor_phase"]["log_conductor"]), -Decimal(row["held_score"]))),
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
    parser.add_argument("--slice-height", type=int, default=50_000)
    parser.add_argument("--conductor-keep", type=int, default=120)
    parser.add_argument("--h5000-keep", type=int, default=36)
    parser.add_argument("--h50000-keep", type=int, default=12)
    parser.add_argument("--h250000-keep", type=int, default=3)
    parser.add_argument("--workers", type=int, default=4)
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
    if not (1 <= args.workers <= 8 and 1 <= args.conductor_keep <= 236):
        raise SystemExit("invalid workers/conductor selection")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the slice artifact")
    configure_engines()
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]

    anchor, _ = anchor_accidental_pivots()
    seed = Seed(
        "icarm245-mestre",
        CANONICAL_ROOTS,
        CANONICAL_PARAMETER,
        50_000,
        "",
        20,
        13,
    )
    slices, candidates = build_slice_records(
        seed,
        {"accidental_pivot_points": anchor["accidental_pivot_points"]},
        height_bound=args.slice_height,
    )
    print(f"slice population closed: slices={len(slices)} candidates={len(candidates)}", flush=True)

    raw_discriminant = CONSTRUCTION.primitive_discriminant_polynomial
    content = engine.polynomial_content(raw_discriminant)
    discriminant = tuple(value.numerator // content for value in raw_discriminant)
    enriched = []
    for candidate in candidates:
        parameter = Q(candidate["parameter"])
        feature = engine.discriminant_feature(
            discriminant, parameter.numerator, parameter.denominator
        )
        if feature["singular"]:
            continue
        record = dict(candidate)
        record.update(
            {
                "numerator": parameter.numerator,
                "denominator": parameter.denominator,
                "discovery_score": engine.score_text(0, parameter, DISCOVERY_PRIMES)[0],
                "held_score": engine.score_text(0, parameter, HELD_PRIMES)[0],
                "discriminant_feature": feature,
            }
        )
        enriched.append(record)
    selected = select_conductors(enriched, args.conductor_keep)
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
    print(f"conductors complete={len(completed)} subtarget={len(subtarget)}", flush=True)

    stages = (
        ("H5000", 5_000, args.h5000_keep, args.h5000_timeout),
        ("H50000", 50_000, args.h50000_keep, args.h50000_timeout),
        ("H250000", 250_000, args.h250000_keep, args.h250000_timeout),
    )
    current = select_points(completed, args.h5000_keep)
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
                    -row["slice_source_count"],
                    -Decimal(row["held_score"]),
                )
            )
            current = current[:keep]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = []
            for row in current:
                forced = tuple(
                    (Q(point["x"]), Q(point["y"]))
                    for point in row["forced_quartic_points"]
                )
                futures.append(
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
                        forced_quartic_points=forced,
                    )
                )
            for row, future in zip(current, futures):
                row.setdefault("point_stages", {})[name] = future.result()
        maximum = max(
            (row["point_stages"][name].get("stable_numerical_rank", -1) for row in current),
            default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={maximum}", flush=True)

    maximum_rank = -1
    finite_attempts = []
    target_hits = []
    rank20_plus = []
    for row in selected:
        for name, stage in row.get("point_stages", {}).items():
            maximum_rank = max(maximum_rank, int(stage.get("stable_numerical_rank", -1)))
            rank = stage.get("finite_reduction_attempt", {}).get(
                "certified_algebraic_rank_lower_bound"
            )
            if rank is None:
                continue
            entry = {"parameter": row["parameter"], "stage": name, "rank": rank}
            finite_attempts.append(entry)
            if rank >= 20:
                rank20_plus.append(
                    {
                        **entry,
                        "conductor": row["conductor_phase"]["conductor"],
                        "log_conductor": row["conductor_phase"]["log_conductor"],
                    }
                )
            if rank >= 30 or (
                rank >= 21
                and row["conductor_phase"]["below_strict_log_conductor_target_numerically"]
            ):
                target_hits.append(rank20_plus[-1] if rank >= 20 else entry)

    artifact = {
        "schema_version": 1,
        "status": "completed bounded exact accidental-slice experiment",
        "mathematical_status": (
            "rank-20 anchor and any reported promoted ranks are exact finite-reduction "
            "certificates; numerical ranks and bounded misses are experiments"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "anchor": anchor,
        "slice_search": {
            "slice_count": len(slices),
            "candidate_count": len(candidates),
            "height_bound": args.slice_height,
            "slices": slices,
        },
        "candidate_screen": {
            "nonsingular_count": len(enriched),
            "conductor_selected": len(selected),
            "conductor_completed": len(completed),
            "below_strict_target": len(subtarget),
        },
        "point_search": {
            "stages": [
                {"name": name, "height": height, "keep": keep, "timeout": timeout}
                for name, height, keep, timeout in stages
            ],
            "maximum_stable_numerical_rank": maximum_rank,
            "finite_reduction_attempts": finite_attempts,
            "rank20_plus": rank20_plus,
            "same_height_retries": 0,
        },
        "records": selected,
        "parameters": {
            key: value for key, value in vars(args).items() if key != "output"
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "reproducing_command": " ".join(sys.argv),
        },
        "software": {"python": platform.python_version()},
        "timings": {"total_wall_seconds": time.monotonic() - started},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "anchor": anchor,
            "slice_count": len(slices),
            "candidate_screen": artifact["candidate_screen"],
            "points": artifact["point_search"],
            "records": selected,
            "target": artifact["target"],
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"complete max_rank={maximum_rank} rank20_plus={len(rank20_plus)} "
        f"target_hits={len(target_hits)} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

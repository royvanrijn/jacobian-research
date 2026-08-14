#!/usr/bin/env python3
"""Bounded alternate-cover search on a specialization of the rank-13 family.

The input parameter is the base-change coordinate ``u`` for the Mestre family
with centers ``(0,23,93,128,133,175)``.  A direct quartic search supplies an
exact finite-reduction basis.  The program then searches deterministic
low-weight and full-coset alternate 2-cover charts in parallel, maps every
returned point exactly to the Jacobian, and certifies the augmented subgroup
dimension over finite reductions.

This is a bounded search, not a complete 2-descent and not a rank upper bound.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any

from alternate_quartic_covers import alternate_cover, short_subset_sum
from mestre_rank13_02393128133175 import (
    CONSTRUCTION,
    base_changed_short_jacobian_coefficients,
    base_parameter,
    known_jacobian_points,
)
from search_mestre_root_tuple_scale import (
    bounded_quartic_points,
    canonical_signless_points,
    point_digest,
    quartic_point_to_jacobian,
    sha256_file,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from search_nagao_u135_alternate_covers import (
    CoverPlan,
    best_cross_ratio_charts,
    enumerate_cover_plans,
    full_coset_identity_frontier,
    run_chart,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_ROOT = Path(
    "artifacts/local/elliptic-curves/mestre-02393128133175-basechange-crt-v1"
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_text(point[0]), "y": rational_text(point[1])}


def candidate_identifier(parameter_u: Fraction) -> str:
    parameter_u = Q(parameter_u)
    return f"u{parameter_u.numerator}_{parameter_u.denominator}"


def direct_pool(
    parameter_u: Fraction, *, height_bound: int, timeout: float, stack_bytes: int
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, str], ...],
    dict[str, Any],
]:
    parameter_t = base_parameter(parameter_u)
    coefficients = base_changed_short_jacobian_coefficients(parameter_u)
    by_x: dict[Fraction, tuple[Fraction, Fraction]] = {
        point[0]: point for point in known_jacobian_points(parameter_u)
    }
    raw = bounded_quartic_points(
        CONSTRUCTION.primitive_quartic_coefficients(parameter_t),
        height_bound=height_bound,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    signless = canonical_signless_points(raw)
    for quartic_point in signless:
        point = quartic_point_to_jacobian(CONSTRUCTION, parameter_t, quartic_point)
        by_x.setdefault(point[0], point)
    pool = tuple(by_x.values())
    certificate = mod3_independence_certificate(coefficients, pool, prime_bound=499)
    return (
        coefficients,
        pool,
        tuple(point_record(point) for point in signless),
        certificate,
    )


def build_plans(
    coefficients: tuple[Fraction, ...],
    basis: tuple[tuple[Fraction, Fraction], ...],
    *,
    low_weight_keep: int,
    full_coset_keep: int,
    charts_per_cover: int,
) -> tuple[CoverPlan, ...]:
    selected: dict[tuple[int, ...], CoverPlan] = {}
    low_weight = enumerate_cover_plans(
        coefficients,
        basis,
        maximum_subset_weight=2,
        charts_per_cover=charts_per_cover,
    )
    for plan in low_weight[:low_weight_keep]:
        selected[plan.subset_indices] = plan
    frontier = full_coset_identity_frontier(
        coefficients, basis, retain_count=full_coset_keep + len(selected)
    )
    added = 0
    for _, subset_indices in frontier:
        if subset_indices in selected:
            continue
        base_point = short_subset_sum(coefficients, basis, subset_indices)
        if base_point is None:
            raise AssertionError("a retained nonzero cover class vanished")
        cover = alternate_cover(coefficients, base_point)
        selected[subset_indices] = CoverPlan(
            subset_indices,
            cover,
            best_cross_ratio_charts(cover, basis, count=charts_per_cover),
        )
        added += 1
        if added == full_coset_keep:
            break
    if added != full_coset_keep:
        raise AssertionError("the full-coset frontier supplied too few new covers")
    return tuple(selected.values())


def run_jobs(
    jobs: list[tuple[CoverPlan, Any]],
    *,
    workers: int,
    stage: str,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> list[tuple[CoverPlan, Any, dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]]:
    def one(job: tuple[CoverPlan, Any]):
        plan, chart = job
        record, points = run_chart(
            plan,
            chart,
            stage=stage,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        return plan, chart, record, points

    results = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(one, job) for job in jobs]
        for index, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            plan, _, record, points = result
            print(
                f"{stage} {index}/{len(jobs)} {plan.identifier} "
                f"status={record['status']} points={len(points)}",
                flush=True,
            )
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--u", type=Q, required=True)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_INPUT_ROOT / "alternate-covers")
    parser.add_argument("--direct-height", type=int, default=200_000)
    parser.add_argument("--pilot-height", type=int, default=100_000)
    parser.add_argument("--deep-height", type=int, default=1_000_000)
    parser.add_argument("--direct-timeout", type=float, default=120.0)
    parser.add_argument("--chart-timeout", type=float, default=30.0)
    parser.add_argument("--deep-timeout", type=float, default=60.0)
    parser.add_argument("--low-weight-keep", type=int, default=16)
    parser.add_argument("--full-coset-keep", type=int, default=12)
    parser.add_argument("--charts-per-cover", type=int, default=3)
    parser.add_argument("--deep-productive-keep", type=int, default=12)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.u == 0:
        raise SystemExit("u=0 is the pole of the base change")
    if not 1 <= args.workers <= 8:
        raise SystemExit("--workers must lie in [1,8]")
    parameter_u = Q(args.u)
    identifier = candidate_identifier(parameter_u)
    input_root = ROOT / args.input_root
    conductor_path = input_root / "conductor-records" / f"{identifier}.json"
    conductor_record = json.loads(conductor_path.read_text()) if conductor_path.exists() else None

    coefficients, pool, quartic_records, baseline = direct_pool(
        parameter_u,
        height_bound=args.direct_height,
        timeout=args.direct_timeout,
        stack_bytes=args.stack_bytes,
    )
    baseline_rank = baseline["combined_exact_rank_over_F3"]
    indices = tuple(index - 1 for index in baseline["independent_subset_indices_one_based"])
    basis = tuple(pool[index] for index in indices)
    if len(basis) != baseline_rank:
        raise AssertionError("the finite-reduction basis length changed")
    plans = build_plans(
        coefficients,
        basis,
        low_weight_keep=args.low_weight_keep,
        full_coset_keep=args.full_coset_keep,
        charts_per_cover=args.charts_per_cover,
    )
    jobs = [(plan, chart) for plan in plans for chart in plan.charts]
    pilot = run_jobs(
        jobs,
        workers=args.workers,
        stage="pilot",
        height_bound=args.pilot_height,
        timeout=args.chart_timeout,
        stack_bytes=args.stack_bytes,
    )

    by_x = {point[0]: point for point in pool}
    pilot_yields = []
    for plan, chart, record, points in pilot:
        before = len(by_x)
        for point in points:
            by_x.setdefault(point[0], point)
        pilot_yields.append((len(by_x) - before, plan, chart, record))
    productive = sorted(
        (item for item in pilot_yields if item[0] > 0),
        key=lambda item: (-item[0], item[1].score, item[2].score),
    )[: args.deep_productive_keep]
    deep = run_jobs(
        [(plan, chart) for _, plan, chart, _ in productive],
        workers=min(args.workers, max(1, len(productive))),
        stage="deep",
        height_bound=args.deep_height,
        timeout=args.deep_timeout,
        stack_bytes=args.stack_bytes,
    ) if productive else []
    for _, _, _, points in deep:
        for point in points:
            by_x.setdefault(point[0], point)
    augmented_pool = tuple(by_x.values())
    augmented = mod3_independence_certificate(coefficients, augmented_pool, prime_bound=799)
    final_rank = augmented["combined_exact_rank_over_F3"]

    def run_record(item: tuple[CoverPlan, Any, dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]):
        plan, chart, record, points = item
        return {
            **record,
            "cover_subset_indices_one_based": [index + 1 for index in plan.subset_indices],
            "normalizing_basis_indices_one_based": [index + 1 for index in chart.basis_indices],
            "matrix_a_b_c_d": list(chart.matrix),
            "exact_curve_points": [point_record(point) for point in points],
        }

    output = ROOT / args.output_root / identifier
    artifact = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed bounded alternate-cover search",
        "candidate": {
            "u": rational_text(parameter_u),
            "T": rational_text(base_parameter(parameter_u)),
            "global_curve": conductor_record.get("global_curve") if conductor_record else None,
        },
        "direct_search": {
            "height_bound": args.direct_height,
            "quartic_points": list(quartic_records),
            "pool_point_count_modulo_inverse": len(pool),
            "pool_point_sha256": point_digest(pool),
            "exact_rank_lower_bound": baseline_rank,
            "finite_reduction_certificate": baseline,
        },
        "basis": {
            "points": [point_record(point) for point in basis],
            "point_sha256": point_digest(basis),
        },
        "budget": {
            "low_weight_cover_count": args.low_weight_keep,
            "full_coset_cover_count": args.full_coset_keep,
            "charts_per_cover": args.charts_per_cover,
            "pilot_chart_count": len(pilot),
            "pilot_height": args.pilot_height,
            "pilot_timeout_seconds_each": args.chart_timeout,
            "productive_deep_chart_count": len(deep),
            "deep_height": args.deep_height,
            "deep_timeout_seconds_each": args.deep_timeout,
            "workers": args.workers,
            "stack_bytes_each": args.stack_bytes,
        },
        "pilot_runs": [run_record(item) for item in pilot],
        "deep_runs": [run_record(item) for item in deep],
        "result": {
            "augmented_pool_point_count_modulo_inverse": len(augmented_pool),
            "augmented_pool_point_sha256": point_digest(augmented_pool),
            "exact_specialization_rank_lower_bound": final_rank,
            "certified_rank_gain": final_rank - baseline_rank,
            "finite_reduction_certificate": augmented,
            "points": [point_record(point) for point in augmented_pool],
        },
        "bounded_search_is_not_a_rank_upper_bound": True,
        "provenance": {
            "script_path": str(Path(__file__).resolve().relative_to(ROOT)),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "conductor_input_path": str(conductor_path.relative_to(ROOT)) if conductor_path.exists() else None,
            "conductor_input_sha256": sha256_file(conductor_path) if conductor_path.exists() else None,
            "actual_command": " ".join(sys.argv),
        },
        "software": {"python": platform.python_version(), "platform": platform.platform()},
    }
    artifact["result_sha256"] = canonical_digest(
        {key: value for key, value in artifact.items() if key != "generated_at_utc"}
    )
    atomic_json(output / "summary.json", artifact)
    print(
        f"complete u={rational_text(parameter_u)} baseline={baseline_rank} "
        f"final={final_rank} pool={len(augmented_pool)} sha={artifact['result_sha256']}",
        flush=True,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run and summarize the five-fibre Kummer-shaped quotient experiment.

The output is deliberately a bounded relation fingerprint.  It compares
identically configured collector lanes for curves 351, 356, 376, 377, and
385, but it never promotes a materialized factor-base dimension to an
S-class, Selmer, or Mordell--Weil rank statement.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
COLLECTOR = ROOT / "elliptic-curves/cas/run_r17_kummer_quotient_sclass_collector.sage"
CURVES = (351, 356, 376, 377, 385)
OBJECTIVES = ("generic", "full-known")
SCHEMA = "elliptic-curves.r17-kummer-quotient-sclass-suite.v1"


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_csv_ints(value: str) -> tuple[int, ...]:
    values = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not values or any(value not in CURVES for value in values):
        raise ValueError(f"curves must be drawn from {CURVES}")
    return values


def parse_csv_objectives(value: str) -> tuple[str, ...]:
    values = tuple(item.strip() for item in value.split(",") if item.strip())
    if not values or any(value not in OBJECTIVES for value in values):
        raise ValueError(f"objectives must be drawn from {OBJECTIVES}")
    return values


def mean(values: Sequence[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median(values: Sequence[int]) -> float | None:
    return statistics.median(values) if values else None


def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for offset in range(start, end):
            ranks[order[offset]] = rank
        start = end
    return ranks


def pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    left_centered = [value - left_mean for value in left]
    right_centered = [value - right_mean for value in right]
    denominator = math.sqrt(
        sum(value * value for value in left_centered)
        * sum(value * value for value in right_centered)
    )
    if denominator == 0:
        return None
    return sum(a * b for a, b in zip(left_centered, right_centered)) / denominator


def correlations(records: Sequence[dict]) -> dict[str, object]:
    gains = [
        float(record["certified_displayed_mw_gain_over_generic_mw17"])
        for record in records
    ]
    metrics = {
        "objective_materialized_quotient_dimension": [
            float(record["objective_materialized_quotient_dimension"])
            for record in records
        ],
        "objective_relation_rank_gain": [
            float(record["objective_relation_rank_gain"]) for record in records
        ],
        "objective_remaining_fraction": [
            float(record["objective_remaining_fraction"]) for record in records
        ],
        "smooth_reductions_per_1000_attempts": [
            float(record["smooth_reductions_per_1000_attempts"])
            for record in records
        ],
    }
    return {
        name: {
            "pearson_with_certified_displayed_mw_gain": pearson(gains, values),
            "spearman_with_certified_displayed_mw_gain": pearson(
                average_ranks(gains), average_ranks(values)
            ),
        }
        for name, values in metrics.items()
    }


def closed_relation_records(data: dict) -> Iterable[dict]:
    yield from data["accepted_quotient_relations"]
    for collision in data["large_prime_collisions"]:
        yield collision


def summarize_checkpoint(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "elliptic-curves.r17-kummer-quotient-sclass-collector.v2":
        raise ValueError(f"unexpected collector schema in {path}")
    objective = data["objective"]
    objective_projection = data["projections"][objective]
    closed = list(closed_relation_records(data))
    target_kinds = Counter(
        relation.get("target", {}).get("kind", "large_prime_cycle")
        for relation in data["accepted_quotient_relations"]
    )
    target_widths = Counter()
    for relation in data["accepted_quotient_relations"]:
        target = relation.get("target", {})
        width = int(target.get("actual_column_count", 1))
        target_widths[width] += 1
    generic_weights = [
        int(record["generic_residual_weight_before"])
        for record in closed
        if "generic_residual_weight_before" in record
    ]
    full_known_weights = [
        int(record["full_known_residual_weight_before"])
        for record in closed
        if "full_known_residual_weight_before" in record
    ]
    base_row_weights = []
    exceptional_support_weights = []
    for record in closed:
        base_hex = record.get(
            "base_row_mask_hex", record.get("combined_base_row_mask_hex")
        )
        exceptional_hex = record.get(
            "exceptional_parity_mask_hex",
            record.get("combined_exceptional_parity_mask_hex"),
        )
        if base_hex is not None:
            base_row_weights.append(int(base_hex, 16).bit_count())
        if exceptional_hex is not None:
            exceptional_support_weights.append(
                int(exceptional_hex, 16).bit_count()
            )
    attempts = int(data["attempts_completed"])
    companion_totals = Counter()
    for stats in data["strategy_statistics"].values():
        for key in (
            "generic_companion_terms",
            "exceptional_companion_terms",
            "negative_exponent_terms",
            "even_exponent_terms",
        ):
            companion_totals[key] += int(stats[key])
    companion_count = (
        companion_totals["generic_companion_terms"]
        + companion_totals["exceptional_companion_terms"]
    )
    return {
        "curve_id": int(data["curve_id"]),
        "objective": objective,
        "status": data["status"],
        "certified_displayed_mw_gain_over_generic_mw17": int(
            data["certified_displayed_mw_gain_over_generic_mw17"]
        ),
        "checkpoint": str(path),
        "checkpoint_sha256": file_sha256(path),
        "attempts_completed": attempts,
        "elapsed_seconds": float(data["elapsed_seconds"]),
        "factor_base_width": len(data["factor_base_columns"]),
        "objective_presentation_width": int(
            objective_projection["presentation_width"]
        ),
        "objective_initial_mod2_rank": int(
            objective_projection["initial_mod2_rank"]
        ),
        "objective_current_mod2_rank": int(
            objective_projection["current_mod2_rank"]
        ),
        "objective_relation_rank_gain": int(
            objective_projection["relation_rank_gain"]
        ),
        "objective_materialized_quotient_dimension": int(
            objective_projection["materialized_quotient_dimension"]
        ),
        "objective_remaining_fraction": (
            objective_projection["materialized_quotient_dimension"]
            / objective_projection["presentation_width"]
        ),
        "generic_projection": data["projections"]["generic"],
        "full_known_projection": data["projections"]["full-known"],
        "fully_factored_reduction_count": int(
            data["fully_factored_reduction_count"]
        ),
        "smooth_reduction_count": int(data["smooth_reduction_count"]),
        "smooth_reductions_per_1000_attempts": (
            1000 * int(data["smooth_reduction_count"]) / attempts
            if attempts
            else 0.0
        ),
        "closed_relation_count": len(closed),
        "unmatched_large_prime_partial_count": int(
            data["unmatched_large_prime_partial_count"]
        ),
        "accepted_target_kind_counts": dict(sorted(target_kinds.items())),
        "accepted_smooth_target_width_counts": {
            str(width): count for width, count in sorted(target_widths.items())
        },
        "base_row_weight_mean": mean(base_row_weights),
        "base_row_weight_median": median(base_row_weights),
        "exceptional_support_weight_mean": mean(exceptional_support_weights),
        "exceptional_support_weight_median": median(exceptional_support_weights),
        "generic_residual_weight_mean": mean(generic_weights),
        "generic_residual_weight_median": median(generic_weights),
        "full_known_residual_weight_mean": mean(full_known_weights),
        "full_known_residual_weight_median": median(full_known_weights),
        "generic_companion_term_count": companion_totals[
            "generic_companion_terms"
        ],
        "exceptional_companion_term_count": companion_totals[
            "exceptional_companion_terms"
        ],
        "exceptional_companion_fraction": (
            companion_totals["exceptional_companion_terms"] / companion_count
            if companion_count
            else None
        ),
        "negative_exponent_term_count": companion_totals[
            "negative_exponent_terms"
        ],
        "even_exponent_term_count": companion_totals["even_exponent_terms"],
        "strategy_statistics": data["strategy_statistics"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curves", default=",".join(map(str, CURVES)))
    parser.add_argument("--objectives", default=",".join(OBJECTIVES))
    parser.add_argument("--factor-base-bound", type=int, default=5000)
    parser.add_argument("--attempts", type=int, default=10000)
    parser.add_argument("--timeout-seconds-per-run", type=float, default=60.0)
    parser.add_argument("--max-target-columns", type=int, default=1)
    parser.add_argument("--companion-strategies", default="single,pair,sparse")
    parser.add_argument("--sparse-min", type=int, default=3)
    parser.add_argument("--sparse-max", type=int, default=6)
    parser.add_argument("--companion-exponent-radius", type=int, default=1)
    parser.add_argument("--signed-companion-exponents", action="store_true")
    parser.add_argument("--exceptional-companion-weight", type=int, default=4)
    parser.add_argument("--max-s-companions", type=int, default=2)
    parser.add_argument("--direction-radius", type=int, default=16)
    parser.add_argument(
        "--reduction-engine",
        choices=("idealred", "idealredmodpower2"),
        default="idealred",
    )
    parser.add_argument("--large-prime-bound", type=int, default=None)
    parser.add_argument("--max-large-primes", type=int, choices=(1, 2), default=1)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--sage", default=None, help="Sage executable (default: PATH)")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    try:
        curves = parse_csv_ints(args.curves)
        objectives = parse_csv_objectives(args.objectives)
    except ValueError as error:
        parser.error(str(error))
    if args.factor_base_bound < 2 or args.attempts <= 0:
        parser.error("factor-base bound and attempts must be positive")
    if args.max_target_columns <= 0:
        parser.error("the maximum target width must be positive")
    if args.timeout_seconds_per_run <= 0:
        parser.error("per-run timeout must be positive")
    if args.summary.exists() and not args.overwrite:
        raise FileExistsError(args.summary)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    checkpoints = [
        args.output_dir
        / f"curve{curve}-{objective}-b{args.factor_base_bound}-a{args.attempts}.json"
        for objective in objectives
        for curve in curves
    ]
    commands = []
    discovered_sage = shutil.which("sage")
    if not args.summarize_only and not (args.sage or discovered_sage):
        raise FileNotFoundError("Sage was not found on PATH; pass --sage")
    sage = args.sage or discovered_sage or "sage"
    for objective in objectives:
        for curve in curves:
            checkpoint = (
                args.output_dir
                / f"curve{curve}-{objective}-b{args.factor_base_bound}-a{args.attempts}.json"
            )
            command = [
                sage,
                "-python",
                str(COLLECTOR),
                "--curve-id",
                str(curve),
                "--objective",
                objective,
                "--factor-base-bound",
                str(args.factor_base_bound),
                "--attempts",
                str(args.attempts),
                "--timeout-seconds",
                str(args.timeout_seconds_per_run),
                "--max-target-columns",
                str(args.max_target_columns),
                "--companion-strategies",
                args.companion_strategies,
                "--sparse-min",
                str(args.sparse_min),
                "--sparse-max",
                str(args.sparse_max),
                "--companion-exponent-radius",
                str(args.companion_exponent_radius),
                "--exceptional-companion-weight",
                str(args.exceptional_companion_weight),
                "--max-s-companions",
                str(args.max_s_companions),
                "--direction-radius",
                str(args.direction_radius),
                "--reduction-engine",
                args.reduction_engine,
                "--max-large-primes",
                str(args.max_large_primes),
                "--seed",
                str(args.seed),
                "--checkpoint",
                str(checkpoint),
            ]
            if args.large_prime_bound is not None:
                command.extend(
                    ["--large-prime-bound", str(args.large_prime_bound)]
                )
            if args.signed_companion_exponents:
                command.append("--signed-companion-exponents")
            if args.overwrite:
                command.append("--overwrite")
            commands.append(command)
            if not args.summarize_only:
                subprocess.run(command, cwd=ROOT, check=True)

    missing = [path for path in checkpoints if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "missing collector checkpoints: " + ", ".join(map(str, missing))
        )
    records = [summarize_checkpoint(path) for path in checkpoints]
    by_objective = {
        objective: {
            "records": [record for record in records if record["objective"] == objective],
            "bounded_correlations": correlations(
                [record for record in records if record["objective"] == objective]
            ),
            "interpretation": (
                "Only relation-rank gain and closed-row structure measure collector "
                "progress. Starting/final raw dimensions are confounded by unequal "
                "factor-base widths; the generic presentation also adjoins one "
                "formal column per displayed exceptional half-ideal."
            ),
        }
        for objective in objectives
    }
    output = {
        "schema": SCHEMA,
        "status": "BOUNDED_COMPARATIVE_RELATION_FINGERPRINT_ONLY",
        "curves": list(curves),
        "objectives": list(objectives),
        "collector": {
            "path": str(COLLECTOR.relative_to(ROOT)),
            "sha256": file_sha256(COLLECTOR),
        },
        "suite_program": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": file_sha256(Path(__file__)),
        },
        "parameters": {
            "factor_base_bound": args.factor_base_bound,
            "attempts": args.attempts,
            "timeout_seconds_per_run": args.timeout_seconds_per_run,
            "max_target_columns": args.max_target_columns,
            "companion_strategies": args.companion_strategies,
            "sparse_min": args.sparse_min,
            "sparse_max": args.sparse_max,
            "companion_exponent_radius": args.companion_exponent_radius,
            "signed_companion_exponents": args.signed_companion_exponents,
            "exceptional_companion_weight": args.exceptional_companion_weight,
            "max_s_companions": args.max_s_companions,
            "direction_radius": args.direction_radius,
            "reduction_engine": args.reduction_engine,
            "large_prime_bound": args.large_prime_bound,
            "max_large_primes": args.max_large_primes,
            "seed": args.seed,
        },
        "commands": commands,
        "comparisons_by_objective": by_objective,
        "claim_boundary": [
            "The displayed-MW-gain column is an exact lower-bound quotient dimension from the pinned public points, not an exact specialization rank.",
            "Every collected closed relation is exact, but the displayed presentation dimensions are bounded factor-base fingerprints without a generation proof.",
            "Five controls cannot establish a predictive law; the reported correlations are descriptive diagnostics only.",
            "No S-class dimension, Selmer dimension, rank upper bound, or causal middle layer is asserted.",
        ],
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.summary.with_suffix(args.summary.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.summary)
    print(
        "R17KUMMERQSUITE|"
        f"runs={len(records)}|status={output['status']}|summary={args.summary}",
        flush=True,
    )


if __name__ == "__main__":
    main()

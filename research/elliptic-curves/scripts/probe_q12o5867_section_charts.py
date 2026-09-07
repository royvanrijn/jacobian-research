#!/usr/bin/env python3
# <!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE f7a8c94736f1b44f -->
"""Gate-protected affine-chart search normalized by baseline sections."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_point_search import (  # noqa: E402
    affine_substitute_polynomial,
    completed_square_coefficients,
    exact_escape_records,
    integral_square_scaled_coefficients,
    novel_points_up_to_sign,
    parse_ratpoints_abscissae,
    point_record,
    points_from_completed_square_abscissa,
    sign_key,
)
from ecsearch.q12o5867_rank_jump_registry import sha256_file  # noqa: E402
from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402
from elkies_residual_selmer_gate import require_gate_for_specialization  # noqa: E402


Q = Fraction
DEFAULT_RATPOINTS = REPOSITORY / "tmp/ratpoints/root/usr/bin/ratpoints"
DEFAULT_LIBRARY = REPOSITORY / "tmp/ratpoints/root/usr/lib/x86_64-linux-gnu"


def parse_relation_primes(text: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("relation primes must be integers") from error
    if not answer:
        raise argparse.ArgumentTypeError("at least one relation prime is required")
    return answer


def fraction_record(value: Fraction) -> str:
    return str(value)


def chart_inventory(
    baseline_x: tuple[Fraction, ...],
    pair_mode: str,
    include_multiplicative: bool,
) -> tuple[tuple[str, Fraction, Fraction], ...]:
    charts: list[tuple[str, Fraction, Fraction]] = []
    if pair_mode == "all":
        pairs = (
            (left, right)
            for left in range(len(baseline_x))
            for right in range(left + 1, len(baseline_x))
        )
    elif pair_mode == "star":
        pairs = ((0, right) for right in range(1, len(baseline_x)))
    elif pair_mode == "adjacent":
        pairs = ((left, left + 1) for left in range(len(baseline_x) - 1))
    else:
        raise ValueError(f"unknown pair mode {pair_mode!r}")
    for left, right in pairs:
        center = baseline_x[left]
        scale = baseline_x[right] - center
        if scale:
            charts.append((f"pair-{left:02d}-{right:02d}", center, scale))
    if include_multiplicative:
        for index, value in enumerate(baseline_x):
            if value:
                charts.append((f"scale-{index:02d}", Q(0), value))
    return tuple(charts)


def run_chart(
    *,
    identifier: str,
    center: Fraction,
    scale: Fraction,
    completed_square: tuple[int, ...],
    executable: Path,
    library: Path,
    height: int,
    denominator_bound: int,
    timeout: float,
    raw_directory: Path,
    overwrite: bool,
) -> tuple[dict[str, Any], tuple[Fraction, ...]]:
    transformed = affine_substitute_polynomial(
        completed_square, center, scale
    )
    integral, ordinate_scale = integral_square_scaled_coefficients(transformed)
    command = [
        str(executable.resolve()),
        " ".join(str(value) for value in integral),
        str(height),
        "-du",
        str(denominator_bound),
        "-q",
        "-y",
    ]
    environment = os.environ.copy()
    old_library = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        str(library.resolve())
        if not old_library
        else str(library.resolve()) + os.pathsep + old_library
    )
    started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            env=environment,
            timeout=timeout,
        )
        stdout = completed.stdout
        stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as error:
        timed_out = True
        stdout = (
            error.stdout.decode()
            if isinstance(error.stdout, bytes)
            else (error.stdout or "")
        )
        stderr = (
            error.stderr.decode()
            if isinstance(error.stderr, bytes)
            else (error.stderr or "")
        )
        returncode = None
    seconds = time.monotonic() - started
    raw_path = raw_directory / f"{identifier}.out"
    mode = "w" if overwrite else "x"
    with raw_path.open(mode) as handle:
        handle.write(stdout)
    transformed_x = ()
    mapped_x = ()
    if not timed_out and returncode == 0:
        transformed_x = parse_ratpoints_abscissae(stdout)
        mapped_x = tuple(center + scale * value for value in transformed_x)
    record = {
        "identifier": identifier,
        "x_change": {
            "formula": "x=center+scale*X",
            "center": fraction_record(center),
            "scale": fraction_record(scale),
        },
        "transformed_coefficients_low_to_high": [
            fraction_record(value) for value in transformed
        ],
        "integral_coefficients_low_to_high": [str(value) for value in integral],
        "ordinate_square_clearing_scale": str(ordinate_scale),
        "command": command,
        "seconds": seconds,
        "timed_out": timed_out,
        "returncode": returncode,
        "raw_output_path": str(raw_path.resolve()),
        "raw_output_sha256": sha256(stdout.encode()).hexdigest(),
        "raw_stderr": stderr,
        "transformed_abscissa_count": len(transformed_x),
        "mapped_abscissae": [fraction_record(value) for value in mapped_x],
    }
    return record, mapped_x


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--residual-selmer-gate", type=Path, required=True)
    parser.add_argument("--pair-mode", choices=("all", "star", "adjacent"), default="all")
    parser.add_argument("--include-multiplicative", action="store_true")
    parser.add_argument("--chart-limit", type=int)
    parser.add_argument("--height", type=int, required=True)
    parser.add_argument("--denominator-bound", type=int, required=True)
    parser.add_argument("--per-chart-timeout", type=float, default=30.0)
    parser.add_argument("--ratpoints", type=Path, default=DEFAULT_RATPOINTS)
    parser.add_argument("--ratpoints-library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--relation-primes", type=parse_relation_primes, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--raw-directory", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.height < 1 or args.denominator_bound < 1:
        parser.error("height and denominator bounds must be positive")
    if args.per_chart_timeout <= 0 or args.per_chart_timeout > 60:
        parser.error("--per-chart-timeout must lie in (0,60]")
    if args.chart_limit is not None and args.chart_limit < 1:
        parser.error("--chart-limit must be positive")
    sys.set_int_max_str_digits(0)

    specialization = json.loads(args.input.read_text())
    if specialization.get("status") != "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND":
        raise SystemExit("input is not an exact certified q12o5867 specialization")
    minimal = specialization["global_minimal_specialization"]
    model = tuple(Q(value) for value in minimal["model"])
    baseline = tuple((Q(x), Q(y)) for x, y in minimal["points"])
    if len(baseline) != 17 or any(
        not is_on_weierstrass_curve(model, point) for point in baseline
    ):
        raise AssertionError("the serialized exact baseline is invalid")
    baseline_x = tuple(point[0] for point in baseline)
    charts = chart_inventory(
        baseline_x, args.pair_mode, args.include_multiplicative
    )
    if args.chart_limit is not None:
        charts = charts[: args.chart_limit]
    gate = require_gate_for_specialization(
        args.residual_selmer_gate,
        specialization,
        requested_search_limits={
            "height": args.height,
            "denominator_bound": args.denominator_bound,
            "chart_count": len(charts),
            "wall_seconds_per_chart": args.per_chart_timeout,
        },
    )
    completed_square = completed_square_coefficients(model)
    raw_directory = args.raw_directory or args.output.with_suffix("")
    raw_directory.mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    records = []
    mapped_abscissae = set()
    for position, (identifier, center, scale) in enumerate(charts, 1):
        record, mapped = run_chart(
            identifier=identifier,
            center=center,
            scale=scale,
            completed_square=completed_square,
            executable=args.ratpoints,
            library=args.ratpoints_library,
            height=args.height,
            denominator_bound=args.denominator_bound,
            timeout=args.per_chart_timeout,
            raw_directory=raw_directory,
            overwrite=args.overwrite,
        )
        records.append(record)
        mapped_abscissae.update(mapped)
        if position % 25 == 0 or position == len(charts):
            print(
                f"charts={position}/{len(charts)} mapped_x={len(mapped_abscissae)}",
                flush=True,
            )
    searched = tuple(
        point
        for x_coordinate in sorted(mapped_abscissae)
        for point in points_from_completed_square_abscissa(model, x_coordinate)
    )
    candidates = novel_points_up_to_sign(model, baseline, searched)
    baseline_keys = {sign_key(model, point) for point in baseline}
    searched_keys = {sign_key(model, point) for point in searched}
    rediscovered = len(baseline_keys & searched_keys)
    escape = exact_escape_records(
        specialization,
        candidates,
        args.relation_primes,
        args.reduction_prime_bound,
    )
    timed_out_count = sum(bool(record["timed_out"]) for record in records)
    error_count = sum(
        not record["timed_out"] and record["returncode"] != 0
        for record in records
    )
    result = {
        "schema": "elliptic-curves.q12o5867-section-normalized-ratpoints-probe.v1",
        "status": (
            "PASS_BOUNDED_SECTION_NORMALIZED_RATPOINTS_PROBE"
            if not timed_out_count and not error_count
            else "PARTIAL_BOUNDED_SECTION_NORMALIZED_RATPOINTS_PROBE"
        ),
        "input_specialization_artifact": str(args.input.resolve()),
        "input_specialization_sha256": sha256_file(args.input),
        "residual_selmer_gate": {
            "path": str(args.residual_selmer_gate.resolve()),
            "sha256": sha256_file(args.residual_selmer_gate),
            "status": gate["status"],
        },
        "parameter": specialization["parameter"],
        "global_minimal_model": [str(value) for value in model],
        "completed_square_coefficients_low_to_high": [
            str(value) for value in completed_square
        ],
        "chart_construction": {
            "pair_mode": args.pair_mode,
            "include_multiplicative": args.include_multiplicative,
            "chart_limit": args.chart_limit,
            "chart_count": len(charts),
            "formula": "x=x_i+(x_j-x_i)*X for pair charts; x=x_i*X for scale charts",
        },
        "bounds": {
            "height": args.height,
            "denominator": args.denominator_bound,
            "per_chart_timeout_seconds": args.per_chart_timeout,
            "relation_primes": list(args.relation_primes),
            "reduction_prime_bound": args.reduction_prime_bound,
        },
        "engine": {
            "ratpoints_path": str(args.ratpoints.resolve()),
            "ratpoints_sha256": sha256_file(args.ratpoints),
            "library_path": str(args.ratpoints_library.resolve()),
        },
        "raw_directory": str(raw_directory.resolve()),
        "chart_records": records,
        "wall_seconds": time.monotonic() - started,
        "timed_out_chart_count": timed_out_count,
        "error_chart_count": error_count,
        "distinct_mapped_abscissa_count": len(mapped_abscissae),
        "exact_returned_affine_point_count": len(searched),
        "rediscovered_generic_baseline_count_up_to_sign": rediscovered,
        "all_returned_affine_points_exactly_verified": True,
        "novel_points_up_to_sign_on_minimal_model": [
            point_record(point) for point in candidates
        ],
        "finite_quotient_escape": escape,
        "promoted": escape["maximum_marginal_dimension"] >= 15,
        "claim_boundary": [
            "This is a finite union of bounded affine-chart searches, not a rank upper bound.",
            "Every mapped point is checked on the original exact minimal model.",
            "No candidate is promoted unless exact quotient gain reaches 15.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={result['status']}")
    print(f"charts={len(charts)}")
    print(f"rediscovered_baseline={rediscovered}")
    print(f"novel_points={len(candidates)}")
    print(f"maximum_marginal_dimension={escape['maximum_marginal_dimension']}")
    print(f"promoted={result['promoted']}")
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()

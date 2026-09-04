#!/usr/bin/env python3
# <!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE 7f8dffe58168acc8 -->
"""Gate-protected direct ratpoints search on one q12o5867 specialization.

The command refuses to start unless a completed unconditional residual
2-Selmer artifact for the identical parameter and minimal curve has passed the
rank-32 dimension gate.
"""

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


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_point_search import (  # noqa: E402
    completed_square_coefficients,
    exact_escape_records,
    novel_points_up_to_sign,
    parse_ratpoints_abscissae,
    point_record,
    points_from_completed_square_abscissa,
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--residual-selmer-gate", type=Path, required=True)
    parser.add_argument("--height", type=int, required=True)
    parser.add_argument("--denominator-bound", type=int, required=True)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--ratpoints", type=Path, default=DEFAULT_RATPOINTS)
    parser.add_argument("--ratpoints-library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--relation-primes", type=parse_relation_primes, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.height < 1 or args.denominator_bound < 1:
        parser.error("height and denominator bounds must be positive")
    if args.timeout <= 0:
        parser.error("timeout must be positive")
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
    gate = require_gate_for_specialization(
        args.residual_selmer_gate,
        specialization,
        requested_search_limits={
            "height": args.height,
            "denominator_bound": args.denominator_bound,
            "wall_seconds": args.timeout,
        },
    )
    coefficients = completed_square_coefficients(model)
    command = [
        str(args.ratpoints.resolve()),
        " ".join(str(value) for value in coefficients),
        str(args.height),
        "-du",
        str(args.denominator_bound),
        "-q",
        "-y",
    ]
    environment = os.environ.copy()
    old_library = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        str(args.ratpoints_library.resolve())
        if not old_library
        else str(args.ratpoints_library.resolve()) + os.pathsep + old_library
    )
    started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            env=environment,
            timeout=args.timeout,
        )
        raw_stdout = completed.stdout
        raw_stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as error:
        timed_out = True
        raw_stdout = (
            error.stdout.decode()
            if isinstance(error.stdout, bytes)
            else (error.stdout or "")
        )
        raw_stderr = (
            error.stderr.decode()
            if isinstance(error.stderr, bytes)
            else (error.stderr or "")
        )
        returncode = None
    search_seconds = time.monotonic() - started
    raw_path = args.raw_output or args.output.with_suffix(".raw.out")
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_mode = "w" if args.overwrite else "x"
    with raw_path.open(raw_mode) as handle:
        handle.write(raw_stdout)

    result: dict[str, object] = {
        "schema": "elliptic-curves.q12o5867-bounded-ratpoints-probe.v1",
        "input_specialization_artifact": str(args.input.resolve()),
        "input_specialization_sha256": sha256_file(args.input),
        "residual_selmer_gate": {
            "path": str(args.residual_selmer_gate.resolve()),
            "sha256": sha256_file(args.residual_selmer_gate),
            "status": gate["status"],
        },
        "parameter": specialization["parameter"],
        "global_minimal_model": [str(value) for value in model],
        "completed_square": {
            "equation": "Y^2=4*x^3+b2*x^2+2*b4*x+b6; Y=2*y+a1*x+a3",
            "coefficients_low_to_high": [str(value) for value in coefficients],
        },
        "bounds": {
            "numerator_absolute_and_default_height": args.height,
            "denominator": args.denominator_bound,
            "wall_timeout_seconds": args.timeout,
            "reduction_prime_bound": args.reduction_prime_bound,
            "relation_primes": list(args.relation_primes),
        },
        "engine": {
            "command": command,
            "ratpoints_path": str(args.ratpoints.resolve()),
            "ratpoints_sha256": sha256_file(args.ratpoints),
            "library_path": str(args.ratpoints_library.resolve()),
        },
        "search_seconds": search_seconds,
        "timed_out": timed_out,
        "returncode": returncode,
        "raw_stdout_path": str(raw_path.resolve()),
        "raw_stdout_sha256": sha256(raw_stdout.encode()).hexdigest(),
        "raw_stderr": raw_stderr,
        "claim_boundary": [
            "This is a bounded rational-point search, not a rank upper bound.",
            "Every reported affine point is reconstructed and checked exactly.",
            "No point is promoted unless its exact quotient gain is certified.",
        ],
    }
    if timed_out:
        result["status"] = "TIMEOUT_BOUNDED_RATPOINTS_PROBE"
    elif returncode != 0:
        result["status"] = "ERROR_BOUNDED_RATPOINTS_PROBE"
    else:
        abscissae = parse_ratpoints_abscissae(raw_stdout)
        searched = tuple(
            point
            for x_coordinate in abscissae
            for point in points_from_completed_square_abscissa(model, x_coordinate)
        )
        candidates = novel_points_up_to_sign(model, baseline, searched)
        escape = exact_escape_records(
            specialization,
            candidates,
            args.relation_primes,
            args.reduction_prime_bound,
        )
        result.update(
            {
                "status": "PASS_BOUNDED_RATPOINTS_PROBE",
                "finite_abscissa_count": len(abscissae),
                "exact_returned_affine_point_count": len(searched),
                "all_returned_affine_points_exactly_verified": True,
                "novel_points_up_to_sign_on_minimal_model": [
                    point_record(point) for point in candidates
                ],
                "finite_quotient_escape": escape,
                "promoted": escape["maximum_marginal_dimension"] >= 15,
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={result['status']}")
    print(f"search_seconds={search_seconds:.6f}")
    print(f"output={args.output.resolve()}")
    if result["status"] == "PASS_BOUNDED_RATPOINTS_PROBE":
        print(f"finite_abscissa_count={result['finite_abscissa_count']}")
        print(f"novel_point_count={len(result['novel_points_up_to_sign_on_minimal_model'])}")
        print(
            "maximum_marginal_dimension="
            f"{result['finite_quotient_escape']['maximum_marginal_dimension']}"
        )
        print(f"promoted={result['promoted']}")
    if result["status"] == "ERROR_BOUNDED_RATPOINTS_PROBE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Single capped H=10^7 extension and saturated descent probe for Nagao u=42.

The script makes exactly one height-10,000,000 ``hyperellratpoints`` attempt,
capped at 120 seconds.  Whether it completes or times out, the already-pinned
height-1,000,000 subset remains available as a fallback.  The selected exact
point basis is saturated at primes below 20, checked exactly, replayed at two
height precisions, and passed once to PARI ``ellrank`` effort zero with a
60-second cap.  There are no retries.

Numerical height ranks and root-number parity remain search evidence.  An
equal PARI rank interval is a software computation, not a stored descent
certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys
from typing import Any, Sequence

from nagao_1994 import PRIMARY_SOURCE, rank13_base_changed_short_jacobian_coefficients
from pari_bridge import pari_version
from search_extra_points import gp_rational, gp_vector, parse_point_vector, parse_precisions, run_gp
from triage_nagao_rank13_finalists import (
    Finalist,
    ellrank_probe,
    exact_candidate_triage,
    height_matrix_replay,
    load_finalists,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
PARAMETER_U = 42
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/extend_nagao_u42_frontier.py"
)


def checkpoint_points(path: Path) -> tuple[tuple[Fraction, Fraction], ...]:
    data = json.loads(path.read_text())
    record = next(
        candidate
        for candidate in data["candidates"]
        if int(candidate["parameter_u"]) == PARAMETER_U
    )
    points = tuple(
        (Q(item["jacobian_x"]), Q(item["jacobian_y"]))
        for item in record["frontier_explicit_numerically_independent_subset"]
    )
    expected = int(record["frontier_stable_pool_numerical_rank"])
    if len(points) != expected:
        raise AssertionError("the u=42 checkpoint subset has the wrong length")
    coefficients = rank13_base_changed_short_jacobian_coefficients(Q(PARAMETER_U))
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a checkpoint point failed exact membership")
    return points


def saturate_exact_basis(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    prime_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    program = "\n".join(
        (
            "default(realprecision,100);",
            f"E=ellinit([{curve}]);",
            f"P=[{point_vector}];",
            "H0=ellheightmatrix(E,P);",
            "gettime();",
            f"S=ellsaturation(E,P,{prime_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("RETURNED_COUNT ",#S);',
            'print("ON_CURVE ",vecsum(vector(#S,i,ellisoncurve(E,S[i]))));',
            "H1=ellheightmatrix(E,S);",
            'print("ORIGINAL_DET ",matdet(H0));',
            'print("SATURATED_DET ",matdet(H1));',
            'print("DET_RATIO ",matdet(H0)/matdet(H1));',
            'print("SATURATED_POINTS_BEGIN");',
            "print(S);",
            'print("SATURATED_POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )

    def value(label: str) -> str:
        match = re.search(rf"^{label} (.+)$", output, re.MULTILINE)
        if match is None:
            raise AssertionError(f"PARI omitted {label}")
        return match.group(1)

    saturated_text = output.split("SATURATED_POINTS_BEGIN\n", 1)[1].split(
        "\nSATURATED_POINTS_END", 1
    )[0]
    saturated = parse_point_vector(saturated_text)
    returned_count = int(value("RETURNED_COUNT"))
    on_curve_count = int(value("ON_CURVE"))
    if len(saturated) != returned_count or on_curve_count != returned_count:
        raise AssertionError("the saturated basis did not parse or check completely")
    if any(not point_on_short_curve(coefficients, point) for point in saturated):
        raise AssertionError("a saturated point failed exact Python membership")
    record = {
        "status": "completed",
        "prime_bound_strict_upper_limit": prime_bound,
        "input_point_count": len(points),
        "returned_point_count": returned_count,
        "exact_returned_points_on_curve": on_curve_count,
        "pari_milliseconds": int(value("PARI_MILLISECONDS")),
        "wall_seconds": wall_seconds,
        "original_height_determinant": value("ORIGINAL_DET"),
        "saturated_height_determinant": value("SATURATED_DET"),
        "height_determinant_ratio": value("DET_RATIO"),
        "saturated_basis_sha256": point_digest(saturated),
        "saturated_basis": [
            {
                "jacobian_x": str(point[0]),
                "jacobian_y": str(point[1]),
                "exact_jacobian_membership_checked": True,
            }
            for point in saturated
        ],
        "scope_warning": (
            "PARI documents ellsaturation under a finite-index hypothesis; "
            "the full rank was unknown when this search was run"
        ),
    }
    return saturated, record


def load_u42_candidate(path: Path) -> Finalist:
    candidates = load_finalists(path, None, ())
    return next(candidate for candidate in candidates if candidate.parameter_u == 42)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scan-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_integer_u.json",
    )
    parser.add_argument(
        "--checkpoint-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_finalist_triage.json",
    )
    parser.add_argument("--height-bound", type=int, default=10_000_000)
    parser.add_argument("--search-timeout", type=float, default=120.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--saturation-bound", type=int, default=20)
    parser.add_argument("--saturation-timeout", type=float, default=10.0)
    parser.add_argument("--ellrank-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--rank-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u42_height_10000000.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height_bound != 10_000_000:
        raise SystemExit("this one-shot extension pins --height-bound=10000000")
    if args.search_timeout <= 0 or args.search_timeout > 120:
        raise SystemExit("--search-timeout must be in (0,120]")
    if args.ellrank_timeout <= 0 or args.ellrank_timeout > 60:
        raise SystemExit("--ellrank-timeout must be in (0,60]")
    if min(args.height_timeout, args.saturation_timeout) <= 0:
        raise SystemExit("height and saturation timeouts must be positive")
    if args.saturation_bound < 3:
        raise SystemExit("--saturation-bound must be at least 3")
    if min(args.stack_bytes, args.rank_stack_bytes) < 8_000_000:
        raise SystemExit("PARI stack bounds are too small")

    candidate = load_u42_candidate(args.scan_input)
    checkpoint = checkpoint_points(args.checkpoint_input)
    coefficients = rank13_base_changed_short_jacobian_coefficients(Q(PARAMETER_U))
    extension: dict[str, Any]
    descent_input = checkpoint
    descent_input_provenance = "pinned H=1000000 checkpoint"
    try:
        extension_result, extension_points = exact_candidate_triage(
            candidate,
            height_bound=args.height_bound,
            precisions=args.precisions,
            search_timeout=args.search_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        extension = {"status": "completed", "result": extension_result}
        descent_input = extension_points
        descent_input_provenance = "completed H=10000000 extension"
        print(
            f"u=42 H={args.height_bound} "
            f"rank={extension_result['frontier_stable_pool_numerical_rank']} "
            f"new={extension_result['bounded_search']['new_distinct_jacobian_images']}",
            flush=True,
        )
    except subprocess.TimeoutExpired:
        extension = {
            "status": "timeout",
            "height_bound": args.height_bound,
            "timeout_seconds": args.search_timeout,
            "interpretation": (
                "the sole declared H=10000000 enumeration attempt timed out; "
                "no bounded-search conclusion is drawn"
            ),
        }
        print("u=42 H=10000000 search timed out; using checkpoint basis", flush=True)
    except RuntimeError as error:
        extension = {
            "status": "pari_error",
            "height_bound": args.height_bound,
            "timeout_seconds": args.search_timeout,
            "error": str(error)[:1000],
        }
        print("u=42 H=10000000 search failed; using checkpoint basis", flush=True)

    saturated, saturation = saturate_exact_basis(
        coefficients,
        descent_input,
        prime_bound=args.saturation_bound,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    saturated_height_runs = height_matrix_replay(
        coefficients,
        saturated,
        precisions=args.precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    saturated_rank = stable_height_rank(saturated_height_runs)
    if saturated_rank != len(saturated):
        raise AssertionError("the saturated basis lost numerical full rank")
    descent = ellrank_probe(
        coefficients,
        saturated,
        timeout=args.ellrank_timeout,
        stack_bytes=args.rank_stack_bytes,
    )
    print(
        f"u=42 saturated={len(saturated)} ellrank={descent['status']} "
        f"bounds={descent.get('lower_bound')},{descent.get('upper_bound')}",
        flush=True,
    )

    script_path = Path(__file__).resolve()
    engine_path = script_path.with_name("triage_nagao_rank13_finalists.py")
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "single capped bounded-search extension plus exact saturated-basis "
            "replay and one capped PARI effort-zero descent probe"
        ),
        "primary_source": PRIMARY_SOURCE,
        "candidate": {
            "parameter_u": PARAMETER_U,
            "parameter_t": str(candidate.parameter_t),
            "log_conductor": candidate.log_conductor,
            "root_number": candidate.root_number,
        },
        "inputs": [
            {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (args.scan_input, args.checkpoint_input)
        ],
        "height_10000000_extension": extension,
        "descent_input_provenance": descent_input_provenance,
        "small_prime_saturation": saturation,
        "saturated_height_matrix_runs": list(saturated_height_runs),
        "stable_saturated_basis_numerical_rank": saturated_rank,
        "pari_ellrank_effort_zero": descent,
        "target_status": {
            "rank21_log_conductor_target_certified": False,
            "rank30_target_certified": False,
        },
        "interpretation": (
            "all stored points pass exact curve membership; height ranks are "
            "numerical, saturation uses PARI's documented finite-index "
            "hypothesis, and the ellrank result is a bounded software computation"
        ),
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "parameters": {
            "height_bound": args.height_bound,
            "height_precisions": list(args.precisions),
            "search_timeout_seconds": args.search_timeout,
            "saturation_prime_bound": args.saturation_bound,
            "saturation_timeout_seconds": args.saturation_timeout,
            "ellrank_effort": 0,
            "ellrank_timeout_seconds": args.ellrank_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "pari_rank_stack_bytes": args.rank_stack_bytes,
            "output": str(args.output),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "triage_engine_sha256": hashlib.sha256(engine_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

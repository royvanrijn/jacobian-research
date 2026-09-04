#!/usr/bin/env python3
# <!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE 7f8dffe58168acc8 -->
"""Run one subprocess-capped eclib search on a q12o5867 specialization.

The parent mode is ordinary Python and enforces a wall-clock timeout around a
Sage worker.  The worker uses eclib's unseeded ``pp=0`` discovery mode, checks
every returned affine point exactly, and measures finite-quotient escape from
the separately certified ordered seventeen-point subgroup.  Avoiding eclib
subgroup processing keeps the bounded discovery probe focused on point search.

Timeout, including timeout during eclib curve initialization, is an expected
bounded experimental result and is preserved in the output artifact.  Parent
and worker modes both require a passing residual-Selmer gate bound to the same
parameter and minimal curve.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, lcm
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any, Sequence


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
)
from finite_quotient_escape import QuotientBlock, analyze_escape  # noqa: E402
from elkies_residual_selmer_gate import require_gate_for_specialization  # noqa: E402


Q = Fraction
RESULT_MARKER = "Q12_MWRANK_RESULT_JSON="


def parse_relation_primes(text: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("relation primes must be integers") from error
    if not answer:
        raise argparse.ArgumentTypeError("at least one relation prime is required")
    return answer


def parse_point(record: Sequence[str]) -> tuple[Fraction, Fraction]:
    if len(record) != 2:
        raise ValueError("an affine point record must have two coordinates")
    return Q(record[0]), Q(record[1])


def point_record(point: tuple[Fraction, Fraction]) -> list[str]:
    return [str(point[0]), str(point[1])]


def negate_point(
    model: Sequence[Fraction | int], point: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    a1, _a2, a3, _a4, _a6 = (Q(value) for value in model)
    x_coordinate, y_coordinate = point
    return x_coordinate, -y_coordinate - a1 * x_coordinate - a3


def sign_key(
    model: Sequence[Fraction | int], point: tuple[Fraction, Fraction]
) -> tuple[tuple[int, int], tuple[int, int]]:
    negative = negate_point(model, point)

    def key(value: tuple[Fraction, Fraction]) -> tuple[tuple[int, int], tuple[int, int]]:
        return tuple((coordinate.numerator, coordinate.denominator) for coordinate in value)  # type: ignore[return-value]

    return min(key(point), key(negative))


def to_mwrank_triple(point: tuple[Fraction, Fraction]) -> list[int]:
    x_coordinate, y_coordinate = point
    projective_denominator = lcm(x_coordinate.denominator, y_coordinate.denominator)
    projective_x = int(x_coordinate * projective_denominator)
    projective_y = int(y_coordinate * projective_denominator)
    common = gcd(
        gcd(abs(projective_x), abs(projective_y)), projective_denominator
    )
    return [
        projective_x // common,
        projective_y // common,
        projective_denominator // common,
    ]


def from_mwrank_triple(triple: Sequence[int]) -> tuple[Fraction, Fraction] | None:
    if len(triple) != 3:
        raise ValueError("an mwrank point must have three projective coordinates")
    projective_x, projective_y, projective_z = (int(value) for value in triple)
    if projective_z == 0:
        return None
    return Q(projective_x, projective_z), Q(projective_y, projective_z)


def exact_escape_records(
    artifact: dict[str, Any],
    candidate_minimal_points: Sequence[tuple[Fraction, Fraction]],
    relation_primes: Sequence[int],
    reduction_prime_bound: int,
) -> dict[str, Any]:
    certificate_block = artifact["finite_quotient_independence"]
    certificate_model = tuple(Q(value) for value in certificate_block["certificate_short_model"])
    baseline_points = tuple(parse_point(point) for point in certificate_block["points"])
    if len(baseline_points) != 17:
        raise ValueError("the specialization artifact does not contain 17 baseline points")
    minimal_to_short = WeierstrassChange.from_values(
        certificate_block["minimal_to_certificate_short_change_u_r_s_t"]
    )
    candidate_points = tuple(
        source_point_to_target(point, minimal_to_short)
        for point in candidate_minimal_points
    )
    if any(
        not is_on_weierstrass_curve(certificate_model, point)
        for point in candidate_points
    ):
        raise AssertionError("an mwrank point missed the exact certificate model")
    all_points = (*baseline_points, *candidate_points)
    labels = tuple(f"mwrank-candidate-{index}" for index in range(len(candidate_points)))
    attempts = []
    maximum_marginal_dimension = 0
    escape_basis_labels: list[str] = []
    if candidate_points:
        for relation_prime in relation_primes:
            certificate = build_finite_quotient_certificate(
                certificate_model,
                all_points,
                relation_prime=int(relation_prime),
                prime_bound=reduction_prime_bound,
            )
            verify_finite_quotient_certificate(certificate_model, all_points, certificate)
            blocks = tuple(
                QuotientBlock.build(
                    modulus=int(relation_prime),
                    rows=signature["rows"],
                    column_count=len(all_points),
                    source=f"good-reduction-p={signature['prime']}",
                )
                for signature in certificate["signatures"]
            )
            profile = analyze_escape(
                blocks, known_column_count=17, candidate_labels=labels
            )
            profile_record = profile.to_record()
            profile_record["finite_quotient_certificate"] = certificate
            attempts.append(profile_record)
            if profile.marginal_dimension > maximum_marginal_dimension:
                maximum_marginal_dimension = profile.marginal_dimension
                escape_basis_labels = list(profile.independent_escape_basis_labels)
    return {
        "known_column_count": 17,
        "candidate_count": len(candidate_points),
        "candidate_points_on_certificate_short_model": [
            point_record(point) for point in candidate_points
        ],
        "relation_prime_profiles": attempts,
        "maximum_marginal_dimension": maximum_marginal_dimension,
        "independent_escape_basis_labels": escape_basis_labels,
        "promotion_threshold": 15,
        "promotion_eligible": maximum_marginal_dimension >= 15,
        "promotion_status": (
            "ELIGIBLE_FOR_32_POINT_CERTIFICATION"
            if maximum_marginal_dimension >= 15
            else "NOT_PROMOTED_QUOTIENT_GAIN_BELOW_15"
        ),
        "claim_boundary": (
            "positive finite-quotient escape is exact; bounded non-escape is not "
            "a dependence proof"
        ),
    }


def run_worker(args: argparse.Namespace) -> int:
    from sage.all import ZZ, version as sage_version
    from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

    started = time.monotonic()
    artifact = json.loads(args.input.read_text())
    if artifact.get("status") != "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND":
        raise ValueError("input is not a certified q12o5867 specialization artifact")
    minimal = artifact["global_minimal_specialization"]
    model = tuple(Q(value) for value in minimal["model"])
    if any(value.denominator != 1 for value in model):
        raise ValueError("the input global minimal model is not integral")
    baseline_points = tuple(parse_point(point) for point in minimal["points"])
    if len(baseline_points) != 17:
        raise ValueError("the input does not have 17 transported generic points")
    if any(not is_on_weierstrass_curve(model, point) for point in baseline_points):
        raise AssertionError("a serialized baseline point missed the minimal model")
    require_gate_for_specialization(args.residual_selmer_gate, artifact)

    print("Q12MW|stage=curve_init|status=start", flush=True)
    stage_started = time.monotonic()
    curve = mwrank_EllipticCurve([ZZ(value.numerator) for value in model])
    subgroup = mwrank_MordellWeil(curve, verbose=False, pp=0, maxr=args.max_rank)
    curve_init_seconds = time.monotonic() - stage_started
    print(
        f"Q12MW|stage=curve_init|status=complete|seconds={curve_init_seconds:.6f}",
        flush=True,
    )

    print(f"Q12MW|stage=search|status=start|height={args.height}", flush=True)
    stage_started = time.monotonic()
    subgroup.search(args.height, verbose=False)
    search_seconds = time.monotonic() - stage_started
    final_projective = tuple(
        tuple(int(value) for value in point) for point in subgroup.points()
    )
    print(
        "Q12MW|stage=search|status=complete|"
        f"rank={len(final_projective)}|seconds={search_seconds:.6f}",
        flush=True,
    )

    final_affine = tuple(
        point
        for triple in final_projective
        if (point := from_mwrank_triple(triple)) is not None
    )
    if any(not is_on_weierstrass_curve(model, point) for point in final_affine):
        raise AssertionError("eclib returned a point outside the exact minimal model")
    known_keys = {sign_key(model, point) for point in baseline_points}
    seen = set(known_keys)
    candidates = []
    for point in final_affine:
        key = sign_key(model, point)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(point)
    print(
        f"Q12MW|stage=exact_escape|status=start|candidates={len(candidates)}",
        flush=True,
    )
    escape = exact_escape_records(
        artifact,
        candidates,
        args.relation_primes,
        args.reduction_prime_bound,
    )
    print(
        "Q12MW|stage=exact_escape|status=complete|"
        f"marginal={escape['maximum_marginal_dimension']}",
        flush=True,
    )
    result = {
        "worker_status": "COMPLETED",
        "sage_version": str(sage_version()),
        "height_limit": args.height,
        "max_rank": args.max_rank,
        "curve_init_seconds": curve_init_seconds,
        "search_seconds": search_seconds,
        "worker_wall_seconds": time.monotonic() - started,
        "global_minimal_model": [str(value) for value in model],
        "baseline_input_count": len(baseline_points),
        "baseline_source": "separately exact-certified specialization artifact",
        "eclib_discovery_mode": "unseeded-pp=0",
        "final_mwrank_projective_points": [list(point) for point in final_projective],
        "exact_nonbaseline_minimal_points": [point_record(point) for point in candidates],
        "all_returned_affine_points_exactly_verified": True,
        "finite_quotient_escape": escape,
    }
    print(RESULT_MARKER + json.dumps(result, sort_keys=True), flush=True)
    return 0


def infer_last_stage(stdout: str) -> str | None:
    stages = []
    for line in stdout.splitlines():
        if not line.startswith("Q12MW|stage="):
            continue
        field = line.split("|", 2)[1]
        if field.startswith("stage="):
            stages.append(field.split("=", 1)[1])
    return stages[-1] if stages else None


def run_parent(args: argparse.Namespace) -> int:
    artifact = json.loads(args.input.read_text())
    require_gate_for_specialization(args.residual_selmer_gate, artifact)
    sage = args.sage or shutil.which("sage")
    if sage is None:
        raise FileNotFoundError("Sage executable not found; pass --sage")
    command = [
        str(sage),
        "-python",
        str(Path(__file__).resolve()),
        "--worker",
        "--input",
        str(args.input.resolve()),
        "--residual-selmer-gate",
        str(args.residual_selmer_gate.resolve()),
        "--height",
        str(args.height),
        "--max-rank",
        str(args.max_rank),
        "--relation-primes",
        ",".join(str(value) for value in args.relation_primes),
        "--reduction-prime-bound",
        str(args.reduction_prime_bound),
    ]
    started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=args.timeout,
        )
        stdout = completed.stdout
        stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as error:
        timed_out = True
        stdout = (error.stdout or b"").decode() if isinstance(error.stdout, bytes) else (error.stdout or "")
        stderr = (error.stderr or b"").decode() if isinstance(error.stderr, bytes) else (error.stderr or "")
        returncode = None
    wall_seconds = time.monotonic() - started
    worker_result = None
    for line in stdout.splitlines():
        if line.startswith(RESULT_MARKER):
            worker_result = json.loads(line[len(RESULT_MARKER) :])
    if timed_out:
        status = "TIMEOUT_BOUNDED_MWRANK_PROBE"
    elif returncode != 0:
        status = "ERROR_BOUNDED_MWRANK_PROBE"
    elif worker_result is None:
        status = "ERROR_MISSING_MWRANK_WORKER_RESULT"
    else:
        status = "PASS_BOUNDED_MWRANK_PROBE"
    result = {
        "schema": "elliptic-curves.q12o5867-bounded-mwrank-probe.v1",
        "status": status,
        "input_specialization_artifact": str(args.input.resolve()),
        "residual_selmer_gate": {
            "path": str(args.residual_selmer_gate.resolve()),
            "sha256": sha256(args.residual_selmer_gate.read_bytes()).hexdigest(),
        },
        "bounds": {
            "mwrank_logarithmic_height_limit": args.height,
            "worker_wall_timeout_seconds": args.timeout,
            "max_rank": args.max_rank,
            "reduction_prime_bound": args.reduction_prime_bound,
            "relation_primes": list(args.relation_primes),
        },
        "command": command,
        "wall_seconds": wall_seconds,
        "returncode": returncode,
        "timed_out": timed_out,
        "last_reported_stage": infer_last_stage(stdout),
        "raw_stdout": stdout,
        "raw_stderr": stderr,
        "worker_result": worker_result,
        "promotion_threshold": 15,
        "promoted": bool(
            worker_result
            and worker_result["finite_quotient_escape"]["maximum_marginal_dimension"]
            >= 15
        ),
        "claim_boundary": [
            "This is a bounded point-search experiment, not a rank upper bound.",
            "Timeout is not evidence that no points exist.",
            "No candidate is promoted unless exact quotient gain reaches 15.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={status}")
    print(f"timed_out={timed_out}")
    print(f"last_reported_stage={result['last_reported_stage']}")
    print(f"output={args.output.resolve()}")
    if worker_result is not None:
        escape = worker_result["finite_quotient_escape"]
        print(f"candidate_count={escape['candidate_count']}")
        print(f"maximum_marginal_dimension={escape['maximum_marginal_dimension']}")
        print(f"promoted={result['promoted']}")
    return 0 if status in {"PASS_BOUNDED_MWRANK_PROBE", "TIMEOUT_BOUNDED_MWRANK_PROBE"} else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--residual-selmer-gate", type=Path, required=True)
    parser.add_argument("--height", type=float, default=10.0)
    parser.add_argument("--max-rank", type=int, default=64)
    parser.add_argument("--relation-primes", type=parse_relation_primes, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--sage", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> None:
    sys.set_int_max_str_digits(0)
    parser = build_parser()
    args = parser.parse_args()
    if args.height <= 0 or args.height > 43.668:
        parser.error("--height must lie in (0,43.668]")
    if args.max_rank < 17:
        parser.error("--max-rank must be at least 17")
    if args.reduction_prime_bound < 3:
        parser.error("--reduction-prime-bound must be at least 3")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.worker:
        raise SystemExit(run_worker(args))
    if args.output is None:
        parser.error("parent mode requires --output")
    raise SystemExit(run_parent(args))


if __name__ == "__main__":
    main()

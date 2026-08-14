#!/usr/bin/env python3
"""Strictly capped PARI 2-descent probe for the 2024 rank-29 curve.

PARI documents the interval returned by ``ellrank`` as unconditional.  The
record curve's irreducible cubic field is much larger than routine examples,
so this runner treats timeout, stack failure, and RSS exhaustion as explicit
computational outcomes.  They are not rank bounds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from typing import Any


REPOSITORY = Path(__file__).resolve().parents[2]
CAS = REPOSITORY / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

from elkies_klagsbrun_rank29 import (  # noqa: E402
    COEFFICIENT_A,
    COEFFICIENT_B,
    PUBLISHED_POINTS,
    point_on_general_curve,
)


DEFAULT_OUTPUT = REPOSITORY / "artifacts/generated-results/elliptic_elkies_klagsbrun_rank29_descent_probe.json"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/tools/probe_elkies_klagsbrun_rank29_descent.py"
)


def gp_rational(value: Any) -> str:
    numerator = value.numerator
    denominator = value.denominator
    return str(numerator) if denominator == 1 else f"({numerator}/{denominator})"


def gp_program(stack_bytes: int) -> str:
    points = ",".join(
        f"[{gp_rational(x_value)},{gp_rational(y_value)}]"
        for x_value, y_value in PUBLISHED_POINTS
    )
    return f"""default(parisizemax,{stack_bytes});
E=ellinit([1,0,0,{COEFFICIENT_A},{COEFFICIENT_B}]);
P=[{points}];
print("PARI_VERSION=",version());
print("EXACT_POINTS_ON_CURVE=",vecsum(vector(#P,i,ellisoncurve(E,P[i]))));
gettime();
R=ellrank(E,0,P);
print("ELLRANK_DONE_MS=",gettime());
print("RANK_LOWER=",R[1]);
print("RANK_UPPER=",R[2]);
print("SHA_PAIRING_RANK=",R[3]);
print("RETURNED_POINTS=",#R[4]);
"""


def observed_rss_bytes(pid: int) -> int:
    result = subprocess.run(
        ["/bin/ps", "-p", str(pid), "-o", "rss="],
        check=True,
        capture_output=True,
        text=True,
    )
    value = result.stdout.strip()
    if not value:
        raise RuntimeError(f"ps returned no RSS for owned pid {pid}")
    return int(value) * 1024


def terminate_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is not None:
        return
    process_group = os.getpgid(process.pid)
    if process_group != process.pid:
        raise RuntimeError(
            f"refusing unexpected process group {process_group} for pid {process.pid}"
        )
    os.killpg(process_group, sig)


def integer_marker(output: str, marker: str) -> int | None:
    match = re.search(rf"^{re.escape(marker)}(\d+)$", output, re.MULTILINE)
    return None if match is None else int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=2_000_000_000)
    parser.add_argument("--rss-limit-bytes", type=int, default=2_000_000_000)
    parser.add_argument("--heartbeat", type=float, default=30.0)
    args = parser.parse_args()
    if args.timeout <= 0 or args.timeout > 600:
        raise SystemExit("--timeout must lie in (0,600]")
    if args.heartbeat <= 0 or args.heartbeat > 60:
        raise SystemExit("--heartbeat must lie in (0,60]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        raise SystemExit("stack and RSS limits must each be at least 64MB")
    if len(PUBLISHED_POINTS) != 29 or not all(
        point_on_general_curve(point) for point in PUBLISHED_POINTS
    ):
        raise AssertionError("the exact public point input changed")

    process = subprocess.Popen(
        ["gp", "-fq"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    assert process.stdin is not None
    process.stdin.write(gp_program(args.stack_bytes))
    process.stdin.close()
    process.stdin = None
    started = time.monotonic()
    next_heartbeat = args.heartbeat
    outcome = "running"
    peak_observed_rss = 0

    try:
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= args.timeout:
                outcome = "strict_wall_timeout"
                terminate_owned(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    terminate_owned(process, signal.SIGKILL)
                    process.wait()
                break
            try:
                rss = observed_rss_bytes(process.pid)
            except (subprocess.CalledProcessError, RuntimeError):
                if process.poll() is None:
                    raise
                break
            peak_observed_rss = max(peak_observed_rss, rss)
            if rss > args.rss_limit_bytes:
                outcome = "strict_rss_limit"
                terminate_owned(process, signal.SIGTERM)
                process.wait(timeout=5)
                break
            if elapsed >= next_heartbeat:
                print(
                    f"descent heartbeat seconds={int(elapsed)} rss_bytes={rss}",
                    flush=True,
                )
                next_heartbeat += args.heartbeat
            time.sleep(1)
    except BaseException:
        if process.poll() is None:
            terminate_owned(process, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                terminate_owned(process, signal.SIGKILL)
                process.wait()
        raise

    stdout = "" if process.stdout is None else process.stdout.read()
    stderr = "" if process.stderr is None else process.stderr.read()
    elapsed = time.monotonic() - started
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 and "***" not in stderr else "pari_failure"

    rank_lower = integer_marker(stdout, "RANK_LOWER=")
    rank_upper = integer_marker(stdout, "RANK_UPPER=")
    ellrank_completed = rank_lower is not None and rank_upper is not None
    if outcome != "completed" and ellrank_completed:
        raise AssertionError("PARI emitted a completed rank interval after a failed outcome")
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "artifact_kind": "strictly_bounded_pari_ellrank_probe",
        "status": (
            "completed_unconditional_rank_interval"
            if ellrank_completed
            else "bounded_probe_no_rank_interval"
        ),
        "claim_scope": {
            "exact_input": "all 29 public points checked exactly before PARI launch",
            "software_probe": (
                "PARI ellrank effort zero with supplied points; timeout or failure "
                "is not a rank bound"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": sys.version.split()[0],
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "resource_contract": {
            "wall_timeout_seconds": args.timeout,
            "pari_stack_max_bytes": args.stack_bytes,
            "observed_rss_limit_bytes": args.rss_limit_bytes,
            "peak_observed_rss_bytes": peak_observed_rss,
        },
        "result": {
            "outcome": outcome,
            "elapsed_wall_seconds": elapsed,
            "process_returncode": process.returncode,
            "exact_points_on_curve_marker": integer_marker(
                stdout, "EXACT_POINTS_ON_CURVE="
            ),
            "ellrank_completed": ellrank_completed,
            "rank_lower": rank_lower,
            "rank_upper": rank_upper,
            "sha_pairing_rank": integer_marker(stdout, "SHA_PAIRING_RANK="),
            "returned_points": integer_marker(stdout, "RETURNED_POINTS="),
            "pari_stdout": stdout,
            "pari_stderr": stderr,
            "rank30_target_hit": bool(rank_lower is not None and rank_lower >= 30),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"outcome={outcome} ellrank_completed={str(ellrank_completed).lower()} "
        f"rank30_target_hit={str(artifact['result']['rank30_target_hit']).lower()}"
    )


if __name__ == "__main__":
    main()

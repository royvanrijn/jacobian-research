#!/usr/bin/env python3
"""Run one strictly capped PARI ``ell2cover`` probe for Nagao ``u=135/2``.

The GP process owns a fresh process group.  This runner observes its resident
set size, enforces wall-clock and RSS limits, terminates only that owned group,
and records success or failure as JSON.  It does not retry, detach, install
software, or submit curve data to an external service.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_rank17_frontier_certificate.json"
)
CERTIFICATE_SHA256 = "7378ce59c72974fe39e0e2a40c740f6a96e8dc555a1361b5aaeef67f4d9e0213"
SHORT_A = Q(-564322920496764715904305097281, 80621568)
SHORT_B = Q(
    376721622053639793561558541510664719787672929,
    1880739938304,
)


def gp_rational(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_program(stack_bytes: int) -> str:
    return f"""default(parisizemax,{stack_bytes});
E=ellinit([0,0,0,{gp_rational(SHORT_A)},{gp_rational(SHORT_B)}]);
print("PARI_VERSION=",version());
print("CERTIFICATE_SHA256={CERTIFICATE_SHA256}");
gettime();
C=ell2cover(E);
print("ELL2COVER_DONE_MS=",gettime());
print("COVER_COUNT=",#C);
for(i=1,#C,
  print("COVER_",i,"_BEGIN");
  print(C[i][1]);
  print(C[i][2][1]);
  print(C[i][2][2]);
  print("COVER_",i,"_END")
);
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


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_artifact(
    path: Path,
    *,
    args: argparse.Namespace,
    status: str,
    elapsed: float,
    maximum_rss: int,
    return_code: int,
    stdout: str,
    stderr: str,
) -> None:
    script_path = Path(__file__).resolve()
    cover_count = None
    for line in stdout.splitlines():
        if line.startswith("COVER_COUNT="):
            cover_count = int(line.split("=", 1)[1])
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": status,
        "candidate": {
            "parameter_u": "135/2",
            "short_weierstrass_coefficients": [
                "0",
                "0",
                "0",
                str(SHORT_A),
                str(SHORT_B),
            ],
            "certified_rank_lower_bound": 17,
        },
        "input": {
            "path": str(CERTIFICATE),
            "expected_sha256": CERTIFICATE_SHA256,
            "actual_sha256": sha256_file(CERTIFICATE),
        },
        "declared_budget": {
            "timeout_seconds": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
            "observed_rss_limit_bytes": args.rss_limit_bytes,
            "heartbeat_seconds": args.heartbeat,
            "one_attempt_no_retry": True,
        },
        "result": {
            "elapsed_seconds": elapsed,
            "maximum_observed_rss_bytes": maximum_rss,
            "return_code": return_code,
            "cover_count": cover_count,
            "raw_stdout": stdout,
            "raw_stderr": stderr,
        },
        "prior_strict_probe": {
            "timeout_seconds": 90,
            "status": "strict_timeout_no_cover_output",
            "maximum_observed_rss_bytes": 177143808,
            "pari_version": "2.17.4",
        },
        "interpretation": (
            "ell2cover returns a basis of everywhere locally soluble 2-covers; "
            "a timeout supplies no Selmer dimension, cover, rank bound, or "
            "evidence of rank"
        ),
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "script_sha256": sha256_file(script_path),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=240.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--rss-limit-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--heartbeat", type=float, default=30.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_u135_ell2cover.json"
        ),
    )
    args = parser.parse_args()
    if args.timeout <= 0 or args.timeout > 600:
        raise SystemExit("--timeout must be in (0,600]")
    if args.heartbeat <= 0 or args.heartbeat > 60:
        raise SystemExit("--heartbeat must be in (0,60]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        raise SystemExit("stack and RSS limits must each be at least 64MB")
    if sha256_file(CERTIFICATE) != CERTIFICATE_SHA256:
        raise SystemExit("refusing a changed rank-17 certificate input")

    with tempfile.TemporaryFile(mode="w+") as stdout_file, tempfile.TemporaryFile(
        mode="w+"
    ) as stderr_file:
        process = subprocess.Popen(
            ["gp", "-fq", "-s", str(args.stack_bytes)],
            stdin=subprocess.PIPE,
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
            start_new_session=True,
        )
        print(f"OWNED_GP_PID={process.pid}", file=sys.stderr, flush=True)
        assert process.stdin is not None
        process.stdin.write(gp_program(args.stack_bytes))
        process.stdin.close()
        started = time.monotonic()
        next_heartbeat = args.heartbeat
        maximum_rss = 0
        status = "completed"
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed >= args.timeout:
                    status = "strict_timeout_no_cover_output"
                    terminate_owned(process, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        terminate_owned(process, signal.SIGKILL)
                        process.wait()
                    break
                if elapsed >= next_heartbeat:
                    rss = observed_rss_bytes(process.pid)
                    maximum_rss = max(maximum_rss, rss)
                    print(
                        f"HEARTBEAT_SECONDS={int(elapsed)} RSS_BYTES={rss}",
                        file=sys.stderr,
                        flush=True,
                    )
                    if rss > args.rss_limit_bytes:
                        status = "strict_rss_limit_no_cover_output"
                        terminate_owned(process, signal.SIGTERM)
                        process.wait(timeout=5)
                        break
                    next_heartbeat += args.heartbeat
                time.sleep(0.5)
        except BaseException:
            if process.poll() is None:
                terminate_owned(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    terminate_owned(process, signal.SIGKILL)
                    process.wait()
            raise
        elapsed = time.monotonic() - started
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read()
        stderr = stderr_file.read()
        if process.returncode not in (0, -signal.SIGTERM, -signal.SIGKILL):
            status = "pari_failure"
        if "COVER_COUNT=" not in stdout and status == "completed":
            status = "pari_failure_without_cover_output"
        write_artifact(
            args.output,
            args=args,
            status=status,
            elapsed=elapsed,
            maximum_rss=maximum_rss,
            return_code=int(process.returncode),
            stdout=stdout,
            stderr=stderr,
        )
        print(
            f"STATUS={status} ELAPSED_SECONDS={elapsed:.3f} "
            f"MAXIMUM_RSS_BYTES={maximum_rss} OUTPUT={args.output}",
            flush=True,
        )


if __name__ == "__main__":
    main()

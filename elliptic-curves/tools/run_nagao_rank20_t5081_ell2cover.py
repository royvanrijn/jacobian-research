#!/usr/bin/env python3
"""Run one strictly capped PARI ``ell2cover`` probe on the section-7 curve.

The input is the independently replayable rank-20 certificate for Nagao's
section-7 specialization at constructor parameter ``T=5081/47``.  The GP
process owns a fresh process group.  This supervisor enforces wall-clock and
resident-memory limits, performs no retry, and records success or failure as
JSON.  A timeout is deliberately interpreted as no mathematical evidence.
"""

from __future__ import annotations

import argparse
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


REPOSITORY = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic-curves/elliptic_nagao_rank20_t5081_rank20_certificate.json"
)
CERTIFICATE_SHA256 = (
    "466946076dc0c3fa02d0c5edd90b947d5ee3d10a4fb8cb16567049ab4380f88d"
)
MINIMAL_MODEL = (
    0,
    -1,
    0,
    -47433564031723813622493745045480,
    124574716166660957649866283198474133374724238272,
)
CONDUCTOR = int(
    "4739512365768104141634183882739432010081578062727282562384233196211945306960"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_certificate() -> dict[str, Any]:
    if sha256_file(CERTIFICATE) != CERTIFICATE_SHA256:
        raise SystemExit("refusing a changed rank-20 certificate input")
    data = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    candidate = data["candidate"]
    rank_certificate = data["exact_rank_certificate"]
    if candidate["constructor_parameter_T"] != "5081/47":
        raise SystemExit("the pinned constructor parameter changed")
    if tuple(candidate["minimal_model"]) != MINIMAL_MODEL:
        raise SystemExit("the pinned minimal model changed")
    if int(candidate["conductor"]) != CONDUCTOR:
        raise SystemExit("the pinned conductor changed")
    if candidate["root_number"] != 1:
        raise SystemExit("the pinned root number changed")
    if rank_certificate["certified_algebraic_rank_lower_bound"] != 20:
        raise SystemExit("the pinned rank lower bound changed")
    return data


def gp_program(stack_bytes: int) -> str:
    model = ",".join(str(value) for value in MINIMAL_MODEL)
    return f"""default(parisizemax,{stack_bytes});
E=ellinit([{model}]);
print("PARI_VERSION=",version());
print("CERTIFICATE_SHA256={CERTIFICATE_SHA256}");
print("CONDUCTOR=",ellglobalred(E)[1]);
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


def parse_marker_integer(stdout: str, marker: str) -> int | None:
    prefix = f"{marker}="
    for line in stdout.splitlines():
        if line.startswith(prefix):
            return int(line.split("=", 1)[1])
    return None


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
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": status,
        "candidate": {
            "constructor_parameter_T": "5081/47",
            "minimal_model": list(MINIMAL_MODEL),
            "conductor": str(CONDUCTOR),
            "root_number": 1,
            "certified_rank_lower_bound": 20,
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
            "replayed_conductor": parse_marker_integer(stdout, "CONDUCTOR"),
            "cover_count": parse_marker_integer(stdout, "COVER_COUNT"),
            "raw_stdout": stdout,
            "raw_stderr": stderr,
        },
        "interpretation": (
            "ell2cover returns a basis of everywhere locally soluble 2-covers; "
            "completion alone is not a rank certificate, while a timeout or "
            "memory stop supplies no Selmer dimension, cover, or rank evidence"
        ),
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "script_sha256": sha256_file(script_path),
        "reproducing_command": (
            ".venv/bin/python elliptic-curves/tools/"
            "run_nagao_rank20_t5081_ell2cover.py --timeout 240 "
            "--stack-bytes 1000000000 --rss-limit-bytes 1000000000"
        ),
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
            / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank20_t5081_ell2cover.json"
        ),
    )
    args = parser.parse_args()
    if args.timeout <= 0 or args.timeout > 600:
        raise SystemExit("--timeout must be in (0,600]")
    if args.heartbeat <= 0 or args.heartbeat > 60:
        raise SystemExit("--heartbeat must be in (0,60]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        raise SystemExit("stack and RSS limits must each be at least 64MB")
    validate_certificate()

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
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            terminate_owned(process, signal.SIGKILL)
                            process.wait()
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
        replayed_conductor = parse_marker_integer(stdout, "CONDUCTOR")
        if replayed_conductor is not None and replayed_conductor != CONDUCTOR:
            status = "pari_conductor_mismatch"
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

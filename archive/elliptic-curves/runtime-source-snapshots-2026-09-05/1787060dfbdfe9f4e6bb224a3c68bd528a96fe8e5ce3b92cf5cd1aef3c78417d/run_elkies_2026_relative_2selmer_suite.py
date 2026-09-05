#!/usr/bin/env python3
"""Supervise generated relative 2-Selmer Magma jobs with hard resource limits."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time


INPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-run.v1"
PROTOCOL = "ELKIESR17REL2"


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_process_group(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def run_one(
    command: str,
    program: Path,
    log: Path,
    *,
    timeout: float,
    rss_limit_bytes: int,
) -> dict[str, object]:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as stdout_file:
        process = subprocess.Popen(
            [command, str(program)],
            stdout=stdout_file,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        started = time.monotonic()
        peak_rss = 0
        outcome = "running"
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= timeout:
                outcome = "strict_wall_timeout"
                stop_process_group(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_process_group(process, signal.SIGKILL)
                    process.wait()
                break
            try:
                peak_rss = max(peak_rss, read_rss_bytes(process.pid))
            except (FileNotFoundError, ProcessLookupError):
                pass
            if peak_rss > rss_limit_bytes:
                outcome = "strict_rss_limit"
                stop_process_group(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_process_group(process, signal.SIGKILL)
                    process.wait()
                break
            time.sleep(0.25)
        wall_seconds = time.monotonic() - started
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 else "backend_failure"
    return {
        "outcome": outcome,
        "returncode": process.returncode,
        "wall_seconds": wall_seconds,
        "peak_observed_rss_bytes": peak_rss,
        "timeout_seconds": timeout,
        "rss_limit_bytes": rss_limit_bytes,
        "log": str(log),
        "log_sha256": file_sha256(log),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--magma-command", default="magma")
    parser.add_argument("--timeout-per-case", type=float, default=86400.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=16_000_000_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout_per_case <= 0 or args.rss_limit_bytes < 64_000_000:
        parser.error("timeout must be positive and RSS limit at least 64 MB")
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit("unexpected relative 2-Selmer input manifest")
    command = shutil.which(args.magma_command)
    runs = []
    if command is None:
        for case in manifest["cases"]:
            runs.append(
                {
                    "case_id": case["case_id"],
                    "outcome": "backend_unavailable",
                    "program": case["program"],
                    "program_sha256": case["program_sha256"],
                }
            )
        status = "INCOMPLETE_MAGMA_BACKEND_UNAVAILABLE"
    else:
        for case in manifest["cases"]:
            program = Path(case["program"])
            if file_sha256(program) != case["program_sha256"]:
                raise SystemExit(f"generated program changed: {program}")
            log = args.log_dir / f"{case['case_id']}.log"
            result = run_one(
                command,
                program,
                log,
                timeout=args.timeout_per_case,
                rss_limit_bytes=args.rss_limit_bytes,
            )
            result.update(
                {
                    "case_id": case["case_id"],
                    "program": str(program),
                    "program_sha256": case["program_sha256"],
                }
            )
            runs.append(result)
            print(
                f"{PROTOCOL}|stage=supervisor|case={case['case_id']}"
                f"|outcome={result['outcome']}|seconds={result['wall_seconds']}",
                flush=True,
            )
        status = (
            "COMPLETE_RAW_TRANSCRIPTS_REQUIRE_PARSER"
            if all(run["outcome"] == "completed" for run in runs)
            else "INCOMPLETE_ONE_OR_MORE_MAGMA_JOBS"
        )
    output = {
        "schema": OUTPUT_SCHEMA,
        "status": status,
        "input_manifest": {"path": str(args.manifest), "sha256": file_sha256(args.manifest)},
        "magma_command_requested": args.magma_command,
        "magma_command_resolved": command,
        "runs": runs,
        "claim_boundary": (
            "Only parsed complete protocol transcripts can establish a Selmer result. "
            "Backend absence, timeout, memory stop, or process failure is missing evidence."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"{PROTOCOL}|stage=supervisor_complete|status={status}|output={args.output}")


if __name__ == "__main__":
    main()

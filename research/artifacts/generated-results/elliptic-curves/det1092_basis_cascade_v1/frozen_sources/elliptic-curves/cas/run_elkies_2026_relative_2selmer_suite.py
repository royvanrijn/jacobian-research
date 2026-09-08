#!/usr/bin/env python3
"""Supervise generated relative 2-Selmer Magma jobs with hard resource limits."""

from __future__ import annotations

from research_runtime.supervisor import Limits, capture, capture_record, captured_run, run as supervised_run

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






def run_one(command, program, log, *, timeout, rss_limit_bytes):
    record = supervised_run([command,str(program)], limits=Limits(timeout,rss_limit_bytes),
        log_path=log,checkpoint_path=log.with_suffix('.supervisor.json'))
    return record


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

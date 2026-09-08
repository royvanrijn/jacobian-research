#!/usr/bin/env python3
"""Run and audit the sealed MW17-only Selmer controls on fibres 356 and 385."""

from __future__ import annotations

from research_runtime.supervisor import Limits, capture, capture_record, captured_run, run as supervised_run

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_mw17_only_selmer_control_inputs_v1.json"
)
LOG_DIR = (
    ROOT
    / "artifacts/local/elliptic-curves"
    / "r17-mw17-only-selmer-control-v1/logs"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_mw17_only_selmer_control_run_v1.json"
)
SOURCE = Path(__file__).resolve()
INPUT_SCHEMA = "elliptic-curves.r17-mw17-only-selmer-control-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.r17-mw17-only-selmer-control-run.v1"
PROTOCOL = "R17MW17SELMER"
PREFIX = f"{PROTOCOL}|"
GENERIC_RANK = 17
EXPECTED_CASES = ("mw17-only-control-356", "mw17-only-control-385")
BITS_RE = re.compile(r"^\[\s*([01](?:\s*,\s*[01])*)?\s*\]$")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parse_record(line: str) -> dict[str, str]:
    if not line.startswith(PREFIX):
        raise ValueError("not an MW17-only replay record")
    record: dict[str, str] = {}
    for field in line[len(PREFIX) :].split("|"):
        if "=" not in field:
            raise ValueError("malformed MW17-only replay field")
        key, value = field.split("=", 1)
        if not key or key in record:
            raise ValueError("empty or duplicate MW17-only replay key")
        record[key] = value
    return record


def protocol_records(text: str) -> list[dict[str, str]]:
    return [
        parse_record(line)
        for line in text.splitlines()
        if line.startswith(PREFIX)
    ]


def parse_bits(value: str, dimension: int) -> list[int]:
    match = BITS_RE.fullmatch(value)
    if match is None:
        raise ValueError("malformed binary vector")
    bits = [] if match.group(1) is None else [
        int(item.strip()) for item in match.group(1).split(",")
    ]
    if len(bits) != dimension:
        raise ValueError(
            f"binary vector has dimension {len(bits)}, expected {dimension}"
        )
    return bits


def f2_rank(rows: Iterable[Sequence[int]]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def one(records: list[dict[str, str]], **filters: str) -> dict[str, str]:
    matches = [
        record
        for record in records
        if all(record.get(key) == value for key, value in filters.items())
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one record for {filters}, found {len(matches)}")
    return matches[0]


def parse_complete_transcript(
    text: str, expected: dict[str, Any], required_jump: int
) -> dict[str, Any]:
    records = protocol_records(text)
    allowed_stages = {
        "input",
        "two_selmer",
        "generic_class",
        "quotient_basis",
        "blind_freeze",
    }
    if not records or any(record.get("stage") not in allowed_stages for record in records):
        raise ValueError("transcript contains a missing or forbidden protocol stage")
    input_record = one(records, stage="input")
    if (
        input_record.get("case") != expected["case_id"]
        or input_record.get("role") != "prospective-mw17-only-selmer-control"
        or input_record.get("known_subgroup") != "MW17"
        or input_record.get("fixture_access") != "false"
        or not input_record.get("magma")
        or int(input_record.get("generic_count", -1)) != GENERIC_RANK
    ):
        raise ValueError("transcript input is not the sealed MW17-only case")

    selmer = one(records, stage="two_selmer", status="complete")
    freeze = one(records, stage="blind_freeze", status="complete")
    if records[-1] is not freeze:
        raise ValueError("blind_freeze is not the final protocol record")
    total_dimension = int(selmer["total_dim"])
    generic_rank = int(selmer["generic_kummer_rank"])
    residual_dimension = int(selmer["residual_dim"])
    if (
        generic_rank != GENERIC_RANK
        or residual_dimension != total_dimension - GENERIC_RANK
        or int(freeze["total_dim"]) != total_dimension
        or int(freeze["generic_kummer_rank"]) != generic_rank
        or int(freeze["residual_dim"]) != residual_dimension
    ):
        raise ValueError("inconsistent MW17-relative Selmer dimensions")

    generic_records = [row for row in records if row.get("stage") == "generic_class"]
    quotient_records = [row for row in records if row.get("stage") == "quotient_basis"]
    if len(generic_records) != GENERIC_RANK:
        raise ValueError("transcript does not contain seventeen generic classes")
    generic_rows = [
        parse_bits(row["selmer_bits"], total_dimension) for row in generic_records
    ]
    if f2_rank(generic_rows) != GENERIC_RANK:
        raise ValueError("the displayed generic Kummer rows lost rank")
    if len(quotient_records) != residual_dimension:
        raise ValueError("transcript omitted an MW17 quotient basis row")
    quotient_rows = [
        parse_bits(row["quotient_bits"], residual_dimension)
        for row in quotient_records
    ]
    if f2_rank(quotient_rows) != residual_dimension:
        raise ValueError("the displayed MW17 quotient rows are not a basis")

    return {
        "case_id": expected["case_id"],
        "curve_id": int(expected["curve_id"]),
        "status": "COMPLETE_SEALED_MW17_ONLY_SELMER_REPLAY",
        "magma_version": input_record["magma"],
        "total_two_selmer_dimension": total_dimension,
        "generic_mw17_kummer_dimension": generic_rank,
        "selmer_modulo_mw17_dimension": residual_dimension,
        "required_public_control_jump": required_jump,
        "blind_control_detected_required_jump": residual_dimension >= required_jump,
        "transcript_ends_at_blind_freeze": True,
    }






def run_one(command, program, log, *, timeout, rss_limit_bytes):
    record = supervised_run([command,str(program)], limits=Limits(timeout,rss_limit_bytes),
        log_path=log,checkpoint_path=log.with_suffix('.supervisor.json'))
    return record


def build_result(
    manifest_path: Path,
    log_dir: Path,
    *,
    execute: bool,
    magma_command: str,
    timeout: float,
    rss_limit_bytes: int,
) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text())
    if (
        manifest.get("schema") != INPUT_SCHEMA
        or tuple(row.get("case_id") for row in manifest.get("cases", []))
        != EXPECTED_CASES
    ):
        raise ArithmeticError("the sealed MW17-only replay manifest changed")
    required_jump = int(
        manifest["unblinding_commitment"]["expected_control_jump_over_mw17"]
    )
    command = shutil.which(magma_command) if execute else None
    runs = []
    for case in manifest["cases"]:
        program = ROOT / case["program"]
        if execute and (
            not program.exists() or digest(program) != case["program_sha256"]
        ):
            raise ArithmeticError(f"sealed program is absent or changed: {program}")
        record: dict[str, Any] = {
            "case_id": case["case_id"],
            "curve_id": int(case["curve_id"]),
            "program": case["program"],
            "program_sha256": case["program_sha256"],
            "calibration": None,
        }
        if not execute:
            record["execution"] = {"outcome": "not_executed"}
        elif command is None:
            record["execution"] = {"outcome": "backend_unavailable"}
        else:
            log = log_dir / f"{case['case_id']}.log"
            record["execution"] = run_one(
                command,
                program,
                log,
                timeout=timeout,
                rss_limit_bytes=rss_limit_bytes,
            )
            if record["execution"]["outcome"] == "completed":
                record["calibration"] = parse_complete_transcript(
                    log.read_text(), case, required_jump
                )
        runs.append(record)

    completed = [row for row in runs if row["calibration"] is not None]
    passed = (
        len(completed) == len(runs)
        and all(
            row["calibration"]["blind_control_detected_required_jump"]
            for row in completed
        )
    )
    if passed:
        status = "PASS_BOTH_SEALED_MW17_ONLY_RECORD_REPLAYS"
    elif not execute:
        status = "INCOMPLETE_BLINDED_REPLAY_NOT_EXECUTED"
    elif command is None:
        status = "INCOMPLETE_MAGMA_BACKEND_UNAVAILABLE"
    else:
        status = "INCOMPLETE_OR_FAILED_MW17_ONLY_RECORD_REPLAY"
    return {
        "schema": OUTPUT_SCHEMA,
        "status": status,
        "input_manifest": {
            "path": relative(manifest_path),
            "sha256": digest(manifest_path),
        },
        "execution": {
            "requested": execute,
            "magma_command_requested": magma_command,
            "magma_command_resolved": command,
            "timeout_seconds_per_case": timeout,
            "rss_limit_bytes": rss_limit_bytes,
        },
        "runs": runs,
        "completed_record_replays": len(completed),
        "required_record_replays": len(runs),
        "selmer_candidate_gate_operationally_calibrated": passed,
        "prospective_sample_stage_authorized": False,
        "claim_boundary": [
            "Only both complete, source-hash-matched blind_freeze transcripts can calibrate the Selmer candidate gate.",
            "The generated executables never receive the twelve post-MW17 control points or their half-ideals.",
            "MW29-relative closure calculations cannot substitute for either MW17-only replay.",
            "This calibration alone does not authorize a prospective sample stage or prove a rank bound for any candidate.",
        ],
        "generation": {
            "script": relative(SOURCE),
            "script_sha256": digest(SOURCE),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--log-dir", type=Path, default=LOG_DIR)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--magma-command", default="magma")
    parser.add_argument("--timeout-per-case", type=float, default=86400.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=16_000_000_000)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.check and args.execute:
        parser.error("--check and --execute are mutually exclusive")
    if args.timeout_per_case <= 0 or args.rss_limit_bytes < 64_000_000:
        parser.error("timeout must be positive and RSS limit at least 64 MB")
    document = build_result(
        args.manifest.resolve(),
        args.log_dir.resolve(),
        execute=args.execute,
        magma_command=args.magma_command,
        timeout=args.timeout_per_case,
        rss_limit_bytes=args.rss_limit_bytes,
    )
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != payload:
            raise ArithmeticError("stored MW17-only replay result differs")
    else:
        if output.exists() and not args.overwrite:
            raise FileExistsError(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload)
    print(
        f"{PROTOCOL}|stage=run|completed={document['completed_record_replays']}|"
        f"operational={str(document['selmer_candidate_gate_operationally_calibrated']).lower()}|"
        f"status={document['status']}"
    )


if __name__ == "__main__":
    main()

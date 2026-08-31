#!/usr/bin/env python3
"""Supervise BNF-free two-cover local audits one cover/place at a time.

The underlying exact audit can spend a long time in a singular Hensel tree.
This wrapper gives every ``(cover, p)`` pair its own Sage worker, wall/RSS
limits, and cache block.  Completed local witnesses and obstructions survive
an unrelated timeout, while every timeout remains explicitly inconclusive.

This is a selected finite-place local audit only.  It never turns norm-one
candidates, partial local witnesses, or timeouts into Selmer classes.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
DEFAULT_SAGE = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")
DEFAULT_CACHE = (
    ROOT
    / "artifacts/local/elliptic-curves/bnf-free-two-cover-local-blocks"
)
PROTOCOL = "BNFFREECOVERLOCALSUP"
BLOCK_SCHEMA = "elliptic-curves.bnf-free-two-cover-local-supervised-block.v1"
INPUT_SCHEMA = "elliptic-curves.bnf-free-2cover-equations.v1"
RESULT_PREFIX = f"{PROTOCOL}|result="


WORKER_TEMPLATE = r'''import json, sys
sys.path.insert(0, __CAS__)
from audit_bnf_free_two_cover_reduction import audit_cover, rational

alpha = [rational(value) for value in __ALPHA__]
coefficients = [rational(value) for value in __COEFFICIENTS__]
result = audit_cover(
    alpha,
    coefficients,
    [__PRIME__],
    __MAX_ENUMERATION_PRIME__,
    __MAX_LIFT_PRECISION__,
    __MAX_LIFT_STATES__,
    __RATIONAL_COVER_WITNESS__,
)
print("BNFFREECOVERLOCALSUP|result=" + json.dumps(result, sort_keys=True), flush=True)
'''


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_primes(text: str) -> list[int]:
    values = sorted({int(part.strip()) for part in text.split(",") if part.strip()})
    if not values or any(value < 2 for value in values):
        raise ValueError("--primes must contain integers at least two")
    return values


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def parse_worker(stdout: str) -> dict | None:
    rows = [
        json.loads(line[len(RESULT_PREFIX) :])
        for line in stdout.splitlines()
        if line.startswith(RESULT_PREFIX)
    ]
    return rows[0] if len(rows) == 1 else None


def worker_source(
    *,
    alpha: list[str],
    coefficients: list[str],
    prime: int,
    max_enumeration_prime: int,
    max_lift_precision: int,
    max_lift_states: int,
    rational_cover_witness: list[str] | None,
) -> str:
    return (
        WORKER_TEMPLATE.replace("__CAS__", repr(str(CAS)))
        .replace("__ALPHA__", repr(alpha))
        .replace("__COEFFICIENTS__", repr(coefficients))
        .replace("__PRIME__", str(prime))
        .replace("__MAX_ENUMERATION_PRIME__", str(max_enumeration_prime))
        .replace("__MAX_LIFT_PRECISION__", str(max_lift_precision))
        .replace("__MAX_LIFT_STATES__", str(max_lift_states))
        .replace("__RATIONAL_COVER_WITNESS__", repr(rational_cover_witness))
    )


def expected_metadata(
    *,
    covers_sha256: str,
    cover_index: int,
    cover: dict,
    prime: int,
    args: argparse.Namespace,
) -> dict:
    metadata = {
        "covers_sha256": covers_sha256,
        "cover_index": cover_index,
        "label": str(cover["label"]),
        "alpha_coefficients": [str(value) for value in cover["alpha_coefficients"]],
        "rational_prime": prime,
        "max_enumeration_prime": args.max_enumeration_prime,
        "max_lift_precision": args.max_lift_precision,
        "max_lift_states": args.max_lift_states,
        "timeout_seconds": args.timeout_per_place,
        "rss_limit_bytes": args.rss_limit_bytes,
    }
    if cover.get("rational_cover_witness") is not None:
        metadata["rational_cover_witness"] = [
            str(value) for value in cover["rational_cover_witness"]
        ]
    return metadata


def validate_cached(block: dict, metadata: dict) -> None:
    if block.get("schema") != BLOCK_SCHEMA or block.get("input") != metadata:
        raise ValueError(
            "stale local-audit cache block; use --overwrite-cache after reviewing it"
        )


def run_block(
    *,
    sage_python: str,
    coefficients: list[str],
    metadata: dict,
    args: argparse.Namespace,
) -> dict:
    source = worker_source(
        alpha=metadata["alpha_coefficients"],
        coefficients=coefficients,
        prime=metadata["rational_prime"],
        max_enumeration_prime=args.max_enumeration_prime,
        max_lift_precision=args.max_lift_precision,
        max_lift_states=args.max_lift_states,
        rational_cover_witness=metadata.get("rational_cover_witness"),
    )
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as handle:
        handle.write(source)
        worker_path = Path(handle.name)
    try:
        with tempfile.TemporaryFile(mode="w+") as stdout_file, tempfile.TemporaryFile(
            mode="w+"
        ) as stderr_file:
            process = subprocess.Popen(
                [sage_python, str(worker_path)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
                start_new_session=True,
            )
            started = time.monotonic()
            peak_rss = 0
            outcome = "running"
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed >= args.timeout_per_place:
                    outcome = "strict_wall_timeout"
                    stop_owned(process, signal.SIGTERM)
                try:
                    peak_rss = max(peak_rss, read_rss_bytes(process.pid))
                except (FileNotFoundError, ProcessLookupError):
                    pass
                if peak_rss > args.rss_limit_bytes and process.poll() is None:
                    outcome = "strict_rss_limit"
                    stop_owned(process, signal.SIGTERM)
                if outcome != "running":
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        stop_owned(process, signal.SIGKILL)
                        process.wait()
                    break
                time.sleep(0.05)
            wall_seconds = time.monotonic() - started
            stdout_file.seek(0)
            stderr_file.seek(0)
            stdout = stdout_file.read()
            stderr = stderr_file.read()
    finally:
        worker_path.unlink(missing_ok=True)

    worker = parse_worker(stdout)
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 and worker else "backend_failure"
    place_result = None
    if outcome == "completed" and worker is not None:
        places = worker.get("finite_places")
        if not isinstance(places, list) or len(places) != 1:
            outcome = "backend_failure"
        else:
            place_result = places[0]
    return {
        "schema": BLOCK_SCHEMA,
        "status": "COMPLETE_LOCAL_PLACE" if place_result else "INCONCLUSIVE_LOCAL_PLACE",
        "input": metadata,
        "worker_source_sha256": sha256(source.encode()).hexdigest(),
        "place_result": place_result,
        "supervisor": {
            "outcome": outcome,
            "returncode": process.returncode,
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "stdout_tail": stdout[-8000:],
            "stderr_tail": stderr[-8000:],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covers", type=Path, required=True)
    parser.add_argument("--primes", required=True)
    parser.add_argument("--max-covers", type=int, default=0)
    parser.add_argument("--max-enumeration-prime", type=int, default=251)
    parser.add_argument("--max-lift-precision", type=int, default=8)
    parser.add_argument("--max-lift-states", type=int, default=10000)
    parser.add_argument("--timeout-per-place", type=float, default=10.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--retry-incomplete", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.timeout_per_place <= 0 or args.rss_limit_bytes < 64_000_000:
        parser.error("worker limits must be positive and RSS at least 64 MB")
    if args.max_lift_precision < 1 or args.max_lift_states < 1:
        parser.error("lift limits must be positive")
    if args.max_covers < 0:
        parser.error("--max-covers cannot be negative")
    primes = parse_primes(args.primes)
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python is unavailable: {args.sage_python}")

    record = json.loads(args.covers.read_text())
    if not isinstance(record, dict) or record.get("schema") != INPUT_SCHEMA:
        raise ValueError("expected BNF-free two-cover equations")
    coefficients = [str(value) for value in record["field_polynomial_ascending"]]
    all_covers = record.get("covers")
    if not isinstance(all_covers, list):
        raise ValueError("cover input lacks a cover list")
    selected = all_covers[: args.max_covers] if args.max_covers else all_covers
    covers_sha256 = file_sha256(args.covers)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    blocks = []
    for cover_index, cover in enumerate(selected):
        for prime in primes:
            metadata = expected_metadata(
                covers_sha256=covers_sha256,
                cover_index=cover_index,
                cover=cover,
                prime=prime,
                args=args,
            )
            cache_path = args.cache_dir / f"cover-{cover_index:04d}-p-{prime}.json"
            block = None
            if cache_path.exists() and not args.overwrite_cache:
                block = json.loads(cache_path.read_text())
                validate_cached(block, metadata)
                if args.retry_incomplete and block.get("status") != "COMPLETE_LOCAL_PLACE":
                    block = None
            if block is None:
                block = run_block(
                    sage_python=sage_python,
                    coefficients=coefficients,
                    metadata=metadata,
                    args=args,
                )
                cache_path.write_text(json.dumps(block, indent=2, sort_keys=True) + "\n")
            blocks.append(block)
            classification = (
                block["place_result"]["classification"]
                if block.get("place_result")
                else f"INCONCLUSIVE_{block['supervisor']['outcome'].upper()}"
            )
            print(
                f"{PROTOCOL}|cover={cover_index + 1}/{len(selected)}|p={prime}"
                f"|classification={classification}|cached={cache_path}",
                flush=True,
            )

    by_cover: dict[int, list[dict]] = {index: [] for index in range(len(selected))}
    for block in blocks:
        cover_index = int(block["input"]["cover_index"])
        if block.get("place_result") is not None:
            place = dict(block["place_result"])
            place["supervisor_outcome"] = block["supervisor"]["outcome"]
        else:
            place = {
                "rational_prime": int(block["input"]["rational_prime"]),
                "classification": (
                    "INCONCLUSIVE_SUPERVISED_"
                    + str(block["supervisor"]["outcome"]).upper()
                ),
                "supervisor": block["supervisor"],
            }
        by_cover[cover_index].append(place)

    output_covers = []
    for index, cover in enumerate(selected):
        output_cover = {
                "label": str(cover["label"]),
                "alpha_coefficients": [str(value) for value in cover["alpha_coefficients"]],
                "finite_places": sorted(
                    by_cover[index], key=lambda item: int(item["rational_prime"])
                ),
            }
        if cover.get("rational_cover_witness") is not None:
            output_cover["rational_cover_witness"] = [
                str(value) for value in cover["rational_cover_witness"]
            ]
        output_covers.append(output_cover)
    worker_complete = sum(block.get("place_result") is not None for block in blocks)
    classifications = [
        place["classification"]
        for cover in output_covers
        for place in cover["finite_places"]
    ]
    certified_points = sum(
        str(classification).startswith("PROVED_QP_POINT")
        for classification in classifications
    )
    certified_obstructions = sum(
        str(classification).startswith("PROVED_NO_QP_POINT")
        for classification in classifications
    )
    mathematically_inconclusive = (
        len(classifications) - certified_points - certified_obstructions
    )
    output = {
        "protocol": "BNFFREECOVERLOCAL-v1",
        "status": "SELECTED_FINITE_LOCAL_REDUCTION_AUDIT_ONLY",
        "execution": "OWNED_PER_COVER_PER_PLACE_WORKERS",
        "covers_input": {
            "path": str(args.covers.resolve()),
            "sha256": covers_sha256,
            "total_cover_count": len(all_covers),
            "selected_cover_count": len(selected),
            "selection_truncated": len(selected) != len(all_covers),
        },
        "tested_rational_primes": primes,
        "expected_place_count": len(blocks),
        "completed_worker_place_count": worker_complete,
        "incomplete_worker_place_count": len(blocks) - worker_complete,
        "certified_local_point_count": certified_points,
        "certified_local_obstruction_count": certified_obstructions,
        "mathematically_inconclusive_place_count": mathematically_inconclusive,
        "max_enumeration_prime": args.max_enumeration_prime,
        "max_lift_precision": args.max_lift_precision,
        "max_lift_states": args.max_lift_states,
        "timeout_per_place_seconds": args.timeout_per_place,
        "rss_limit_bytes": args.rss_limit_bytes,
        "covers": output_covers,
        "claim_boundary": [
            "A PROVED_QP_POINT classification is a local witness at that one finite place.",
            "A PROVED_NO_QP_POINT classification is a certified local obstruction.",
            "Timeouts, backend failures, untested places, the real place, and global class-group completeness remain unresolved.",
            "This artifact is not a Selmer upper bound and never authorizes point search.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|covers={len(selected)}|places={len(blocks)}"
        f"|worker_completed={worker_complete}"
        f"|local_points={certified_points}|local_obstructions={certified_obstructions}"
        f"|inconclusive={mathematically_inconclusive}"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()

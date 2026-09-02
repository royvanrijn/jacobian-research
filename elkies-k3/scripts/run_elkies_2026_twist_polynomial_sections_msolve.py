#!/usr/bin/env python3
"""Run and summarize msolve on an exported modular twist-section scheme.

Systems with the same SHA-256 digest are solved once.  For the present short
Weierstrass exports this identifies the two signs of Y, whose X-equations are
identical.  Each subprocess has a hard wall-clock timeout, and checkpoint
records make the batch resumable.

The result is an exhaustive statement about the displayed finite-field
polynomial-section scheme only.  It is not a characteristic-zero rank bound.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
from time import monotonic


ROOT = Path(__file__).resolve().parents[2]
SUMMARY_SCHEMA = "elkies-k3.elkies-2026-twist-polynomial-section-msolve.v1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def parse_metrics(path: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    if not path.exists():
        return result
    for line in path.read_text(errors="replace").splitlines():
        key, separator, value = line.strip().partition(": ")
        if not separator:
            continue
        if key == "Maximum resident set size (kbytes)":
            result["maximum_resident_kbytes"] = int(value)
        elif key == "User time (seconds)":
            result["user_seconds"] = float(value)
        elif key == "System time (seconds)":
            result["system_seconds"] = float(value)
    return result


def classify_solution(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"classification": "missing_output"}
    text = path.read_text(errors="replace").strip()
    if text == "[-1]:":
        return {"classification": "empty_over_algebraic_closure"}
    if text.startswith("[1,"):
        return {"classification": "positive_dimensional"}
    match = re.match(r"\[0,\s*\[(\d+),\s*(\d+),\s*(\d+),", text)
    if match:
        return {
            "classification": "zero_dimensional",
            "field_characteristic": int(match.group(1)),
            "variable_count": int(match.group(2)),
            "quotient_dimension": int(match.group(3)),
        }
    return {"classification": "unparsed_output", "output_prefix": text[:120]}


def process_tree_high_water_kbytes(root_pid: int) -> int | None:
    """Sample VmHWM for a Linux process tree immediately before timeout."""

    processes: dict[int, tuple[int, int]] = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text().split()
            parent_pid = int(stat[3])
            high_water = 0
            for line in (entry / "status").read_text().splitlines():
                if line.startswith("VmHWM:"):
                    high_water = int(line.split()[1])
                    break
            processes[int(entry.name)] = (parent_pid, high_water)
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
            continue
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, (parent_pid, unused_high_water) in processes.items():
            if parent_pid in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    values = [processes[pid][1] for pid in descendants if pid in processes]
    return sum(values) if values else None


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--export", type=Path, required=True)
parser.add_argument("--threads", type=int, default=4, help="threads per msolve process")
parser.add_argument("--jobs", type=int, default=2, help="concurrent msolve processes")
parser.add_argument("--timeout", type=float, default=180.0, help="seconds per distinct system")
parser.add_argument("--max-groups", type=int, help="optional pilot limit")
parser.add_argument(
    "--block",
    type=int,
    action="append",
    help="solve only the distinct-system group containing this block (repeatable)",
)
parser.add_argument("--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve"))
parser.add_argument("--output-dir", type=Path)
parser.add_argument("--summary", type=Path)
args = parser.parse_args()
if args.threads <= 0 or args.jobs <= 0 or args.timeout <= 0:
    parser.error("threads, jobs, and timeout must be positive")
if args.max_groups is not None and args.max_groups <= 0:
    parser.error("--max-groups must be positive")

export_path = args.export.resolve()
export = json.loads(export_path.read_text())
if export.get("schema") != "elkies-k3.elkies-2026-twist-polynomial-section-msolve-export.v1":
    raise ValueError("unexpected modular section export schema")
candidate = export["candidate"]
tag = (
    f"singleton-{candidate['key']}"
    if candidate["kind"] == "singleton"
    else f"product-{candidate['key'].replace(':', '-')}"
    if candidate["kind"] == "product"
    else f"genus-one-{candidate['key']}"
)
prime = int(export["prime"])
if args.output_dir is None:
    args.output_dir = (
        ROOT / "artifacts/local/elkies-k3/twist-polynomial-sections" / tag / f"p{prime}/msolve"
    )
output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
if args.summary is None:
    args.summary = (
        ROOT
        / "artifacts/generated-results"
        / f"elkies-2026-twist-polynomial-sections-{tag}-p{prime}-msolve.json"
    )
summary_path = args.summary.resolve()
summary_path.parent.mkdir(parents=True, exist_ok=True)

groups_by_digest: dict[str, list[dict[str, object]]] = {}
for system in export["systems"]:
    system_path = ROOT / system["path"]
    actual_digest = digest(system_path)
    if actual_digest != system["sha256"]:
        raise ArithmeticError(f"system digest mismatch: {system_path}")
    groups_by_digest.setdefault(actual_digest, []).append(system)
groups = sorted(groups_by_digest.values(), key=lambda group: int(group[0]["block_index"]))
if args.block:
    requested_blocks = set(args.block)
    groups = [
        group
        for group in groups
        if requested_blocks.intersection(int(item["block_index"]) for item in group)
    ]
    covered_requested = {
        int(item["block_index"])
        for group in groups
        for item in group
    }.intersection(requested_blocks)
    if covered_requested != requested_blocks:
        raise ValueError(f"unknown requested blocks: {sorted(requested_blocks - covered_requested)}")
if args.max_groups is not None:
    groups = groups[: args.max_groups]

msolve_path = args.msolve.resolve()
msolve_digest = digest(msolve_path)


def run_group(group: list[dict[str, object]]) -> dict[str, object]:
    representative = group[0]
    block_index = int(representative["block_index"])
    system_path = (ROOT / str(representative["path"])).resolve()
    stem = f"block-{block_index:03d}"
    solution_path = output_dir / f"{stem}.solve"
    log_path = output_dir / f"{stem}.log"
    metrics_path = output_dir / f"{stem}.time"
    checkpoint_path = output_dir / f"{stem}.run.json"
    run_key = {
        "system_sha256": representative["sha256"],
        "msolve_sha256": msolve_digest,
        "threads": args.threads,
        "timeout_seconds": args.timeout,
    }
    if checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text())
        if checkpoint.get("run_key") == run_key and checkpoint.get("status") == "completed":
            return checkpoint

    command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(metrics_path),
        str(msolve_path),
        "-f",
        str(system_path),
        "-o",
        str(solution_path),
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    for path in (solution_path, log_path, metrics_path):
        path.unlink(missing_ok=True)
    started = monotonic()
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    timed_out = False
    timeout_process_tree_high_water_kbytes = None
    try:
        stdout, unused_stderr = process.communicate(timeout=args.timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        timeout_process_tree_high_water_kbytes = process_tree_high_water_kbytes(
            process.pid
        )
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, unused_stderr = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, unused_stderr = process.communicate()
    elapsed = monotonic() - started
    log_path.write_text(stdout)
    classification = (
        {"classification": "timeout"}
        if timed_out
        else classify_solution(solution_path)
    )
    status = "completed" if not timed_out and process.returncode == 0 else "timeout" if timed_out else "error"
    result = {
        "run_key": run_key,
        "status": status,
        "representative_block": block_index,
        "equivalent_blocks": [int(item["block_index"]) for item in group],
        "leading_x_y": [item["leading_x_y"] for item in group],
        "system": relative(system_path),
        "solution": relative(solution_path) if solution_path.exists() else None,
        "log": relative(log_path),
        "metrics": relative(metrics_path) if metrics_path.exists() else None,
        "returncode": process.returncode,
        "elapsed_seconds": elapsed,
        "timeout_process_tree_high_water_kbytes": (
            timeout_process_tree_high_water_kbytes
        ),
        **parse_metrics(metrics_path),
        **classification,
    }
    checkpoint_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


results = []
with ThreadPoolExecutor(max_workers=args.jobs) as executor:
    future_map = {executor.submit(run_group, group): group for group in groups}
    for future in as_completed(future_map):
        result = future.result()
        results.append(result)
        print(
            f"MSOLVE|block={result['representative_block']}|equivalent={result['equivalent_blocks']}"
            f"|classification={result['classification']}|elapsed={result['elapsed_seconds']:.2f}",
            flush=True,
        )

results.sort(key=lambda result: int(result["representative_block"]))
for result in results:
    for label in ("solution", "log", "metrics"):
        stored_path = result.get(label)
        if stored_path:
            materialized_path = ROOT / str(stored_path)
            if materialized_path.exists():
                result[f"{label}_sha256"] = digest(materialized_path)
counts: dict[str, int] = {}
for result in results:
    classification = str(result["classification"])
    counts[classification] = counts.get(classification, 0) + 1
covered_blocks = sum(len(result["equivalent_blocks"]) for result in results)
complete = not args.block and len(results) == len(groups_by_digest) and all(
    result["status"] == "completed" and result["classification"] in {
        "empty_over_algebraic_closure",
        "zero_dimensional",
        "positive_dimensional",
    }
    for result in results
)
payload = {
    "schema": SUMMARY_SCHEMA,
    "status": (
        "PASS_COMPLETE_MODP_POLYNOMIAL_SECTION_SCHEME_SOLVED"
        if complete
        else "INCOMPLETE_BOUNDED_MSOLVE_RUN"
    ),
    "candidate": candidate,
    "prime": prime,
    "export": relative(export_path),
    "export_sha256": digest(export_path),
    "total_blocks": len(export["systems"]),
    "distinct_systems": len(groups_by_digest),
    "attempted_distinct_systems": len(results),
    "covered_blocks": covered_blocks,
    "classification_counts_by_distinct_system": counts,
    "run_parameters": {
        "jobs": args.jobs,
        "threads_per_job": args.threads,
        "timeout_seconds_per_distinct_system": args.timeout,
        "requested_blocks": args.block,
        "msolve": str(msolve_path),
        "msolve_sha256": msolve_digest,
    },
    "results": results,
    "proof_boundary": (
        "A complete result describes only polynomial P.O=0 sections within the chi degree "
        "bounds over the displayed finite field and chart. It is not a characteristic-zero "
        "Mordell-Weil rank upper bound; sections can have bad reduction at this prime."
    ),
}
temporary_summary = summary_path.with_suffix(summary_path.suffix + ".tmp")
temporary_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
temporary_summary.replace(summary_path)
print(
    f"ELKIES2026TWISTMSOLVE|kind={candidate['kind']}|key={candidate['key']}|p={prime}"
    f"|distinct={len(results)}/{len(groups_by_digest)}|blocks={covered_blocks}/{len(export['systems'])}"
    f"|status={payload['status']}|output={summary_path}",
    flush=True,
)

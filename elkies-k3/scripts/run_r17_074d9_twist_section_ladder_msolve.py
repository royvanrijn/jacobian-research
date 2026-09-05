#!/usr/bin/env python3
"""Run checkpointed msolve leading-ideal tests for section-ladder exports."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
from research_runtime.supervisor import Limits, capture, WorkerFailure, preserve_previous
from research_runtime.section_gate import guard_export


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("exports", nargs="+", type=Path)
parser.add_argument("--workers", type=int, default=8)
parser.add_argument("--threads-per-job", type=int, default=1)
parser.add_argument("--timeout", type=int, default=120)
parser.add_argument("--rss-limit-bytes", type=int, default=2_000_000_000)
parser.add_argument("--reduction-only", action="store_true")
parser.add_argument("--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve"))
args = parser.parse_args()

if args.workers < 1 or args.threads_per_job < 1 or args.timeout < 1:
    raise ValueError("worker, thread, and timeout parameters must be positive")
msolve = args.msolve.resolve()
if not msolve.is_file():
    raise FileNotFoundError(msolve)

tasks = []
export_inputs = []
for export_argument in args.exports:
    export_path = export_argument.resolve()
    export = json.loads(export_path.read_text())
    if export.get("schema") != "elkies-k3.r17-074d9-twist-section-ladder-msolve-export.v1":
        raise ValueError(f"{export_path}: unexpected export schema")
    guard_export(export, ROOT, reduction_only=args.reduction_only,
                 limits={"wall_seconds": args.timeout, "rss_bytes": args.rss_limit_bytes})
    export_inputs.append({"path": str(export_path.relative_to(ROOT)), "sha256": digest(export_path)})
    for system in export["systems"]:
        input_path = ROOT / system["path"]
        if digest(input_path) != system["sha256"]:
            raise ArithmeticError(f"{input_path}: stale system hash")
        output_path = input_path.with_suffix(".leading-ideal")
        tasks.append(
            {
                "label": export["label"],
                "prime": int(export["prime"]),
                "P_dot_O": int(export["intersection_P_dot_O"]),
                "chart_index": int(system["chart_index"]),
                "block_index": int(system["block_index"]),
                "leading_x_y": system["leading_x_y"],
                "input_path": input_path,
                "input_sha256": system["sha256"],
                "output_path": output_path,
            }
        )


def run(task):
    command = [
        str(msolve),
        "-t",
        str(args.threads_per_job),
        "-g",
        "1",
        "-f",
        str(task["input_path"]),
        "-o",
        str(task["output_path"]),
    ]
    try:
        preserve_previous(task["output_path"])
        completed = capture(command, limits=Limits(args.timeout, args.rss_limit_bytes),
                            log_path=task["output_path"].with_suffix(".log"), check=False)
    except subprocess.TimeoutExpired as error:
        task["outcome"] = "TIMEOUT"
        task["driver_output"] = (error.stdout or "").decode() if isinstance(error.stdout, bytes) else (error.stdout or "")
        return task
    except WorkerFailure as error:
        task["outcome"] = "ERROR"
        task["supervision"] = error.record
        return task
    task["returncode"] = completed.returncode
    task["driver_output"] = completed.stdout
    if completed.returncode != 0 or not task["output_path"].is_file():
        task["outcome"] = "ERROR"
        return task
    output_text = task["output_path"].read_text()
    task["output_sha256"] = digest(task["output_path"])
    if re.search(r"(?:^|\n)\[\]:\s*$", output_text):
        task["outcome"] = "UNIT_IDEAL"
    else:
        task["outcome"] = "NONUNIT_LEADING_IDEAL"
    return task


rows = []
with ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = [executor.submit(run, task) for task in tasks]
    for completed_count, future in enumerate(as_completed(futures), start=1):
        row = future.result()
        rows.append(row)
        print(
            f"R17074D9MSOLVE|done={completed_count}/{len(tasks)}|label={row['label']}"
            f"|P.O={row['P_dot_O']}|chart={row['chart_index']}"
            f"|block={row['block_index']}|outcome={row['outcome']}",
            flush=True,
        )

rows.sort(key=lambda row: (row["label"], row["P_dot_O"], row["chart_index"], row["block_index"]))
for row in rows:
    row["input_path"] = str(row["input_path"].relative_to(ROOT))
    row["output_path"] = str(row["output_path"].relative_to(ROOT))

labels = sorted({row["label"] for row in rows})
po_values = sorted({row["P_dot_O"] for row in rows})
tag = "--".join(label.removeprefix("074d9-orbit-") for label in labels)
po_tag = "-".join(map(str, po_values))
summary_path = (
    ROOT
    / "artifacts/local/elkies-k3/r17-074d9-twist-section-ladder"
    / f"msolve-{tag}-po-{po_tag}.json"
)
summary_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-074d9-twist-section-ladder-msolve-run.v1",
    "status": (
        "PASS_ALL_LEADING_IDEALS_CLASSIFIED"
        if all(row["outcome"] in {"UNIT_IDEAL", "NONUNIT_LEADING_IDEAL"} for row in rows)
        else "INCOMPLETE_BOUNDED_MSOLVE_RUN"
    ),
    "proof_boundary": (
        "UNIT_IDEAL excludes geometric mod-p solutions in that block. A nonunit "
        "leading ideal only retains a candidate block; rational-point extraction "
        "and characteristic-zero lifting are separate."
    ),
    "msolve": {"path": str(msolve), "sha256": digest(msolve)},
    "workers": args.workers,
    "threads_per_job": args.threads_per_job,
    "timeout_seconds_per_block": args.timeout,
    "exports": export_inputs,
    "counts": {
        outcome: sum(row["outcome"] == outcome for row in rows)
        for outcome in ("UNIT_IDEAL", "NONUNIT_LEADING_IDEAL", "TIMEOUT", "ERROR")
    },
    "blocks": rows,
}
summary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17074D9MSOLVE|blocks={len(rows)}|counts={payload['counts']}"
    f"|summary={summary_path}|status={payload['status']}",
    flush=True,
)

#!/usr/bin/env python3
"""Run bounded exact PARI rank intervals on shortlisted V4 base Jacobians.

status: ACTIVE_SEARCH
claim: bounded exact rank-interval screen for declared genus-one base Jacobians
inputs: exact 64-pair alternate-Q80 V4 shortlist
outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json
supersedes: none; refines pair priority without asserting unresolved ranks

Each curve is isolated in a Sage subprocess with a hard timeout.  Completed
``ellrank`` calls return exact lower and upper bounds.  A timeout leaves the
rank ``UNKNOWN`` and is never interpreted as a rank bound.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
SHORTLIST = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
)
CHECKPOINTS = ROOT / "artifacts/local/elkies-k3/r17-norm12-11952-v4-base-ranks"
SCHEMA = "elkies-k3.r17-norm12-11952-v4-base-rank-screen.v1"
WORKER_PREFIX = "ALTV4RANKWORKER|"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def worker(shortlist_path: Path, index: int, effort: int) -> None:
    from sage.all import EllipticCurve, QQ, pari  # type: ignore

    payload = json.loads(shortlist_path.read_text())
    pair = payload["pairs"][index]
    curve = EllipticCurve(
        QQ,
        [QQ(value) for value in pair["base_jacobian_integral_a1_a2_a3_a4_a6"]],
    )
    answer = pari(curve).ellrank(effort)
    result = {
        "pair_key": pair["pair_key"],
        "shortlist_rank": int(pair["shortlist_rank"]),
        "rank_lower_bound": int(answer[0]),
        "rank_upper_bound": int(answer[1]),
        "sha_information": int(answer[2]),
        "independent_points_found": [str(point) for point in answer[3]],
        "pari_effort": effort,
    }
    print(WORKER_PREFIX + json.dumps(result, sort_keys=True), flush=True)


def verify(path: Path) -> None:
    payload = json.loads(path.read_text())
    if payload.get("schema") != SCHEMA:
        raise ValueError("unexpected V4 base-rank screen schema")
    for name, expected in payload["inputs"].items():
        if digest(ROOT / name) != expected:
            raise ArithmeticError(f"rank-screen input digest changed: {name}")
    print(
        "ALTV4RANKCHECK|"
        f"completed={payload['summary']['completed']}|"
        f"timeouts={payload['summary']['timeouts']}|output={display_path(path)}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, default=SHORTLIST)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--limit", type=int, default=64)
    parser.add_argument("--timeout", type=float, default=15.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--effort", type=int, default=0)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--worker-index", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_index is not None:
        worker(args.shortlist, args.worker_index, args.effort)
        return
    if args.check:
        verify(args.output)
        return
    if args.limit < 1 or args.timeout <= 0 or args.jobs < 1 or args.effort < 0:
        parser.error("limit, timeout, jobs, and effort must be nonnegative/positive")

    shortlist = json.loads(args.shortlist.read_text())
    if shortlist.get("status") != "PASS_EXACT_BOUNDED_RATIONAL_V4_PAIR_SHORTLIST":
        raise ValueError("input is not the exact rational V4 shortlist")
    pairs = shortlist["pairs"][: args.limit]
    if len(pairs) != args.limit:
        raise ValueError("--limit exceeds the exact shortlist")
    sage = shutil.which("sage")
    if sage is None:
        raise FileNotFoundError("the Sage launcher is required")
    sage_version = subprocess.run(
        [sage, "--version"], text=True, capture_output=True, check=True
    ).stdout.strip()
    pari_version = subprocess.run(
        [sage, "-python", "-c", "from sage.all import pari; print(pari.version())"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    def run(index: int) -> dict[str, object]:
        pair = pairs[index]
        checkpoint_path = CHECKPOINTS / f"pair-{int(pair['shortlist_rank']):03d}.json"
        run_key = {
            "script_sha256": digest(Path(__file__).resolve()),
            "shortlist_sha256": digest(args.shortlist),
            "pair_key": pair["pair_key"],
            "timeout_seconds": args.timeout,
            "pari_effort": args.effort,
        }
        if checkpoint_path.exists():
            old = json.loads(checkpoint_path.read_text())
            if old.get("run_key") == run_key:
                return old
        command = [
            sage,
            "-python",
            str(Path(__file__).resolve()),
            "--shortlist",
            str(args.shortlist),
            "--worker-index",
            str(index),
            "--effort",
            str(args.effort),
        ]
        try:
            completed = subprocess.run(
                command,
                cwd=ROOT,
                text=True,
                capture_output=True,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            result = {
                "run_key": run_key,
                "status": "timeout",
                "pair_key": pair["pair_key"],
                "shortlist_rank": int(pair["shortlist_rank"]),
                "rank_status": "UNKNOWN",
            }
        else:
            marker = next(
                (
                    line[len(WORKER_PREFIX) :]
                    for line in completed.stdout.splitlines()
                    if line.startswith(WORKER_PREFIX)
                ),
                None,
            )
            if completed.returncode or marker is None:
                result = {
                    "run_key": run_key,
                    "status": "error",
                    "pair_key": pair["pair_key"],
                    "shortlist_rank": int(pair["shortlist_rank"]),
                    "rank_status": "UNKNOWN",
                    "returncode": completed.returncode,
                    "stderr_tail": completed.stderr[-1000:],
                }
            else:
                rank = json.loads(marker)
                result = {
                    "run_key": run_key,
                    "status": "completed",
                    "rank_status": "EXACT_INTERVAL",
                    **rank,
                }
        checkpoint_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return result

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run, index): index for index in range(len(pairs))}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                "ALTV4RANK|"
                f"pair={result['pair_key']}|status={result['status']}|"
                f"interval={result.get('rank_lower_bound', '?')}:{result.get('rank_upper_bound', '?')}",
                flush=True,
            )
    results.sort(key=lambda row: int(row["shortlist_rank"]))
    completed = [row for row in results if row["status"] == "completed"]
    ranked = sorted(
        completed,
        key=lambda row: (
            int(row["rank_upper_bound"]),
            int(row["rank_lower_bound"]),
            int(row["shortlist_rank"]),
        ),
    )
    result = {
        "schema": SCHEMA,
        "status": "PASS_BOUNDED_EXACT_BASE_JACOBIAN_RANK_INTERVAL_SCREEN",
        "inputs": {
            display_path(path): digest(path)
            for path in (Path(__file__).resolve(), args.shortlist)
        },
        "limits": {
            "shortlist_prefix": args.limit,
            "timeout_seconds_per_pair": args.timeout,
            "concurrent_workers": args.jobs,
            "pari_effort": args.effort,
        },
        "summary": {
            "completed": len(completed),
            "timeouts": sum(row["status"] == "timeout" for row in results),
            "errors": sum(row["status"] == "error" for row in results),
            "exact_rank_count": sum(
                row["rank_lower_bound"] == row["rank_upper_bound"] for row in completed
            ),
            "minimum_completed_upper_bound": (
                min(int(row["rank_upper_bound"]) for row in completed)
                if completed
                else None
            ),
        },
        "results": results,
        "completed_pairs_ranked_by_upper_then_lower_bound": [
            row["pair_key"] for row in ranked
        ],
        "software_assumptions": {
            "sage": sage_version,
            "pari_via_sage": pari_version,
        },
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Completed PARI ellrank calls give exact displayed rank intervals. A timeout "
            "or error leaves the rank UNKNOWN. The screen does not alter the exact "
            "64-base shortlist and says nothing about product-twist sections."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ALTV4RANKSCREEN|"
        f"completed={len(completed)}/{len(results)}|"
        f"exact={result['summary']['exact_rank_count']}|"
        f"status={result['status']}|output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

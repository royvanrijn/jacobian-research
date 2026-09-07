#!/usr/bin/env python3
"""Run independent exact alternate-fibration searches over a frame beam."""

from __future__ import annotations

import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


BASE = Path(__file__).resolve().parents[2]
SEARCH = BASE / "elkies-k3" / "scripts" / "search_alternate_fibrations.sage"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--qmin", type=int, default=4)
    ap.add_argument("--qmax", type=int, default=16)
    ap.add_argument("--enum-baseline-cap", type=int, default=200)
    ap.add_argument("--enum-restarts", type=int, default=4)
    ap.add_argument("--enum-cap", type=int, default=200)
    ap.add_argument("--enum-seed", type=int, default=20260820)
    ap.add_argument("--per-root-data-cap", type=int, default=2)
    ap.add_argument("--report", type=int, default=80)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rank-order", choices=("low", "high"), default="low")
    return ap.parse_args()


def read_nodes(path: Path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(not row.get("node_id") or not row.get("frame") for row in rows):
        raise SystemExit("input TSV needs nonempty node_id and frame columns")
    if len({row["node_id"] for row in rows}) != len(rows):
        raise SystemExit("input TSV node_id values must be unique")
    for row in rows:
        frame = (BASE / row["frame"]).resolve() if not Path(row["frame"]).is_absolute() else Path(row["frame"]).resolve()
        if not frame.is_file():
            raise SystemExit(f"missing frame: {frame}")
        row["frame_path"] = frame
    return rows


def run_node(args, row, out_dir: Path):
    node_id = row["node_id"]
    node_dir = out_dir / node_id
    node_dir.mkdir(parents=True, exist_ok=False)
    frames_dir = node_dir / "frames"
    result_path = node_dir / "search.txt"
    log_path = node_dir / "search.log"
    command = [
        "sage",
        str(SEARCH),
        "--frame",
        str(row["frame_path"]),
        "--min-qnorm",
        str(args.qmin),
        "--max-qnorm",
        str(args.qmax),
        "--enum-baseline-cap",
        str(args.enum_baseline_cap),
        "--enum-restarts",
        str(args.enum_restarts),
        "--enum-cap",
        str(args.enum_cap),
        "--enum-seed",
        str(args.enum_seed),
        "--proper-factors-only",
        "--one-factor-order",
        "--per-root-data-cap",
        str(args.per_root_data_cap),
        "--quiet-candidates",
        "--root-method",
        "pari",
        "--rank-order",
        args.rank_order,
        "--report",
        str(args.report),
        "--out",
        str(result_path),
        "--frames-dir",
        str(frames_dir),
    ]
    with log_path.open("w") as log:
        proc = subprocess.run(command, cwd=BASE, stdout=log, stderr=subprocess.STDOUT)
    return node_id, proc.returncode, frames_dir / "hits.tsv", command


def main():
    args = parse_args()
    nodes = read_nodes(args.input.resolve())
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=False)

    completed = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_node, args, row, out_dir): row for row in nodes}
        for future in as_completed(futures):
            node_id, returncode, hits_path, command = future.result()
            print(f"BEAM|node={node_id}|returncode={returncode}", flush=True)
            completed.append((node_id, returncode, hits_path, command))

    failed = [node_id for node_id, code, _, _ in completed if code]
    if failed:
        raise SystemExit("failed beam nodes: " + ",".join(sorted(failed)))

    rows = []
    for node_id, _, hits_path, command in completed:
        with hits_path.open() as handle:
            for hit in csv.DictReader(handle, delimiter="\t"):
                rows.append({
                    "parent_node": node_id,
                    "MW": hit["MW"],
                    "root_rank": hit["root_rank"],
                    "roots": hit["roots"],
                    "rootdet": hit["rootdet"],
                    "q": hit["q"],
                    "a": hit["a"],
                    "b": hit["b"],
                    "v": hit["v"],
                    "frame": str((hits_path.parent / hit["frame"]).relative_to(BASE)),
                    "parent_frame": hit["parent_frame"],
                })
    if args.rank_order == "low":
        rows.sort(key=lambda row: (
            int(row["MW"]), -int(row["roots"]), int(row["rootdet"]),
            int(row["q"]), row["parent_node"],
        ))
    else:
        rows.sort(key=lambda row: (
            -int(row["MW"]), int(row["roots"]), int(row["rootdet"]),
            int(row["q"]), row["parent_node"],
        ))
    fields = [
        "parent_node", "MW", "root_rank", "roots", "rootdet", "q",
        "a", "b", "v", "frame", "parent_frame",
    ]
    with (out_dir / "beam-hits.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (out_dir / "commands.txt").open("w") as handle:
        for node_id, _, _, command in sorted(completed):
            handle.write(node_id + "\t" + " ".join(command) + "\n")

    best = rows[0]
    print(
        "BEAM|status=done"
        f"|nodes={len(nodes)}|retained={len(rows)}|best_MW={best['MW']}"
        f"|best_parent={best['parent_node']}|best_frame={best['frame']}",
        flush=True,
    )


if __name__ == "__main__":
    main()

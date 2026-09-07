#!/usr/bin/env python3
from pathlib import Path
import argparse
import concurrent.futures
import json
import subprocess


ap = argparse.ArgumentParser(description="Run complete canonical P3 audits on target-P2 surfaces.")
ap.add_argument("--input", required=True, help="JSONL from search_mw3_a10_p1p2_multihit.sage")
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--dir", default="artifacts/local/elkies-k3/mw3-a10-p1/p3-targets")
args = ap.parse_args()

records = [json.loads(line) for line in Path(args.input).read_text().splitlines() if line.strip()]
surfaces = {}
for record in records:
    surfaces.setdefault((record["seed"], record["hit"]), record)
root = Path(args.dir)
root.mkdir(parents=True, exist_ok=True)


def run_one(record):
    seed = record["seed"]
    hit = record["hit"]
    command = [
        "sage",
        "elkies-k3/scripts/search_mw3_a10_p3_sliced.sage",
        "--p", str(record["p"]),
        "--A", ",".join(map(str, record["A"])),
        "--B", ",".join(map(str, record["B"])),
        "--lam", str(record["point"]["lam"]),
        "--nodes", ",".join(map(str, record["nodes"])),
        "--sinf", str(record["sinf"]),
        "--max-hits", "100",
    ]
    process = subprocess.run(command, text=True, capture_output=True)
    log = root / f"seed{seed}-hit{hit}.p3.log"
    log.write_text(process.stdout + process.stderr)
    if process.returncode:
        return seed, hit, "error", process.stderr[-1000:]
    hits = [line for line in process.stdout.splitlines() if line.startswith("MW3A10P3_HIT|")]
    return seed, hit, "P3-HIT" if hits else "no-p3", "\n".join(hits)


counts = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = [executor.submit(run_one, record) for record in surfaces.values()]
    for future in concurrent.futures.as_completed(futures):
        seed, hit, status, detail = future.result()
        counts[status] = counts.get(status, 0) + 1
        print(f"MW3A10P3BATCH|seed={seed}|hit={hit}|status={status}", flush=True)
        if detail:
            print(detail, flush=True)

print(
    f"MW3A10P3BATCH|done=1|surfaces={len(surfaces)}|p2_records={len(records)}"
    f"|counts={json.dumps(counts, sort_keys=True)}",
    flush=True,
)

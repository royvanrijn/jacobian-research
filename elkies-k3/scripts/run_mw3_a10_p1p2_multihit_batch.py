#!/usr/bin/env python3
from pathlib import Path
import argparse
import concurrent.futures
import json
import subprocess


ap = argparse.ArgumentParser(description="Partition the persistent A10/P1+P2 reconstruction gate.")
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--dir", required=True)
ap.add_argument("--formula-meta", required=True)
ap.add_argument("--target-jsonl", required=True)
args = ap.parse_args()

count = args.seed_end - args.seed_start + 1
if count <= 0:
    raise SystemExit("empty seed range")
workers = min(args.workers, count)
base, extra = divmod(count, workers)
ranges = []
lower = args.seed_start
for worker in range(workers):
    size = base + (1 if worker < extra else 0)
    upper = lower + size - 1
    ranges.append((worker, lower, upper))
    lower = upper + 1

target_path = Path(args.target_jsonl)
target_path.parent.mkdir(parents=True, exist_ok=True)


def run_range(job):
    worker, lower, upper = job
    part = target_path.with_name(target_path.name + f".part{worker}")
    command = [
        "sage", "elkies-k3/scripts/search_mw3_a10_p1p2_multihit.sage",
        "--seed-start", str(lower), "--seed-end", str(upper),
        "--dir", args.dir, "--formula-meta", args.formula_meta,
        "--target-jsonl", str(part),
    ]
    process = subprocess.run(command, text=True, capture_output=True)
    final = next(
        (line for line in reversed(process.stdout.splitlines()) if "|done=1|" in line),
        "",
    )
    return worker, lower, upper, part, process.returncode, final, process.stderr[-1000:]


parts = []
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
    futures = [executor.submit(run_range, job) for job in ranges]
    for future in concurrent.futures.as_completed(futures):
        worker, lower, upper, part, returncode, final, error = future.result()
        status = "ok" if returncode == 0 else "error"
        print(
            f"MW3A10MULTIP2BATCH|start={lower}|end={upper}|status={status}|summary={final}",
            flush=True,
        )
        if error:
            print(error, flush=True)
        if returncode:
            raise SystemExit(f"P2 worker {worker} failed")
        parts.append(part)

records = []
for part in parts:
    records.extend(
        json.loads(line) for line in part.read_text().splitlines() if line.strip()
    )
records.sort(key=lambda record: (record["seed"], record["hit"], record["P2"]["r"]))
target_path.write_text("".join(json.dumps(record, sort_keys=True) + "\n" for record in records))
print(
    f"MW3A10MULTIP2BATCH|done=1|seeds={count}|workers={workers}|targets={len(records)}",
    flush=True,
)

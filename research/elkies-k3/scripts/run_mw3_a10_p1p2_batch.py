#!/usr/bin/env python3
from pathlib import Path
import argparse
import ast
import concurrent.futures
import json
import re
import subprocess


ap = argparse.ArgumentParser(
    description="Reconstruct scanned A10/P1 hits and run the fixed-surface P2 class-6 test."
)
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument("--workers", type=int, default=4)
ap.add_argument(
    "--dir", default="artifacts/local/elkies-k3/mw3-a10-p1"
)
ap.add_argument(
    "--formula-meta",
    default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.meta.txt",
)
args = ap.parse_args()

root = Path(args.dir)
hit_re = re.compile(r"MW3A10SCAN_HIT\|(.+)")
jobs = []
for seed in range(args.seed_start, args.seed_end + 1):
    meta = root / f"component2-valid-seed{seed}.meta.txt"
    scan = root / f"component2-valid-seed{seed}.scan.log"
    if not meta.exists() or not scan.exists():
        continue
    values_line = next(
        line.split("=", 1)[1]
        for line in meta.read_text().splitlines()
        if line.startswith("values=")
    )
    fixed = ast.literal_eval(values_line)
    hit_index = 0
    for line in scan.read_text().splitlines():
        match = hit_re.fullmatch(line)
        if not match:
            continue
        hit_index += 1
        hit = {}
        for item in re.split(r"[|,]", match.group(1)):
            name, value = item.split("=", 1)
            hit[name] = int(value)
        point = {**fixed, **hit}
        jobs.append((seed, hit_index, point))


def run_one(job):
    seed, hit_index, point = job
    point_arg = ",".join(
        f"{name}={point[name]}"
        for name in ("rho", "r1", "s1", "lam", "x2", "x3", "y4")
    )
    reconstruction = subprocess.run(
        [
            "sage", "elkies-k3/scripts/reconstruct_mw3_a10_p1_hit.sage",
            "--meta", args.formula_meta, "--point", point_arg,
        ],
        text=True,
        capture_output=True,
    )
    if reconstruction.returncode:
        return seed, hit_index, "reconstruct-error", reconstruction.stderr[-1000:]
    line = next(
        line for line in reconstruction.stdout.splitlines() if line.startswith("MW3A10FAST|")
    )
    record = json.loads(line.split("|", 1)[1])
    if not record["valid_semistable"]:
        return seed, hit_index, "boundary", json.dumps(record, sort_keys=True)

    out = root / f"p2c6-seed{seed}-hit{hit_index}.ms"
    command = [
        "sage", "elkies-k3/scripts/build_mw3_a10_p2_component6.sage",
        "--out", str(out),
        "--A", ",".join(map(str, record["A"])),
        "--B", ",".join(map(str, record["B"])),
        "--lam", str(record["point"]["lam"]),
        "--nodes", ",".join(map(str, record["nodes"])),
        "--sinf", str(record["sinf"]),
        "--X1", ",".join(map(str, record["X1"])),
        "--Y1", ",".join(map(str, record["Y1"])),
        "--solve",
    ]
    p2 = subprocess.run(command, text=True, capture_output=True)
    if p2.returncode:
        return seed, hit_index, "p2-error", p2.stderr[-1000:]
    interesting = [
        line for line in p2.stdout.splitlines()
        if line.startswith("MW3A10P2C6_HIT|") or line.startswith("MW3A10P2C6_SOLVE|")
    ]
    if any("target_intersection=1" in line for line in interesting):
        status = "TARGET-P2-HIT"
    elif any("_HIT|" in line for line in interesting):
        status = "P2-HIT"
    else:
        status = "no-p2"
    return seed, hit_index, status, "\n".join(interesting)


counts = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = [executor.submit(run_one, job) for job in jobs]
    for future in concurrent.futures.as_completed(futures):
        seed, hit_index, status, detail = future.result()
        counts[status] = counts.get(status, 0) + 1
        print(
            f"MW3A10BATCH|seed={seed}|hit={hit_index}|status={status}",
            flush=True,
        )
        if status in ("P2-HIT", "TARGET-P2-HIT"):
            print(detail, flush=True)

print(
    "MW3A10BATCH|done=1|jobs=" + str(len(jobs))
    + "|counts=" + json.dumps(counts, sort_keys=True),
    flush=True,
)

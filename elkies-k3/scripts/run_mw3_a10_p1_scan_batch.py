#!/usr/bin/env python3
from pathlib import Path
import argparse
import concurrent.futures
import subprocess


ap = argparse.ArgumentParser(description="Run complete normalized A10/P1 coordinate slices in parallel.")
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument("--workers", type=int, default=3)
ap.add_argument("--dir", default="artifacts/local/elkies-k3/mw3-a10-p1")
ap.add_argument(
    "--input", default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.ms"
)
ap.add_argument(
    "--open-input", default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.open.ms"
)
args = ap.parse_args()

root = Path(args.dir)
root.mkdir(parents=True, exist_ok=True)


def run_seed(seed):
    sliced = root / f"component2-valid-seed{seed}.ms"
    build = subprocess.run(
        [
            "sage", "elkies-k3/scripts/build_mw3_a10_p1_slice.sage",
            "--input", args.input, "--open-input", args.open_input,
            "--out", str(sliced), "--seed", str(seed),
            "--kill", "rho,r1,lam", "--saturate", "",
        ],
        text=True,
        capture_output=True,
    )
    if build.returncode:
        return seed, "build-error", build.stderr[-1000:]
    scan = subprocess.run(
        [
            "sage", "elkies-k3/scripts/search_mw3_a10_p1_slice.sage",
            "--input", str(sliced),
            "--open-input", str(sliced.with_suffix(".open.ms")),
            "--nonzero", "s1,y4", "--max-hits", "20",
        ],
        text=True,
        capture_output=True,
    )
    if scan.returncode:
        return seed, "scan-error", scan.stderr[-1000:]
    (root / f"component2-valid-seed{seed}.scan.log").write_text(scan.stdout)
    hits = [line for line in scan.stdout.splitlines() if line.startswith("MW3A10SCAN_HIT|")]
    return seed, "ok", hits


total_hits = 0
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = [executor.submit(run_seed, seed) for seed in range(args.seed_start, args.seed_end + 1)]
    for future in concurrent.futures.as_completed(futures):
        seed, status, detail = future.result()
        if status != "ok":
            print(f"MW3A10SCANBATCH|seed={seed}|status={status}|detail={detail}", flush=True)
            continue
        total_hits += len(detail)
        print(f"MW3A10SCANBATCH|seed={seed}|status=ok|hits={len(detail)}", flush=True)

print(
    f"MW3A10SCANBATCH|done=1|seeds={args.seed_end-args.seed_start+1}|hits={total_hits}",
    flush=True,
)

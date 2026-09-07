#!/usr/bin/env python3
import argparse
import concurrent.futures
import subprocess


ap = argparse.ArgumentParser(description="Partition the persistent A10/P1 multislice scan.")
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--dir", default="artifacts/local/elkies-k3/mw3-a10-p1")
ap.add_argument("--input", default=None)
ap.add_argument("--open-input", default=None)
ap.add_argument("--deterministic-fixed", action="store_true")
ap.add_argument("--fixed-names", default="rho,r1,lam")
ap.add_argument("--nonzero-keep", default="s1,y4")
ap.add_argument("--prefix", default="component2-valid")
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
    ranges.append((lower, upper))
    lower = upper + 1


def run_range(bounds):
    lower, upper = bounds
    command = [
        "sage", "elkies-k3/scripts/search_mw3_a10_p1_multislice.sage",
        "--seed-start", str(lower), "--seed-end", str(upper), "--dir", args.dir,
        "--fixed-names", args.fixed_names,
        "--nonzero-keep", args.nonzero_keep,
        "--prefix", args.prefix,
    ]
    if args.input:
        command += ["--input", args.input]
    if args.open_input:
        command += ["--open-input", args.open_input]
    if args.deterministic_fixed:
        command.append("--deterministic-fixed")
    process = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )
    final = next(
        (line for line in reversed(process.stdout.splitlines()) if "|done=1|" in line),
        "",
    )
    return lower, upper, process.returncode, final, process.stderr[-1000:]


with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
    futures = [executor.submit(run_range, bounds) for bounds in ranges]
    for future in concurrent.futures.as_completed(futures):
        lower, upper, returncode, final, error = future.result()
        status = "ok" if returncode == 0 else "error"
        print(
            f"MW3A10MULTIBATCH|start={lower}|end={upper}|status={status}|summary={final}",
            flush=True,
        )
        if error:
            print(error, flush=True)

print(f"MW3A10MULTIBATCH|done=1|seeds={count}|workers={workers}", flush=True)

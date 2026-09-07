#!/usr/bin/env python3
from pathlib import Path
import argparse
import subprocess
import concurrent.futures
import time
import shutil

ap = argparse.ArgumentParser()
ap.add_argument("--p", type=int, default=101)
ap.add_argument("--seeds", default="1,2,3,4")
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--workers", type=int, default=2)
ap.add_argument("--timeout", type=int, default=300)
ap.add_argument(
    "--outdir",
    default="artifacts/local/elkies-k3/e6-square",
)
args = ap.parse_args()

if not shutil.which("sage"):
    raise SystemExit("sage missing")

if not shutil.which("msolve"):
    raise SystemExit("msolve missing")

O = Path(args.outdir)
O.mkdir(parents=True, exist_ok=True)

seeds = [
    int(x)
    for x in args.seeds.split(",")
    if x.strip()
]

def one(seed):
    inp = O / f"p{args.p}-seed{seed}-square.ms"
    sol = O / f"p{args.p}-seed{seed}.solve"
    log = O / f"p{args.p}-seed{seed}.log"

    build = subprocess.run(
        [
            "python3",
            "elkies-k3/scripts/export_e6_p1_square.py",
            "--p", str(args.p),
            "--seed", str(seed),
            "--slices", "3",
            "--out", str(inp),
        ],
        text=True,
        capture_output=True,
    )

    if build.returncode:
        return (
            seed,
            "build-error",
            0.0,
            build.stdout + "\n" + build.stderr,
        )

    t0 = time.time()

    try:
        with log.open("w") as h:
            solve = subprocess.run(
                [
                    "msolve",
                    "-t", str(args.threads),
                    "-v", "2",
                    "-f", str(inp),
                    "-o", str(sol),
                ],
                stdout=h,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=args.timeout,
            )

        dt = time.time() - t0

        logtail = (
            log.read_text(errors="replace")[-5000:]
            if log.exists()
            else ""
        )

        soltail = (
            sol.read_text(errors="replace")[-5000:]
            if sol.exists()
            else ""
        )

        return (
            seed,
            f"exit-{solve.returncode}",
            dt,
            build.stdout
            + "\nSOLVE\n"
            + soltail
            + "\nLOG\n"
            + logtail,
        )

    except subprocess.TimeoutExpired:
        dt = time.time() - t0

        logtail = (
            log.read_text(errors="replace")[-5000:]
            if log.exists()
            else ""
        )

        return (
            seed,
            "timeout",
            dt,
            build.stdout + "\nLOG\n" + logtail,
        )

print(
    f"E6SQBATCH|stage=start|p={args.p}|seeds={seeds}",
    flush=True,
)

with concurrent.futures.ThreadPoolExecutor(
    max_workers=args.workers
) as ex:

    futures = [
        ex.submit(one, seed)
        for seed in seeds
    ]

    for fut in concurrent.futures.as_completed(futures):
        seed, status, dt, text = fut.result()

        print(
            f"E6SQBATCH|seed={seed}|status={status}|seconds={dt:.1f}",
            flush=True,
        )

        for line in text.splitlines():
            if (
                line.startswith("E6SLICE|")
                or line.startswith("E6SQUARE|")
            ):
                print(line, flush=True)

        if status == "build-error":
            print(
                f"E6SQBATCH_ERROR|seed={seed}|"
                + text[-4000:].replace("\n", " | "),
                flush=True,
            )

        elif status.startswith("exit--"):
            print(
                f"E6SQBATCH_CRASH|seed={seed}|"
                + text[-4000:].replace("\n", " | "),
                flush=True,
            )

        elif status == "exit-0":
            print(
                f"E6SQBATCH_SOLVED|seed={seed}|"
                + text[-3000:].replace("\n", " | "),
                flush=True,
            )

print("E6SQBATCH|stage=done", flush=True)

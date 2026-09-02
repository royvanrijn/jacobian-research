#!/usr/bin/env python3
"""Run the contiguous 24A1 positive seven-octad completion frontier.

Invoke this with the repository's Sage Python.  It shards all exact five-prefix
orbits into disjoint half-open intervals, runs the existing shard enumerator,
and fails if any subprocess fails.  The manifest builder remains the authority
for gap/overlap checks and artifact hashes.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREFIX = (
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json"
)
ENUMERATOR = HERE / "enumerate_24a1_octad_rank7_completion_shard.sage"


def intervals(start: int, stop: int, size: int) -> list[tuple[int, int]]:
    return [(left, min(left + size, stop)) for left in range(start, stop, size)]


def run_interval(bounds: tuple[int, int], check: bool) -> tuple[int, int, str]:
    start, stop = bounds
    command = [
        sys.executable,
        str(ENUMERATOR),
        "--prefix-start",
        str(start),
        "--prefix-stop",
        str(stop),
    ]
    if check:
        command.append("--check")
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"24A1 completion shard {start}:{stop} failed\n{completed.stdout}"
        )
    return start, stop, completed.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix-start", type=int, default=0)
    parser.add_argument("--prefix-stop", type=int)
    parser.add_argument("--shard-size", type=int, default=250)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = json.loads(PREFIX.read_text())
    assert payload["schema"] == "elkies-k3.24a1-octad-prefix-orbits.v1"
    assert payload["status"] == "PASS_EXACT_M24_OCTAD_SUBSET_ORBITS_THROUGH_SIZE_5"
    assert payload["layers"][-1]["size"] == 5
    frontier_stop = len(payload["layers"][-1]["orbits"])
    stop = arguments.prefix_stop if arguments.prefix_stop is not None else frontier_stop
    assert 0 <= arguments.prefix_start < stop <= frontier_stop
    assert arguments.shard_size > 0 and arguments.jobs > 0
    work = intervals(arguments.prefix_start, stop, arguments.shard_size)
    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = {
            executor.submit(run_interval, bounds, arguments.check): bounds
            for bounds in work
        }
        for future in as_completed(futures):
            start, shard_stop, output = future.result()
            final_line = output.splitlines()[-1]
            print(
                f"OCTADFRONTIER|completed={start}:{shard_stop}|{final_line}",
                flush=True,
            )
    mode = "CHECK" if arguments.check else "WRITE"
    print(
        f"OCTADFRONTIER|mode={mode}|shards={len(work)}|"
        f"prefixes={arguments.prefix_start}:{stop}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

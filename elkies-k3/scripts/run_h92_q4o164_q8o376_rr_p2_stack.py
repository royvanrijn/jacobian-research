#!/usr/bin/env python3
"""Run modular q8/orbit376 scans, exact QQ reconstruction and pointing.

Invoke this with the repository's Sage Python. Independent prime scans are
run in parallel, complete v2 outputs are reused unless ``--force`` is set,
and the final prime is held out by default. After exact reconstruction the
known-section pointing probe runs unless explicitly disabled.
"""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
PROBE = ROOT / "elkies-k3/scripts/probe_h92_q4o164_q8o376_rr_p2_modp.sage"
RECONSTRUCT = ROOT / "elkies-k3/scripts/reconstruct_h92_q4o164_q8o376_rr_p2_qq.sage"
POINT = ROOT / "elkies-k3/scripts/point_h92_q4o164_q8o376_from_known_sections_qq.sage"
EXPECTED_SCHEMA = "elkies-k3.q4o164-q8o376-rr-p2-scan-modp.v2"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--primes",
    default="41,43,47,53,59",
    help="comma-separated distinct good primes; the final one is held out by default",
)
parser.add_argument(
    "--held-out",
    type=int,
    help="held-out prime; default is the final value from --primes",
)
parser.add_argument("--workers", type=int, default=2)
parser.add_argument("--python", default=sys.executable)
parser.add_argument("--output-dir", type=Path, default=LOCAL)
parser.add_argument(
    "--reconstruction-output",
    type=Path,
    default=LOCAL / "q4o164-q8o376-rr-p2-qq.json",
)
parser.add_argument(
    "--pointing-output",
    type=Path,
    default=LOCAL / "q4o164-q8o376-known-section-pointing-qq.json",
)
parser.add_argument("--skip-pointing", action="store_true")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

primes = [int(value.strip()) for value in args.primes.split(",") if value.strip()]
if len(primes) < 3 or len(set(primes)) != len(primes):
    raise SystemExit("provide at least three distinct primes")
held_out = args.held_out if args.held_out is not None else primes[-1]
if held_out not in primes:
    raise SystemExit("--held-out must occur in --primes")
if args.workers <= 0:
    raise SystemExit("--workers must be positive")

python = Path(args.python)
output_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
reconstruction_output = (
    args.reconstruction_output
    if args.reconstruction_output.is_absolute()
    else ROOT / args.reconstruction_output
)
pointing_output = (
    args.pointing_output
    if args.pointing_output.is_absolute()
    else ROOT / args.pointing_output
)
output_dir.mkdir(parents=True, exist_ok=True)


def scan_path(prime):
    return output_dir / f"q4o164-q8o376-rr-p2-scan-mod{prime}.json"


def reusable(path, prime):
    if args.force or not path.exists():
        return False
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if data.get("schema") != EXPECTED_SCHEMA:
        return False
    search = data.get("search_space", {})
    processed = int(search.get("processed", 0))
    complete = bool(
        search.get("complete", False)
        and search.get("stop_after") is None
        and processed == int(search.get("expected_size", -1))
    )
    return bool(data.get("prime") == prime and complete)


def run_scan(prime):
    output = scan_path(prime)
    if reusable(output, prime):
        return prime, output, "REUSED"
    command = [
        str(python),
        str(PROBE),
        "--prime",
        str(prime),
        "--output",
        str(output),
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    if not reusable(output, prime):
        raise RuntimeError(f"prime {prime} scan did not produce a complete reusable v2 artifact")
    return prime, output, "COMPUTED"


results = {}
with ThreadPoolExecutor(max_workers=min(args.workers, len(primes))) as executor:
    futures = {executor.submit(run_scan, prime): prime for prime in primes}
    for future in as_completed(futures):
        prime, path, status = future.result()
        results[prime] = path
        print(f"Q8O376STACK|prime={prime}|scan={status}|output={path}", flush=True)

construction = [prime for prime in primes if prime != held_out]
command = [str(python), str(RECONSTRUCT)]
command.extend(str(results[prime]) for prime in construction)
command.extend(["--held-out", str(results[held_out])])
command.extend(["--output", str(reconstruction_output)])
subprocess.run(command, cwd=ROOT, check=True)
print(
    "Q8O376STACK|construction={}|held_out={}|reconstruction=PASS|output={}".format(
        ",".join(map(str, construction)), held_out, reconstruction_output
    ),
    flush=True,
)

if not args.skip_pointing:
    # The pointing script consumes the reconstruction at its canonical local
    # path. A non-default reconstruction output is useful for experiments but
    # cannot be pointed implicitly without changing the proof input.
    canonical_reconstruction = LOCAL / "q4o164-q8o376-rr-p2-qq.json"
    if reconstruction_output.resolve() != canonical_reconstruction.resolve():
        raise RuntimeError(
            "known-section pointing requires the canonical reconstruction output; "
            "use --skip-pointing with a non-default --reconstruction-output"
        )
    subprocess.run([str(python), str(POINT)], cwd=ROOT, check=True)
    if pointing_output.resolve() != (
        LOCAL / "q4o164-q8o376-known-section-pointing-qq.json"
    ).resolve():
        pointing_output.parent.mkdir(parents=True, exist_ok=True)
        pointing_output.write_bytes(
            (LOCAL / "q4o164-q8o376-known-section-pointing-qq.json").read_bytes()
        )
    pointing = json.loads(pointing_output.read_text())
    print(
        "Q8O376STACK|pointing_status={}|degree1={}|output={}".format(
            pointing["status"], pointing["degree_one_count"], pointing_output
        ),
        flush=True,
    )

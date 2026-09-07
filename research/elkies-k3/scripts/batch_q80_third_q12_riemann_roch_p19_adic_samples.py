#!/usr/bin/env python3
"""Compile a residue-distinct batch of high-precision p=19 RR samples."""

import argparse
import concurrent.futures
import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
WORKER = ROOT / "elkies-k3/scripts/compile_q80_third_q12_riemann_roch_p19_adic_sample.sage"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--source",
    type=Path,
    default=RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision1028.json",
)
parser.add_argument(
    "--lift",
    type=Path,
    default=RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision1024.json",
)
parser.add_argument(
    "--basis",
    type=Path,
    default=RESULTS / "q80-third-q12-integral-basis-mod19-power-precision1024.json",
)
parser.add_argument(
    "--seed-pattern",
    default=str(RESULTS / "q80-third-q12-p19-adic-U*.json"),
)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=LOCAL / "q80-third-q12-p19-adic-precision1024-samples",
)
parser.add_argument(
    "--manifest",
    type=Path,
    default=RESULTS / "q80-third-q12-p19-adic-precision1024-sample-manifest.json",
)
parser.add_argument("--workers", type=int, default=8)
parser.add_argument(
    "--limit",
    type=int,
    default=20,
    help=(
        "number of residue-distinct seeds to compile; defaults to the 17 "
        "support-determining samples plus three p-adic holdouts"
    ),
)
args = parser.parse_args()
args.source = args.source.resolve()
args.lift = args.lift.resolve()
args.basis = args.basis.resolve()
args.output_dir = args.output_dir.resolve()
args.manifest = args.manifest.resolve()
if args.workers < 1:
    raise ValueError("worker count must be positive")
if args.limit < 1:
    raise ValueError("sample limit must be positive")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


lift = json.loads(args.lift.read_text())
digits = int(lift["specialization"]["digits"])
if lift.get("status") != "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19":
    raise ValueError("factor lift is not certified")

seeds = []
for path in sorted(Path().glob(args.seed_pattern) if not Path(args.seed_pattern).is_absolute() else []):
    seeds.append(path.resolve())
if not seeds:
    # pathlib does not support absolute glob patterns.
    import glob

    seeds = [Path(path).resolve() for path in sorted(glob.glob(args.seed_pattern))]
seeds = seeds[: args.limit]
if not seeds:
    raise ValueError("no residue-distinct seed samples found")

records = []
seen = set()
for seed in seeds:
    payload = json.loads(seed.read_text())
    base = tuple(int(value) for value in payload["specialization"]["base_U_coefficients_1_omega"])
    if base in seen:
        raise ArithmeticError(f"duplicate seed base {base}")
    seen.add(base)
    records.append((seed, base, args.output_dir / seed.name))

args.output_dir.mkdir(parents=True, exist_ok=True)


def output_is_current(path, base):
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return False
    return (
        payload.get("status") == "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE"
        and int(payload["specialization"]["digits"]) == digits
        and tuple(payload["specialization"]["base_U_coefficients_1_omega"]) == base
        and payload["inputs"]["source"]["sha256"] == sha256(args.source)
        and payload["inputs"]["lift"]["sha256"] == sha256(args.lift)
        and payload["inputs"]["basis"]["sha256"] == sha256(args.basis)
    )


def compile_one(record):
    seed, base, output = record
    if output_is_current(output, base):
        return output, "cached"
    command = [
        "sage",
        "-python",
        str(WORKER.relative_to(ROOT)),
        "--source",
        str(args.source.relative_to(ROOT)),
        "--lift",
        str(args.lift.relative_to(ROOT)),
        "--basis",
        str(args.basis.relative_to(ROOT)),
        "--base-constant",
        str(base[0]),
        "--base-anti",
        str(base[1]),
        "--output",
        str(output.relative_to(ROOT)),
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode:
        raise RuntimeError(
            f"sample U={base} failed with exit {completed.returncode}:\n{completed.stdout}"
        )
    if not output_is_current(output, base):
        raise ArithmeticError(f"sample U={base} did not produce a current certificate")
    return output, completed.stdout.strip().splitlines()[-1]


completed_paths = []
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = {executor.submit(compile_one, record): record for record in records}
    for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
        output, message = future.result()
        completed_paths.append(output)
        print(f"Q80THIRDQ12PADICBATCH|completed={index}/{len(records)}|{output.name}|{message}", flush=True)

completed_paths.sort()
manifest = {
    "schema": "elkies-k3.q80-third-q12-p19-adic-sample-manifest.v1",
    "status": "PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE_BATCH",
    "specialization": {
        "u": "-2",
        "prime": 19,
        "digits": digits,
        "sample_count": len(completed_paths),
        "residue_distinct_bases": True,
    },
    "inputs": {
        "source": {"path": str(args.source.relative_to(ROOT)), "sha256": sha256(args.source)},
        "lift": {"path": str(args.lift.relative_to(ROOT)), "sha256": sha256(args.lift)},
        "basis": {"path": str(args.basis.relative_to(ROOT)), "sha256": sha256(args.basis)},
    },
    "samples": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in completed_paths
    ],
    "worker": {
        "path": str(WORKER.relative_to(ROOT)),
        "sha256": sha256(WORKER),
    },
    "batch_worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "every listed sample independently passes the p-adic Riemann--Roch, long Weierstrass, and two-way map identities",
            "the listed sample bases are residue-distinct",
        ],
        "not_proved": [
            "a characteristic-zero Jacobian or map reconstruction",
            "generic interpolation from fewer than the support-determined number of samples",
        ],
    },
    "reproduce": shlex.join(
        [
            "python3",
            "elkies-k3/scripts/batch_q80_third_q12_riemann_roch_p19_adic_samples.py",
            "--source",
            str(args.source.relative_to(ROOT)),
            "--lift",
            str(args.lift.relative_to(ROOT)),
            "--basis",
            str(args.basis.relative_to(ROOT)),
            "--seed-pattern",
            args.seed_pattern,
            "--output-dir",
            str(args.output_dir.relative_to(ROOT)),
            "--manifest",
            str(args.manifest.relative_to(ROOT)),
            "--workers",
            str(args.workers),
            "--limit",
            str(args.limit),
        ]
    ),
}
serialized = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
args.manifest.parent.mkdir(parents=True, exist_ok=True)
args.manifest.write_text(serialized)
print(
    f"Q80THIRDQ12PADICBATCH|digits={digits}|samples={len(completed_paths)}|"
    "status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_P19_ADIC_SAMPLE_BATCH"
)

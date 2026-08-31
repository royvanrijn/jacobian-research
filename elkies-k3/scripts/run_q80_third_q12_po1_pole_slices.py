#!/usr/bin/env python3
"""Solve the Q80 pole-one marking scheme one finite pole at a time.

The direct pole-one export has a pole variable ``z``.  Eliminating ``z`` in
the full system is much more expensive than fixing its value in the prime
field.  This runner appends each equation ``z-c`` to the exported msolve
system, executes all prime-field slices with bounded parallelism, and writes
a resumable machine-readable summary.

This is an exhaustive computation for the displayed finite-field scheme
only.  A nonempty slice still has to be parsed and substituted into the
Weierstrass equation before it becomes a marked-section certificate.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
from time import monotonic


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "elkies-k3.q80-third-q12-po1-pole-slices.v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def parse_system(path: Path) -> tuple[list[str], int, list[str]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError("truncated msolve input")
    variables = [item.strip() for item in lines[0].split(",")]
    if not variables or variables[0] != "z" or len(set(variables)) != len(variables):
        raise ValueError("expected distinct variables with pole variable z first")
    prime = int(lines[1])
    equations = []
    for index, line in enumerate(lines[2:]):
        has_comma = line.endswith(",")
        if index + 1 < len(lines[2:]) and not has_comma:
            raise ValueError("nonfinal msolve equation lacks a comma")
        equations.append(line[:-1].strip() if has_comma else line)
    if any(not equation for equation in equations):
        raise ValueError("empty msolve equation")
    return variables, prime, equations


def classify_solution(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"classification": "missing_output"}
    text = path.read_text(errors="replace").strip()
    if text == "[-1]:":
        return {"classification": "empty_over_algebraic_closure"}
    if text.startswith("[1,"):
        return {"classification": "positive_dimensional"}
    match = re.match(r"\[0,\s*\[(\d+),\s*(\d+),\s*(\d+),", text)
    if match:
        return {
            "classification": "zero_dimensional",
            "field_characteristic": int(match.group(1)),
            "variable_count": int(match.group(2)),
            "quotient_dimension": int(match.group(3)),
            "output_prefix": text[:240],
        }
    return {"classification": "unparsed_output", "output_prefix": text[:240]}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--system", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--summary", type=Path, required=True)
parser.add_argument("--jobs", type=int, default=4)
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--timeout", type=float, default=120.0)
parser.add_argument(
    "--pole", action="append", type=int, default=[],
    help="pole value to include; repeat (default: every prime-field value)",
)
parser.add_argument(
    "--slice-leading", action="store_true",
    help="also split by every nonzero value of the leading coefficient l",
)
parser.add_argument(
    "--leading", action="append", type=int, default=[],
    help="nonzero leading value to include; repeat and use with --slice-leading",
)
parser.add_argument(
    "--substitute-slices", action="store_true",
    help=(
        "evaluate sliced variables and combine coefficients before msolve; "
        "requires execution with `sage -python`"
    ),
)
parser.add_argument("--msolve", type=Path, default=Path(shutil.which("msolve") or "msolve"))
args = parser.parse_args()
if args.jobs < 1 or args.threads < 1 or args.timeout <= 0:
    parser.error("jobs, threads, and timeout must be positive")

system_path = args.system.resolve()
variables, prime, equations = parse_system(system_path)
msolve_path = args.msolve.resolve()
output_dir = args.output_dir.resolve()
summary_path = args.summary.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
summary_path.parent.mkdir(parents=True, exist_ok=True)

run_key = {
    "source_sha256": digest(system_path),
    "msolve_sha256": digest(msolve_path),
    "jobs": args.jobs,
    "threads_per_process": args.threads,
    "timeout_seconds_per_slice": args.timeout,
    "slice_leading": args.slice_leading,
    "substitute_slices": args.substitute_slices,
}

poles = sorted(set(value % prime for value in (args.pole or range(prime))))
if args.leading and not args.slice_leading:
    parser.error("--leading requires --slice-leading")
leading_values = (
    sorted(set(value % prime for value in (args.leading or range(1, prime))))
    if args.slice_leading else [None]
)
if args.slice_leading and (not leading_values or 0 in leading_values):
    parser.error("leading slices must be nonzero modulo the prime")


def slice_paths(value: int, leading: int | None) -> tuple[Path, Path, Path, Path]:
    stem = f"z-{value:02d}" + ("" if leading is None else f"-l-{leading:02d}")
    return tuple(output_dir / f"{stem}.{suffix}" for suffix in ("ms", "solve", "log", "json"))


def write_slice(value: int, leading: int | None, path: Path) -> None:
    if args.substitute_slices:
        sliced_names = {"z"}
        sliced_values = {"z": value}
        if leading is not None:
            sliced_names.add("l")
            sliced_values["l"] = leading
        remaining = [name for name in variables if name not in sliced_names]
        target_ring = PolynomialRing(finite_field, names=remaining, order="degrevlex")
        target_by_name = target_ring.gens_dict()
        images = [
            finite_field(sliced_values[name])
            if name in sliced_values else target_by_name[name]
            for name in variables
        ]
        specialization = source_ring.hom(images, target_ring)
        specialized = []
        seen = set()
        for equation in source_equations:
            image = target_ring(specialization(equation))
            if not image:
                continue
            text = str(image).replace("**", "^")
            if text not in seen:
                specialized.append(text)
                seen.add(text)
        output_variables = remaining
        sliced = specialized
    else:
        output_variables = variables
        sliced = equations + [f"z-{value}"]
        if leading is not None:
            sliced.append(f"l-{leading}")
    with path.open("w", encoding="utf-8") as stream:
        stream.write(",".join(output_variables) + "\n")
        stream.write(str(prime) + "\n")
        for index, equation in enumerate(sliced):
            stream.write(equation)
            stream.write(",\n" if index + 1 < len(sliced) else "\n")


def run_slice(value: int, leading: int | None) -> dict[str, object]:
    input_path, solution_path, log_path, checkpoint_path = slice_paths(value, leading)
    key = {**run_key, "pole": value, "leading": leading}
    if checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text())
        if checkpoint.get("run_key") == key and checkpoint.get("status") in {
            "completed", "timeout", "error"
        }:
            return checkpoint
    solution_path.unlink(missing_ok=True)
    command = [
        str(msolve_path), "-f", str(input_path), "-o", str(solution_path),
        "-t", str(args.threads), "-v", "1",
    ]
    started = monotonic()
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )
    timed_out = False
    try:
        stdout, unused_stderr = process.communicate(timeout=args.timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, unused_stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, unused_stderr = process.communicate()
    log_path.write_text(stdout, encoding="utf-8")
    classification = (
        {"classification": "timeout"}
        if timed_out
        else classify_solution(solution_path)
    )
    status = (
        "timeout" if timed_out else "completed" if process.returncode == 0 else "error"
    )
    result = {
        "run_key": key,
        "status": status,
        "pole": value,
        "leading": leading,
        "elapsed_seconds": monotonic() - started,
        "returncode": process.returncode,
        "system": relative(input_path),
        "system_sha256": digest(input_path),
        "solution": relative(solution_path) if solution_path.exists() else None,
        "log": relative(log_path),
        **classification,
    }
    checkpoint_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


started = monotonic()
results = []
if args.substitute_slices:
    try:
        from sage.all import GF, PolynomialRing
    except ImportError as error:
        raise RuntimeError(
            "--substitute-slices requires `sage -python SCRIPT ...`"
        ) from error
    finite_field = GF(prime)
    source_ring = PolynomialRing(finite_field, names=variables, order="degrevlex")
    source_equations = tuple(
        source_ring(equation.replace("^", "**")) for equation in equations
    )
for value in poles:
    for leading in leading_values:
        input_path, unused_solution, unused_log, unused_checkpoint = slice_paths(
            value, leading
        )
        write_slice(value, leading, input_path)
with ThreadPoolExecutor(max_workers=args.jobs) as executor:
    futures = {
        executor.submit(run_slice, value, leading): (value, leading)
        for value in poles for leading in leading_values
    }
    for future in as_completed(futures):
        result = future.result()
        results.append(result)
        print(
            "Q80PO1SLICE|pole={}|leading={}|status={}|classification={}|seconds={:.3f}".format(
                result["pole"], result["leading"], result["status"], result["classification"],
                result["elapsed_seconds"],
            ),
            flush=True,
        )
results.sort(key=lambda item: (int(item["pole"]), -1 if item["leading"] is None else int(item["leading"])))

counts = {}
for result in results:
    key = str(result["classification"])
    counts[key] = counts.get(key, 0) + 1
status = (
    "PASS_EXHAUSTIVE_POLE_SLICES"
    if all(result["status"] == "completed" for result in results)
    else "INCOMPLETE_BOUNDED_POLE_SLICES"
)
output = {
    "schema": SCHEMA,
    "status": status,
    "source": {"path": relative(system_path), "sha256": run_key["source_sha256"]},
    "prime": prime,
    "selected_poles": poles,
    "selected_leading_values": leading_values,
    "run_key": run_key,
    "classification_counts": counts,
    "slices": results,
    "runtime_seconds": monotonic() - started,
    "claim_boundary": {
        "proved": [
            (
                "exhaustive selected pole/leading partition of the displayed "
                "finite-field scheme"
            )
        ] if status.startswith("PASS") else [],
        "not_proved": [
            "parsed finite-field section coordinates",
            "literal section substitution",
            "characteristic-zero lifting",
        ],
    },
}
summary_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80PO1SLICE|prime={prime}|counts={counts}|output={relative(summary_path)}|"
    f"status={status}",
    flush=True,
)

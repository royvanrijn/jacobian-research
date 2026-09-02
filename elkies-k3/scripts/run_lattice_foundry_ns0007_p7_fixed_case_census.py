#!/usr/bin/env python3
"""Exhaust base-field leading cases for the reduced NS0007 pole-zero system.

The reduced system has six especially useful base-field coordinates:

``a2_4,a2_3,a2_2,a2_1,si_0,sl_0``.

Fixing them by six linear equations changes a difficult monolithic Gröbner
problem into small independent cases.  This runner enumerates all ``p^6``
cases, appends the linear equations to the pinned common system, and runs
msolve with bounded parallelism.  Temporary inputs use a copy-on-write clone
when the filesystem supports ``FICLONE``; the fallback is an ordinary copy.

In the default ``leading-ideal`` mode, the result is an exact unit-versus-
nonunit-ideal census of the declared fixed-lambda finite-field chart when no
case times out or errors.  Exact index intervals make long runs resumable;
each worker has a private current and temporary directory.  A nonunit ideal
still needs a clean expanded-system replay, solution decoding, base-field
rationality, and the independent geometric gates in the source certificate.
"""

from __future__ import annotations

import argparse
import base64
import fcntl
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import tempfile
import threading
import time
import zlib
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/p7-lambda2.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-census-mod7.json"
)
FIXED_VARIABLES = ["a2_4", "a2_3", "a2_2", "a2_1", "si_0", "sl_0"]
FICLONE = 0x40049409
STATUS_CODE = {
    "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE": 0,
    "NONUNIT_IDEAL": 1,
    "FINITE_SOLUTION_SET": 2,
    "POSITIVE_DIMENSIONAL_SOLUTION_SET": 3,
    "TIMEOUT": 4,
    "SOLVER_ERROR": 5,
}


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clone_file(source: Path, target: Path) -> str:
    with source.open("rb") as source_handle, target.open("w+b") as target_handle:
        try:
            fcntl.ioctl(target_handle.fileno(), FICLONE, source_handle.fileno())
            return "FICLONE_COPY_ON_WRITE"
        except OSError:
            source_handle.seek(0)
            shutil.copyfileobj(source_handle, target_handle)
            return "ORDINARY_COPY_FALLBACK"


def classify_output(text: str, solver_mode: str) -> str:
    stripped = text.strip()
    if solver_mode == "leading-ideal":
        if "#length of basis:" not in stripped:
            return "SOLVER_ERROR"
        if stripped.splitlines()[-1].strip() == "[]:":
            return "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE"
        return "NONUNIT_IDEAL"
    if stripped == "[-1]:":
        return "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE"
    if stripped.startswith("[1,") and ",-1,[]]:" in stripped.replace(" ", ""):
        return "POSITIVE_DIMENSIONAL_SOLUTION_SET"
    if stripped:
        return "FINITE_SOLUTION_SET"
    return "SOLVER_ERROR"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--msolve", default="msolve")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--solver-threads", type=int, default=1)
    parser.add_argument(
        "--solver-mode",
        choices=("leading-ideal", "solve"),
        default="leading-ideal",
        help=(
            "leading-ideal classifies unit versus nonunit ideals without the "
            "concurrency-sensitive rational-parametrization stage"
        ),
    )
    parser.add_argument(
        "--allow-factored-input-experiment",
        action="store_true",
        help=(
            "allow a non-expanded msolve input for diagnostics only; factored "
            "syntax has produced false nonunit cases and is not certifying"
        ),
    )
    parser.add_argument("--timeout-seconds", type=float, default=5.0)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--progress-every", type=int, default=5000)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.solver_threads < 1:
        raise SystemExit("worker and solver-thread counts must be positive")
    if arguments.timeout_seconds <= 0:
        raise SystemExit("--timeout-seconds must be positive")
    if arguments.start_index < 0:
        raise SystemExit("--start-index must be nonnegative")
    if arguments.limit is not None and arguments.limit < 1:
        raise SystemExit("--limit must be positive")

    metadata_path = arguments.input.resolve()
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("schema") != (
        "elkies-k3.lattice-foundry-ns0007-pole0-reduced-modp-system.v1"
    ):
        raise ValueError("unexpected reduced NS0007 system schema")
    prime = int(metadata["prime"])
    system_path = (ROOT / metadata["system"]["msolve_input"]).resolve()
    polynomial_encoding = metadata["system"].get(
        "msolve_polynomial_encoding", "FULLY_EXPANDED"
    )
    if (
        polynomial_encoding != "FULLY_EXPANDED"
        and not arguments.allow_factored_input_experiment
    ):
        raise SystemExit(
            "refusing non-expanded msolve syntax: a GF(7) cross-check produced "
            "false nonunit cases; pass --allow-factored-input-experiment only "
            "for diagnostics"
        )
    if digest(system_path) != metadata["system"]["msolve_sha256"]:
        raise ArithmeticError("reduced NS0007 system digest mismatch")
    system_lines = system_path.read_text().splitlines()
    names = system_lines[0].split(",")
    if any(name not in names for name in FIXED_VARIABLES):
        raise ValueError("a fixed-case variable is absent from the source system")

    total_cases = prime ** len(FIXED_VARIABLES)
    if arguments.start_index >= total_cases:
        raise SystemExit("--start-index lies beyond the fixed-case census")
    available_cases = total_cases - arguments.start_index
    selected_cases = available_cases if arguments.limit is None else min(
        available_cases, arguments.limit
    )
    case_iterator = itertools.islice(
        itertools.product(range(prime), repeat=len(FIXED_VARIABLES)),
        arguments.start_index,
        arguments.start_index + selected_cases,
    )
    statuses = bytearray([255]) * selected_cases
    exceptional = []
    completed = 0
    clone_methods = Counter()
    started = time.monotonic()
    thread_state = threading.local()

    with tempfile.TemporaryDirectory(prefix="ns0007-fixed-census-") as temp_name:
        temp_dir = Path(temp_name)
        common_path = temp_dir / "common.ms"
        common_path.write_text(system_path.read_text().rstrip() + ",\n")

        def run_case(
            local_index: int,
            global_index: int,
            values: tuple[int, ...],
        ) -> dict:
            if not hasattr(thread_state, "directory"):
                thread_state.directory = temp_dir / f"worker-{threading.get_ident()}"
                thread_state.directory.mkdir()
            worker_dir = thread_state.directory
            input_path = worker_dir / "case.ms"
            output_path = worker_dir / "case.solve"
            clone_method = clone_file(common_path, input_path)
            fixed_equations = [
                f"{name}-{value}"
                for name, value in zip(FIXED_VARIABLES, values)
            ]
            with input_path.open("a") as handle:
                handle.write(",\n".join(fixed_equations) + "\n")
            case_started = time.monotonic()
            try:
                command = [
                        arguments.msolve,
                        "-f",
                        str(input_path),
                        "-o",
                        str(output_path),
                        "-t",
                        str(arguments.solver_threads),
                        "-v",
                        "0",
                    ]
                if arguments.solver_mode == "leading-ideal":
                    command += ["-g", "1"]
                environment = dict(os.environ)
                environment["TMPDIR"] = str(worker_dir)
                process = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=arguments.timeout_seconds,
                    check=False,
                    cwd=worker_dir,
                    env=environment,
                )
                elapsed = time.monotonic() - case_started
                output_text = output_path.read_text() if output_path.is_file() else ""
                if process.returncode:
                    status = "SOLVER_ERROR"
                else:
                    status = classify_output(output_text, arguments.solver_mode)
                detail = {
                    "index": global_index,
                    "local_index": local_index,
                    "values": list(values),
                    "status": status,
                    "returncode": process.returncode,
                    "solution_output": output_text if status != "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE" else None,
                    "stdout_tail": process.stdout[-2000:] if process.stdout else "",
                    "stderr_tail": process.stderr[-2000:] if process.stderr else "",
                    "clone_method": clone_method,
                }
            except subprocess.TimeoutExpired as error:
                elapsed = time.monotonic() - case_started
                detail = {
                    "index": global_index,
                    "local_index": local_index,
                    "values": list(values),
                    "status": "TIMEOUT",
                    "returncode": None,
                    "solution_output": None,
                    "stdout_tail": (error.stdout or "")[-2000:] if isinstance(error.stdout, str) else "",
                    "stderr_tail": (error.stderr or "")[-2000:] if isinstance(error.stderr, str) else "",
                    "clone_method": clone_method,
                }
            finally:
                input_path.unlink(missing_ok=True)
                output_path.unlink(missing_ok=True)
            return detail

        pending = {}
        with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
            for local_index, values in enumerate(case_iterator):
                global_index = arguments.start_index + local_index
                future = executor.submit(
                    run_case, local_index, global_index, values
                )
                pending[future] = local_index
                if len(pending) < 2 * arguments.workers:
                    continue
                done, unused = wait(pending, return_when=FIRST_COMPLETED)
                for finished in done:
                    detail = finished.result()
                    pending.pop(finished)
                    status = detail["status"]
                    statuses[detail["local_index"]] = STATUS_CODE[status]
                    clone_methods[detail["clone_method"]] += 1
                    if status != "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE":
                        exceptional.append(detail)
                    completed += 1
                    if arguments.progress_every and completed % arguments.progress_every == 0:
                        print(
                            "FOUNDRYNS0007CENSUS|"
                            f"completed={completed}/{selected_cases}|"
                            f"exceptional={len(exceptional)}|"
                            f"elapsed={time.monotonic()-started:.1f}|status=RUNNING",
                            flush=True,
                        )
            while pending:
                done, unused = wait(pending, return_when=FIRST_COMPLETED)
                for finished in done:
                    detail = finished.result()
                    pending.pop(finished)
                    status = detail["status"]
                    statuses[detail["local_index"]] = STATUS_CODE[status]
                    clone_methods[detail["clone_method"]] += 1
                    if status != "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE":
                        exceptional.append(detail)
                    completed += 1

    if any(code == 255 for code in statuses):
        raise AssertionError("fixed-case status vector is incomplete")
    status_histogram = Counter(statuses)
    inverse_status = {code: name for name, code in STATUS_CODE.items()}
    named_histogram = {
        inverse_status[code]: count
        for code, count in sorted(status_histogram.items())
    }
    exhaustive = arguments.start_index == 0 and selected_cases == total_cases
    complete_without_open_cases = not any(
        named_histogram.get(name, 0)
        for name in ("TIMEOUT", "SOLVER_ERROR")
    )
    nonunit_count = named_histogram.get("NONUNIT_IDEAL", 0)
    if (
        exhaustive
        and complete_without_open_cases
        and arguments.solver_mode == "leading-ideal"
    ):
        status = (
            "PASS_EXHAUSTIVE_FIXED_CASE_UNIT_IDEAL_CENSUS"
            if not nonunit_count
            else "PASS_EXHAUSTIVE_FIXED_CASE_LEADING_IDEAL_CENSUS_WITH_NONUNIT_CASES"
        )
    elif exhaustive and complete_without_open_cases:
        status = "PASS_EXHAUSTIVE_FIXED_CASE_ALGEBRAIC_CLOSURE_CENSUS"
    elif complete_without_open_cases and arguments.solver_mode == "leading-ideal":
        status = "PASS_BOUNDED_FIXED_CASE_RANGE_LEADING_IDEAL_CENSUS"
    elif complete_without_open_cases:
        status = "PASS_BOUNDED_FIXED_CASE_RANGE_ALGEBRAIC_CLOSURE_CENSUS"
    else:
        status = "OPEN_TIMEOUT_OR_SOLVER_ERROR_CASES_REMAIN"
    compressed_statuses = zlib.compress(bytes(statuses), level=9)
    output = {
        "schema": "elkies-k3.lattice-foundry-ns0007-pole0-fixed-case-census-modp.v1",
        "status": status,
        "input": {
            "metadata": relative(metadata_path),
            "metadata_sha256": digest(metadata_path),
            "system": relative(system_path),
            "system_sha256": digest(system_path),
        },
        "prime": prime,
        "lambda": metadata["lambda"],
        "fixed_variables_in_enumeration_order": FIXED_VARIABLES,
        "search": {
            "total_cases": total_cases,
            "selected_start_index": arguments.start_index,
            "selected_stop_index_exclusive": arguments.start_index + selected_cases,
            "selected_cases": selected_cases,
            "exhaustive": exhaustive,
            "workers": arguments.workers,
            "solver_threads_per_case": arguments.solver_threads,
            "solver_mode": arguments.solver_mode,
            "msolve_polynomial_encoding": polynomial_encoding,
            "certifying_input_syntax": polynomial_encoding == "FULLY_EXPANDED",
            "timeout_seconds_per_case": arguments.timeout_seconds,
            "msolve": arguments.msolve,
            "clone_methods": dict(sorted(clone_methods.items())),
        },
        "accounting": {
            "status_histogram": named_histogram,
            "exceptional_case_count": len(exceptional),
        },
        "status_vector": {
            "encoding": (
                "zlib-compressed bytes for the selected half-open index range "
                "in lexicographic product order; codes "
                + json.dumps(STATUS_CODE, sort_keys=True)
            ),
            "base64": base64.b64encode(compressed_statuses).decode("ascii"),
            "uncompressed_sha256": hashlib.sha256(bytes(statuses)).hexdigest(),
        },
        "exceptional_cases": sorted(exceptional, key=lambda row: row["index"]),
        "proof_boundary": {
            "proved": (
                "Each completed case adds exact linear equations for the six "
                "displayed base-field values. Leading-ideal mode proves unit "
                "versus nonunit ideal by an exact finite-field Groebner basis; "
                "solve mode additionally requests dimension and parametrization."
            ),
            "complete_if": (
                "This artifact exactly covers its declared half-open index range "
                "when all selected cases complete without timeout or solver error; "
                "the fixed-lambda chart is exhaustive only after contiguous coverage "
                "of all p^6 cases."
            ),
            "open": (
                "Finite solutions require decoding, base-field rationality, exact "
                "Kodaira orders, residual squarefreeness, NS marking, and a "
                "characteristic-zero lift."
            ),
        },
        "reproduce": " ".join(
            [
                "python3",
                "elkies-k3/scripts/run_lattice_foundry_ns0007_p7_fixed_case_census.py",
                "--input",
                relative(metadata_path),
                "--output",
                relative(arguments.output.resolve()),
                "--msolve",
                arguments.msolve,
                "--workers",
                str(arguments.workers),
                "--solver-threads",
                str(arguments.solver_threads),
                "--solver-mode",
                arguments.solver_mode,
                "--timeout-seconds",
                str(arguments.timeout_seconds),
                "--start-index",
                str(arguments.start_index),
                "--limit",
                str(selected_cases),
                *(
                    ["--allow-factored-input-experiment"]
                    if arguments.allow_factored_input_experiment
                    else []
                ),
            ]
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("NS0007 fixed-case census artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007CENSUS|"
        f"cases={selected_cases}/{total_cases}|"
        f"exceptional={len(exceptional)}|status={status}",
        flush=True,
    )


if __name__ == "__main__":
    main()

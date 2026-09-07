#!/usr/bin/env sage
"""Search target-fitted quartics from distinct norm-eight R17 traces.

status: ACTIVE_SEARCH
claim: exhaustive mixed-trace simultaneous-splitting search in a compact box
inputs: published R17 data and the exact genus-one bisection constructor
outputs: generated mixed-trace splitting artifact

The seven equation-cheapest finite-pole norm-eight traces each give eleven
quartics fitted through the public rank-28 complement.  This script constructs
all 77 exactly, then searches every primitive ``t=a/b`` in the requested box.
The modular hot loop retains a parameter only while quartics from at least two
distinct trace pencils can still split; all survivors receive exact square
tests.  The result is bounded and is not a global rational-point theorem.
"""

from __future__ import annotations

import argparse
from array import array
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
CONSTRUCTOR = SCRIPTS / "search_elkies_2026_rank28_genus_one_bisections.sage"
SCANNER = SCRIPTS / "scan_elkies_2026_rank28_mixed_trace_splitting.cpp"
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-mixed-trace-splitting-h10000-v1.json"
)
LOCAL = ROOT / "artifacts/local/elkies-k3/r17-rank28-mixed-trace-splitting"
MAGIC = 0x4D475331
PRIMES = (59, 73, 97, 101, 107, 131, 137, 149, 163, 167, 173, 179, 181, 199, 211, 223)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rational_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def normalized_integer_quartic(record: dict) -> tuple[int, tuple[int, ...]]:
    coefficients = [
        Fraction(value)
        for value in record["branch_polynomial_q_coefficients_low_to_high"]
    ]
    denominator = 1
    for value in coefficients:
        denominator = math.lcm(denominator, value.denominator)
    integral = tuple(
        value.numerator * (denominator // value.denominator)
        for value in coefficients
    )
    if len(integral) != 5:
        raise ArithmeticError("expected five quartic coefficients")
    return denominator, integral


def homogeneous_value(curve, numerator: int, denominator: int) -> int:
    _scalar, coefficients = curve
    b2 = denominator * denominator
    b3 = b2 * denominator
    b4 = b2 * b2
    return (
        (((coefficients[4] * numerator + coefficients[3] * denominator) * numerator
          + coefficients[2] * b2) * numerator + coefficients[1] * b3) * numerator
        + coefficients[0] * b4
    )


def exact_square_root(curve, numerator: int, denominator: int) -> int | None:
    scalar, _coefficients = curve
    radicand = scalar * homogeneous_value(curve, numerator, denominator)
    if radicand < 0:
        return None
    root = math.isqrt(radicand)
    return root if root * root == radicand else None


def construct_quartics(path: Path, trace_count: int) -> tuple[list[str], dict]:
    command = [
        sys.executable,
        str(CONSTRUCTOR),
        "--trace-limit",
        str(trace_count),
        "--output",
        str(path),
    ]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    document = json.loads(path.read_text())
    if document.get("selected_trace_count") != trace_count:
        raise ArithmeticError("constructor returned the wrong trace count")
    if document.get("successful_trace_target_pairs") != 11 * trace_count:
        raise ArithmeticError("not every mixed trace-target construction succeeded")
    return command, document


def export_tables(path: Path, curves, trace_count: int, targets_per_trace: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(
            struct.pack(
                "<IIIII",
                MAGIC,
                len(curves),
                trace_count,
                targets_per_trace,
                len(PRIMES),
            )
        )
        for prime in PRIMES:
            if any(scalar % prime == 0 for scalar, _ in curves):
                raise ArithmeticError(f"bad mixed-trace residue prime {prime}")
            squares = {value * value % prime for value in range(prime)}
            masks = array("Q")
            reduced = [
                (scalar % prime, tuple(value % prime for value in coefficients))
                for scalar, coefficients in curves
            ]
            for numerator in range(prime):
                for denominator in range(prime):
                    low = 0
                    high = 0
                    for index, curve in enumerate(reduced):
                        value = curve[0] * homogeneous_value(
                            (1, curve[1]), numerator, denominator
                        )
                        if value % prime not in squares:
                            continue
                        if index < 64:
                            low |= 1 << index
                        else:
                            high |= 1 << (index - 64)
                    masks.extend((low, high))
            stream.write(struct.pack("<I", prime))
            stream.write(masks.tobytes())


def compile_and_scan(table: Path, candidates: Path, binary: Path, height: int):
    compiler = shutil.which("g++")
    if compiler is None:
        raise SystemExit("g++ is required for the exhaustive mixed-trace scan")
    compile_command = [
        compiler,
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(SCANNER),
        "-o",
        str(binary),
    ]
    subprocess.run(compile_command, check=True)
    scan_command = [str(binary), str(table), str(height), str(candidates)]
    completed = subprocess.run(scan_command, check=True, text=True, capture_output=True)
    summary = completed.stdout.strip()
    if "|status=PASS" not in summary:
        raise ArithmeticError("mixed-trace scanner did not return PASS")
    return compile_command, scan_command, summary


def set_indices(low: int, high: int, count: int):
    for index in range(count):
        if index < 64:
            present = (low >> index) & 1
        else:
            present = (high >> (index - 64)) & 1
        if present:
            yield index


def exact_results(candidate_path: Path, curves, records, targets_per_trace: int):
    exact_tests = 0
    simultaneous = []
    maximum_exact_trace_count = 0
    maximum_exact_quartic_count = 0
    for line in candidate_path.read_text().splitlines():
        numerator, denominator, low, high = map(int, line.split())
        split = []
        for index in set_indices(low, high, len(curves)):
            exact_tests += 1
            root = exact_square_root(curves[index], numerator, denominator)
            if root is not None:
                trace_index = index // targets_per_trace
                target_index = index % targets_per_trace
                split.append(
                    {
                        "flat_index_zero_based": index,
                        "trace_index_one_based": trace_index + 1,
                        "target_index_one_based": target_index + 1,
                        "target_label": records[index]["target_label"],
                        "integral_square_root": str(root),
                    }
                )
        trace_indices = sorted({row["trace_index_one_based"] for row in split})
        maximum_exact_trace_count = max(maximum_exact_trace_count, len(trace_indices))
        maximum_exact_quartic_count = max(maximum_exact_quartic_count, len(split))
        if len(trace_indices) >= 2:
            simultaneous.append(
                {
                    "t": rational_text(Fraction(numerator, denominator)),
                    "projective_pair": [numerator, denominator],
                    "split_trace_indices_one_based": trace_indices,
                    "split_quartic_count": len(split),
                    "split_quartics": split,
                }
            )
    return {
        "modular_candidate_count": sum(1 for _ in candidate_path.open()),
        "exact_square_tests": exact_tests,
        "maximum_exact_split_trace_count": maximum_exact_trace_count,
        "maximum_exact_split_quartic_count": maximum_exact_quartic_count,
        "mixed_trace_simultaneous_splits": simultaneous,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=10_000)
    parser.add_argument("--trace-count", type=int, default=7)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--local-directory", type=Path, default=LOCAL)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.height < 1 or not 2 <= args.trace_count <= 11:
        parser.error("require positive height and 2 <= trace-count <= 11")

    args.local_directory.mkdir(parents=True, exist_ok=True)
    construction_path = args.local_directory / "fitted-quartics.json"
    table_path = args.local_directory / "mixed-square-residue-masks.bin"
    candidate_path = args.local_directory / "mixed-compact-candidates.txt"
    binary_path = args.local_directory / "scan-mixed-trace-splitting"
    construction_command, document = construct_quartics(
        construction_path, args.trace_count
    )
    traces = document["traces"]
    targets_per_trace = len(traces[0]["targets"])
    if targets_per_trace != 11 or any(
        len(trace["targets"]) != targets_per_trace for trace in traces
    ):
        raise ArithmeticError("mixed traces do not share the eleven target labels")
    records = [record for trace in traces for record in trace["targets"]]
    curves = [normalized_integer_quartic(record) for record in records]
    export_tables(table_path, curves, args.trace_count, targets_per_trace)
    compile_command, scan_command, scanner_summary = compile_and_scan(
        table_path, candidate_path, binary_path, args.height
    )
    exact = exact_results(candidate_path, curves, records, targets_per_trace)
    t0_rows = [
        row
        for row in exact["mixed_trace_simultaneous_splits"]
        if row["t"] == document["parameter"]
    ]
    expected_count = args.trace_count * targets_per_trace
    if len(t0_rows) != 1 or t0_rows[0]["split_quartic_count"] != expected_count:
        raise ArithmeticError("the mixed-trace t0 positive control was not recovered")
    new_rows = [
        row
        for row in exact["mixed_trace_simultaneous_splits"]
        if row["t"] != document["parameter"]
    ]
    status = (
        "PASS_NEW_MIXED_TRACE_SPLITS_REQUIRE_INDEPENDENCE"
        if new_rows
        else "PASS_EXHAUSTIVE_BOUNDED_NO_NEW_MIXED_TRACE_SPLIT"
    )
    result = {
        "schema": "elkies-k3.r17-rank28-mixed-trace-splitting.v1",
        "status": status,
        "construction": {
            "trace_count": args.trace_count,
            "targets_per_trace": targets_per_trace,
            "quartic_count": len(curves),
            "trace_published_basis_vectors": [
                trace["published_basis_w"] for trace in traces
            ],
            "all_quartics_irreducible_squarefree": all(
                record["branch_polynomial_irreducible_over_Q"]
                and record["branch_polynomial_squarefree"]
                for record in records
            ),
            "constructor_command": construction_command,
            "constructor_output_sha256": file_sha256(construction_path),
        },
        "compact_projective_scan": {
            "region": (
                f"primitive t=a/b with |a| <= {args.height}, "
                f"1 <= b <= {args.height}"
            ),
            "height": args.height,
            "residue_primes": list(PRIMES),
            "scanner_summary": scanner_summary,
            **exact,
            "new_mixed_trace_splits_away_from_t0": new_rows,
        },
        "specialization_and_quotient_independence": (
            "required for every new row before promotion"
            if new_rows
            else "vacuous: no new exact mixed-trace split survived"
        ),
        "proof_boundary": (
            "The t-line enumeration is exhaustive only in the displayed primitive box and "
            "only for the selected equation-cheapest trace prefix. An empty result is not a "
            "global rational-point obstruction and does not test the remaining norm-eight "
            "trace classes or trisections."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_elkies_2026_rank28_mixed_trace_splitting.sage"
        ),
        "implementation": {
            "script": str(Path(__file__).relative_to(ROOT)),
            "script_sha256": file_sha256(Path(__file__)),
            "constructor": str(CONSTRUCTOR.relative_to(ROOT)),
            "constructor_sha256": file_sha256(CONSTRUCTOR),
            "scanner": str(SCANNER.relative_to(ROOT)),
            "scanner_sha256": file_sha256(SCANNER),
            "compile_command": compile_command,
            "scan_command": scan_command,
            "residue_table_sha256": file_sha256(table_path),
        },
    }
    if args.check:
        if json.loads(args.output.read_text()) != result:
            raise ArithmeticError("stored mixed-trace search artifact changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"R17MIXEDTRACEG1|traces={args.trace_count}|quartics={len(curves)}|"
        f"height={args.height}|modular={exact['modular_candidate_count']}|"
        f"new={len(new_rows)}|status={status}|output={args.output}"
    )


if __name__ == "__main__":
    main()

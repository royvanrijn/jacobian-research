#!/usr/bin/env sage
"""Exhaust deterministic (pole root, leading coefficient) slices of a Q80 P.O=1 scheme.

The input is the exact msolve system exported by
``certify_q80_fixed_u_marked_third_q12.sage``.  This worker does not rebuild or
alter the surface: it specializes ``z`` and ``l`` over the declared prime,
runs independent bounded msolve jobs, and persists enough data to distinguish
proved empty slices, solution-bearing zero-dimensional slices, and unresolved
timeouts.
"""

from sage.all import GF, PolynomialRing

import argparse
import ast
import concurrent.futures
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--pole-roots",
    default=None,
    help="comma-separated residues (default: every element of the input field)",
)
parser.add_argument(
    "--leading-coefficients",
    default=None,
    help="comma-separated nonzero residues (default: every nonzero element)",
)
parser.add_argument("--workers", type=int, default=4)
parser.add_argument("--msolve-threads", type=int, default=1)
parser.add_argument("--timeout", type=int, default=30)
parser.add_argument("--max-slices", type=int, default=None)
args = parser.parse_args()

if min(args.workers, args.msolve_threads, args.timeout) < 1:
    raise ValueError("workers, msolve threads, and timeout must be positive")
if args.max_slices is not None and args.max_slices < 1:
    raise ValueError("max-slices must be positive")
msolve = shutil.which("msolve")
if msolve is None:
    raise RuntimeError("msolve is unavailable")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def residues(text, default, prime, nonzero):
    values = default if text is None else [int(value) % prime for value in text.split(",")]
    values = sorted(set(values))
    if nonzero and 0 in values:
        raise ValueError("leading coefficient slices must be nonzero")
    return values


lines = [line.strip() for line in args.input.read_text().splitlines() if line.strip()]
names = [name.strip() for name in lines[0].split(",") if name.strip()]
prime = int(lines[1])
if "l" not in names:
    raise ValueError("input is not a denominator-one section scheme")
finite_pole_chart = "z" in names
if not finite_pole_chart and args.pole_roots is not None:
    raise ValueError("--pole-roots is invalid for the pole-at-infinity chart")
field = GF(prime)
ring = PolynomialRing(field, names=names, order="degrevlex")
variables = ring.gens_dict()
equations = [
    ring((line[:-1] if line.endswith(",") else line).replace("^", "**"))
    for line in lines[2:]
]
pole_roots = (
    residues(args.pole_roots, range(prime), prime, False)
    if finite_pole_chart
    else [None]
)
leading_coefficients = residues(
    args.leading_coefficients, range(1, prime), prime, True
)
slices = [(z, leading) for z in pole_roots for leading in leading_coefficients]
if args.max_slices is not None:
    slices = slices[: args.max_slices]

keep_names = [name for name in names if name != "l" and (name != "z" or not finite_pole_chart)]
slice_ring = PolynomialRing(field, names=keep_names, order="degrevlex")
slice_variables = slice_ring.gens_dict()
args.output_dir.mkdir(parents=True, exist_ok=True)


def export_slice(z, leading):
    homomorphism = ring.hom(
        [
            slice_ring(z)
            if name == "z" and finite_pole_chart
            else slice_ring(leading)
            if name == "l"
            else slice_variables[name]
            for name in names
        ],
        slice_ring,
    )
    specialized = [homomorphism(equation) for equation in equations]
    specialized = [equation for equation in specialized if equation]
    stem = (
        f"z{z:02d}-l{leading:02d}"
        if finite_pole_chart
        else f"infinity-l{leading:02d}"
    )
    input_path = args.output_dir / f"{stem}.ms"
    with input_path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(keep_names) + "\n")
        handle.write(str(prime) + "\n")
        for index, equation in enumerate(specialized):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if index + 1 < len(specialized) else "\n")
    return {
        "pole_root": int(z) if finite_pole_chart else None,
        "pole_location": "finite" if finite_pole_chart else "infinity",
        "leading_coefficient": int(leading),
        "input_path": str(input_path),
        "input_sha256": sha256(input_path),
        "variables": len(keep_names),
        "equations": len(specialized),
        "maximum_equation_terms": max(len(equation.dict()) for equation in specialized),
    }


def evaluate_coefficient_list(coefficients, value):
    answer = field.zero()
    for coefficient in reversed(coefficients):
        answer = answer * value + field(coefficient)
    return answer


def decode_degree_one_solution(solution_text, record):
    """Decode and replay the simple degree-one RUR emitted by msolve."""
    text = solution_text.strip()
    parsed = ast.literal_eval(text[:-1] if text.endswith(":") else text)
    if parsed[0] != 0:
        raise ArithmeticError("msolve output is not a zero-dimensional parametrization")
    payload = parsed[1]
    output_prime, variable_count, quotient_degree, output_names = payload[:4]
    separating_vector = payload[4]
    if output_prime != prime or variable_count != len(keep_names) or quotient_degree != 1:
        raise ArithmeticError("msolve parametrization header is not degree one")
    if output_names != keep_names:
        raise ArithmeticError("msolve variable order changed")
    if separating_vector.count(1) != 1 or any(value not in (0, 1) for value in separating_vector):
        raise ArithmeticError("unsupported separating linear form")
    separating_index = separating_vector.index(1)
    components = payload[5]
    if components[0] != 1:
        raise ArithmeticError("degree-one output has multiple parametrization blocks")
    elimination, denominator, coordinate_data = components[1]
    elimination_coefficients = elimination[1]
    if len(elimination_coefficients) != 2 or not elimination_coefficients[1] % prime:
        raise ArithmeticError("degree-one elimination polynomial is malformed")
    root = -field(elimination_coefficients[0]) / field(elimination_coefficients[1])
    if denominator != [0, [1]]:
        raise ArithmeticError("unsupported degree-one parametrization denominator")
    nonseparating_names = [
        name for index, name in enumerate(output_names) if index != separating_index
    ]
    if len(coordinate_data) != len(nonseparating_names):
        raise ArithmeticError("degree-one coordinate count mismatch")
    assignments = {output_names[separating_index]: root}
    for name, coordinate in zip(nonseparating_names, coordinate_data):
        if len(coordinate) != 1 or coordinate[0][0] != 0:
            raise ArithmeticError("unsupported nonconstant degree-one coordinate")
        # msolve stores the coordinate relation as x_i*den + numerator = 0.
        assignments[name] = -evaluate_coefficient_list(coordinate[0][1], root)
    assignments["l"] = field(record["leading_coefficient"])
    if finite_pole_chart:
        assignments["z"] = field(record["pole_root"])
    substitution = {variables[name]: value for name, value in assignments.items()}
    if any(equation.subs(substitution) for equation in equations):
        raise ArithmeticError("decoded degree-one solution fails the original scheme")

    leading = assignments["l"]
    n_coefficients = [assignments[f"n{index}"] for index in range(6)] + [leading**2]
    m_coefficients = [assignments[f"m{index}"] for index in range(9)] + [leading**3]
    if finite_pole_chart:
        z = assignments["z"]
        x_denominator = [z**2, -2 * z, field.one()]
        y_denominator = [-z**3, 3 * z**2, -3 * z, field.one()]
    else:
        x_denominator = [field.one()]
        y_denominator = [field.one()]
    return {
        "assignments": {name: int(assignments[name]) for name in names},
        "separating_variable": output_names[separating_index],
        "separating_value": int(root),
        "scheme_equations_replayed": True,
        "x": {
            "numerator_coefficients_low_to_high": list(map(int, n_coefficients)),
            "denominator_coefficients_low_to_high": list(map(int, x_denominator)),
        },
        "y": {
            "numerator_coefficients_low_to_high": list(map(int, m_coefficients)),
            "denominator_coefficients_low_to_high": list(map(int, y_denominator)),
        },
    }


def solve_slice(record):
    input_path = Path(record["input_path"])
    solution_path = input_path.with_suffix(".solve")
    log_path = input_path.with_suffix(".log")
    command = [
        msolve,
        "-t",
        str(args.msolve_threads),
        "-v",
        "1",
        "-f",
        str(input_path),
        "-o",
        str(solution_path),
    ]
    started = time.monotonic()
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=args.timeout,
            check=False,
        )
        log_path.write_text(completed.stdout + completed.stderr, encoding="utf-8")
        solution_text = (
            solution_path.read_text(encoding="utf-8", errors="replace")
            if solution_path.exists()
            else ""
        )
        if completed.returncode != 0:
            status = "MSOLVE_FAILED"
        elif solution_text.lstrip().startswith("[-1]"):
            status = "PROVED_EMPTY_SLICE"
        elif solution_text.lstrip().startswith("[0"):
            status = "SOLUTION_PARAMETRIZATION"
        else:
            status = "UNCLASSIFIED_MSOLVE_OUTPUT"
        record.update(
            {
                "status": status,
                "msolve_returncode": int(completed.returncode),
                "solution_path": str(solution_path),
                "solution_sha256": sha256(solution_path) if solution_path.exists() else None,
                "solution_size_bytes": solution_path.stat().st_size if solution_path.exists() else 0,
                "log_path": str(log_path),
                "runtime_seconds": time.monotonic() - started,
            }
        )
        if status == "SOLUTION_PARAMETRIZATION":
            record["decoded_degree_one_solution"] = decode_degree_one_solution(
                solution_text, record
            )
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout or ""
        stderr = error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        log_path.write_text(stdout + stderr, encoding="utf-8")
        record.update(
            {
                "status": "MSOLVE_TIMEOUT",
                "log_path": str(log_path),
                "runtime_seconds": time.monotonic() - started,
            }
        )
    return record


started = time.monotonic()
records = [export_slice(z, leading) for z, leading in slices]
print(
    f"Q80PO1SLICE|stage=start|prime={prime}|slices={len(records)}"
    f"|workers={args.workers}|threads_each={args.msolve_threads}",
    flush=True,
)
completed_records = []
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = {executor.submit(solve_slice, record): record for record in records}
    for future in concurrent.futures.as_completed(futures):
        record = future.result()
        completed_records.append(record)
        print(
            "Q80PO1SLICE|z={}|l={}|status={}|seconds={:.3f}".format(
                record["pole_root"] if record["pole_root"] is not None else "infinity",
                record["leading_coefficient"],
                record["status"],
                record["runtime_seconds"],
            ),
            flush=True,
        )

completed_records.sort(
    key=lambda record: (
        -1 if record["pole_root"] is None else record["pole_root"],
        record["leading_coefficient"],
    )
)
counts = {
    status: sum(record["status"] == status for record in completed_records)
    for status in sorted({record["status"] for record in completed_records})
}
candidates = [
    {
        "pole_root": record["pole_root"],
        "leading_coefficient": record["leading_coefficient"],
        "solution_path": record.get("solution_path"),
        "solution_sha256": record.get("solution_sha256"),
        "decoded_degree_one_solution": record.get("decoded_degree_one_solution"),
    }
    for record in completed_records
    if record["status"] == "SOLUTION_PARAMETRIZATION"
]
unresolved = sum(
    record["status"] not in ("PROVED_EMPTY_SLICE", "SOLUTION_PARAMETRIZATION")
    for record in completed_records
)
complete_finite_pole_chart = (
    finite_pole_chart
    and pole_roots == list(range(prime))
    and leading_coefficients == list(range(1, prime))
    and args.max_slices is None
)
complete_infinity_pole_chart = (
    not finite_pole_chart
    and leading_coefficients == list(range(1, prime))
    and args.max_slices is None
)
output = {
    "schema": "elkies-k3.q80-fixed-u-po1-msolve-slices.v1",
    "status": (
        "PASS_COMPLETE_FINITE_POLE_SLICE_SCAN_WITH_CANDIDATES"
        if complete_finite_pole_chart and not unresolved and candidates
        else "PASS_COMPLETE_FINITE_POLE_SLICE_SCAN_EMPTY"
        if complete_finite_pole_chart and not unresolved
        else "PASS_COMPLETE_INFINITY_POLE_SLICE_SCAN_WITH_CANDIDATES"
        if complete_infinity_pole_chart and not unresolved and candidates
        else "PASS_COMPLETE_INFINITY_POLE_SLICE_SCAN_EMPTY"
        if complete_infinity_pole_chart and not unresolved
        else "BOUNDED_FINITE_POLE_SLICE_SCAN_WITH_CANDIDATES"
        if candidates
        else "BOUNDED_FINITE_POLE_SLICE_SCAN_EMPTY"
        if not unresolved
        else "INCOMPLETE_FINITE_POLE_SLICE_SCAN"
    ),
    "input": {"path": str(args.input), "sha256": sha256(args.input)},
    "prime": prime,
    "pole_roots": pole_roots,
    "leading_coefficients": leading_coefficients,
    "declared_slice_count": len(slices),
    "complete_finite_pole_chart": complete_finite_pole_chart,
    "complete_infinity_pole_chart": complete_infinity_pole_chart,
    "pole_location": "finite" if finite_pole_chart else "infinity",
    "solver": {
        "workers": args.workers,
        "threads_each": args.msolve_threads,
        "timeout_seconds_per_slice": args.timeout,
    },
    "counts": counts,
    "candidates": candidates,
    "slices": completed_records,
    "claim_boundary": (
        "A solution parametrization is only a modular marking-generator candidate. "
        "It must be decoded, substituted into the Q80 equation, assigned exact "
        "Shioda/component data, aligned across primes, and reconstructed over QQ. "
        "The finite-pole chart does not include a section whose sole intersection "
        "with O lies over base infinity."
    ),
    "runtime_seconds": time.monotonic() - started,
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80PO1SLICE|stage=done|counts={counts}|candidates={len(candidates)}"
    f"|output={args.output}|status={output['status']}",
    flush=True,
)

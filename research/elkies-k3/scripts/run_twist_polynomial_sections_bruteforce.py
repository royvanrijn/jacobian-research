#!/usr/bin/env python3
"""Compile and run the bounded finite-field polynomial-section enumerator."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
from research_runtime.supervisor import Limits, capture
from research_runtime.section_gate import guard_export
SOURCE = Path(__file__).with_name("bruteforce_twist_polynomial_sections_modp.cpp")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def convolution(left, right, prime):
    result = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] = (result[i + j] + left_value * right_value) % prime
    return result


def matrix_rank(rows, prime):
    rows = [[value % prime for value in row] for row in rows]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(value * inverse) % prime for value in rows[rank]]
        for index in range(len(rows)):
            if index == rank or not rows[index][column]:
                continue
            scale = rows[index][column]
            rows[index] = [
                (value - scale * pivot_value) % prime
                for value, pivot_value in zip(rows[index], rows[rank])
            ]
        rank += 1
    return rank


def full_shell_tangent_rank(solution, coefficient_a, coefficient_b, prime):
    """Rank with the leading X and Y coefficients restored as variables."""

    X = solution["X_coefficients_low_to_high"]
    Y = solution["Y_coefficients_low_to_high"]
    x_degree = len(X) - 1
    y_degree = len(Y) - 1
    if 2 * y_degree != 3 * x_degree:
        raise ArithmeticError("incompatible polynomial-section degree bounds")
    variable_count = x_degree + 2
    residual_degree_count = 2 * y_degree + 1
    x_square = convolution(X, X, prime)
    derivative_rhs = []
    kernel = [0] * max(len(x_square), len(coefficient_a))
    for index in range(len(kernel)):
        kernel[index] = (
            3 * (x_square[index] if index < len(x_square) else 0)
            + (coefficient_a[index] if index < len(coefficient_a) else 0)
        ) % prime
    for variable in range(x_degree + 1):
        derivative_rhs.append([0] * variable + kernel)
    derivative_rhs.append([0] * residual_degree_count)

    derivative_y = [[0] * variable_count for unused in range(y_degree + 1)]
    derivative_y[y_degree][-1] = 1
    denominator = 2 * Y[y_degree] % prime
    inverse_denominator = pow(denominator, -1, prime)
    for degree in range(2 * y_degree - 1, y_degree - 1, -1):
        index = degree - y_degree
        known = sum(
            Y[left] * Y[degree - left]
            for left in range(index + 1, y_degree + 1)
            if index < degree - left <= y_degree
        ) % prime
        # The stored Y already gives the quotient value.  Differentiate its
        # defining recursion, including the variable leading-Y denominator.
        for variable in range(variable_count):
            derivative_known = sum(
                derivative_y[left][variable] * Y[degree - left]
                + Y[left] * derivative_y[degree - left][variable]
                for left in range(index + 1, y_degree + 1)
                if index < degree - left <= y_degree
            ) % prime
            rhs_derivative = (
                derivative_rhs[variable][degree]
                if degree < len(derivative_rhs[variable])
                else 0
            )
            denominator_derivative = 2 if variable == variable_count - 1 else 0
            derivative_y[index][variable] = (
                (rhs_derivative - derivative_known)
                - Y[index] * denominator_derivative
            ) * inverse_denominator % prime

    rows = []
    for degree in [2 * y_degree, *range(y_degree)]:
        row = []
        for variable in range(variable_count):
            y_derivative = sum(
                derivative_y[left][variable] * Y[degree - left]
                + Y[left] * derivative_y[degree - left][variable]
                for left in range(y_degree + 1)
                if 0 <= degree - left <= y_degree
            ) % prime
            rhs_derivative = (
                derivative_rhs[variable][degree]
                if degree < len(derivative_rhs[variable])
                else 0
            )
            row.append((y_derivative - rhs_derivative) % prime)
        rows.append(row)
    return matrix_rank(rows, prime)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--export", type=Path, required=True)
parser.add_argument("--known-section", type=Path)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--wall-seconds", type=float, default=60)
parser.add_argument("--rss-limit-bytes", type=int, default=2_000_000_000)
parser.add_argument("--reduction-only", action="store_true")
args = parser.parse_args()

export_path = args.export.resolve()
export = json.loads(export_path.read_text())
section_gate = guard_export(export, ROOT, reduction_only=args.reduction_only,
    limits={"wall_seconds": args.wall_seconds, "rss_bytes": args.rss_limit_bytes})
prime = int(export["prime"])
chi = int(export["candidate"]["chi"])
systems = export["systems"]
work_dir = export_path.parent
input_path = work_dir / "bruteforce-input.txt"
binary_path = work_dir / "bruteforce-twist-polynomial-sections"

lines = [f"{prime} {chi} {len(systems)}"]
for key in ("twist_A_coefficients_low_to_high", "twist_B_coefficients_low_to_high"):
    values = [int(value) for value in export[key]]
    lines.append(f"{len(values)} " + " ".join(map(str, values)))
for system in systems:
    leading_x, leading_y = system["leading_x_y"]
    lines.append(f"{system['block_index']} {leading_x} {leading_y}")
input_path.write_text("\n".join(lines) + "\n")

capture(["g++", "-O3", "-std=c++17", str(SOURCE), "-o", str(binary_path)],
        limits=Limits(min(args.wall_seconds, 60), args.rss_limit_bytes), log_path=work_dir/"compile.log")
completed = capture([str(binary_path), str(input_path)],
    limits=Limits(args.wall_seconds, args.rss_limit_bytes), log_path=work_dir/"enumeration.log")

solutions = []
summary = None
for line in completed.stdout.splitlines():
    pieces = line.split()
    if pieces[0] == "SOLUTION":
        values = list(map(int, pieces[1:]))
        x_count = 2 * chi + 1
        y_count = 3 * chi + 1
        if len(values) != 2 + x_count + y_count:
            raise ArithmeticError("enumerator emitted a malformed solution row")
        solutions.append(
            {
                "representative_leading_x_y": values[:2],
                "X_coefficients_low_to_high": values[2 : 2 + x_count],
                "Y_coefficients_low_to_high": values[
                    2 + x_count : 2 + x_count + y_count
                ],
            }
        )
    elif pieces[0] == "SUMMARY":
        summary = {
            "x_polynomials_tested": int(pieces[1]),
            "passed_value_sieve": int(pieces[2]),
            "representative_sign_solutions": int(pieces[3]),
            "enumeration_strategy": (
                pieces[4] if len(pieces) > 4 else "coefficient_odometer"
            ),
        }
if summary is None:
    raise ArithmeticError("enumerator did not emit a summary")

coefficient_a = [int(value) for value in export["twist_A_coefficients_low_to_high"]]
coefficient_b = [int(value) for value in export["twist_B_coefficients_low_to_high"]]
for solution in solutions:
    solution["full_shell_tangent_rank"] = full_shell_tangent_rank(
        solution, coefficient_a, coefficient_b, prime
    )

known_match_indices = []
known_path = None
if args.known_section is not None:
    known_path = args.known_section.resolve()
    known = json.loads(known_path.read_text())["modular_identification"]
    known_x = [int(value) for value in known["section_X_coefficients_low_to_high"]]
    known_match_indices = [
        index
        for index, solution in enumerate(solutions)
        if solution["X_coefficients_low_to_high"] == known_x
    ]

payload = {
    "schema": "elkies-k3.twist-polynomial-section-bruteforce.v1",
    "status": "PASS_EXHAUSTIVE_FINITE_FIELD_ENUMERATION_OF_EXPORTED_BLOCKS",
    "prime": prime,
    "candidate": export["candidate"],
    "export_status": export["status"],
    "enumeration": summary,
    "unique_leading_x_values": sorted(
        {int(system["leading_x_y"][0]) for system in systems}
    ),
    "solutions": solutions,
    "full_shell_tangent_rank_histogram": {
        str(rank): sum(
            solution["full_shell_tangent_rank"] == rank for solution in solutions
        )
        for rank in sorted({solution["full_shell_tangent_rank"] for solution in solutions})
    },
    "known_section_match_indices": known_match_indices,
    "proof_boundary": (
        "The enumeration is exhaustive only for the blocks in the input export. "
        "It neither supplies a characteristic-zero lift nor bounds the full "
        "Mordell-Weil rank. The export proof boundary remains in force."
    ),
    "inputs": {
        relative(export_path): digest(export_path),
        relative(SOURCE): digest(SOURCE),
        **(
            {relative(known_path): digest(known_path)}
            if known_path is not None
            else {}
        ),
    },
}
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"TWISTPOLYBRUTE|p={prime}|tested={summary['x_polynomials_tested']}"
    f"|solutions={len(solutions)}|known_matches={len(known_match_indices)}"
    f"|output={relative(output_path)}|status=PASS_EXHAUSTIVE_BLOCKS",
    flush=True,
)

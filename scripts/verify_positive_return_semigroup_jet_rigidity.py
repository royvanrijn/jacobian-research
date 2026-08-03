#!/usr/bin/env python3
"""Exact replays for positive return-semigroup full-jet rigidity.

The unbounded proofs are in
``extended-geometry/POSITIVE_RETURN_SEMIGROUP_JET_RIGIDITY.md``.
This checker verifies bounded lattice-kernel instances without dependencies
and uses Singular, when available, to reconstruct the centered (2,2)
rectangle's exact fourth-jet standard-basis certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
from math import comb, factorial
from pathlib import Path
import shutil
import subprocess


SCRIPT_DIRECTORY = Path(__file__).resolve().parent


def rational_rank(rows: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in rows]
    if not work:
        return 0
    height = len(work)
    width = len(work[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (
                row
                for row in range(pivot_row, height)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(height):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(
                    work[row],
                    work[pivot_row],
                    strict=True,
                )
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def integer_kernel_rank(matrix: list[list[int]]) -> int:
    return len(matrix[0]) - rational_rank(matrix)


def return_vectors(
    columns: list[tuple[int, ...]],
    total: int,
) -> list[tuple[int, ...]]:
    answer = []
    width = len(columns)

    def visit(prefix: list[int], remaining: int) -> None:
        if len(prefix) == width - 1:
            candidate = tuple(prefix + [remaining])
            if all(
                sum(
                    multiplicity * columns[index][row]
                    for index, multiplicity in enumerate(candidate)
                ) == 0
                for row in range(len(columns[0]))
            ):
                answer.append(candidate)
            return
        for multiplicity in range(remaining + 1):
            visit(prefix + [multiplicity], remaining - multiplicity)

    visit([], total)
    return answer


def verify_bounded_group_completion(max_radius: int) -> tuple[int, int]:
    configurations = 0
    returns = 0
    for left in range(1, max_radius + 1):
        for middle in range(1, max_radius + 1):
            for right in range(1, max_radius + 1):
                # The Cartesian product has the explicit strictly positive
                # return obtained from the two one-dimensional means.
                x_columns = (-left, 0, right)
                y_columns = (-middle, 0, middle)
                columns = [
                    (x_value, y_value)
                    for x_value, y_value in product(x_columns, y_columns)
                ]
                # Symmetry kills the y mean.  These x weights have mean zero,
                # including a positive central weight.
                x_weights = {
                    -left: right,
                    0: left + right,
                    right: left,
                }
                positive_return = tuple(
                    x_weights[x_value]
                    for x_value, _ in columns
                )
                assert all(entry > 0 for entry in positive_return)
                assert all(
                    sum(
                        positive_return[index] * columns[index][row]
                        for index in range(len(columns))
                    ) == 0
                    for row in range(2)
                )
                configurations += 1

                sampled = []
                positive_total = sum(positive_return)
                for total in range(1, min(positive_total, 2 * max_radius + 4) + 1):
                    sampled.extend(return_vectors(columns, total))
                sampled.append(positive_return)
                returns += len(sampled)

                # Every sampled return annihilates the two row directions.
                assert all(
                    sum(
                        vector[index] * columns[index][row]
                        for index in range(len(columns))
                    ) == 0
                    for vector in sampled
                    for row in range(2)
                )

                # Adding/subtracting the positive return realizes sampled
                # lattice differences without changing the A-image.
                differences = []
                for vector in sampled:
                    difference = [
                        vector[index] - positive_return[index]
                        for index in range(len(columns))
                    ]
                    differences.append(difference)
                    assert all(
                        sum(
                            difference[index] * columns[index][row]
                            for index in range(len(columns))
                        ) == 0
                        for row in range(2)
                    )
                kernel_rank = integer_kernel_rank(
                    [[column[row] for column in columns] for row in range(2)]
                )
                assert kernel_rank == len(columns) - 2
                assert rational_rank(differences) == kernel_rank
    return configurations, returns


VARIABLE_CELLS = (
    (0, 0),
    (0, 2),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
)


def centered_square_returns(total: int) -> list[tuple[int, ...]]:
    cells = [(x, y) for x in range(3) for y in range(3)]
    answer = []

    def visit(prefix: list[int], remaining: int) -> None:
        if len(prefix) == 8:
            candidate = tuple(prefix + [remaining])
            if (
                sum(candidate[index] * cells[index][0] for index in range(9))
                == total
                and sum(
                    candidate[index] * cells[index][1]
                    for index in range(9)
                )
                == total
            ):
                answer.append(candidate)
            return
        for multiplicity in range(remaining + 1):
            visit(prefix + [multiplicity], remaining - multiplicity)

    visit([], total)
    return answer


def derivative_expression(total: int, order: int) -> str:
    cells = [(x, y) for x in range(3) for y in range(3)]
    terms = []
    for vector in centered_square_returns(total):
        weight = factorial(total)
        labels = []
        for index, multiplicity in enumerate(vector):
            weight //= factorial(multiplicity)
            x_value, y_value = cells[index]
            weight *= (
                comb(2, x_value) * comb(2, y_value)
            ) ** multiplicity
        for variable, cell in zip("abcdefg", VARIABLE_CELLS, strict=True):
            multiplicity = vector[cells.index(cell)]
            if multiplicity:
                labels.append(
                    variable if multiplicity == 1 else f"{multiplicity}*{variable}"
                )
        linear_form = "+".join(labels) if labels else "0"
        if linear_form != "0":
            terms.append(f"{weight}*({linear_form})^{order}")
    return "+".join(terms) if terms else "0"


EXPECTED_STAGES = {
    "LINEAR": (4, -1, 3),
    "QUADRATIC": (2, -1, 6),
    "CUBIC": (1, -1, 9),
    "QUARTIC": (0, 40, 11),
}

REQUIRED_BASIS_MEMBERS = (
    "d+f",
    "c",
    "a+b+3d+e+3f+g",
    "g4",
    "4f4+2f2g2-bg3-eg3",
    "e3-3e2f-4bf2-ef2-f3+e2g-2efg-7f2g+2bg2+3eg2-fg2+g3",
    "2b2+3be+2e2+bf-ef+5f2+bg+eg+g2",
)


def singular_program() -> str:
    equations_by_order = {
        order: [
            derivative_expression(total, order)
            for total in range(1, 4)
        ]
        for order in range(1, 5)
    }
    lines = ["ring R=0,(a,b,c,d,e,f,g),dp;"]
    accumulated = []
    for order, name in enumerate(EXPECTED_STAGES, start=1):
        accumulated.extend(equations_by_order[order])
        generators = ",".join(accumulated)
        lines.extend(
            [
                f"ideal I{order}={generators};",
                f"ideal G{order}=std(I{order});",
                f'print("@@{name}");',
                f"print(dim(G{order}));",
                f"print(vdim(G{order}));",
                f"print(size(G{order}));",
            ]
        )
    lines.extend(['print("@@BASIS");', "G4;"])
    return "\n".join(lines) + "\n"


def compact_polynomial(value: str) -> str:
    polynomial = value.split("=", 1)[-1]
    return polynomial.replace(" ", "").replace("*", "").replace("^", "")


def verify_centered_square_certificate(require_singular: bool) -> bool:
    singular = shutil.which("Singular")
    if singular is None:
        if require_singular:
            raise SystemExit("Singular is required but was not found")
        print("SKIP centered-square fourth-jet certificate: Singular not found")
        return False
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_program(),
        text=True,
        capture_output=True,
        check=True,
        cwd=SCRIPT_DIRECTORY.parent,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    for name, expected in EXPECTED_STAGES.items():
        marker = f"@@{name}"
        index = lines.index(marker)
        observed = tuple(int(lines[index + offset]) for offset in range(1, 4))
        assert observed == expected, (name, observed, expected)
    basis_index = lines.index("@@BASIS")
    basis = {
        compact_polynomial(line.rstrip(","))
        for line in lines[basis_index + 1 :]
        if not line.startswith("//")
    }
    for required in REQUIRED_BASIS_MEMBERS:
        assert required in basis, (required, sorted(basis))
    print(
        "PASS centered-square exact jet ladder: dimensions 4,2,1,0; "
        "quartic quotient dimension 40 and standard-basis size 11"
    )
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-radius", type=int, default=3)
    parser.add_argument("--require-singular", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.max_radius < 1:
        raise SystemExit("max-radius must be positive")
    configurations, returns = verify_bounded_group_completion(arguments.max_radius)
    print(
        "PASS positive-return group-completion and polarized paired-point "
        "regressions: "
        f"{configurations} Cartesian configurations and {returns} sampled returns"
    )
    verify_centered_square_certificate(arguments.require_singular)
    print(
        "STATUS: identity twists and independently marked paired points "
        "proved rigid; common-quotient Hall marking remains unproved in "
        "the parked route and is bypassed by Hall-envelope separation"
    )


if __name__ == "__main__":
    main()

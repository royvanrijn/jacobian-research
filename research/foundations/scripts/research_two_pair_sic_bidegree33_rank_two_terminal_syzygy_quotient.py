#!/usr/bin/env python3
"""Test the terminal syzygy in the 14-dimensional critical quotient.

The expanded triangular reduction of the terminal target has severe fill-in.
The first invariant obstruction is much smaller: reduce

    T/Q  and  H(R)=C(D_u+3)R-A R_t

in the interior critical algebra B=F_p[u^{+-1},t]/(A,C).  This script
constructs the exact 14-by-14 multiplication and H matrices, checks that Q
is invertible on B, and decides whether the terminal class lies in im(H).

Membership is necessary for a last-step polynomial Koszul correction.
Failure is therefore an exact modular no-go for that repair; success only
removes the quotient obstruction and still requires a polynomial lift.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from flint import nmod_mat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_two_pair_sic_bidegree33_rank_two_relative_divergence import (  # noqa: E402
    OPERATOR_INPUT,
)
from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


TERMINAL = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m0.poly"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_quotient_research.json"
)


def singular_code(prime: int, terminal: str) -> str:
    return f"""
ring ambient={prime},(u,U,t),dp;
ideal relation=std(ideal(u*U-1));
qring r=relation;
poly Q={q_expression(0)};
poly A=u*diff(Q,u)-3*Q;
poly C=diff(Q,t);
ideal critical=std(ideal(A,C));
if(vdim(critical)!=14){{ERROR("interior quotient length");}}
ideal basis=kbase(critical);
if(size(basis)!=14){{ERROR("interior basis size");}}

proc dump_poly(poly f)
{{
  intvec exponent;
  f=reduce(f,critical);
  print("POLY_BEGIN");
  while(f!=0)
  {{
    exponent=leadexp(f);
    print(leadcoef(f));
    print(exponent[1]);
    print(exponent[2]);
    print(exponent[3]);
    f=f-lead(f);
  }}
  print("POLY_END");
}}

int index;
poly R;
poly H;
for(index=1;index<=14;index++)
{{
  dump_poly(basis[index]);
}}
poly terminal={terminal};
dump_poly(terminal);
for(index=1;index<=14;index++)
{{
  dump_poly(Q*basis[index]);
}}
for(index=1;index<=14;index++)
{{
  R=basis[index];
  H=C*(u*diff(R,u)-U*diff(R,U)+3*R)-A*diff(R,t);
  dump_poly(H);
}}
print("PASS terminal syzygy quotient data");
"""


def parse_polynomials(
    output: str,
    prime: int,
) -> list[dict[tuple[int, int, int], int]]:
    blocks = re.findall(
        r"POLY_BEGIN\s*(.*?)\s*POLY_END",
        output,
        flags=re.DOTALL,
    )
    polynomials = []
    for block in blocks:
        values = [int(value) for value in block.split()]
        if len(values) % 4:
            raise RuntimeError(block)
        polynomial = {}
        for offset in range(0, len(values), 4):
            coefficient, u_degree, inverse_degree, t_degree = values[
                offset : offset + 4
            ]
            exponent = (u_degree, inverse_degree, t_degree)
            polynomial[exponent] = coefficient % prime
        polynomials.append(polynomial)
    return polynomials


def vector(
    polynomial: dict[tuple[int, int, int], int],
    basis_index: dict[tuple[int, int, int], int],
) -> list[int]:
    answer = [0] * len(basis_index)
    for exponent, coefficient in polynomial.items():
        if exponent not in basis_index:
            raise RuntimeError(f"unexpected quotient monomial {exponent}")
        answer[basis_index[exponent]] = coefficient
    return answer


def solve_modular(
    rows: list[list[int]],
    target: list[int],
    prime: int,
) -> list[int] | None:
    augmented = [
        [value % prime for value in row] + [right % prime]
        for row, right in zip(rows, target, strict=True)
    ]
    row_count = len(augmented)
    column_count = len(rows[0])
    pivot_row = 0
    pivots = []
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        augmented[pivot_row], augmented[pivot] = (
            augmented[pivot],
            augmented[pivot_row],
        )
        inverse = pow(augmented[pivot_row][column], -1, prime)
        augmented[pivot_row] = [
            value * inverse % prime for value in augmented[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            scalar = augmented[row][column]
            if scalar:
                augmented[row] = [
                    (left - scalar * right) % prime
                    for left, right in zip(
                        augmented[row],
                        augmented[pivot_row],
                        strict=True,
                    )
                ]
        pivots.append(column)
        pivot_row += 1
    for row in range(pivot_row, row_count):
        if not any(augmented[row][:-1]) and augmented[row][-1]:
            return None
    solution = [0] * column_count
    for row, column in enumerate(pivots):
        solution[column] = augmented[row][-1]
    return solution


def centered_expression(
    coefficients: list[int],
    basis: list[tuple[int, int, int]],
    prime: int,
) -> str:
    terms = []
    for coefficient, (u_degree, inverse_degree, t_degree) in zip(
        coefficients,
        basis,
        strict=True,
    ):
        if not coefficient:
            continue
        centered = (
            coefficient
            if coefficient <= prime // 2
            else coefficient - prime
        )
        factors = [str(centered)]
        if u_degree:
            factors.append(
                "u" if u_degree == 1 else f"u^{u_degree}"
            )
        if inverse_degree:
            factors.append(
                "U" if inverse_degree == 1 else f"U^{inverse_degree}"
            )
        if t_degree:
            factors.append(
                "t" if t_degree == 1 else f"t^{t_degree}"
            )
        terms.append("*".join(factors))
    return "+".join(terms).replace("+-", "-") if terms else "0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument("--terminal", type=Path, default=TERMINAL)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    operator = json.loads(arguments.operator.read_text())
    prime = int(operator["modular_operator"]["prime"])
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(
            prime,
            arguments.terminal.read_text().strip(),
        ),
        text=True,
        capture_output=True,
        timeout=arguments.timeout,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "?" in combined or "error occurred" in combined:
        raise RuntimeError(combined)
    if "PASS terminal syzygy quotient data" not in combined:
        raise RuntimeError(combined)
    polynomials = parse_polynomials(combined, prime)
    if len(polynomials) != 43:
        raise RuntimeError(f"expected 43 polynomials, got {len(polynomials)}")
    basis_polynomials = polynomials[:14]
    basis = []
    for polynomial in basis_polynomials:
        if len(polynomial) != 1:
            raise RuntimeError("kbase element is not a monomial")
        exponent, coefficient = next(iter(polynomial.items()))
        if coefficient != 1:
            raise RuntimeError("kbase element is not monic")
        basis.append(exponent)
    basis_index = {
        exponent: index for index, exponent in enumerate(basis)
    }
    if len(basis_index) != 14:
        raise RuntimeError("duplicate quotient basis monomial")
    terminal_vector = vector(polynomials[14], basis_index)
    q_images = [
        vector(polynomial, basis_index)
        for polynomial in polynomials[15:29]
    ]
    h_images = [
        vector(polynomial, basis_index)
        for polynomial in polynomials[29:43]
    ]
    q_rows = [
        [q_images[column][row] for column in range(14)]
        for row in range(14)
    ]
    h_rows = [
        [h_images[column][row] for column in range(14)]
        for row in range(14)
    ]
    q_matrix = nmod_mat(q_rows, prime)
    h_matrix = nmod_mat(h_rows, prime)
    q_rank = q_matrix.rank()
    if q_rank != 14:
        raise RuntimeError(f"Q multiplication rank is only {q_rank}")
    divided_target = solve_modular(q_rows, terminal_vector, prime)
    if divided_target is None:
        raise RuntimeError("Q multiplication unexpectedly not invertible")
    solution = solve_modular(h_rows, divided_target, prime)
    augmented_rows = [
        row + [right]
        for row, right in zip(h_rows, divided_target, strict=True)
    ]
    augmented_rank = nmod_mat(augmented_rows, prime).rank()
    h_rank = h_matrix.rank()
    in_image = solution is not None
    if in_image != (augmented_rank == h_rank):
        raise RuntimeError("rank and solve tests disagree")

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "terminal-syzygy-quotient-research-v1"
        ),
        "status": (
            "exact modular quotient solution for the final-syzygy repair"
            if in_image
            else (
                "exact modular quotient obstruction to the "
                "final-syzygy repair"
            )
        ),
        "prime": prime,
        "point": operator["point"],
        "critical_length": 14,
        "basis_exponents_u_U_t": [list(exponent) for exponent in basis],
        "Q_multiplication_rank": q_rank,
        "H_rank": h_rank,
        "H_augmented_rank": augmented_rank,
        "terminal_class_nonzero": any(divided_target),
        "terminal_class_in_H_image": in_image,
        "quotient_solution": (
            centered_expression(solution, basis, prime)
            if solution is not None
            else None
        ),
        "conclusion": (
            (
                "the terminal class has a quotient-level Koszul correction; "
                "lift its ideal remainder before claiming a polynomial "
                "certificate"
            )
            if in_image
            else (
                "no last-degree polynomial Koszul syzygy can remove the "
                "terminal residual; higher-degree syzygy freedom or the "
                "endpoint-extended connection is necessary"
            )
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS multiplication by Q has rank {q_rank}")
    print(f"PASS terminal syzygy map has rank {h_rank}")
    print(
        "PASS terminal class is "
        f"{'in' if in_image else 'not in'} the syzygy-map image"
    )
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()

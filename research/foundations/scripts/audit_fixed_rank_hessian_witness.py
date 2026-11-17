#!/usr/bin/env python3
"""Audit the generic Hessian rank of the 42-variable quartic HN witness.

The expanded quartic is unnecessary.  If ``V=x+H`` is the 21-variable
cubic-homogeneous source, its cotangent potential is

    p(x,y) = y.H(x),

and, before the orthogonal change to the standard Laplacian, its Hessian is

                 [ A(x,y)  JH(x)^t ]
    Hess(p)  =   [                  ],
                 [ JH(x)       0   ]

where A=sum_i y_i Hess(H_i).  Congruence by the invertible orthogonal-change
matrix preserves generic rank.  The script proves rank 38 in two ways:

* deterministic good-prime specializations give rank at least 38; and
* Singular computes four generically independent polynomial syzygies, giving
  rank at most 42-4=38 over Q(x,y).

It also records the rank ladder through the existing 24/22/21-variable cubic
artifacts and a finite-field nilpotency-index specialization after the
orthogonal change.  The latter is a diagnostic, not the all-points HN proof.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import random
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = ROOT / "artifacts" / "generated-results"
SOURCE = ARTIFACT_DIR / "essential_bcw_21_counterexample.json"
PRIME = 1_000_003
COMPLEX_PRIME = 1_000_033  # 1 mod 4; 350504^2 = -1.
SQRT_MINUS_ONE = 350_504

Exponent = tuple[int, ...]
SparsePoly = dict[Exponent, Fraction]
SparseMatrix = list[list[SparsePoly]]


def add_term(poly: SparsePoly, exponent: Exponent, coefficient: Fraction) -> None:
    value = poly.get(exponent, Fraction(0)) + coefficient
    if value:
        poly[exponent] = value
    else:
        poly.pop(exponent, None)


def derivative(poly: SparsePoly, variable: int) -> SparsePoly:
    answer: SparsePoly = {}
    for exponent, coefficient in poly.items():
        power = exponent[variable]
        if power:
            reduced = list(exponent)
            reduced[variable] -= 1
            add_term(answer, tuple(reduced), coefficient * power)
    return answer


def decode_h(stored: dict[str, object]) -> list[SparsePoly]:
    dimension = int(stored["dimension"])
    answer: list[SparsePoly] = []
    for component in stored["H"]:
        poly: SparsePoly = {}
        for term in component:
            exponent = [0] * dimension
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_term(poly, tuple(exponent), Fraction(term["coefficient"]))
        answer.append(poly)
    return answer


def cotangent_blocks(h: list[SparsePoly]) -> tuple[SparseMatrix, SparseMatrix]:
    """Return JH and sum y_i Hess(H_i), both in 2n variables."""

    dimension = len(h)
    zero = (0,) * dimension
    jacobian_small = [
        [derivative(h[output], variable) for variable in range(dimension)]
        for output in range(dimension)
    ]
    jacobian = [
        [
            {exponent + zero: coefficient for exponent, coefficient in entry.items()}
            for entry in row
        ]
        for row in jacobian_small
    ]
    upper_left: SparseMatrix = []
    for first in range(dimension):
        row: list[SparsePoly] = []
        for second in range(dimension):
            entry: SparsePoly = {}
            for output in range(dimension):
                for exponent, coefficient in derivative(
                    jacobian_small[output][first], second
                ).items():
                    lifted = list(exponent) + [0] * dimension
                    lifted[dimension + output] = 1
                    add_term(entry, tuple(lifted), coefficient)
            row.append(entry)
        upper_left.append(row)
    return jacobian, upper_left


def cotangent_hessian(h: list[SparsePoly]) -> tuple[SparseMatrix, SparseMatrix, SparseMatrix]:
    jacobian, upper_left = cotangent_blocks(h)
    dimension = len(h)
    zero: SparsePoly = {}
    matrix = [
        upper_left[row] + [jacobian[column][row] for column in range(dimension)]
        for row in range(dimension)
    ]
    matrix.extend(jacobian[row] + [zero] * dimension for row in range(dimension))
    return jacobian, upper_left, matrix


def residue(value: Fraction, prime: int) -> int:
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def evaluate(poly: SparsePoly, point: list[int], prime: int) -> int:
    answer = 0
    for exponent, coefficient in poly.items():
        value = residue(coefficient, prime)
        for variable, power in zip(point, exponent):
            value = value * pow(variable, power, prime) % prime
        answer = (answer + value) % prime
    return answer


def evaluate_matrix(matrix: SparseMatrix, point: list[int], prime: int) -> list[list[int]]:
    return [[evaluate(entry, point, prime) for entry in row] for row in matrix]


def modular_rank(matrix: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [value * inverse % prime for value in work[pivot_row]]
        for row in range(pivot_row + 1, len(work)):
            if work[row][column]:
                factor = work[row][column]
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def multiply(left: list[list[int]], right: list[list[int]], prime: int) -> list[list[int]]:
    columns = list(zip(*right))
    return [
        [sum(a * b for a, b in zip(row, column)) % prime for column in columns]
        for row in left
    ]


def deterministic_point(dimension: int, prime: int, seed: int) -> list[int]:
    generator = random.Random(seed)
    return [generator.randrange(1, prime) for _ in range(dimension)]


def specialization_profile(
    jacobian: SparseMatrix,
    upper_left: SparseMatrix,
    hessian: SparseMatrix,
    seed: int,
    prime: int = PRIME,
) -> tuple[int, int, int]:
    point = deterministic_point(len(hessian), prime, seed)
    return (
        modular_rank(evaluate_matrix(jacobian, point, prime), prime),
        modular_rank(evaluate_matrix(upper_left, point, prime), prime),
        modular_rank(evaluate_matrix(hessian, point, prime), prime),
    )


def constant_kernel_coefficient_rank(matrix: SparseMatrix, prime: int = PRIME) -> int:
    """Rank of all coefficient equations for a constant kernel vector."""

    monomials = sorted({exponent for row in matrix for entry in row for exponent in entry})
    equations = [
        [residue(matrix[row][column].get(monomial, Fraction(0)), prime) for column in range(len(row_values))]
        for row, row_values in enumerate(matrix)
        for monomial in monomials
    ]
    return modular_rank(equations, prime)


def singular_polynomial(poly: SparsePoly) -> str:
    terms: list[str] = []
    for exponent, coefficient in sorted(poly.items(), reverse=True):
        if not coefficient:
            continue
        scalar = (
            str(coefficient.numerator)
            if coefficient.denominator == 1
            else f"({coefficient.numerator}/{coefficient.denominator})"
        )
        factors = [scalar]
        factors.extend(
            f"x{index}" + (f"^{power}" if power > 1 else "")
            for index, power in enumerate(exponent)
            if power
        )
        terms.append("*".join(factors))
    return "+".join(terms).replace("+-", "-") or "0"


def exact_singular_certificate(
    matrix: SparseMatrix,
    point: list[int],
    *,
    expected_kernel_rank: int = 4,
) -> tuple[int, int, bool]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the exact rank-upper-bound certificate")
    dimension = len(matrix)
    variables = ",".join(f"x{index}" for index in range(dimension))
    entries = ",".join(singular_polynomial(entry) for row in matrix for entry in row)
    commands = [
        f"ring r=0,({variables}),dp;",
        f"matrix M[{dimension}][{dimension}]={entries};",
        "option(redSB);",
        "module K=syz(M);",
        'print("SYZYGY_COLUMNS");',
        "ncols(K);",
        "matrix Z=M*matrix(K);",
        "int product_zero=1;",
        "int ii,jj;",
        f"for (ii=1; ii<={dimension}; ii++) {{ for (jj=1; jj<=ncols(K); jj++) {{ if (Z[ii,jj]!=0) {{ product_zero=0; }} }} }}",
        'if (product_zero==1) { print("PRODUCT_ZERO"); } else { print("PRODUCT_NONZERO"); }',
        "matrix Ms=M;",
        "module Ks=K;",
    ]
    for index, value in enumerate(point):
        commands.extend(
            [f"Ms=subst(Ms,x{index},{value});", f"Ks=subst(Ks,x{index},{value});"]
        )
    commands.extend(
        [
            'print("SPECIALIZED_RANKS");',
            "rank(Ms);",
            "rank(Ks);",
            "quit;",
        ]
    )
    result = subprocess.run(
        [executable, "-q"],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    syzygy_marker = lines.index("SYZYGY_COLUMNS")
    ranks_marker = lines.index("SPECIALIZED_RANKS")
    return (
        int(lines[syzygy_marker + 1]),
        int(lines[ranks_marker + 1]),
        "PRODUCT_ZERO" in lines
        and int(lines[ranks_marker + 2]) == expected_kernel_rank,
    )


def transformed_hessian(
    jacobian_value: list[list[int]], upper_left_value: list[list[int]], prime: int
) -> list[list[int]]:
    """Congruence to (u,v), omitting one irrelevant nonzero scalar."""

    dimension = len(jacobian_value)
    imaginary = SQRT_MINUS_ONE
    assert prime == COMPLEX_PRIME and imaginary * imaginary % prime == prime - 1
    transpose = [list(column) for column in zip(*jacobian_value)]
    top_left = [
        [
            (upper_left_value[i][j] + transpose[i][j] + jacobian_value[i][j]) % prime
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    top_right = [
        [
            imaginary
            * (upper_left_value[i][j] - transpose[i][j] + jacobian_value[i][j])
            % prime
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    bottom_left = [
        [
            imaginary
            * (upper_left_value[i][j] + transpose[i][j] - jacobian_value[i][j])
            % prime
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    bottom_right = [
        [
            (-upper_left_value[i][j] + transpose[i][j] + jacobian_value[i][j]) % prime
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    return [
        top_left[row] + top_right[row] for row in range(dimension)
    ] + [bottom_left[row] + bottom_right[row] for row in range(dimension)]


def nilpotency_profile(
    jacobian: SparseMatrix, upper_left: SparseMatrix, seed: int
) -> list[int]:
    point = deterministic_point(2 * len(jacobian), COMPLEX_PRIME, seed)
    j_value = evaluate_matrix(jacobian, point, COMPLEX_PRIME)
    a_value = evaluate_matrix(upper_left, point, COMPLEX_PRIME)
    matrix = transformed_hessian(j_value, a_value, COMPLEX_PRIME)
    ranks: list[int] = []
    power = matrix
    for _ in range(42):
        rank = modular_rank(power, COMPLEX_PRIME)
        ranks.append(rank)
        if rank == 0:
            break
        power = multiply(power, matrix, COMPLEX_PRIME)
    return ranks


def artifact_profile(filename: str) -> tuple[int, int, int]:
    stored = json.loads((ARTIFACT_DIR / filename).read_text())
    jacobian, upper_left, hessian = cotangent_hessian(decode_h(stored))
    return specialization_profile(jacobian, upper_left, hessian, 20_260_723)


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    assert stored["dimension"] == 21 and stored["H"][20] == []
    jacobian, upper_left, hessian = cotangent_hessian(decode_h(stored))

    profiles = [
        specialization_profile(jacobian, upper_left, hessian, 20_260_723 + offset)
        for offset in range(3)
    ]
    assert profiles == [(18, 15, 38)] * 3

    exact_point = deterministic_point(42, PRIME, 20_260_723)
    syzygies, specialized_rank, exact_kernel_check = exact_singular_certificate(
        hessian, exact_point
    )
    assert syzygies == 4 and specialized_rank == 38 and exact_kernel_check

    # Column 41 is zero because H_20=0, and the coefficient equations have
    # rank 41.  Thus the polynomial Hessian has exactly one constant kernel
    # direction even though its generic kernel has dimension four.
    assert all(not row[41] for row in hessian)
    assert constant_kernel_coefficient_rank(hessian) == 41

    artifact_profiles = {
        filename: artifact_profile(filename)
        for filename in (
            "rank_compressed_bcw_24_counterexample.json",
            "constant_kernel_bcw_22_counterexample.json",
            "essential_bcw_21_counterexample.json",
        )
    }
    assert set(artifact_profiles.values()) == {(18, 15, 38)}

    power_ranks = nilpotency_profile(jacobian, upper_left, 20_260_723)
    assert power_ranks[-2:] == [1, 0] and len(power_ranks) == 35

    print("PASS fixed-rank HN: three deterministic specializations have ranks JH=18, A=15, Hess=38")
    print("PASS fixed-rank HN: four independent Singular syzygies certify generic Hessian rank 38")
    print("PASS fixed-rank HN: the 24D, 22D, and 21D cubic sources all retain cotangent-Hessian rank 38")
    print("PASS fixed-rank HN: generic Hessian corank is 4 but constant kernel dimension is only 1")
    print("PASS fixed-rank HN: transformed finite-field specialization has nilpotency index 35")


if __name__ == "__main__":
    main()

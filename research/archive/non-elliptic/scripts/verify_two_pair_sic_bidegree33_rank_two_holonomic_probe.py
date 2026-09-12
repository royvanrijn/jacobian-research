#!/usr/bin/env python3
"""Exact setup and modular recurrence probe for rank-two bidegree (3,3).

The exact part verifies the beta/constant-term period formula and its
rank-two factorization.  The modular part compiles the adjacent C++ helper,
computes 501 normalized moments at two exact rank-two factor points over
three prime fields, and fits an order-27, degree-11 scalar recurrence with
139 unused equations in every case.

This is deliberately not called a creative-telescoping certificate.  The
fitted recurrence coefficients have not been reconstructed in the universal
rank-two parameter ring.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "scripts"
    / "explore_two_pair_sic_bidegree33_rank_two_recurrence.cpp"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_holonomic_probe.json"
)
PRIMES = (1_000_003, 1_000_033, 1_000_037)
POINTS = (
    (
        ((1, 0), (0, 1), (2, 3), (5, 7)),
        ((11, 13, 17, 19), (23, 29, 31, 37)),
    ),
    (
        ((1, 2), (3, 5), (7, 11), (13, 17)),
        ((19, 23, 29, 31), (37, 41, 43, 47)),
    ),
)
LINEAR_SHIFTS = (71, 73, 74, 76, 77, 79, 80, 82)
RECURRENCE_ORDER = 27
RECURRENCE_DEGREE = 11
MAXIMUM_MOMENT = 500

Polynomial = dict[tuple[int, int], Fraction]
Matrix = list[list[Fraction]]


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for (i, j), left_value in left.items():
        for (a, b), right_value in right.items():
            exponent = (i + a, j + b)
            answer[exponent] = (
                answer.get(exponent, Fraction(0))
                + left_value * right_value
            )
    return {exponent: value for exponent, value in answer.items() if value}


def power(polynomial: Polynomial, order: int) -> Polynomial:
    answer = {(0, 0): Fraction(1)}
    for _ in range(order):
        answer = multiply(answer, polynomial)
    return answer


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(
                (
                    left[row][inner] * right[inner][column]
                    for inner in range(2)
                ),
                Fraction(0),
            )
            for column in range(4)
        ]
        for row in range(4)
    ]


def matrix_rank(matrix: Matrix) -> int:
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                index
                for index in range(row, len(work))
                if work[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_value = work[row][column]
        work[row] = [value / pivot_value for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            scale = work[index][column]
            work[index] = [
                left - scale * right
                for left, right in zip(
                    work[index],
                    work[row],
                    strict=True,
                )
            ]
        row += 1
    return row


def substituted_polynomial(matrix: Matrix) -> Polynomial:
    """Return Phi_C(1,u,t,(1-t)/u), with exponents (u,t)."""

    answer: Polynomial = {}
    for i, j in product(range(4), repeat=2):
        for extra in range(4 - j):
            exponent = (j - i, j + extra)
            answer[exponent] = (
                answer.get(exponent, Fraction(0))
                + matrix[i][j]
                * (-1) ** extra
                * comb(3 - j, extra)
            )
    return {exponent: value for exponent, value in answer.items() if value}


def factor_substituted_polynomial(u: Matrix, w: Matrix) -> Polynomial:
    answer: Polynomial = {}
    for inner in range(2):
        dual = {
            (3 - i, 0): u[i][inner]
            for i in range(4)
            if u[i][inner]
        }
        coordinate: Polynomial = {}
        for j in range(4):
            for extra in range(4 - j):
                exponent = (j - 3, j + extra)
                coordinate[exponent] = (
                    coordinate.get(exponent, Fraction(0))
                    + w[inner][j]
                    * (-1) ** extra
                    * comb(3 - j, extra)
                )
        for exponent, value in multiply(dual, coordinate).items():
            answer[exponent] = answer.get(exponent, Fraction(0)) + value
    return {exponent: value for exponent, value in answer.items() if value}


def constant_term_integral(polynomial: Polynomial) -> Fraction:
    return sum(
        (
            coefficient / Fraction(t_degree + 1)
            for (u_degree, t_degree), coefficient in polynomial.items()
            if u_degree == 0
        ),
        Fraction(0),
    )


def pure_moment(matrix: Matrix, order: int) -> Fraction:
    coefficient_polynomial = {
        (i, j): matrix[i][j]
        for i, j in product(range(4), repeat=2)
        if matrix[i][j]
    }
    polynomial_power = power(coefficient_polynomial, order)
    degree = 3 * order
    return sum(
        (
            factorial(index)
            * factorial(degree - index)
            * polynomial_power.get((index, index), Fraction(0))
            for index in range(degree + 1)
        ),
        Fraction(0),
    )


def polynomial_multiply(
    left: list[Fraction],
    right: list[Fraction],
) -> list[Fraction]:
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return answer


def expected_monic_leading_polynomial() -> list[Fraction]:
    answer = [Fraction(1)]
    for shift in LINEAR_SHIFTS:
        answer = polynomial_multiply(
            answer,
            [Fraction(shift), Fraction(3)],
        )
    answer = polynomial_multiply(
        answer,
        [
            Fraction(-5_973_460, 3),
            Fraction(210_730, 9),
            Fraction(-204),
            Fraction(1),
        ],
    )
    leading = answer[-1]
    return [coefficient / leading for coefficient in answer]


def reduce_fraction(value: Fraction, prime: int) -> int:
    return (
        value.numerator
        * pow(value.denominator % prime, -1, prime)
        % prime
    )


def parse_probe(output: str) -> list[int]:
    lines = output.strip().splitlines()
    assert lines[0] == "FOUND order=27 degree=11"
    leading = lines[-1].split()
    assert leading[0] == "27"
    coefficients = [int(value) for value in leading[1:]]
    assert len(coefficients) == RECURRENCE_DEGREE + 1
    return coefficients


def run_probe(executable: Path, prime: int, point: int, degree: int) -> str:
    completed = subprocess.run(
        [
            str(executable),
            str(prime),
            str(MAXIMUM_MOMENT),
            str(RECURRENCE_ORDER),
            str(degree),
            str(point),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    return completed.stdout.strip()


def main() -> None:
    exact_records = []
    for raw_u, raw_w in POINTS:
        u = [[Fraction(value) for value in row] for row in raw_u]
        w = [[Fraction(value) for value in row] for row in raw_w]
        matrix = matrix_product(u, w)
        assert matrix_rank(matrix) == 2
        substituted = substituted_polynomial(matrix)
        assert substituted == factor_substituted_polynomial(u, w)
        moments = []
        for order in range(1, 5):
            moment = pure_moment(matrix, order)
            period = constant_term_integral(power(substituted, order))
            assert moment == factorial(3 * order + 1) * period
            moments.append(int(moment))
        exact_records.append(
            {
                "U": raw_u,
                "W": raw_w,
                "rank": 2,
                "pure_moments_1_through_4": moments,
            }
        )

    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required for the modular recurrence probe")
    expected = expected_monic_leading_polynomial()
    modular_records = []
    with tempfile.TemporaryDirectory(prefix="sic33-holonomic-probe-") as path:
        executable = Path(path) / "recurrence-probe"
        subprocess.run(
            [
                compiler,
                "-O3",
                "-std=c++17",
                str(SOURCE),
                "-o",
                str(executable),
            ],
            check=True,
            timeout=30,
        )
        for point, prime in product(range(len(POINTS)), PRIMES):
            output = run_probe(
                executable,
                prime,
                point,
                RECURRENCE_DEGREE,
            )
            leading = parse_probe(output)
            expected_modular = [
                reduce_fraction(value, prime) for value in expected
            ]
            assert leading == expected_modular
            modular_records.append(
                {
                    "point": point,
                    "prime": prime,
                    "fit": {
                        "order": RECURRENCE_ORDER,
                        "coefficient_degree": RECURRENCE_DEGREE,
                        "moments_computed": MAXIMUM_MOMENT + 1,
                        "linear_fit_equations": (
                            (RECURRENCE_ORDER + 1)
                            * (RECURRENCE_DEGREE + 1)
                            - 1
                        ),
                        "unused_verification_equations": (
                            MAXIMUM_MOMENT
                            - RECURRENCE_ORDER
                            + 1
                            - (
                                (RECURRENCE_ORDER + 1)
                                * (RECURRENCE_DEGREE + 1)
                                - 1
                            )
                        ),
                    },
                    "monic_leading_coefficients_low_to_high": leading,
                }
            )

        # Degree ten fails with a large holdout at both exact points over the
        # first prime.  This is a bounded minimal-degree observation only.
        for point in range(len(POINTS)):
            assert run_probe(executable, PRIMES[0], point, 10) == "NONE"

    integer_cubic = lambda value: (
        9 * value**3
        - 1836 * value**2
        + 210_730 * value
        - 17_920_380
    )
    assert all(integer_cubic(value) % 29 for value in range(29))

    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-holonomic-probe-v1",
        "status": (
            "exact period setup plus modular recurrence evidence; "
            "not a universal creative-telescoping certificate"
        ),
        "factor_chart": "C=U*W with U in Mat(4,2), W in Mat(2,4)",
        "period_identity": (
            "nu_m=mu_m/(3m+1)!="
            "CT_u integral_0^1 Phi_C(1,u,t,(1-t)/u)^m dt"
        ),
        "period_identity_checked_orders": 4,
        "exact_points": exact_records,
        "generic_laurent_support": {
            "newton_polygon_vertices": [[-3, 0], [0, 0], [3, 3], [-3, 3]],
            "normalized_volume": 27,
            "forced_endpoint_face": (
                "the u=-3 face is c_30*u^-3*(1-t)^3, so the constrained "
                "beta family lies on the unconstrained face discriminant"
            ),
        },
        "modular_probe": {
            "sequence": "nu_m=mu_m/(3m+1)!",
            "records": modular_records,
            "common_monic_forward_coefficient": {
                "factorization": (
                    "product_(k in {71,73,74,76,77,79,80,82})(3m+k)"
                    "*(9m^3-1836m^2+210730m-17920380)/(9*3^8)"
                ),
                "no_nonnegative_integer_roots": True,
                "integer_root_check": (
                    "the cubic has no root modulo 29; all eight linear "
                    "factors are positive for m>=0"
                ),
            },
            "degree_10_fit_at_order_27": (
                "none at either exact point modulo 1000003"
            ),
        },
        "unresolved_universal_steps": [
            (
                "construct and independently verify telescoping "
                "certificates over the rank-two parameter ring"
            ),
            (
                "recover the parameter denominator/discriminant D(U,W) "
                "multiplying the common forward coefficient"
            ),
            (
                "stratify D=0 and reduce the required bridge moments on "
                "mu_1,...,mu_12,mu_14"
            ),
        ],
        "corrected_system_status": (
            "not decided: an order-27 recurrence needs additional bridge "
            "values unless the corrected ideal supplies them"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS exact rank-two beta/constant-term identity at two points")
    print("PASS six modular order-27 degree-11 recurrence fits")
    print("PASS 139 unused recurrence equations in every modular fit")
    print("PASS common forward coefficient has no integer singular step")
    print("PASS order-27 degree-10 ansatz fails at both sampled points")
    print("PASS result remains modular evidence, not a recurrence certificate")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

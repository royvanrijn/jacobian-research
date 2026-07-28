#!/usr/bin/env python3
"""Modular full-anchor test for two-pair SIC in bidegree (3,3).

Normalize the non-null Sym^2 component to 2*X*T.  Its stabilizer is a
one-dimensional torus.  The locus with only torus-weight-zero higher
components is already closed by the diagonal-slice theorem.  The remaining
locus is covered, up to Weyl reflection, by five charts obtained by setting
one of s0,s1,s2,t0,t1 to one.

This script generates the restricted moments by their torus weights and
tests the five full Sym^6+Sym^4+Sym^2 charts over a finite field.  Its output
is evidence only unless a separate characteristic-zero certificate is made.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from math import factorial
from pathlib import Path
import re
import shutil
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = (
    "s0", "s1", "s2", "s3", "s4", "s5", "s6",
    "t0", "t1", "t2", "t3", "t4",
)
WEIGHTS = (3, 2, 1, 0, -1, -2, -3, 2, 1, 0, -1, -2)

# After removing the unmatched X or Y power prescribed by WEIGHTS, each
# irreducible-basis coefficient is a polynomial in q=X*Y.
Q_POLYNOMIALS = (
    (1,),
    (-3, 3),
    (3, -9, 3),
    (-1, 9, -9, 1),
    (-3, 9, -3),
    (-3, 3),
    (-1,),
    (1, 1),
    (-2, 0, 2),
    (1, -3, -3, 1),
    (2, 0, -2),
    (1, 1),
)

# The normalized quadratic 2*X*T has biform coefficient polynomial
# -1-q+q^2+q^3 in the chosen divided-power basis.
QUADRATIC_Q = (-1, -1, 1, 1)
REPRESENTATIVE_CHARTS = (0, 1, 2, 7, 8)


def convolve(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if left_coefficient % prime == 0:
            continue
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] = (
                answer[left_index + right_index]
                + left_coefficient * right_coefficient
            ) % prime
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_powers(
    polynomial: tuple[int, ...], maximum: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(convolve(powers[-1], polynomial, prime))
    return tuple(powers)


def moment_terms(order: int, prime: int) -> dict[tuple[int, ...], int]:
    """Return mu_order after setting the quadratic scale c to one."""

    factorials = [factorial(index) % prime for index in range(3 * order + 1)]
    inverse_factorials = [
        pow(factorial(index) % prime, -1, prime)
        for index in range(order + 1)
    ]
    basis_powers = tuple(
        polynomial_powers(polynomial, order, prime)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = polynomial_powers(QUADRATIC_Q, order, prime)
    parameter_order = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], int] = defaultdict(int)
    order_factorial = factorials[order]

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(position: int, degree_left: int) -> tuple[int, int]:
        remaining_weights = [WEIGHTS[index] for index in parameter_order[position:]]
        if not remaining_weights or degree_left == 0:
            return 0, 0
        return (
            min(0, degree_left * min(remaining_weights)),
            max(0, degree_left * max(remaining_weights)),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        inverse_denominator: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(parameter_order):
            if weight != 0:
                return
            quadratic_exponent = order - used_degree
            product = convolve(
                q_polynomial,
                quadratic_powers[quadratic_exponent],
                prime,
            )
            scalar = (
                order_factorial
                * inverse_denominator
                * inverse_factorials[quadratic_exponent]
            ) % prime
            contraction = 0
            for q_degree, coefficient in enumerate(product):
                diagonal = shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = scalar * contraction % prime
            if coefficient:
                exponent_tuple = tuple(exponents)
                answer[exponent_tuple] = (
                    answer[exponent_tuple] + coefficient
                ) % prime
            return

        parameter_index = parameter_order[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(position + 1, degree_left)
            if not minimum <= -new_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                inverse_denominator * inverse_factorials[exponent] % prime,
                convolve(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                    prime,
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    return {
        exponents: coefficient
        for exponents, coefficient in answer.items()
        if coefficient % prime
    }


def chart_expression(
    terms: dict[tuple[int, ...], int],
    fixed_index: int,
    prime: int,
) -> str:
    combined: dict[tuple[int, ...], int] = defaultdict(int)
    for exponents, coefficient in terms.items():
        reduced = exponents[:fixed_index] + exponents[fixed_index + 1 :]
        combined[reduced] = (combined[reduced] + coefficient) % prime

    serialized: list[str] = []
    variable_names = PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    for exponents, coefficient in sorted(combined.items()):
        coefficient %= prime
        if not coefficient:
            continue
        factors: list[str] = []
        for variable, exponent in zip(variable_names, exponents):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if not monomial:
            serialized.append(str(coefficient))
        elif coefficient == 1:
            serialized.append(monomial)
        else:
            serialized.append(f"{coefficient}*{monomial}")
    return "+".join(serialized) or "0"


def run_chart(
    singular: str,
    fixed_index: int,
    expressions: list[str],
    prime: int,
    timeout: int,
) -> tuple[int, int, bool, float]:
    variables = PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring anchor={prime},({",".join(variables)}),dp;
option(redSB);
ideal I={",".join(expressions)};
ideal G=std(I);
print("ANCHOR "+string(dim(G))+" "+string(size(G))+" "+string(G[1]==1));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - started
    marker = re.search(r"(?m)^ANCHOR (-?\\d+) (\\d+) ([01])$", completed.stdout)
    if marker is None:
        raise AssertionError(completed.stdout[-2000:])
    dimension, basis_size, unit = marker.groups()
    return int(dimension), int(basis_size), unit == "1", elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument(
        "--charts",
        default="s0,s1,s2,t0,t1",
        help="comma-separated representative charts",
    )
    arguments = parser.parse_args()
    assert arguments.max_order < arguments.prime
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"

    requested = tuple(
        PARAMETERS.index(name)
        for name in arguments.charts.split(",")
        if name
    )
    assert set(requested) <= set(REPRESENTATIVE_CHARTS)

    all_terms: dict[int, dict[tuple[int, ...], int]] = {}
    for order in range(2, arguments.max_order + 1):
        started = time.monotonic()
        terms = moment_terms(order, arguments.prime)
        all_terms[order] = terms
        print(
            f"MOMENT {order} terms={len(terms)} "
            f"seconds={time.monotonic()-started:.2f}",
            flush=True,
        )

    for fixed_index in requested:
        expressions = [
            chart_expression(all_terms[order], fixed_index, arguments.prime)
            for order in range(2, arguments.max_order + 1)
        ]
        try:
            dimension, basis_size, unit, elapsed = run_chart(
                singular,
                fixed_index,
                expressions,
                arguments.prime,
                arguments.timeout,
            )
        except subprocess.TimeoutExpired:
            print(
                f"CHART {PARAMETERS[fixed_index]} TIMEOUT "
                f"seconds={arguments.timeout}",
                flush=True,
            )
            continue
        print(
            f"CHART {PARAMETERS[fixed_index]} dimension={dimension} "
            f"basis={basis_size} unit={int(unit)} seconds={elapsed:.2f}",
            flush=True,
        )

    print(
        "EVIDENCE ONLY: modular full-anchor charts; "
        "no characteristic-zero theorem is promoted"
    )


if __name__ == "__main__":
    main()

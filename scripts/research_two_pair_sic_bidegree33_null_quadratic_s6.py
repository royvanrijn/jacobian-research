#!/usr/bin/env python3
"""Exact null-quadratic, nonzero-s6 chart for bidegree (3,3).

In Clebsch--Gordan coordinates normalize the null quadratic to ``r0=-1``.
On the semistable chart ``s6 != 0``, the residual torus and overall scaling
normalize ``s6=-1``; the unipotent stabilizer of ``r0`` then sets ``s5=0``.
The remaining variables are

    s0,s1,s2,s3,s4,t0,t1,t2,t3,t4.

This script constructs the exact restricted contraction moments by their
torus weights and sends a chosen finite prefix to msolve.  A nonempty output
is only a finite-prefix calculation until its points/components are decoded
and checked by an all-order argument.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
import json
from math import factorial, gcd
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import sympy as sp

from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    Q_POLYNOMIALS,
    WEIGHTS,
)
from research_two_pair_sic_bidegree33_anti_weyl import (
    convolve_integer,
    polynomial_powers_integer,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_null_quadratic_s6.json"
)
VARIABLE_PARAMETER_INDICES = (0, 1, 2, 3, 4, 7, 8, 9, 10, 11)
VARIABLES = tuple(PARAMETERS[index] for index in VARIABLE_PARAMETER_INDICES)
S6_INDEX = 6
R0_Q = (1, 2, 1)
R0_WEIGHT = 1


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orders",
        default="2,3,4,5,6,7,8,9,10,11,12,14",
    )
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--skip-solver", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def parse_orders(specification: str) -> tuple[int, ...]:
    orders = tuple(
        sorted({
            int(piece)
            for piece in specification.split(",")
            if piece.strip()
        })
    )
    if not orders or orders[0] < 1:
        raise ValueError("--orders must contain positive integers")
    return orders


@lru_cache(maxsize=None)
def restricted_moment(order: int) -> dict[tuple[int, ...], int]:
    """Return the primitive moment on r0=s6=-1, s5=0."""

    # Ten variable higher-form terms, followed by the fixed s6 and r0 terms.
    term_names = [
        PARAMETERS[index]
        for index in VARIABLE_PARAMETER_INDICES
    ] + ["s6_fixed", "r0_fixed"]
    term_weights = [
        WEIGHTS[index]
        for index in VARIABLE_PARAMETER_INDICES
    ] + [WEIGHTS[S6_INDEX], R0_WEIGHT]
    term_q = [
        Q_POLYNOMIALS[index]
        for index in VARIABLE_PARAMETER_INDICES
    ] + [Q_POLYNOMIALS[S6_INDEX], R0_Q]
    term_fixed_sign = [1] * len(VARIABLES) + [-1, -1]
    term_powers = [
        polynomial_powers_integer(polynomial, order)
        for polynomial in term_q
    ]
    factorials = [factorial(index) for index in range(3 * order + 1)]
    exponents = [0] * len(VARIABLES)
    answer: dict[tuple[int, ...], int] = defaultdict(int)

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(
        position: int,
        degree_left: int,
    ) -> tuple[int, int]:
        remaining = term_weights[position:]
        if not remaining or degree_left == 0:
            return 0, 0
        return (
            degree_left * min(remaining),
            degree_left * max(remaining),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        exponent_factorials: int,
        sign: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(term_names):
            if used_degree != order or weight:
                return
            contraction = 0
            for q_degree, coefficient in enumerate(q_polynomial):
                diagonal = shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = (
                sign
                * factorial(order)
                // exponent_factorials
                * contraction
            )
            if coefficient:
                answer[tuple(exponents)] += coefficient
            return

        available = order - used_degree
        term_weight = term_weights[position]
        for exponent in range(available + 1):
            new_weight = weight + exponent * term_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(
                position + 1,
                degree_left,
            )
            if not minimum <= -new_weight <= maximum:
                continue
            if position < len(VARIABLES):
                exponents[position] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(term_weight, 0) * exponent,
                exponent_factorials * factorial(exponent),
                sign * term_fixed_sign[position] ** exponent,
                convolve_integer(
                    q_polynomial,
                    term_powers[position][exponent],
                ),
            )
        if position < len(VARIABLES):
            exponents[position] = 0

    visit(0, 0, 0, 0, 1, 1, (1,))
    answer = {
        exponent: coefficient
        for exponent, coefficient in answer.items()
        if coefficient
    }
    content = 0
    for coefficient in answer.values():
        content = gcd(content, abs(coefficient))
    if content:
        answer = {
            exponent: coefficient // content
            for exponent, coefficient in answer.items()
        }
    return answer


def evaluate(
    polynomial: dict[tuple[int, ...], int],
    values: tuple[int | tuple[int, int], ...],
) -> int | tuple[int, int]:
    """Evaluate using pairs (numerator,denominator) without SymPy."""

    from fractions import Fraction

    result = Fraction(0)
    converted = [
        Fraction(value[0], value[1])
        if isinstance(value, tuple)
        else Fraction(value)
        for value in values
    ]
    for exponents, coefficient in polynomial.items():
        term = Fraction(coefficient)
        for value, exponent in zip(converted, exponents, strict=True):
            term *= value**exponent
        result += term
    return (
        result.numerator
        if result.denominator == 1
        else (result.numerator, result.denominator)
    )


def local_survivor_certificate() -> dict[str, object]:
    """Prove exact formal isolation of the Rodrigues point through mu_11."""

    truncation = 5
    free_index = VARIABLES.index("t4")
    pivot_indices = tuple(
        index
        for index in range(len(VARIABLES))
        if index != free_index
    )
    point = [sp.Rational(0)] * len(VARIABLES)
    point[VARIABLES.index("t3")] = sp.Rational(1, 2)
    equations = [restricted_moment(order) for order in range(2, 11)]

    def derivative_at(
        polynomial: dict[tuple[int, ...], int],
        variable_index: int,
    ) -> sp.Rational:
        answer = sp.Rational(0)
        for exponents, coefficient in polynomial.items():
            exponent = exponents[variable_index]
            if not exponent:
                continue
            term = sp.Integer(coefficient * exponent)
            for index, value in enumerate(point):
                power = exponents[index] - (
                    1 if index == variable_index else 0
                )
                term *= value**power
            answer += term
        return sp.Rational(answer)

    jacobian = sp.Matrix([
        [
            derivative_at(polynomial, variable_index)
            for variable_index in pivot_indices
        ]
        for polynomial in equations
    ])
    determinant = sp.factor(jacobian.det())
    expected_determinant = sp.Integer(
        3445505947738252325099075904000
    )
    if determinant != expected_determinant:
        raise AssertionError("unexpected local Jacobian determinant")
    inverse_jacobian = jacobian.inv()

    series = [
        [sp.Rational(0)] * (truncation + 1)
        for _ in VARIABLES
    ]
    for index, value in enumerate(point):
        series[index][0] = value
    series[free_index][1] = sp.Rational(1)

    def convolve(
        left: list[sp.Rational],
        right: list[sp.Rational],
    ) -> list[sp.Rational]:
        answer = [sp.Rational(0)] * (truncation + 1)
        for left_degree, left_value in enumerate(left):
            if not left_value:
                continue
            for right_degree in range(
                truncation - left_degree + 1
            ):
                right_value = right[right_degree]
                if right_value:
                    answer[left_degree + right_degree] += (
                        left_value * right_value
                    )
        return answer

    def power_cache(
        polynomials: list[dict[tuple[int, ...], int]],
    ) -> dict[tuple[int, int], list[sp.Rational]]:
        cache: dict[tuple[int, int], list[sp.Rational]] = {}
        for variable_index in range(len(VARIABLES)):
            cache[(variable_index, 0)] = [
                sp.Rational(1),
                *([sp.Rational(0)] * truncation),
            ]
            maximum = max(
                (
                    exponents[variable_index]
                    for polynomial in polynomials
                    for exponents in polynomial
                ),
                default=0,
            )
            for exponent in range(1, maximum + 1):
                cache[(variable_index, exponent)] = convolve(
                    cache[(variable_index, exponent - 1)],
                    series[variable_index],
                )
        return cache

    def series_coefficient(
        polynomial: dict[tuple[int, ...], int],
        degree: int,
        cache: dict[tuple[int, int], list[sp.Rational]],
    ) -> sp.Rational:
        answer = sp.Rational(0)
        for exponents, coefficient in polynomial.items():
            term = [
                sp.Rational(coefficient),
                *([sp.Rational(0)] * truncation),
            ]
            for variable_index, exponent in enumerate(exponents):
                if exponent:
                    term = convolve(
                        term,
                        cache[(variable_index, exponent)],
                    )
            answer += term[degree]
        return sp.factor(answer)

    formal_coefficients: dict[str, dict[str, str]] = {}
    for degree in range(1, truncation + 1):
        cache = power_cache(equations)
        residual = sp.Matrix([
            -series_coefficient(
                polynomial,
                degree,
                cache,
            )
            for polynomial in equations
        ])
        solution = inverse_jacobian * residual
        nonzero = {}
        for variable_index, value in zip(
            pivot_indices,
            solution,
            strict=True,
        ):
            value = sp.factor(value)
            series[variable_index][degree] = value
            if value:
                nonzero[VARIABLES[variable_index]] = str(value)
        formal_coefficients[str(degree)] = nonzero

    mu11 = restricted_moment(11)
    cache = power_cache([mu11])
    mu11_coefficients = [
        sp.factor(series_coefficient(mu11, degree, cache))
        for degree in range(truncation + 1)
    ]
    expected_leading = sp.Rational(
        558209902860000,
        44871740771,
    )
    if (
        any(mu11_coefficients[:truncation])
        or mu11_coefficients[truncation] != expected_leading
    ):
        raise AssertionError("unexpected formal mu_11 leading term")

    return {
        "point": {
            variable: (
                "1/2" if variable == "t3" else "0"
            )
            for variable in VARIABLES
        },
        "implicit_equations": "mu_2,...,mu_10",
        "free_parameter": "epsilon=t4",
        "pivot_variables": [
            VARIABLES[index]
            for index in pivot_indices
        ],
        "jacobian_determinant": str(determinant),
        "formal_coefficients_through_order_5": formal_coefficients,
        "mu11_restriction_orders_0_through_5": [
            str(value)
            for value in mu11_coefficients
        ],
        "local_intersection_multiplicity": 5,
        "conclusion": (
            "the point is an isolated length-five local component of "
            "(mu_2,...,mu_11) on this normalized chart"
        ),
    }


def five_variable_slice_certificate() -> dict[str, object]:
    """Compute the exact slice containing the formal tangent direction."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the slice certificate")
    retained = ("s0", "s3", "t0", "t3", "t4")
    retained_indices = tuple(VARIABLES.index(name) for name in retained)
    equations = []
    for order in range(2, 11):
        restricted: dict[tuple[int, ...], int] = {}
        for exponents, coefficient in restricted_moment(order).items():
            if any(
                exponent
                for index, exponent in enumerate(exponents)
                if index not in retained_indices
            ):
                continue
            retained_exponents = tuple(
                exponents[index]
                for index in retained_indices
            )
            restricted[retained_exponents] = (
                restricted.get(retained_exponents, 0)
                + coefficient
            )
        content = 0
        for coefficient in restricted.values():
            content = gcd(content, abs(coefficient))
        if content:
            restricted = {
                exponents: coefficient // content
                for exponents, coefficient in restricted.items()
            }
        terms = []
        for exponents in sorted(restricted, reverse=True):
            coefficient = restricted[exponents]
            monomial = "*".join(
                (
                    variable
                    if exponent == 1
                    else f"{variable}^{exponent}"
                )
                for variable, exponent in zip(
                    retained,
                    exponents,
                    strict=True,
                )
                if exponent
            )
            terms.append(
                f"{coefficient}*{monomial}"
                if monomial
                else str(coefficient)
            )
        equations.append("+".join(terms).replace("+-", "-"))

    expected_basis = (
        "17*t0-70*t4,\n"
        "s0,\n"
        "t4^2,\n"
        "7*t3*t4+17*s3,\n"
        "s3*t4,\n"
        "4*t3^2-1,\n"
        "68*s3*t3+7*t4,\n"
        "s3^2"
    )
    script = f"""
ring r=0,({','.join(retained)}),dp;
ideal I={','.join(equations)};
option(redSB);
ideal G=std(I);
print("BASIS_BEGIN");
print(G);
print("BASIS_END");
print("DIM_BEGIN");
print(dim(G));
print("DIM_END");
print("VDIM_BEGIN");
print(vdim(G));
print("VDIM_END");
"""
    with tempfile.TemporaryDirectory(
        prefix="sic33-null-s6-slice-"
    ) as directory:
        input_path = Path(directory) / "slice.sing"
        input_path.write_text(script, encoding="utf-8")
        completed = subprocess.run(
            [singular, "-q", str(input_path)],
            text=True,
            capture_output=True,
            timeout=60,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-2000:])

    def between(start: str, stop: str) -> str:
        return completed.stdout.split(start, 1)[1].split(stop, 1)[0].strip()

    basis = between("BASIS_BEGIN", "BASIS_END")
    dimension = between("DIM_BEGIN", "DIM_END")
    quotient_length = between("VDIM_BEGIN", "VDIM_END")
    if basis != expected_basis or dimension != "0" or quotient_length != "4":
        raise AssertionError("unexpected five-variable slice certificate")
    return {
        "retained_variables": list(retained),
        "zeroed_variables": [
            variable
            for variable in VARIABLES
            if variable not in retained
        ],
        "moments": [2, 10],
        "reduced_groebner_basis": basis.splitlines(),
        "dimension": 0,
        "quotient_length": 4,
        "reduced_support": [
            {
                "s0": "0",
                "s3": "0",
                "t0": "0",
                "t3": "1/2",
                "t4": "0",
            },
            {
                "s0": "0",
                "s3": "0",
                "t0": "0",
                "t3": "-1/2",
                "t4": "0",
            },
        ],
        "conclusion": (
            "the slice contains exactly the two normalized Weyl/torus "
            "copies of the Rodrigues orbit; each has scheme length two"
        ),
    }


def msolve_expression(polynomial: dict[tuple[int, ...], int]) -> str:
    terms = []
    for exponents in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponents]
        monomial = "*".join(
            (
                variable
                if exponent == 1
                else f"{variable}^{exponent}"
            )
            for variable, exponent in zip(
                VARIABLES,
                exponents,
                strict=True,
            )
            if exponent
        )
        if monomial:
            terms.append(f"{coefficient}*{monomial}")
        else:
            terms.append(str(coefficient))
    return "+".join(terms).replace("+-", "-")


def decoded_subprocess_output(value: str | bytes | None) -> str:
    """Normalize TimeoutExpired output across Python subprocess versions."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def main() -> None:
    arguments = parse_arguments()
    orders = parse_orders(arguments.orders)
    msolve = None if arguments.skip_solver else shutil.which("msolve")
    if not arguments.skip_solver and msolve is None:
        raise RuntimeError("msolve is required")

    started = time.monotonic()
    local_certificate = local_survivor_certificate()
    slice_certificate = five_variable_slice_certificate()
    moments = {}
    profiles = {}
    survivor_values = (0, 0, 0, 0, 0, 0, 0, 0, (1, 2), 0)
    for order in orders:
        build_started = time.monotonic()
        moment = restricted_moment(order)
        moments[order] = moment
        survivor_value = evaluate(moment, survivor_values)
        if survivor_value != 0:
            raise AssertionError(
                f"Rodrigues survivor fails restricted mu_{order}: "
                f"{survivor_value}"
            )
        profiles[order] = {
            "terms": len(moment),
            "total_degree": max(
                (sum(exponents) for exponents in moment),
                default=-1,
            ),
            "seconds": round(time.monotonic() - build_started, 6),
        }
        print(
            f"NULL_S6_BUILD order={order} terms={len(moment)} "
            f"seconds={profiles[order]['seconds']:.3f}",
            flush=True,
        )

    if arguments.skip_solver:
        completed = None
        timed_out = False
        stdout = ""
        stderr = ""
        solve_seconds = 0.0
        result = ""
    else:
        with tempfile.TemporaryDirectory(
            prefix="sic33-null-s6-msolve-"
        ) as directory:
            input_path = Path(directory) / "system.ms"
            output_path = Path(directory) / "result.ms"
            input_path.write_text(
                ",".join(VARIABLES)
                + "\n"
                + str(arguments.prime)
                + "\n"
                + ",\n".join(
                    msolve_expression(moments[order])
                    for order in orders
                    if moments[order]
                )
                + "\n",
                encoding="utf-8",
            )
            solve_started = time.monotonic()
            try:
                completed = subprocess.run(
                    [
                        msolve,
                        "-f",
                        str(input_path),
                        "-o",
                        str(output_path),
                        "-t",
                        str(arguments.threads),
                        "-l",
                        "2",
                        "-v",
                        "1",
                    ],
                    text=True,
                    capture_output=True,
                    timeout=(
                        arguments.timeout
                        if arguments.timeout
                        else None
                    ),
                    check=False,
                )
                timed_out = False
            except subprocess.TimeoutExpired as error:
                completed = None
                timed_out = True
                stdout = decoded_subprocess_output(error.stdout)
                stderr = decoded_subprocess_output(error.stderr)
            else:
                stdout = completed.stdout
                stderr = completed.stderr
            solve_seconds = time.monotonic() - solve_started
            result = (
                output_path.read_text(encoding="utf-8").strip()
                if output_path.exists()
                else ""
            )

    payload = {
        "calculation": "two_pair_sic_bidegree33_null_quadratic_s6",
        "scope": (
            "r0=-1, s6=-1, s5=0 normalized null-quadratic chart; "
            "finite moment prefix only"
        ),
        "orders": list(orders),
        "prime": arguments.prime,
        "variables": list(VARIABLES),
        "normalization": {
            "r0": "-1",
            "s6": "-1",
            "s5": "0",
        },
        "known_survivor": {
            "t3": "1/2",
            "all_other_variables": "0",
            "checked_on_every_exported_moment": True,
        },
        "local_survivor_certificate": local_certificate,
        "five_variable_slice_certificate": slice_certificate,
        "profiles": profiles,
        "solver": {
            "skipped": arguments.skip_solver,
            "timed_out": timed_out,
            "returncode": (
                None if completed is None else completed.returncode
            ),
            "seconds": round(solve_seconds, 6),
            "result_head": result[:4000],
            "result_tail": result[-4000:],
            "stdout_tail": stdout[-4000:],
            "stderr_tail": stderr[-4000:],
        },
        "status": (
            "proved_local_certificate"
            if arguments.skip_solver
            else (
                "experiment"
                if timed_out or arguments.prime
                else "exact_finite_prefix_calculation"
            )
        ),
        "seconds": round(time.monotonic() - started, 6),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "orders": payload["orders"],
        "prime": payload["prime"],
        "profiles": profiles,
        "solver": {
            "timed_out": timed_out,
            "returncode": payload["solver"]["returncode"],
            "seconds": payload["solver"]["seconds"],
            "result_head": result[:500],
        },
        "seconds": payload["seconds"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

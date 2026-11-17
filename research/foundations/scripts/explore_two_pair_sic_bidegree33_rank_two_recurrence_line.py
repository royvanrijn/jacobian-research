#!/usr/bin/env python3
"""Reconstruct recurrence coefficients on a rank-two parameter pencil.

This exploratory script samples the modular order-27 recurrence produced by
``explore_two_pair_sic_bidegree33_rank_two_recurrence.cpp`` on

    (U(s), W(s)) = (U_0 + s U_1, W_0 + s W_1)

and searches for low-degree rational interpolants in ``s`` for selected
recurrence coefficients.  It estimates the size of the common parameter
denominator that a universal reconstruction must recover.  A fitted rational
function is modular evidence, not a creative-telescoping certificate.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import sympy as sp


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
    / "two_pair_sic_bidegree33_rank_two_recurrence_line.json"
)
ORDER = 27
DEGREE = 11
MAXIMUM_MOMENT = 500
SELECTED = (
    (0, 0),
    (0, 11),
    (1, 11),
    (7, 11),
    (13, 11),
    (20, 11),
    (26, 0),
    (26, 9),
)


def inverse(value: int, prime: int) -> int:
    return pow(value % prime, prime - 2, prime)


def solve_square(
    matrix: list[list[int]],
    right: list[int],
    prime: int,
) -> list[int] | None:
    size = len(matrix)
    work = [
        [entry % prime for entry in row] + [value % prime]
        for row, value in zip(matrix, right, strict=True)
    ]
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        work[column], work[pivot] = work[pivot], work[column]
        scale = inverse(work[column][column], prime)
        for index in range(column, size + 1):
            work[column][index] = work[column][index] * scale % prime
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            for index in range(column, size + 1):
                work[row][index] = (
                    work[row][index]
                    - scale * work[column][index]
                ) % prime
    return [work[row][-1] for row in range(size)]


def evaluate(coefficients: list[int], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % prime
    return answer


def polynomial_trim(polynomial: list[int]) -> list[int]:
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def polynomial_add(
    left: list[int],
    right: list[int],
    prime: int,
    scale: int = 1,
) -> list[int]:
    answer = left[:] + [0] * max(0, len(right) - len(left))
    for index, coefficient in enumerate(right):
        answer[index] = (
            answer[index] + scale * coefficient
        ) % prime
    return polynomial_trim(answer)


def polynomial_multiply(
    left: list[int],
    right: list[int],
    prime: int,
) -> list[int]:
    if not left or not right:
        return []
    answer = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] = (
                answer[i + j] + left_value * right_value
            ) % prime
    return polynomial_trim(answer)


def polynomial_divmod(
    numerator: list[int],
    denominator: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    remainder = numerator[:]
    quotient = [0] * max(0, len(numerator) - len(denominator) + 1)
    denominator_inverse = inverse(denominator[-1], prime)
    while remainder and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] * denominator_inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % prime
        polynomial_trim(remainder)
    return polynomial_trim(quotient), remainder


def interpolate_consecutive(values: list[int], prime: int) -> list[int]:
    """Return the degree-<N interpolant through (0,y_0),...,(N-1,y_N-1)."""

    differences = [value % prime for value in values]
    answer: list[int] = []
    falling = [1]
    factorial = 1
    for order in range(len(values)):
        coefficient = (
            differences[0] * inverse(factorial, prime) % prime
        )
        answer = polynomial_add(
            answer,
            [
                coefficient * value % prime
                for value in falling
            ],
            prime,
        )
        differences = [
            (right - left) % prime
            for left, right in zip(differences, differences[1:])
        ]
        falling = polynomial_multiply(
            falling,
            [-order % prime, 1],
            prime,
        )
        factorial = factorial * (order + 1) % prime
    return answer


def rational_interpolate(
    samples: list[tuple[int, int]],
    prime: int,
    holdout: int,
) -> dict[str, object] | None:
    fit_samples = samples[: len(samples) - holdout]
    assert [
        parameter for parameter, _ in fit_samples
    ] == list(range(len(fit_samples)))
    fit_values = [value for _, value in fit_samples]
    interpolation = interpolate_consecutive(fit_values, prime)
    modulus = [1]
    for parameter, _ in fit_samples:
        modulus = polynomial_multiply(
            modulus,
            [-parameter % prime, 1],
            prime,
        )

    remainders = (modulus, interpolation)
    cofactors = ([], [1])
    candidates = []
    while remainders[1]:
        numerator = remainders[1]
        denominator = cofactors[1]
        if denominator and all(
            (
                evaluate(numerator, parameter, prime)
                - value * evaluate(denominator, parameter, prime)
            )
            % prime
            == 0
            for parameter, value in samples
        ):
            candidates.append((numerator, denominator))
        quotient, remainder = polynomial_divmod(
            remainders[0],
            remainders[1],
            prime,
        )
        next_cofactor = polynomial_add(
            cofactors[0],
            polynomial_multiply(quotient, cofactors[1], prime),
            prime,
            scale=-1,
        )
        remainders = (remainders[1], remainder)
        cofactors = (cofactors[1], next_cofactor)

    if candidates:
        numerator, denominator = min(
            candidates,
            key=lambda pair: len(pair[0]) + len(pair[1]),
        )
        scale = inverse(denominator[-1], prime)
        numerator = [
            coefficient * scale % prime
            for coefficient in numerator
        ]
        denominator = [
            coefficient * scale % prime
            for coefficient in denominator
        ]
        if denominator:
            variable = sp.symbols("s")
            denominator_expression = sum(
                coefficient * variable**index
                for index, coefficient in enumerate(denominator)
            )
            factorization = sp.factor_list(
                denominator_expression,
                modulus=prime,
            )
            return {
                "numerator_degree": len(numerator) - 1,
                "denominator_degree": len(denominator) - 1,
                "numerator_coefficients_low_to_high": numerator,
                "denominator_coefficients_low_to_high": denominator,
                "denominator_factorization_mod_prime": {
                    "unit": int(factorization[0]) % prime,
                    "factors": [
                        {
                            "coefficients_high_to_low": [
                                int(value) % prime
                                for value in sp.Poly(
                                    factor,
                                    variable,
                                    modulus=prime,
                                ).all_coeffs()
                            ],
                            "multiplicity": multiplicity,
                        }
                        for factor, multiplicity in factorization[1]
                    ],
                },
                "samples_used_for_fit": len(fit_samples),
                "holdout_samples": holdout,
            }
    return None


def parse_recurrence(output: str) -> list[list[int]]:
    lines = output.strip().splitlines()
    if not lines or lines[0] != "FOUND order=27 degree=11":
        raise RuntimeError(output or "empty recurrence-probe output")
    answer = [[0] * (DEGREE + 1) for _ in range(ORDER + 1)]
    for line in lines[1:]:
        values = [int(value) for value in line.split()]
        shift = values[0]
        assert 0 <= shift <= ORDER
        assert len(values) == DEGREE + 2
        answer[shift] = values[1:]
    return answer


def sample_point(
    executable: Path,
    prime: int,
    parameter: int,
    maximum_moment: int,
    family: str,
) -> tuple[int, list[list[int]]]:
    point_code = (
        parameter + 2
        if family == "generic"
        else parameter + 10000
    )
    completed = subprocess.run(
        [
            str(executable),
            str(prime),
            str(maximum_moment),
            str(ORDER),
            str(DEGREE),
            str(point_code),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode:
        raise RuntimeError(
            f"parameter {parameter}: {completed.stdout}\n{completed.stderr}"
        )
    return parameter, parse_recurrence(completed.stdout)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1_000_003)
    parser.add_argument("--samples", type=int, default=72)
    parser.add_argument("--holdout", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument(
        "--family",
        choices=("generic", "scaling"),
        default="generic",
    )
    parser.add_argument(
        "--maximum-moment",
        type=int,
        default=MAXIMUM_MOMENT,
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")
    with tempfile.TemporaryDirectory(prefix="sic33-line-recurrence-") as path:
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
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            records = list(
                executor.map(
                    lambda parameter: sample_point(
                        executable,
                        args.prime,
                        parameter,
                        args.maximum_moment,
                        args.family,
                    ),
                    range(args.samples),
                )
            )

    records.sort()
    reconstructed = {}
    selected_samples = {}
    for shift, exponent in SELECTED:
        samples = [
            (parameter, recurrence[shift][exponent])
            for parameter, recurrence in records
        ]
        label = f"shift_{shift}_m_degree_{exponent}"
        selected_samples[label] = [
            {"parameter": parameter, "value": value}
            for parameter, value in samples
        ]
        reconstructed[label] = (
            rational_interpolate(
                samples,
                args.prime,
                args.holdout,
            )
        )

    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-recurrence-line-v1",
        "status": "modular rational-interpolation experiment",
        "prime": args.prime,
        "parameter_pencil": (
            "(U0+s*U1,W0+s*W1)"
            if args.family == "generic"
            else "((1+s)*U0,W0)"
        ),
        "samples": args.samples,
        "recurrence_order": ORDER,
        "recurrence_m_degree": DEGREE,
        "maximum_moment": args.maximum_moment,
        "sample_values": selected_samples,
        "selected_coefficients": reconstructed,
        "warning": (
            "Rational interpolation of guessed scalar recurrences is not "
            "a creative-telescoping certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")
    for label, result in reconstructed.items():
        if result is None:
            print(f"OPEN {label}: no interpolant within sampled bound")
        else:
            print(
                f"PASS {label}: numerator degree "
                f"{result['numerator_degree']}, denominator degree "
                f"{result['denominator_degree']}, "
                f"{result['holdout_samples']} holdouts"
            )
    try:
        displayed_output = args.output.relative_to(ROOT)
    except ValueError:
        displayed_output = args.output
    print(f"PASS wrote {displayed_output}")


if __name__ == "__main__":
    main()

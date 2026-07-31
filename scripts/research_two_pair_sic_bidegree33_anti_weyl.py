#!/usr/bin/env python3
"""Exact anti-Weyl reduction for the bidegree-(3,3) moment problem.

On the normalized non-null quadratic chart, let

    w_q(X,T) = (a*T, -a^(-1)*X),  q = a^2.

The condition ``w_q(F) = -F`` is

    s4=-q*s2, s5=q^2*s1, s6=-q^3, t2=0,
    t3=q*t1, t4=-q^2*t0

after fixing ``s0=1``.  Haar invariance then kills every odd moment.
By default the residual torus conjugates the involution to ``a=1`` and
the six retained variables are ``s0,s1,s2,s3,t0,t1``.  The original
``s0=1`` ratio chart is also available.  The script constructs the
remaining even moments exactly and optionally sends them, incrementally,
to Singular.  A modular unit is evidence only; a characteristic-zero
unit computation is an exact exclusion.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from math import factorial, gcd
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    Q_POLYNOMIALS,
    QUADRATIC_Q,
    WEIGHTS,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_anti_weyl_research.json"
)
PARAMETERIZATIONS = {
    "normalized": {
        "variables": ("s0", "s1", "s2", "s3", "t0", "t1"),
        "substitution": {
            "s4": "-s2",
            "s5": "s1",
            "s6": "-s0",
            "t2": "0",
            "t3": "t1",
            "t4": "-t0",
        },
    },
    "s0-chart": {
        "variables": ("q", "s1", "s2", "s3", "t0", "t1"),
        "substitution": {
            "s0": "1",
            "s4": "-q*s2",
            "s5": "q^2*s1",
            "s6": "-q^3",
            "t2": "0",
            "t3": "q*t1",
            "t4": "-q^2*t0",
        },
    },
}
PARAMETER_ORDER = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)


def convolve_integer(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if not left_coefficient:
            continue
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_powers_integer(
    polynomial: tuple[int, ...],
    maximum: int,
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(convolve_integer(powers[-1], polynomial))
    return tuple(powers)


def integer_moment_terms(order: int) -> dict[tuple[int, ...], int]:
    """Return the exact raw moment in the twelve higher-form parameters."""

    factorials = [factorial(index) for index in range(3 * order + 1)]
    basis_powers = tuple(
        polynomial_powers_integer(polynomial, order)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = polynomial_powers_integer(QUADRATIC_Q, order)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], int] = defaultdict(int)

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(
        position: int,
        degree_left: int,
    ) -> tuple[int, int]:
        remaining_weights = [
            WEIGHTS[index]
            for index in PARAMETER_ORDER[position:]
        ]
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
        exponent_factorials: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(PARAMETER_ORDER):
            if weight:
                return
            quadratic_exponent = order - used_degree
            product = convolve_integer(
                q_polynomial,
                quadratic_powers[quadratic_exponent],
            )
            multinomial = factorial(order) // (
                exponent_factorials * factorial(quadratic_exponent)
            )
            contraction = 0
            for q_degree, coefficient in enumerate(product):
                diagonal = shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = multinomial * contraction
            if coefficient:
                answer[tuple(exponents)] += coefficient
            return

        parameter_index = PARAMETER_ORDER[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(
                position + 1,
                degree_left,
            )
            if not minimum <= -new_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                exponent_factorials * factorial(exponent),
                convolve_integer(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    return dict(answer)


def anti_weyl_reduce(
    terms: dict[tuple[int, ...], int],
    parameterization: str,
) -> dict[tuple[int, ...], int]:
    """Apply the six-variable anti-Weyl substitution exactly."""

    answer: dict[tuple[int, ...], int] = defaultdict(int)
    for exponents, coefficient in terms.items():
        if exponents[9]:
            continue
        if parameterization == "normalized":
            reduced_exponents = (
                exponents[0] + exponents[6],
                exponents[1] + exponents[5],
                exponents[2] + exponents[4],
                exponents[3],
                exponents[7] + exponents[11],
                exponents[8] + exponents[10],
            )
        else:
            if parameterization != "s0-chart":
                raise ValueError(parameterization)
            reduced_exponents = (
                exponents[4]
                + 2 * exponents[5]
                + 3 * exponents[6]
                + exponents[10]
                + 2 * exponents[11],
                exponents[1] + exponents[5],
                exponents[2] + exponents[4],
                exponents[3],
                exponents[7] + exponents[11],
                exponents[8] + exponents[10],
            )
        if (exponents[4] + exponents[6] + exponents[11]) % 2:
            coefficient = -coefficient
        answer[reduced_exponents] += coefficient
    return {
        exponents: coefficient
        for exponents, coefficient in answer.items()
        if coefficient
    }


def primitive_polynomial(
    polynomial: dict[tuple[int, ...], int],
) -> tuple[dict[tuple[int, ...], int], int]:
    if not polynomial:
        return {}, 0
    content = 0
    for coefficient in polynomial.values():
        content = gcd(content, abs(coefficient))
    primitive = {
        exponents: coefficient // content
        for exponents, coefficient in polynomial.items()
    }
    first_exponents = max(primitive)
    if primitive[first_exponents] < 0:
        primitive = {
            exponents: -coefficient
            for exponents, coefficient in primitive.items()
        }
        content = -content
    return primitive, content


def polynomial_profile(
    polynomial: dict[tuple[int, ...], int],
    content: int,
    variables: tuple[str, ...],
) -> dict[str, object]:
    return {
        "terms": len(polynomial),
        "total_degree": max(
            (sum(exponents) for exponents in polynomial),
            default=-1,
        ),
        "variable_degrees": {
            variable: max(
                (exponents[index] for exponents in polynomial),
                default=-1,
            )
            for index, variable in enumerate(variables)
        },
        "content": str(content),
        "maximum_coefficient_bits": max(
            (abs(coefficient).bit_length() for coefficient in polynomial.values()),
            default=0,
        ),
    }


def singular_expression(
    polynomial: dict[tuple[int, ...], int],
    prime: int,
    variables: tuple[str, ...],
) -> str:
    pieces: list[str] = []
    for exponents, raw_coefficient in sorted(
        polynomial.items(),
        key=lambda item: (sum(item[0]), item[0]),
        reverse=True,
    ):
        coefficient = raw_coefficient % prime if prime else raw_coefficient
        if not coefficient:
            continue
        factors = [
            variable if exponent == 1 else f"{variable}^{exponent}"
            for variable, exponent in zip(variables, exponents)
            if exponent
        ]
        monomial = "*".join(factors)
        if prime:
            pieces.append(
                str(coefficient)
                if not monomial
                else f"{coefficient}*{monomial}"
            )
        elif not monomial:
            pieces.append(str(coefficient))
        elif coefficient == 1:
            pieces.append(monomial)
        elif coefficient == -1:
            pieces.append(f"-{monomial}")
        else:
            pieces.append(f"{coefficient}*{monomial}")
    if prime:
        return "+".join(pieces) or "0"
    return "+".join(pieces).replace("+-", "-") or "0"


def run_singular(
    polynomials: dict[int, dict[tuple[int, ...], int]],
    variables: tuple[str, ...],
    prime: int,
    algorithm: str,
    timeout: int,
    full_system: bool,
) -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        return {"status": "missing-singular"}
    characteristic = str(prime) if prime else "0"
    declarations = "\n".join(
        f"poly mu{order}={singular_expression(polynomial, prime, variables)};"
        for order, polynomial in polynomials.items()
    )
    updates: list[str] = []
    if full_system:
        moment_names = ",".join(
            f"mu{order}"
            for order in reversed(polynomials)
        )
        updates.append(f"ideal I={moment_names};")
        updates.append(f"ideal G={algorithm}(I);")
        updates.append(
            f'print("ANTI_WEYL {max(polynomials)} "+string(dim(G))+" "'
            f'+string(size(G))+" "+string(G[1]==1));'
        )
    else:
        for position, order in enumerate(polynomials):
            if position == 0:
                updates.append(f"ideal G={algorithm}(mu{order});")
            else:
                updates.append(f"G={algorithm}(G+mu{order});")
            updates.append(
                f'print("ANTI_WEYL {order} "+string(dim(G))+" "'
                f'+string(size(G))+" "+string(G[1]==1));'
            )
    program = f"""
ring antiweyl={characteristic},({",".join(variables)}),dp;
option(redSB);
{declarations}
{"".join(updates)}
int finalDimension=dim(G);
int finalLength=-1;
if (finalDimension==0) {{ finalLength=vdim(G); }}
print("ANTI_WEYL_FINAL "+string(finalDimension)+" "+string(finalLength)
      +" "+string(size(G))+" "+string(G[1]==1));
"""
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout if timeout else None,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        return {
            "status": "timeout",
            "seconds": round(time.monotonic() - started, 6),
            "stdout_tail": (error.stdout or "")[-12000:],
            "stderr_tail": (error.stderr or "")[-12000:],
        }
    markers = re.findall(
        r"(?m)^ANTI_WEYL (\d+) (-?\d+) (\d+) ([01])$",
        completed.stdout,
    )
    final = re.search(
        r"(?m)^ANTI_WEYL_FINAL (-?\d+) (-?\d+) (\d+) ([01])$",
        completed.stdout,
    )
    return {
        "status": (
            "completed"
            if completed.returncode == 0 and final is not None
            else "failed"
        ),
        "returncode": completed.returncode,
        "seconds": round(time.monotonic() - started, 6),
        "stages": [
            {
                "order": int(order),
                "dimension": int(dimension),
                "basis_size": int(size),
                "unit": unit == "1",
            }
            for order, dimension, size, unit in markers
        ],
        "final": (
            {
                "dimension": int(final.group(1)),
                "length": int(final.group(2)),
                "basis_size": int(final.group(3)),
                "unit": final.group(4) == "1",
            }
            if final is not None
            else None
        ),
        "stdout_tail": completed.stdout[-12000:],
        "stderr_tail": completed.stderr[-12000:],
        "scope": (
            "exact characteristic zero"
            if not prime
            else f"finite-field scout over GF({prime})"
        ),
    }


def run_msolve(
    polynomials: dict[int, dict[tuple[int, ...], int]],
    variables: tuple[str, ...],
    prime: int,
    timeout: int,
    threads: int,
    linear_algebra: int,
) -> dict[str, object]:
    """Classify the complete modular system with msolve."""

    msolve = shutil.which("msolve")
    if msolve is None:
        return {"status": "missing-msolve"}
    expressions = [
        singular_expression(polynomial, prime, variables)
        for polynomial in polynomials.values()
    ]
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="sic33-anti-weyl-msolve-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(variables)
            + "\n"
            + str(prime)
            + "\n"
            + ",\n".join(expressions)
            + "\n",
            encoding="utf-8",
        )
        try:
            completed = subprocess.run(
                [
                    msolve,
                    "-f",
                    str(input_path),
                    "-o",
                    str(output_path),
                    "-t",
                    str(threads),
                    "-v",
                    "1",
                    "-l",
                    str(linear_algebra),
                ],
                text=True,
                capture_output=True,
                timeout=timeout if timeout else None,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "status": "timeout",
                "seconds": round(time.monotonic() - started, 6),
                "stdout_tail": (error.stdout or "")[-12000:],
                "stderr_tail": (error.stderr or "")[-12000:],
            }
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )
    compact = result.replace(" ", "")
    if result in ("[-1]:", "[-1]"):
        status = "unit"
    elif result.startswith("[1,") and ",-1,[]" in compact:
        status = "positive-dimensional"
    elif completed.returncode == 0 and result:
        status = "finite-nonempty"
    else:
        status = "failed"
    return {
        "status": status,
        "returncode": completed.returncode,
        "seconds": round(time.monotonic() - started, 6),
        "stdout_tail": completed.stdout[-12000:],
        "stderr_tail": completed.stderr[-12000:],
        "result_head": result[:2000],
        "result_tail": result[-2000:],
        "scope": (
            "exact characteristic zero"
            if not prime
            else f"finite-field scout over GF({prime})"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=14)
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument(
        "--backend",
        choices=("singular", "msolve"),
        default="singular",
    )
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--linear-algebra", type=int, default=2)
    parser.add_argument("--full-system", action="store_true")
    parser.add_argument(
        "--parameterization",
        choices=tuple(PARAMETERIZATIONS),
        default="normalized",
    )
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    started = time.monotonic()
    parameterization = PARAMETERIZATIONS[arguments.parameterization]
    variables = parameterization["variables"]
    reduced: dict[int, dict[tuple[int, ...], int]] = {}
    profiles: dict[str, object] = {}
    odd_zero_orders: list[int] = []
    for order in range(1, arguments.through + 1):
        order_started = time.monotonic()
        polynomial = anti_weyl_reduce(
            integer_moment_terms(order),
            arguments.parameterization,
        )
        primitive, content = primitive_polynomial(polynomial)
        if order % 2:
            if primitive:
                raise AssertionError(f"odd moment mu{order} did not vanish")
            odd_zero_orders.append(order)
        else:
            reduced[order] = primitive
        profiles[f"mu{order}"] = {
            **polynomial_profile(primitive, content, variables),
            "seconds": round(time.monotonic() - order_started, 6),
        }
        print(
            f"ANTI_WEYL_BUILD order={order} terms={len(primitive)} "
            f"seconds={time.monotonic() - order_started:.3f}",
            flush=True,
        )
    payload: dict[str, object] = {
        "calculation": "two_pair_sic_bidegree33_anti_weyl",
        "status": "research",
        "parameterization": arguments.parameterization,
        "variables": variables,
        "substitution": parameterization["substitution"],
        "odd_moments_identically_zero": odd_zero_orders,
        "profiles": profiles,
        "generation_seconds": round(time.monotonic() - started, 6),
    }
    if not arguments.generate_only:
        if arguments.backend == "singular":
            payload["singular"] = run_singular(
                reduced,
                variables,
                arguments.prime,
                arguments.algorithm,
                arguments.timeout,
                arguments.full_system,
            )
        else:
            payload["msolve"] = run_msolve(
                reduced,
                variables,
                arguments.prime,
                arguments.timeout,
                arguments.threads,
                arguments.linear_algebra,
            )
    payload["seconds"] = round(time.monotonic() - started, 6)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

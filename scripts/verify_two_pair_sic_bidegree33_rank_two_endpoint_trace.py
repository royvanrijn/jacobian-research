#!/usr/bin/env python3
"""Audit the endpoint trace of a complete modular interior certificate.

At the fixed rank-two point,

    P(u,0) = 11 + 23/u + 91/u^2 + 216/u^3,
    P(u,1) = 354 + 149u + 37u^2 + 19u^3.

Consequently every constant term in an endpoint expression is an
exponential-polynomial in m.  This verifier extracts the Laurent
coefficients of all stored Y_r, constructs those two polynomials exactly
over the certificate prime, and tests the full boundary trace

    CT_u(Y(m,u,1)P(u,1)^m - Y(m,u,0)P(u,0)^m).

The result is a finite-field certificate only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_two_pair_sic_bidegree33_rank_two_relative_divergence import (  # noqa: E402
    EXPECTED_M_DEGREE,
    OPERATOR_INPUT,
    singular_string,
)
from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_endpoint_trace_research.json"
)
TERMINAL_SYZYGY = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing"
)


def trim(polynomial: list[int]) -> list[int]:
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def add(
    left: list[int],
    right: list[int],
    prime: int,
) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] = value
    for index, value in enumerate(right):
        result[index] = (result[index] + value) % prime
    return trim(result)


def multiply(
    left: list[int],
    right: list[int],
    prime: int,
) -> list[int]:
    if not left or not right:
        return []
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            result[left_degree + right_degree] = (
                result[left_degree + right_degree]
                + left_value * right_value
            ) % prime
    return trim(result)


def scale(
    polynomial: list[int],
    scalar: int,
    prime: int,
) -> list[int]:
    return trim([(scalar * value) % prime for value in polynomial])


def shift(polynomial: list[int], degree: int) -> list[int]:
    return ([0] * degree + polynomial) if polynomial else []


def binomial_polynomials(max_degree: int, prime: int) -> list[list[int]]:
    values = [[1]]
    current = [1]
    for degree in range(1, max_degree + 1):
        current = multiply(current, [-(degree - 1) % prime, 1], prime)
        current = scale(current, pow(degree, -1, prime), prime)
        values.append(current)
    return values


def normalized_power_coefficients(
    endpoint: dict[int, int],
    max_exponent: int,
    prime: int,
) -> tuple[int, list[list[int]]]:
    base = endpoint[0] % prime
    inverse_base = pow(base, -1, prime)
    tail = {
        exponent: coefficient * inverse_base % prime
        for exponent, coefficient in endpoint.items()
        if exponent
    }
    binomials = binomial_polynomials(max_exponent, prime)
    powers = [[0] * (max_exponent + 1) for _ in range(max_exponent + 1)]
    powers[0][0] = 1
    for power in range(1, max_exponent + 1):
        for previous_exponent, previous in enumerate(powers[power - 1]):
            if previous == 0:
                continue
            for exponent, coefficient in tail.items():
                total = previous_exponent + exponent
                if total <= max_exponent:
                    powers[power][total] = (
                        powers[power][total] + previous * coefficient
                    ) % prime
    coefficients: list[list[int]] = []
    for exponent in range(max_exponent + 1):
        polynomial: list[int] = []
        for power in range(exponent + 1):
            coefficient = powers[power][exponent]
            if coefficient:
                polynomial = add(
                    polynomial,
                    scale(binomials[power], coefficient, prime),
                    prime,
                )
        coefficients.append(polynomial)
    return base, coefficients


def certificate_degrees(paths: list[Path]) -> list[int]:
    degrees = []
    for path in paths:
        with path.open() as source:
            for line in source:
                match = re.match(r"poly Y(\d+)=", line)
                if match is not None:
                    degrees.append(int(match.group(1)))
    if len(degrees) != len(set(degrees)):
        raise ValueError("duplicate Y certificate degrees")
    return sorted(degrees, reverse=True)


def singular_code(
    prime: int,
    certificates: list[Path],
    degrees: list[int],
    terminal_syzygy: Path | None,
) -> str:
    loads = "\n".join(
        f"execute(read({singular_string(path)}));"
        for path in certificates
    )
    correction = ""
    if terminal_syzygy is not None:
        correction = f"""
execute(read({singular_string(terminal_syzygy)}));
Y0=reduce(Y0-A*R,relation);
"""
    dumps = []
    for degree in degrees:
        dumps.append(
            f"dump_endpoint(subst(Y{degree},t,0),{degree},0);"
        )
        dumps.append(
            f"dump_endpoint(subst(Y{degree},t,1),{degree},1);"
        )
    return f"""
ring r={prime},(u,U,t),dp;
poly Q={q_expression(0)};
poly P=Q*U^3;
poly A=u*diff(Q,u)-3*Q;
poly rel=u*U-1;
ideal relation=std(ideal(rel));
{loads}
{correction}

proc dump_endpoint(poly f,int certificate_degree,int side)
{{
  intvec exponent;
  f=reduce(f,relation);
  while(f!=0)
  {{
    exponent=leadexp(f);
    print("ENDPOINT_TERM");
    print(certificate_degree);
    print(side);
    print(leadcoef(f));
    print(exponent[1]);
    print(exponent[2]);
    f=f-lead(f);
  }}
}}

proc dump_P(poly f,int side)
{{
  intvec exponent;
  f=reduce(f,relation);
  while(f!=0)
  {{
    exponent=leadexp(f);
    print("P_TERM");
    print(side);
    print(leadcoef(f));
    print(exponent[1]);
    print(exponent[2]);
    f=f-lead(f);
  }}
}}

dump_P(subst(P,t,0),0);
dump_P(subst(P,t,1),1);
{"".join(dumps)}
print("PASS dumped endpoint trace");
"""


def evaluate(polynomial: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * value + coefficient) % prime
    return result


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument(
        "--certificate",
        type=Path,
        action="append",
        required=True,
        help="certificate chunk; repeat for every chunk",
    )
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument(
        "--terminal-syzygy",
        type=Path,
        help="apply the final correction Y0 <- Y0-A*R",
    )
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    operator = json.loads(arguments.operator.read_text())
    prime = int(operator["modular_operator"]["prime"])
    degrees = certificate_degrees(arguments.certificate)
    expected_degrees = list(range(EXPECTED_M_DEGREE))
    complete = sorted(degrees) == expected_degrees
    if not complete and not arguments.allow_partial:
        missing = sorted(set(expected_degrees) - set(degrees))
        raise ValueError(f"incomplete certificate; missing Y degrees {missing}")

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(
            prime,
            arguments.certificate,
            degrees,
            arguments.terminal_syzygy,
        ),
        text=True,
        capture_output=True,
        timeout=arguments.timeout,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "?" in combined or "error occurred" in combined:
        raise RuntimeError(combined)
    if "PASS dumped endpoint trace" not in combined:
        raise RuntimeError(combined)

    endpoint_terms = [
        tuple(int(value) for value in values)
        for values in re.findall(
            r"ENDPOINT_TERM\s+(\d+)\s+(\d+)\s+(-?\d+)\s+"
            r"(\d+)\s+(\d+)",
            combined,
        )
    ]
    p_terms = [
        tuple(int(value) for value in values)
        for values in re.findall(
            r"P_TERM\s+(\d+)\s+(-?\d+)\s+(\d+)\s+(\d+)",
            combined,
        )
    ]
    endpoints: dict[int, dict[int, int]] = {0: {}, 1: {}}
    for side, coefficient, u_degree, inverse_degree in p_terms:
        signed_exponent = u_degree - inverse_degree
        directed_exponent = (
            -signed_exponent if side == 0 else signed_exponent
        )
        if directed_exponent < 0:
            raise RuntimeError("endpoint P support points inward")
        endpoints[side][directed_exponent] = coefficient % prime

    relevant: dict[int, list[tuple[int, int, int]]] = {0: [], 1: []}
    discarded_terms = {0: 0, 1: 0}
    for degree, side, coefficient, u_degree, inverse_degree in endpoint_terms:
        signed_exponent = u_degree - inverse_degree
        directed_exponent = signed_exponent if side == 0 else -signed_exponent
        if directed_exponent < 0:
            discarded_terms[side] += 1
            continue
        relevant[side].append(
            (degree, directed_exponent, coefficient % prime)
        )

    component_polynomials: dict[int, list[int]] = {}
    bases: dict[int, int] = {}
    for side in (0, 1):
        max_exponent = max(
            (exponent for _, exponent, _ in relevant[side]),
            default=0,
        )
        base, power_coefficients = normalized_power_coefficients(
            endpoints[side],
            max_exponent,
            prime,
        )
        component: list[int] = []
        for degree, exponent, coefficient in relevant[side]:
            component = add(
                component,
                shift(
                    scale(
                        power_coefficients[exponent],
                        coefficient,
                        prime,
                    ),
                    degree,
                ),
                prime,
            )
        bases[side] = base
        component_polynomials[side] = component

    degrees_of_components = {
        side: len(component_polynomials[side]) - 1
        for side in (0, 1)
    }
    order_bound = sum(
        max(0, degree + 1)
        for degree in degrees_of_components.values()
    )
    if order_bound >= prime:
        raise RuntimeError("endpoint confluent order reaches characteristic")
    trace_prefix = []
    for moment in range(order_bound):
        t0 = (
            pow(bases[0], moment, prime)
            * evaluate(component_polynomials[0], moment, prime)
        ) % prime
        t1 = (
            pow(bases[1], moment, prime)
            * evaluate(component_polynomials[1], moment, prime)
        ) % prime
        trace_prefix.append((t1 - t0) % prime)
    components_zero = all(
        not polynomial for polynomial in component_polynomials.values()
    )
    trace_zero = not any(trace_prefix)
    if complete and trace_zero != components_zero:
        raise RuntimeError(
            "confluent endpoint representations disagree on zero trace"
        )

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "endpoint-trace-research-v1"
        ),
        "status": (
            "exact modular zero endpoint trace"
            if complete and trace_zero
            else "exact modular nonzero endpoint trace"
            if complete
            else "partial modular endpoint trace; not conclusive"
        ),
        "prime": prime,
        "point": operator["point"],
        "complete_certificate": complete,
        "terminal_syzygy": (
            display_path(arguments.terminal_syzygy)
            if arguments.terminal_syzygy is not None
            else None
        ),
        "certificate_degrees": degrees,
        "endpoint_P_directed_coefficients": endpoints,
        "endpoint_exponential_bases": bases,
        "relevant_Y_terms": {
            side: len(relevant[side]) for side in (0, 1)
        },
        "discarded_Y_terms": discarded_terms,
        "component_degrees": degrees_of_components,
        "component_nonzero_terms": {
            side: len(component_polynomials[side]) for side in (0, 1)
        },
        "confluent_order_bound": order_bound,
        "trace_prefix_nonzero_count": sum(
            value != 0 for value in trace_prefix
        ),
        "trace_prefix_first_nonzero": next(
            (
                {"m": moment, "value": value}
                for moment, value in enumerate(trace_prefix)
                if value
            ),
            None,
        ),
        "conclusion": (
            (
                "both endpoint exponential-polynomial components vanish, "
                "so the modular interior divergence has zero integrated "
                "boundary for every nonnegative m"
            )
            if complete and trace_zero
            else (
                "the stored interior divergence has a nonzero endpoint "
                "trace and must be extended by endpoint states"
            )
            if complete
            else (
                "run again with all 58 certificate degrees before drawing "
                "an endpoint conclusion"
            )
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    if complete:
        print(
            "PASS exact modular endpoint trace is "
            f"{'zero' if trace_zero else 'nonzero'}"
        )
    else:
        print(f"PASS audited {len(degrees)} partial certificate degrees")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()

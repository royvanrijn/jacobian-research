#!/usr/bin/env python3
"""Try to remove the terminal residual by the final Koszul syzygy.

For the interior generators

    A = u Q_u - 3Q,  C = Q_t,

all polynomial syzygies are generated (up to the Laurent relation) by
(C,-A).  Replacing the last certificate pair by

    (X_0,Y_0) + R(C,-A)

changes its divergence by

    H(R) = C(D_u+3)R - A R_t.

This script reconstructs the exact terminal target from the stored X_0,Y_0
certificate, verifies it against the independently checked checkpoint, and
performs a sparse triangular image reduction for H over the certificate
field.  Success gives a final-syzygy certificate; failure gives an exact
leading obstruction for this last-step repair.  Neither outcome is a
characteristic-zero statement.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import heapq
import json
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_sic_bidegree33_rank_two_holonomic_probe import (  # noqa: E402
    POINTS,
    matrix_product,
    substituted_polynomial,
)


OPERATOR = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json"
)
CERTIFICATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m18_m0.sing"
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
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_research.json"
)
R_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing"
)

Exponent = tuple[int, int]
Polynomial = dict[Exponent, int]


def add_term(
    polynomial: Polynomial,
    exponent: Exponent,
    coefficient: int,
    prime: int,
) -> None:
    value = (polynomial.get(exponent, 0) + coefficient) % prime
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def parse_singular_polynomial(expression: str, prime: int) -> Polynomial:
    expression = expression.strip().rstrip(";")
    if expression == "0":
        return {}
    answer: Polynomial = {}
    for raw_term in re.findall(r"[+-]?[^+-]+", expression):
        sign = -1 if raw_term.startswith("-") else 1
        term = raw_term.lstrip("+-")
        coefficient_match = re.match(r"\d+", term)
        if coefficient_match is None:
            coefficient = 1
            variables = term
        else:
            coefficient = int(coefficient_match.group())
            variables = term[coefficient_match.end() :]
        u_exponent = 0
        t_exponent = 0
        consumed = ""
        for variable_match in re.finditer(r"(u|U|t)(\d*)", variables):
            variable = variable_match.group(1)
            exponent = (
                int(variable_match.group(2))
                if variable_match.group(2)
                else 1
            )
            if variable == "u":
                u_exponent += exponent
            elif variable == "U":
                u_exponent -= exponent
            else:
                t_exponent += exponent
            consumed += variable_match.group()
        if consumed != variables:
            raise ValueError(f"cannot parse Singular term {raw_term!r}")
        add_term(
            answer,
            (u_exponent, t_exponent),
            sign * coefficient,
            prime,
        )
    return answer


def named_certificate_polynomial(
    path: Path,
    name: str,
    prime: int,
) -> Polynomial:
    prefix = f"poly {name}="
    with path.open() as source:
        for line in source:
            if line.startswith(prefix):
                return parse_singular_polynomial(
                    line[len(prefix) :],
                    prime,
                )
    raise ValueError(f"{name} not found in {path}")


def fixed_point_polynomials(prime: int) -> tuple[Polynomial, Polynomial]:
    raw_u, raw_w = POINTS[0]
    u = [[Fraction(value) for value in row] for row in raw_u]
    w = [[Fraction(value) for value in row] for row in raw_w]
    raw_p = substituted_polynomial(matrix_product(u, w))
    p: Polynomial = {}
    q: Polynomial = {}
    for exponent, coefficient in raw_p.items():
        modular = (
            coefficient.numerator
            * pow(coefficient.denominator % prime, -1, prime)
        ) % prime
        p[exponent] = modular
        q[(exponent[0] + 3, exponent[1])] = modular
    return p, q


def add(left: Polynomial, right: Polynomial, prime: int) -> Polynomial:
    answer = dict(left)
    for exponent, coefficient in right.items():
        add_term(answer, exponent, coefficient, prime)
    return answer


def scale(polynomial: Polynomial, scalar: int, prime: int) -> Polynomial:
    scalar %= prime
    if scalar == 0:
        return {}
    return {
        exponent: coefficient * scalar % prime
        for exponent, coefficient in polynomial.items()
    }


def multiply(left: Polynomial, right: Polynomial, prime: int) -> Polynomial:
    answer: Polynomial = {}
    for (left_u, left_t), left_coefficient in left.items():
        for (right_u, right_t), right_coefficient in right.items():
            add_term(
                answer,
                (left_u + right_u, left_t + right_t),
                left_coefficient * right_coefficient,
                prime,
            )
    return answer


def divergence(
    x: Polynomial,
    y: Polynomial,
    prime: int,
) -> Polynomial:
    answer: Polynomial = {}
    for (u_degree, t_degree), coefficient in x.items():
        add_term(
            answer,
            (u_degree, t_degree),
            coefficient * u_degree,
            prime,
        )
    for (u_degree, t_degree), coefficient in y.items():
        if t_degree:
            add_term(
                answer,
                (u_degree, t_degree - 1),
                coefficient * t_degree,
                prime,
            )
    return answer


def target_polynomial(
    p: Polynomial,
    operator_coefficients: list[list[int]],
    x0: Polynomial,
    y0: Polynomial,
    prime: int,
) -> Polynomial:
    powers = [{(0, 0): 1}]
    for _ in range(1, len(operator_coefficients)):
        powers.append(multiply(powers[-1], p, prime))
    target: Polynomial = {}
    for shift, coefficients in enumerate(operator_coefficients):
        target = add(
            target,
            scale(powers[shift], coefficients[0], prime),
            prime,
        )
    return add(
        target,
        scale(divergence(x0, y0, prime), -1, prime),
        prime,
    )


def h_monomial(
    u_degree: int,
    t_degree: int,
    q: Polynomial,
    prime: int,
) -> Polynomial:
    answer: Polynomial = {}
    for (q_u, q_t), q_coefficient in q.items():
        coefficient = q_coefficient * (
            q_t * (u_degree + 3) - t_degree * (q_u - 3)
        )
        if coefficient:
            add_term(
                answer,
                (u_degree + q_u, t_degree + q_t - 1),
                coefficient,
                prime,
            )
    return answer


def leading_exponent(polynomial: Polynomial) -> Exponent:
    return max(polynomial, key=lambda exponent: (exponent[1], exponent[0]))


def pivot_for(
    exponent: Exponent,
    q_terms: list[tuple[Exponent, int]],
    q: Polynomial,
    prime: int,
) -> tuple[int, int, int, Polynomial] | None:
    target_u, target_t = exponent
    for (q_u, q_t), _ in q_terms:
        r_u = target_u - q_u
        r_t = target_t - q_t + 1
        if r_t < 0:
            continue
        image = h_monomial(r_u, r_t, q, prime)
        if image and leading_exponent(image) == exponent:
            return r_u, r_t, image[exponent], image
    return None


def singular_expression(polynomial: Polynomial, prime: int) -> str:
    terms = []
    for (u_degree, t_degree), coefficient in sorted(
        polynomial.items(),
        key=lambda item: (item[0][1], item[0][0]),
        reverse=True,
    ):
        centered = (
            coefficient
            if coefficient <= prime // 2
            else coefficient - prime
        )
        factors = [str(centered)]
        if u_degree > 0:
            factors.append("u" if u_degree == 1 else f"u{u_degree}")
        elif u_degree < 0:
            inverse_degree = -u_degree
            factors.append(
                "U" if inverse_degree == 1 else f"U{inverse_degree}"
            )
        if t_degree:
            factors.append("t" if t_degree == 1 else f"t{t_degree}")
        terms.append("".join(factors))
    return "+".join(terms).replace("+-", "-") if terms else "0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--terminal", type=Path, default=TERMINAL)
    parser.add_argument(
        "--max-steps",
        type=int,
        default=1000,
        help=(
            "bounded monomial pivots; the unbounded expansion has severe "
            "fill-in, so use the adjacent t-block solver for completion"
        ),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--R-output", type=Path, default=R_OUTPUT)
    arguments = parser.parse_args()

    operator = json.loads(arguments.operator.read_text())
    modular = operator["modular_operator"]
    prime = int(modular["prime"])
    coefficients = modular["coefficients"]
    p, q = fixed_point_polynomials(prime)
    x0 = named_certificate_polynomial(
        arguments.certificate,
        "X0",
        prime,
    )
    y0 = named_certificate_polynomial(
        arguments.certificate,
        "Y0",
        prime,
    )
    target = target_polynomial(p, coefficients, x0, y0, prime)
    terminal = parse_singular_polynomial(
        arguments.terminal.read_text(),
        prime,
    )
    if multiply(q, target, prime) != terminal:
        raise RuntimeError(
            "reconstructed terminal target does not match checkpoint"
        )

    remainder = dict(target)
    correction: Polynomial = {}
    heap = [
        (-t_degree, -u_degree, u_degree, t_degree)
        for u_degree, t_degree in remainder
    ]
    heapq.heapify(heap)
    q_terms = sorted(
        q.items(),
        key=lambda item: (item[0][1], item[0][0]),
        reverse=True,
    )
    steps = 0
    obstruction = None
    while remainder:
        if arguments.max_steps is not None and steps >= arguments.max_steps:
            break
        while heap:
            _, _, u_degree, t_degree = heapq.heappop(heap)
            exponent = (u_degree, t_degree)
            if exponent in remainder:
                break
        else:
            raise RuntimeError("remainder heap exhausted")
        pivot = pivot_for(exponent, q_terms, q, prime)
        if pivot is None:
            obstruction = {
                "u_degree": exponent[0],
                "t_degree": exponent[1],
                "coefficient": remainder[exponent],
            }
            break
        r_u, r_t, pivot_coefficient, image = pivot
        scalar = (
            remainder[exponent] * pow(pivot_coefficient, -1, prime)
        ) % prime
        add_term(correction, (r_u, r_t), scalar, prime)
        for image_exponent, image_coefficient in image.items():
            was_present = image_exponent in remainder
            add_term(
                remainder,
                image_exponent,
                -scalar * image_coefficient,
                prime,
            )
            if not was_present and image_exponent in remainder:
                heapq.heappush(
                    heap,
                    (
                        -image_exponent[1],
                        -image_exponent[0],
                        image_exponent[0],
                        image_exponent[1],
                    ),
                )
        steps += 1

    complete = not remainder
    if complete:
        if h_monomial(0, 0, {}, prime):
            raise AssertionError
        reconstructed: Polynomial = {}
        for (u_degree, t_degree), coefficient in correction.items():
            reconstructed = add(
                reconstructed,
                scale(
                    h_monomial(u_degree, t_degree, q, prime),
                    coefficient,
                    prime,
                ),
                prime,
            )
        if reconstructed != target:
            raise RuntimeError("final syzygy reconstruction failed")
        arguments.R_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.R_output.write_text(
            "// Exact modular final-syzygy correction\n"
            f"// prime={prime}\n"
            f"poly R={singular_expression(correction, prime)};\n"
        )

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "terminal-syzygy-research-v1"
        ),
        "status": (
            "exact modular final-syzygy certificate"
            if complete
            else "bounded triangular terminal-syzygy reduction"
            if arguments.max_steps is not None and obstruction is None
            else "exact modular leading obstruction to the final-syzygy repair"
        ),
        "prime": prime,
        "point": operator["point"],
        "terminal_terms": len(terminal),
        "target_after_dividing_Q_terms": len(target),
        "X0_terms": len(x0),
        "Y0_terms": len(y0),
        "triangular_order": "descending (t_degree,u_degree)",
        "reduction_steps": steps,
        "remaining_terms": len(remainder),
        "leading_obstruction": obstruction,
        "R_terms": len(correction),
        "R_certificate": str(arguments.R_output) if complete else None,
        "conclusion": (
            "the terminal residual is removed by the last Koszul syzygy"
            if complete
            else (
                "the triangular monomial images of "
                "H=C(D_u+3)-A*d_t do not span the terminal target; "
                "higher-degree syzygy freedom or endpoint states remain "
                "necessary"
            )
            if obstruction is not None
            else (
                "increase --max-steps to continue the exact sparse "
                "triangular reduction"
            )
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    if complete:
        print(
            f"PASS final syzygy removes all {len(target)} target terms"
        )
        print(f"PASS wrote {arguments.R_output}")
    elif obstruction is not None:
        print(
            "PASS exact final-syzygy leading obstruction at "
            f"u^{obstruction['u_degree']}t^{obstruction['t_degree']}"
        )
    else:
        print(
            f"PASS completed {steps} bounded terminal-syzygy reductions"
        )
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()

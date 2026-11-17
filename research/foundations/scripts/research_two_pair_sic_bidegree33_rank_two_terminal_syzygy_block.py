#!/usr/bin/env python3
"""Solve the final-syzygy equation by descending t-degree blocks.

The monomial triangular reduction has severe fill-in.  Since Q has
t-degree three, the top t-coefficient of

    H(R)=Q_t(D_u+3)R-(uQ_u-3Q)R_t

depends on one Laurent polynomial R_b(u) at a time.  Each block is a
banded seven-diagonal system controlled by Q_3(u).  This script solves
those systems exactly over the certificate field, descends in t, and
either writes R or records the first exact compatibility obstruction.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_two_pair_sic_bidegree33_rank_two_terminal_syzygy import (  # noqa: E402
    CERTIFICATE,
    OPERATOR,
    R_OUTPUT,
    TERMINAL,
    Polynomial,
    add_term,
    fixed_point_polynomials,
    h_monomial,
    multiply,
    named_certificate_polynomial,
    parse_singular_polynomial,
    singular_expression,
    target_polynomial,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_block_research.json"
)

Affine = tuple[int, int]


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def affine_add(left: Affine, right: Affine, prime: int) -> Affine:
    return (
        (left[0] + right[0]) % prime,
        (left[1] + right[1]) % prime,
    )


def affine_scale(value: Affine, scalar: int, prime: int) -> Affine:
    return (
        value[0] * scalar % prime,
        value[1] * scalar % prime,
    )


def solve_block(
    target: dict[int, int],
    t_degree: int,
    q3: dict[int, int],
    prime: int,
    lower_margin: int,
) -> tuple[dict[int, int] | None, dict | None]:
    if not target:
        return {}, None
    target_min = min(target)
    target_max = max(target)
    r_min = target_min - lower_margin
    r_max = target_max
    values: dict[int, Affine] = {}
    resonance = t_degree - 3

    def contribution(output_degree: int, include_i6: bool) -> Affine:
        answer = (target.get(output_degree, 0), 0)
        maximum_i = 6 if include_i6 else 5
        for q_degree, q_coefficient in q3.items():
            if q_degree > maximum_i:
                continue
            r_degree = output_degree - q_degree
            if r_degree not in values:
                continue
            coefficient = q_coefficient * (
                3 * (r_degree + 3)
                - t_degree * (q_degree - 3)
            )
            answer = affine_add(
                answer,
                affine_scale(values[r_degree], -coefficient, prime),
                prime,
            )
        return answer

    for output_degree in range(r_max + 6, r_min + 5, -1):
        r_degree = output_degree - 6
        right = contribution(output_degree, include_i6=False)
        diagonal = q3[6] * 3 * (r_degree + 3 - t_degree) % prime
        if r_degree == resonance:
            if diagonal != 0:
                raise AssertionError
            if right != (0, 0):
                return None, {
                    "kind": "resonant_top_block",
                    "R_u_degree": r_degree,
                    "output_u_degree": output_degree,
                    "compatibility": list(right),
                }
            values[r_degree] = (0, 1)
        else:
            if diagonal == 0:
                raise AssertionError
            inverse = pow(diagonal, -1, prime)
            values[r_degree] = affine_scale(right, inverse, prime)

    constraints: list[tuple[int, Affine]] = []
    for output_degree in range(r_min + 5, r_min - 1, -1):
        right = contribution(output_degree, include_i6=True)
        if right != (0, 0):
            constraints.append((output_degree, right))
    parameter = None
    for output_degree, (constant, linear) in constraints:
        if linear:
            candidate = -constant * pow(linear, -1, prime) % prime
            if parameter is None:
                parameter = candidate
            elif parameter != candidate:
                return None, {
                    "kind": "inconsistent_lower_tail",
                    "output_u_degree": output_degree,
                    "first_parameter": parameter,
                    "second_parameter": candidate,
                }
        elif constant:
            return None, {
                "kind": "nonzero_lower_tail",
                "output_u_degree": output_degree,
                "compatibility": [constant, linear],
            }
    if parameter is None:
        parameter = 0
    solution = {}
    for degree, (constant, linear) in values.items():
        value = (constant + linear * parameter) % prime
        if value:
            solution[degree] = value

    reconstructed: dict[int, int] = {}
    for r_degree, r_coefficient in solution.items():
        for q_degree, q_coefficient in q3.items():
            coefficient = q_coefficient * (
                3 * (r_degree + 3)
                - t_degree * (q_degree - 3)
            )
            exponent = r_degree + q_degree
            value = (
                reconstructed.get(exponent, 0)
                + r_coefficient * coefficient
            ) % prime
            if value:
                reconstructed[exponent] = value
            else:
                reconstructed.pop(exponent, None)
    if reconstructed != target:
        return None, {
            "kind": "block_reconstruction_mismatch",
            "target_terms": len(target),
            "reconstructed_terms": len(reconstructed),
        }
    return solution, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--terminal", type=Path, default=TERMINAL)
    parser.add_argument("--lower-margin", type=int, default=12)
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--R-output", type=Path, default=R_OUTPUT)
    arguments = parser.parse_args()
    if arguments.lower_margin < 6:
        raise ValueError("lower margin must be at least six")

    operator = json.loads(arguments.operator.read_text())
    modular = operator["modular_operator"]
    prime = int(modular["prime"])
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
    target = target_polynomial(
        p,
        modular["coefficients"],
        x0,
        y0,
        prime,
    )
    terminal = parse_singular_polynomial(
        arguments.terminal.read_text(),
        prime,
    )
    if multiply(q, target, prime) != terminal:
        raise RuntimeError("terminal checkpoint mismatch")
    q3 = {
        u_degree: coefficient
        for (u_degree, t_degree), coefficient in q.items()
        if t_degree == 3
    }
    if set(q3) != set(range(7)) or q3[6] == 0:
        raise RuntimeError("unexpected top t-coefficient support")

    remainder = dict(target)
    correction: Polynomial = {}
    maximum_t = max(t_degree for _, t_degree in remainder)
    maximum_r_t = maximum_t - 2
    blocks_completed = 0
    obstruction = None
    block_statistics = []
    for r_t in range(maximum_r_t, -1, -1):
        if (
            arguments.max_blocks is not None
            and blocks_completed >= arguments.max_blocks
        ):
            break
        output_t = r_t + 2
        block_target = {
            u_degree: coefficient
            for (u_degree, t_degree), coefficient in remainder.items()
            if t_degree == output_t
        }
        solution, block_obstruction = solve_block(
            block_target,
            r_t,
            q3,
            prime,
            arguments.lower_margin,
        )
        if block_obstruction is not None:
            obstruction = {
                "R_t_degree": r_t,
                "output_t_degree": output_t,
                **block_obstruction,
            }
            break
        assert solution is not None
        for u_degree, coefficient in solution.items():
            add_term(
                correction,
                (u_degree, r_t),
                coefficient,
                prime,
            )
            image = h_monomial(u_degree, r_t, q, prime)
            for exponent, image_coefficient in image.items():
                add_term(
                    remainder,
                    exponent,
                    -coefficient * image_coefficient,
                    prime,
                )
        block_statistics.append(
            {
                "R_t_degree": r_t,
                "target_terms": len(block_target),
                "R_terms": len(solution),
                "remaining_terms": len(remainder),
            }
        )
        blocks_completed += 1

    complete = not remainder
    bounded = (
        arguments.max_blocks is not None
        and blocks_completed >= arguments.max_blocks
        and obstruction is None
        and not complete
    )
    if complete:
        arguments.R_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.R_output.write_text(
            "// Exact modular final-syzygy correction\n"
            f"// prime={prime}, block solver\n"
            f"poly R={singular_expression(correction, prime)};\n"
        )
        serialized = named_certificate_polynomial(
            arguments.R_output,
            "R",
            prime,
        )
        if serialized != correction:
            raise RuntimeError("serialized R does not match sparse R")
        reconstructed: Polynomial = {}
        for (u_degree, t_degree), coefficient in serialized.items():
            image = h_monomial(u_degree, t_degree, q, prime)
            for exponent, image_coefficient in image.items():
                add_term(
                    reconstructed,
                    exponent,
                    coefficient * image_coefficient,
                    prime,
                )
        if reconstructed != target:
            raise RuntimeError("serialized R does not reconstruct target")
    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "terminal-syzygy-block-research-v1"
        ),
        "status": (
            "exact modular final-syzygy certificate"
            if complete
            else "bounded exact t-block reduction"
            if bounded
            else "exact modular t-block obstruction to final syzygy"
        ),
        "prime": prime,
        "point": operator["point"],
        "terminal_terms": len(terminal),
        "target_after_dividing_Q_terms": len(target),
        "maximum_target_t_degree": maximum_t,
        "maximum_R_t_degree": maximum_r_t,
        "lower_u_margin": arguments.lower_margin,
        "blocks_completed": blocks_completed,
        "remaining_terms": len(remainder),
        "R_terms": len(correction),
        "obstruction": obstruction,
        "last_blocks": block_statistics[-10:],
        "R_certificate": (
            display_path(arguments.R_output) if complete else None
        ),
        "conclusion": (
            "the terminal residual is removed by the last Koszul syzygy"
            if complete
            else (
                "the top-t banded compatibility fails, so the final "
                "Koszul correction is impossible in the tested Laurent "
                "support class"
            )
            if obstruction is not None
            else "increase --max-blocks to continue the exact reduction"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    if complete:
        print(
            f"PASS final block syzygy removes all {len(target)} terms"
        )
        print(f"PASS wrote {arguments.R_output}")
    elif obstruction is not None:
        print(
            "PASS exact block obstruction at "
            f"R t-degree {obstruction['R_t_degree']}"
        )
    else:
        print(f"PASS completed {blocks_completed} bounded t-blocks")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()

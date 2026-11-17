#!/usr/bin/env python3
"""Explore the uniform diagonal core at the first degree-eight threshold.

For the normalized mixed family

    V_(d,d) + V_(1,d+1) + V_(d+1,1),

the d=2 and d=3 proofs reduce every pure contraction on the nonzero
positive branch to a diagonal balanced block plus one opposite corner.
This script constructs that core directly for d=4 and checks moments
two through six over two finite fields.

The modular Groebner calculations are exact in their stated
characteristics.  They are evidence only for the characteristic-zero
radical and are deliberately not registered as a theorem.
"""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_mixed_diagonal_core.json"
)
DEGREE = 4
PRIMES = (101, 1009)


def compositions(
    total: int, parts: int, prefix: tuple[int, ...] = ()
):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def diagonal_core_moment(
    order: int,
    diagonal: tuple[sp.Expr, ...],
    corner: sp.Symbol,
) -> sp.Expr:
    """Return the exact scalar contraction of the order-th core power."""

    result = sp.Integer(0)
    for corner_count in range(order // 2 + 1):
        diagonal_count = order - 2 * corner_count
        for counts in compositions(diagonal_count, DEGREE + 1):
            xi1_degree = (DEGREE + 1) * corner_count + sum(
                index * counts[index] for index in range(DEGREE + 1)
            )
            xi2_degree = corner_count + sum(
                (DEGREE - index) * counts[index]
                for index in range(DEGREE + 1)
            )
            multinomial = factorial(order) // factorial(corner_count) ** 2
            for count in counts:
                multinomial //= factorial(count)
            term = (
                sp.Integer(multinomial)
                * factorial(xi1_degree)
                * factorial(xi2_degree)
                * corner**corner_count
            )
            for coefficient, count in zip(diagonal, counts):
                term *= coefficient**count
            result += term
    return sp.expand(result)


def singular_modular_result(
    singular: str,
    prime: int,
    polynomials: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, object]:
    expressions = []
    for polynomial in polynomials:
        integral = sp.Poly(
            polynomial, *variables, domain=sp.QQ
        ).clear_denoms(convert=True)[1].as_expr()
        expressions.append(str(integral).replace("**", "^"))

    variable_names = ",".join(str(variable) for variable in variables)
    commands = f"""
ring R={prime},({variable_names}),dp;
ideal I={",".join(expressions)};
option(redSB);
ideal G=std(I);
size(G);
vdim(G);
reduce(c1^40,G);
reduce(c2^40,G);
reduce(c3^40,G);
reduce(c4^40,G);
reduce(h^20,G);
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=commands,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    numeric_lines = [
        int(line)
        for line in completed.stdout.splitlines()
        if re.fullmatch(r"-?[0-9]+", line.strip())
    ]
    assert numeric_lines == [132, 360, 0, 0, 0, 0, 0]
    return {
        "prime": prime,
        "groebner_basis_size": numeric_lines[0],
        "quotient_dimension": numeric_lines[1],
        "power_remainders": {
            "c1^40": numeric_lines[2],
            "c2^40": numeric_lines[3],
            "c3^40": numeric_lines[4],
            "c4^40": numeric_lines[5],
            "h^20": numeric_lines[6],
        },
    }


def main() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the modular core audit")

    coefficients = sp.symbols("c0:5")
    corner = sp.symbols("h")
    leading = -sum(
        coefficients[index]
        * factorial(index)
        * factorial(DEGREE - index)
        for index in range(1, DEGREE + 1)
    ) / factorial(DEGREE)
    diagonal = (leading,) + coefficients[1:]
    variables = (*coefficients[1:], corner)

    normalized_moments: dict[int, sp.Expr] = {}
    contents: dict[int, sp.Rational] = {}
    term_counts: dict[int, int] = {}
    first = diagonal_core_moment(1, diagonal, corner)
    assert first == 0
    for order in range(2, 7):
        moment = diagonal_core_moment(order, diagonal, corner)
        content, primitive = sp.Poly(moment, *variables).primitive()
        normalized_moments[order] = primitive.as_expr()
        contents[order] = content
        term_counts[order] = len(primitive.terms())

    jacobian = sp.Matrix(
        [
            [
                sp.diff(normalized_moments[order], variable)
                for variable in variables
            ]
            for order in range(2, 7)
        ]
    )
    point = {
        variable: index + 1 for index, variable in enumerate(variables)
    }
    assert jacobian.subs(point).rank() == 5

    modular_results = [
        singular_modular_result(
            singular,
            prime,
            list(normalized_moments.values()),
            variables,
        )
        for prime in PRIMES
    ]

    payload = {
        "claim_boundary": (
            "Exact finite-field evidence only; no characteristic-zero "
            "radical or full d=4 triangularization theorem is claimed"
        ),
        "degree": DEGREE,
        "core": (
            "sum_(i=0)^4 c_i (xi1*z1)^i (xi2*z2)^(4-i) "
            "+ xi2*z1^5 + h*xi1^5*z2, with mu1=0"
        ),
        "moment_orders": list(range(2, 7)),
        "contents": {
            str(order): str(value) for order, value in contents.items()
        },
        "term_counts": {
            str(order): value for order, value in term_counts.items()
        },
        "normalized_moments": {
            str(order): str(polynomial)
            for order, polynomial in normalized_moments.items()
        },
        "exact_jacobian_rank_at_1_2_3_4_5": 5,
        "modular_results": modular_results,
        "interpretation": (
            "The diagonal core has origin radical over both tested finite "
            "fields.  Rational Groebner reconstruction and the preceding "
            "d=4 upper-triangularization gate remain open."
        ),
        "status": "finite-field experiment",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(
        "PASS d=4 diagonal core experiment: exact rank five; "
        "origin radicals over GF(101) and GF(1009)"
    )


if __name__ == "__main__":
    main()

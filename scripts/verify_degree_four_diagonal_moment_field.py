#!/usr/bin/env python3
"""Exact degree-two moment-field theorem on the diagonal quartic slice.

For C=diag(a_0,...,a_4), the apolar adjoint reverses the five coordinates.
This checker proves over QQ, using Singular's modular standard-basis
algorithm, that:

* mu_1,...,mu_5 form a homogeneous system of parameters of the polynomial
  ring QQ[a_0,...,a_4], with quotient length 1*2*3*4*5=120;
* the fiber of mu_1,...,mu_6 through (2,3,5,7,11) is exactly the two
  reduced reversal-related points (2,3,5,7,11) and (11,7,5,3,2).

Finiteness over the first five moments makes the coordinate ring finite
over the full moment algebra.  The two-point fiber bounds the generic
degree by two, while the nontrivial reversal fixing every moment bounds it
below by two.  Thus the diagonal full moment field is exactly the
reversal-fixed field.

This is a slice theorem, not a proof of the corresponding statement for
the full 22-dimensional invariant quotient.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_diagonal_moment_field.json"
)
POINT = (2, 3, 5, 7, 11)


def diagonal_moments(cutoff: int) -> tuple[tuple[sp.Symbol, ...], list[sp.Expr]]:
    variables = sp.symbols("a0:5")
    power: dict[int, sp.Expr] = {0: sp.Integer(1)}
    moments = []
    for order in range(1, cutoff + 1):
        current: dict[int, sp.Expr] = {}
        for degree, coefficient in power.items():
            for index, variable in enumerate(variables):
                current[degree + index] = (
                    current.get(degree + index, sp.Integer(0))
                    + coefficient * variable
                )
        power = current
        moments.append(
            sp.expand(
                sum(
                    factorial(index)
                    * factorial(4 * order - index)
                    * coefficient
                    for index, coefficient in power.items()
                )
            )
        )
    return variables, moments


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


def run_singular(code: str) -> list[str]:
    executable = shutil.which("Singular")
    if executable is None:
        raise SystemExit("Singular is required on PATH")
    completed = subprocess.run(
        [executable, "-q"],
        input=code,
        text=True,
        capture_output=True,
        check=True,
        timeout=300,
    )
    if completed.stderr.strip():
        raise AssertionError(completed.stderr)
    return [
        line.strip()
        for line in completed.stdout.splitlines()
        if line.strip()
    ]


def main() -> None:
    variables, moments = diagonal_moments(6)
    substitution = dict(zip(variables, POINT, strict=True))
    targets = [int(moment.subs(substitution)) for moment in moments]
    reversed_substitution = dict(
        zip(variables, reversed(POINT), strict=True)
    )
    assert [
        int(moment.subs(reversed_substitution)) for moment in moments
    ] == targets

    parameter_code = (
        'LIB "modstd.lib";\n'
        "ring r=0,(a0,a1,a2,a3,a4),dp;\n"
        + "ideal I="
        + ",".join(singular_expression(moment) for moment in moments[:5])
        + ";\n"
        + "ideal G=modStd(I);\n"
        + "vdim(G);\n"
        + "size(G);\n"
    )
    parameter_output = run_singular(parameter_code)
    parameter_length = int(parameter_output[0])
    parameter_basis_size = int(parameter_output[1])
    assert parameter_length == factorial(5) == 120

    fiber_generators = [
        "a0+a4-13",
        "9*a1+4*a4-71",
        "a2-5",
        "9*a3-4*a4-19",
        "(a4-2)*(a4-11)",
    ]
    fiber_code = (
        'LIB "modstd.lib";\n'
        "ring r=0,(a0,a1,a2,a3,a4),dp;\n"
        + "ideal I="
        + ",".join(
            f"({singular_expression(moment)})-({target})"
            for moment, target in zip(moments, targets, strict=True)
        )
        + ";\n"
        + "ideal J="
        + ",".join(fiber_generators)
        + ";\n"
        + "ideal GI=modStd(I);\n"
        + "ideal GJ=std(J);\n"
        + "vdim(GI);\n"
        + "size(GI);\n"
        + "reduce(I,GJ);\n"
        + "reduce(J,GI);\n"
    )
    fiber_output = run_singular(fiber_code)
    fiber_length = int(fiber_output[0])
    fiber_basis_size = int(fiber_output[1])
    assert fiber_length == 2
    assert fiber_basis_size == 5
    reduction_lines = fiber_output[2:]
    assert len(reduction_lines) == 11
    assert all(line.endswith("=0") for line in reduction_lines)

    payload = {
        "format": "degree-four-diagonal-moment-field-v1",
        "slice": "C=diag(a_0,a_1,a_2,a_3,a_4)",
        "apolar_involution": (
            "(a_0,a_1,a_2,a_3,a_4)"
            " -> (a_4,a_3,a_2,a_1,a_0)"
        ),
        "first_five_parameter_quotient_length": parameter_length,
        "first_five_standard_basis_size": parameter_basis_size,
        "test_point": list(POINT),
        "first_six_moment_values": [str(value) for value in targets],
        "fiber_ideal_generators": fiber_generators,
        "fiber_length": fiber_length,
        "fiber_standard_basis_size": fiber_basis_size,
        "fiber_points": [list(POINT), list(reversed(POINT))],
        "mutual_ideal_reductions_zero": True,
        "generic_full_moment_degree_on_slice": 2,
        "full_moment_field_on_slice": "reversal-fixed field",
        "scope_warning": (
            "slice theorem only; does not determine the degree on "
            "Frac(R_4)"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    print(
        "PASS diagonal quartic: first-five parameter quotient length",
        parameter_length,
    )
    print(
        "PASS diagonal quartic: first-six fiber is exactly the two "
        "reversal-related reduced points"
    )
    print(
        "PASS diagonal quartic: full moment field has exact generic degree two"
    )


if __name__ == "__main__":
    main()

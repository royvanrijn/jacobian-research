#!/usr/bin/env python3
"""Exact diagonal moment-field theorems for d=3,4,5.

On C=diag(a_0,...,a_d), apolar adjunction reverses the coordinates.  For
each requested degree this checker proves:

* the first d+1 moments form a homogeneous parameter system; and
* the first d+2 moments through one integral point cut out exactly that
  point and its coordinate reversal.

Parameter finiteness uses an exact finite-field homogeneous quotient of
length (d+1)!.  Since the corresponding projective family is proper over
Z, an empty projective fiber implies an empty characteristic-zero
projective fiber.  The d+1 homogeneous forms are then a regular sequence,
so their characteristic-zero quotient also has length (d+1)!.

The two-point fiber is checked directly over QQ.  An invertible linear
change sends the point and its reversal to s=+1 and s=-1 on the line
y_0=...=y_(d-1)=0.  Mutual ideal reductions prove that the transformed
moment ideal is exactly (y_0,...,y_(d-1),s^2-1).

These are slice theorems.  They do not determine the generic moment-field
degree on the full invariant quotients R_d.
"""

from __future__ import annotations

import argparse
import json
from math import factorial
from pathlib import Path

import sympy as sp

from verify_degree_four_diagonal_moment_field import (
    run_singular,
    singular_expression,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "completed_moment_diagonal_fields.json"
)
PARAMETER_PRIME = 32003
POINTS = {
    3: (2, 3, 5, 7),
    4: (2, 3, 5, 7, 11),
    5: (2, 3, 5, 7, 11, 13),
}


def diagonal_moments(
    d: int,
    cutoff: int,
) -> tuple[tuple[sp.Symbol, ...], list[sp.Expr]]:
    variables = sp.symbols(f"a0:{d + 1}")
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
                    * factorial(d * order - index)
                    * coefficient
                    for index, coefficient in power.items()
                )
            )
        )
    return variables, moments


def parameter_certificate(
    d: int,
    variables: tuple[sp.Symbol, ...],
    moments: list[sp.Expr],
) -> tuple[int, int]:
    names = ",".join(str(variable) for variable in variables)
    code = (
        f"ring r={PARAMETER_PRIME},({names}),dp;\n"
        + "ideal I="
        + ",".join(
            singular_expression(moment) for moment in moments[: d + 1]
        )
        + ";\n"
        + "ideal G=std(I);\n"
        + "vdim(G);\n"
        + "size(G);\n"
    )
    output = run_singular(code)
    quotient_length = int(output[0])
    basis_size = int(output[1])
    assert quotient_length == factorial(d + 1)
    return quotient_length, basis_size


def adapted_fiber_certificate(
    d: int,
    variables: tuple[sp.Symbol, ...],
    moments: list[sp.Expr],
    point: tuple[int, ...],
) -> dict[str, object]:
    substitution = dict(zip(variables, point, strict=True))
    targets = [int(moment.subs(substitution)) for moment in moments]
    reversed_point = tuple(reversed(point))
    reversed_substitution = dict(
        zip(variables, reversed_point, strict=True)
    )
    assert [
        int(moment.subs(reversed_substitution)) for moment in moments
    ] == targets

    midpoint = [
        sp.Rational(left + right, 2)
        for left, right in zip(point, reversed_point, strict=True)
    ]
    direction = [
        sp.Rational(left - right, 2)
        for left, right in zip(point, reversed_point, strict=True)
    ]
    assert direction[-1] != 0

    source_names = ",".join(str(variable) for variable in variables)
    y_names = [f"y{index}" for index in range(d)]
    target_names = ",".join(["s"] + y_names)
    s = sp.Symbol("s")
    map_entries = [
        singular_expression(midpoint[index] + direction[index] * s)
        + f"+{y_names[index]}"
        for index in range(d)
    ]
    map_entries.append(
        singular_expression(midpoint[d] + direction[d] * s)
    )
    source_ideal = ",".join(
        f"({singular_expression(moment)})-({target})"
        for moment, target in zip(moments, targets, strict=True)
    )
    expected_ideal = ",".join(y_names + ["s^2-1"])
    code = (
        'LIB "modstd.lib";\n'
        + f"ring ra=0,({source_names}),dp;\n"
        + f"ideal IA={source_ideal};\n"
        + f"ring rb=0,({target_names}),dp;\n"
        + f"map phi=ra,{','.join(map_entries)};\n"
        + "ideal I=phi(IA);\n"
        + f"ideal J={expected_ideal};\n"
        + "ideal GI=modStd(I);\n"
        + "ideal GJ=std(J);\n"
        + "vdim(GI);\n"
        + "size(GI);\n"
        + "size(GJ);\n"
        + "reduce(I,GJ);\n"
        + "reduce(J,GI);\n"
    )
    output = run_singular(code)
    fiber_length = int(output[0])
    fiber_basis_size = int(output[1])
    expected_basis_size = int(output[2])
    reduction_lines = output[3:]
    assert fiber_length == 2
    assert expected_basis_size == d + 1
    assert len(reduction_lines) == (d + 2) + (d + 1)
    assert all(line.endswith("=0") for line in reduction_lines)
    return {
        "test_point": list(point),
        "reversed_point": list(reversed_point),
        "first_d_plus_2_moment_values": [
            str(target) for target in targets
        ],
        "adapted_coordinates": {
            "midpoint": [str(value) for value in midpoint],
            "direction": [str(value) for value in direction],
            "map": map_entries,
            "expected_fiber_ideal": y_names + ["s^2-1"],
        },
        "fiber_length": fiber_length,
        "fiber_standard_basis_size": fiber_basis_size,
        "expected_ideal_standard_basis_size": expected_basis_size,
        "mutual_ideal_reductions_zero": True,
    }


def verify_degree(d: int) -> dict[str, object]:
    point = POINTS[d]
    variables, moments = diagonal_moments(d, d + 2)
    parameter_length, parameter_basis_size = parameter_certificate(
        d, variables, moments
    )
    fiber = adapted_fiber_certificate(
        d, variables, moments, point
    )
    print(
        f"PASS diagonal d={d}: first {d + 1} moments give "
        f"parameter length {parameter_length}"
    )
    print(
        f"PASS diagonal d={d}: first {d + 2} fiber is exactly "
        "the two reversal-related reduced points"
    )
    return {
        "slice": f"C=diag(a_0,...,a_{d})",
        "apolar_involution": "coordinate reversal",
        "parameter_prime": PARAMETER_PRIME,
        "finite_field_parameter_quotient_length": parameter_length,
        "finite_field_parameter_standard_basis_size": (
            parameter_basis_size
        ),
        "characteristic_zero_parameter_quotient_length": (
            parameter_length
        ),
        "parameter_lift_argument": (
            "homogeneous finite special fiber makes the projective "
            "characteristic-zero fiber empty by properness; regular "
            "sequence Hilbert series gives length (d+1)!"
        ),
        **fiber,
        "generic_full_moment_degree_on_slice": 2,
        "full_moment_field_on_slice": "reversal-fixed field",
        "scope_warning": (
            "slice theorem only; does not determine the degree on "
            f"Frac(R_{d})"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--degrees",
        nargs="+",
        type=int,
        choices=(3, 4, 5),
        default=(3, 4, 5),
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    payload = {
        "format": "completed-moment-diagonal-fields-v1",
        "degrees": {},
    }
    for d in arguments.degrees:
        payload["degrees"][str(d)] = verify_degree(d)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact degree-two moment fields on single-phase quartic slices.

For each nonzero phase h=1,2,3,4, choose a tau-even coefficient direction
B of matrix phase +h and a tau-even direction C of phase -h.  On

    diag(a_0,...,a_4) + b B + c C

the residual diagonal torus scales b and c inversely, so every moment
depends on them only through z=bc.  The quotient slice is therefore
Spec QQ[a_0,...,a_4,z], with tau reversing the a_i and fixing z.

There are ten coordinate pairs of tau-even directions: four for h=1,
four for h=2, and one each for h=3,4.  For every pair this checker proves
over QQ that:

* mu_1,...,mu_6 form a weighted homogeneous system of parameters, with
  quotient length (1*2*3*4*5*6)/2=360;
* the first-seven-moment fiber through
  (a_0,...,a_4,z)=(2,3,5,7,11,221) is exactly the two reduced points
  related by reversal.

Thus the full moment field on every displayed quotient slice is the
tau-fixed field and has generic degree two.  Since both finite-parameter
and reduced-fiber conditions are open in the direction coefficients, this
also proves the same statement for a nonempty Zariski-open family of
direction pairs in every phase.

These slice theorems do not determine the degree on the full
22-dimensional invariant quotient.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from functools import reduce
from math import factorial
from operator import mul
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_single_phase_moment_fields.json"
)
POINT = (2, 3, 5, 7, 11, 221)
VARIABLES = ("a0", "a1", "a2", "a3", "a4", "z")
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]
Direction = dict[tuple[int, int], int]


def tau_position(position: tuple[int, int]) -> tuple[int, int]:
    row, column = position
    return 4 - column, 4 - row


def tau_even_directions(phase: int, positive: bool) -> list[Direction]:
    positions = [
        (
            (index + phase, index)
            if positive
            else (index, index + phase)
        )
        for index in range(5 - phase)
    ]
    seen: set[int] = set()
    result = []
    for index, position in enumerate(positions):
        if index in seen:
            continue
        partner = 4 - phase - index
        seen.update((index, partner))
        direction = {position: 1}
        if partner != index:
            direction[positions[partner]] = (-1) ** phase
        for target, coefficient in direction.items():
            source = tau_position(target)
            expected = (-1) ** sum(target) * direction.get(source, 0)
            assert coefficient == expected
        result.append(direction)
    return result


def restricted_moments(
    positive: Direction,
    negative: Direction,
    cutoff: int,
) -> list[Polynomial]:
    # Intermediate variables are a_0,...,a_4,b,c.  After phase balance,
    # equal b- and c-exponents become the exponent of z=bc.
    terms = [(index, index, index, 1) for index in range(5)]
    terms.extend(
        (row, column, 5, coefficient)
        for (row, column), coefficient in positive.items()
    )
    terms.extend(
        (row, column, 6, coefficient)
        for (row, column), coefficient in negative.items()
    )
    state: dict[tuple[int, int, tuple[int, ...]], int] = {
        (0, 0, (0,) * 7): 1
    }
    moments = []
    for order in range(1, cutoff + 1):
        updated: dict[tuple[int, int, tuple[int, ...]], int] = {}
        for (left, right, exponents), value in state.items():
            for delta_left, delta_right, variable, coefficient in terms:
                new_exponents = list(exponents)
                new_exponents[variable] += 1
                key = (
                    left + delta_left,
                    right + delta_right,
                    tuple(new_exponents),
                )
                updated[key] = updated.get(key, 0) + value * coefficient
        state = {key: value for key, value in updated.items() if value}

        moment: Polynomial = {}
        for (left, right, exponents), value in state.items():
            if left != right or exponents[5] != exponents[6]:
                continue
            quotient_exponents = exponents[:5] + (exponents[5],)
            moment[quotient_exponents] = (
                moment.get(quotient_exponents, 0)
                + factorial(left)
                * factorial(4 * order - left)
                * value
            )
        moments.append(
            {
                exponents: value
                for exponents, value in moment.items()
                if value
            }
        )
    return moments


def polynomial_string(polynomial: Polynomial) -> str:
    terms = []
    for exponents, coefficient in polynomial.items():
        sign = "+" if coefficient >= 0 and terms else ""
        absolute = abs(coefficient)
        factors = []
        if absolute != 1 or not any(exponents):
            factors.append(str(absolute))
        for variable, exponent in zip(VARIABLES, exponents, strict=True):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors) if factors else "0"
        terms.append(sign + ("" if coefficient >= 0 else "-") + monomial)
    return "".join(terms) if terms else "0"


def evaluate(polynomial: Polynomial) -> int:
    return sum(
        coefficient
        * reduce(
            mul,
            (
                value**exponent
                for value, exponent in zip(POINT, exponents, strict=True)
            ),
            1,
        )
        for exponents, coefficient in polynomial.items()
    )


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


def direction_payload(direction: Direction) -> list[dict[str, int]]:
    return [
        {"row": row, "column": column, "coefficient": coefficient}
        for (row, column), coefficient in sorted(direction.items())
    ]


def verify_slice(
    phase: int,
    positive_index: int,
    positive: Direction,
    negative_index: int,
    negative: Direction,
) -> dict[str, object]:
    moments = restricted_moments(positive, negative, 7)
    moment_strings = [polynomial_string(moment) for moment in moments]
    targets = [evaluate(moment) for moment in moments]

    parameter_code = (
        'LIB "modstd.lib";\n'
        "ring r=0,(a0,a1,a2,a3,a4,z),dp;\n"
        + "ideal I="
        + ",".join(moment_strings[:6])
        + ";\n"
        + "ideal G=modStd(I);\n"
        + "dim(G);\n"
        + "vdim(G);\n"
        + "size(G);\n"
    )
    parameter_output = run_singular(parameter_code)
    assert int(parameter_output[0]) == 0
    parameter_length = int(parameter_output[1])
    parameter_basis_size = int(parameter_output[2])
    assert parameter_length == 360

    fiber_generators = [
        "a0+a4-13",
        "9*a1+4*a4-71",
        "a2-5",
        "9*a3-4*a4-19",
        "z-221",
        "(a4-2)*(a4-11)",
    ]
    fiber_code = (
        'LIB "modstd.lib";\n'
        "ring r=0,(a0,a1,a2,a3,a4,z),dp;\n"
        + "ideal I="
        + ",".join(
            f"({moment})-({target})"
            for moment, target in zip(
                moment_strings,
                targets,
                strict=True,
            )
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
    assert fiber_basis_size == 6
    reduction_lines = fiber_output[2:]
    assert len(reduction_lines) == 13
    assert all(line.endswith("=0") for line in reduction_lines)

    print(
        "PASS single-phase slice:",
        f"h={phase}",
        f"positive={positive_index}",
        f"negative={negative_index}",
        "parameter length 360, fiber length 2",
    )
    return {
        "phase": phase,
        "positive_direction_index": positive_index,
        "positive_direction": direction_payload(positive),
        "negative_direction_index": negative_index,
        "negative_direction": direction_payload(negative),
        "moment_term_counts_orders_1_through_7": [
            len(moment) for moment in moments
        ],
        "first_seven_moment_values": [str(value) for value in targets],
        "first_six_parameter_quotient_length": parameter_length,
        "first_six_standard_basis_size": parameter_basis_size,
        "first_seven_fiber_length": fiber_length,
        "first_seven_fiber_standard_basis_size": fiber_basis_size,
        "mutual_fiber_ideal_reductions_zero": True,
    }


def main() -> None:
    slices = []
    for phase in range(1, 5):
        positive_directions = tau_even_directions(phase, True)
        negative_directions = tau_even_directions(phase, False)
        for positive_index, positive in enumerate(positive_directions):
            for negative_index, negative in enumerate(negative_directions):
                slices.append(
                    verify_slice(
                        phase,
                        positive_index,
                        positive,
                        negative_index,
                        negative,
                    )
                )
    assert len(slices) == 10

    payload = {
        "format": "degree-four-single-phase-moment-fields-v1",
        "quotient_coordinates": ["a0", "a1", "a2", "a3", "a4", "z=bc"],
        "apolar_involution": (
            "(a0,a1,a2,a3,a4,z)"
            " -> (a4,a3,a2,a1,a0,z)"
        ),
        "test_point": list(POINT),
        "slice_count": len(slices),
        "slices": slices,
        "generic_full_moment_degree_on_each_slice": 2,
        "full_moment_field_on_each_slice": "apolar-fixed field",
        "open_family_consequence": (
            "for every phase h=1,2,3,4, a nonempty Zariski-open family "
            "of tau-even positive/negative direction pairs has the same "
            "degree-two fixed-field property"
        ),
        "scope_warning": (
            "single-phase quotient slices only; no determination of the "
            "degree on Frac(R_4)"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()

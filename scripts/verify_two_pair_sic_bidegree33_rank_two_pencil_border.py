#!/usr/bin/env python3
"""Relative logarithmic-Jacobian border basis on a rank-two pencil.

This computes over F_p(s), so the result is exact along the chosen
one-parameter family but remains modular evidence rather than a universal
parameter-space calculation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import warnings

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

from sympy import Poly, factor_list, gcd, lcm, symbols
from sympy.utilities.exceptions import SymPyDeprecationWarning

from verify_two_pair_sic_bidegree33_rank_two_holonomic_probe import POINTS


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_pencil_border.json"
)
PRIMES = (1_000_003, 1_000_033, 1_000_037)
EXPECTED_BORDER = ("u5", "u4t", "u3t2", "u2t4", "ut5", "t6")
EXPECTED_DENOMINATOR_DEGREES = (74, 74, 88, 94, 74, 74)
EXPECTED_COMPONENT_CLASSIFICATIONS = {
    1_000_003: {
        "extra_degree_14:relative_length_change": 4,
    },
    1_000_033: {
        "common_degree_74:stable_length_profile_chart_failure": 1,
        "extra_degree_14:relative_length_change": 4,
        "extra_degree_20:relative_length_change": 1,
    },
    1_000_037: {
        "common_degree_74:stable_length_profile_chart_failure": 3,
        "extra_degree_14:relative_length_change": 2,
    },
}
S = symbols("s")


def linear(left: int, right: int, value: int | None = None) -> str:
    parameter = "s" if value is None else str(value)
    return f"({left}+({right})*{parameter})"


def q_expression(value: int | None = None) -> str:
    """Return Q=u^3 Phi(1,u,t,(1-t)/u) on U0+sU1,W0+sW1."""

    (u0, w0), (u1, w1) = POINTS
    summands = []
    for inner in range(2):
        dual = "+".join(
            f"{linear(u0[i][inner], u1[i][inner], value)}*u^{3-i}"
            for i in range(4)
        )
        coordinate = "+".join(
            (
                f"{linear(w0[inner][j], w1[inner][j], value)}"
                f"*u^{j}*t^{j}*(1-t)^{3-j}"
            )
            for j in range(4)
        )
        summands.append(f"({dual})*({coordinate})")
    return "+".join(summands)


def run_singular(code: str) -> str:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=code,
        text=True,
        capture_output=True,
        timeout=180,
        check=True,
    )
    if "?" in completed.stdout or completed.stderr:
        raise RuntimeError(completed.stdout + completed.stderr)
    return completed.stdout


def pencil_border_basis(prime: int) -> str:
    return run_singular(
        f"""
LIB "elim.lib";
ring r=({prime},s),(u,t),dp;
poly Q={q_expression()};
poly A=u*diff(Q,u)-3*Q;
poly C=t*(1-t)*diff(Q,t);
list saturation=sat_with_exp(ideal(A,C),ideal(u));
ideal relative=saturation[1];
if(saturation[2]!=6){{ERROR("saturation exponent");}}
if(vdim(relative)!=18){{ERROR("relative length");}}
option(redSB);
ideal border=std(relative);
print("SATURATION_EXPONENT");
print(saturation[2]);
print("RELATIVE_LENGTH");
print(vdim(relative));
print("BASIS_SIZE");
print(size(border));
int i;
poly normalized;
poly cleared;
for(i=1;i<=size(border);i++){{
  normalized=border[i]/leadcoef(border[i]);
  cleared=cleardenom(normalized);
  print("RELATION");
  print(i);
  print("BORDER");
  print(leadmonom(normalized));
  print("DENOMINATOR");
  print(leadcoef(cleared));
  print("TERMS");
  print(size(normalized));
}}
"""
    )


def modular_matrix_rank(matrix: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                index
                for index in range(row, len(work))
                if work[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [value * inverse % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or work[index][column] == 0:
                continue
            scale = work[index][column]
            work[index] = [
                (left - scale * right) % prime
                for left, right in zip(work[index], work[row], strict=True)
            ]
        row += 1
    return row


def factor_chart_rank(prime: int, value: int) -> int:
    (u0, w0), (u1, w1) = POINTS
    u = [
        [
            (u0[row][column] + value * u1[row][column]) % prime
            for column in range(2)
        ]
        for row in range(4)
    ]
    w = [
        [
            (w0[row][column] + value * w1[row][column]) % prime
            for column in range(4)
        ]
        for row in range(2)
    ]
    product = [
        [
            sum(u[row][inner] * w[inner][column] for inner in range(2))
            % prime
            for column in range(4)
        ]
        for row in range(4)
    ]
    return modular_matrix_rank(product, prime)


def specialized_relative_fiber(prime: int, value: int) -> dict[str, object]:
    output = run_singular(
        f"""
LIB "elim.lib";
ring r={prime},(u,t),dp;
poly Q={q_expression(value)};
poly A=u*diff(Q,u)-3*Q;
poly C=t*(1-t)*diff(Q,t);
list saturation=sat_with_exp(ideal(A,C),ideal(u));
ideal relative=std(saturation[1]);
ideal endpoint0=std(sat(ideal(A,t),ideal(u)));
ideal endpoint1=std(sat(ideal(A,t-1),ideal(u)));
ideal interior=std(
  sat(sat(ideal(A,diff(Q,t)),ideal(u)),ideal(t*(1-t)))
);
print("SATURATION_EXPONENT");
print(saturation[2]);
print("RELATIVE_LENGTH");
print(vdim(relative));
print("ENDPOINT_ZERO_LENGTH");
print(vdim(endpoint0));
print("ENDPOINT_ONE_LENGTH");
print(vdim(endpoint1));
print("INTERIOR_LENGTH");
print(vdim(interior));
"""
    )
    lines = [
        line.strip()
        for line in output.splitlines()
        if line.strip() and not line.startswith("//")
    ]
    expected_labels = (
        "SATURATION_EXPONENT",
        "RELATIVE_LENGTH",
        "ENDPOINT_ZERO_LENGTH",
        "ENDPOINT_ONE_LENGTH",
        "INTERIOR_LENGTH",
    )
    values = {}
    cursor = 0
    for label in expected_labels:
        assert lines[cursor] == label
        values[label] = int(lines[cursor + 1])
        cursor += 2
    assert cursor == len(lines)

    rank = factor_chart_rank(prime, value)
    profile = (
        values["ENDPOINT_ZERO_LENGTH"],
        values["ENDPOINT_ONE_LENGTH"],
        values["INTERIOR_LENGTH"],
    )
    relative_length = values["RELATIVE_LENGTH"]
    if rank < 2:
        classification = "coefficient_rank_drop"
    elif relative_length != 18:
        classification = "relative_length_change"
    elif profile != (2, 2, 14):
        classification = "boundary_interior_redistribution"
    else:
        classification = "stable_length_profile_chart_failure"
    return {
        "parameter_value": value,
        "coefficient_matrix_rank": rank,
        "saturation_exponent": values["SATURATION_EXPONENT"],
        "relative_length": relative_length,
        "endpoint_t0_length": profile[0],
        "endpoint_t1_length": profile[1],
        "interior_length": profile[2],
        "classification": classification,
    }


def parse_polynomial(text: str, prime: int) -> Poly:
    expression = text.strip().removeprefix("(").removesuffix(")")
    expression = re.sub(r"s([0-9]+)", r"s**\1", expression)
    expression = re.sub(r"(?<=[0-9])s", "*s", expression)
    return Poly(expression, S, modulus=prime).monic()


def parse_output(output: str, prime: int) -> dict[str, object]:
    lines = [
        line.strip()
        for line in output.splitlines()
        if line.strip() and not line.startswith("//")
    ]
    assert lines[:6] == [
        "SATURATION_EXPONENT",
        "6",
        "RELATIVE_LENGTH",
        "18",
        "BASIS_SIZE",
        "6",
    ]
    relations = []
    denominators = []
    cursor = 6
    for index in range(1, 7):
        assert lines[cursor : cursor + 2] == ["RELATION", str(index)]
        assert lines[cursor + 2] == "BORDER"
        border = lines[cursor + 3].replace("^", "").replace("*", "")
        assert lines[cursor + 4] == "DENOMINATOR"
        denominator = parse_polynomial(lines[cursor + 5], prime)
        assert lines[cursor + 6] == "TERMS"
        terms = int(lines[cursor + 7])
        relations.append(
            {
                "border_monomial": border,
                "terms": terms,
                "denominator_degree": denominator.degree(),
            }
        )
        denominators.append(denominator)
        cursor += 8
    assert cursor == len(lines)
    assert {record["border_monomial"] for record in relations} == set(
        EXPECTED_BORDER
    )
    assert [record["terms"] for record in relations] == [19] * 6
    assert tuple(
        record["denominator_degree"] for record in relations
    ) == EXPECTED_DENOMINATOR_DEGREES

    unique: list[Poly] = []
    for denominator in denominators:
        if all(denominator != previous for previous in unique):
            unique.append(denominator)
    assert [item.degree() for item in unique] == [74, 88, 94]
    gcd_degrees = [
        gcd(unique[left], unique[right]).degree()
        for left in range(3)
        for right in range(left)
    ]
    assert gcd_degrees == [74, 74, 74]
    extra_14, remainder_14 = unique[1].div(unique[0])
    extra_20, remainder_20 = unique[2].div(unique[0])
    assert remainder_14.is_zero and remainder_20.is_zero
    assert (extra_14.degree(), extra_20.degree()) == (14, 20)
    assert gcd(extra_14, extra_20).degree() == 0
    exceptional = lcm(lcm(unique[0], unique[1]), unique[2]).monic()
    assert exceptional.degree() == 108

    factor_degrees = []
    for denominator in unique:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SymPyDeprecationWarning)
            factors = factor_list(denominator, modulus=prime)[1]
        degrees = [
            [factor.degree(), multiplicity]
            for factor, multiplicity in factors
        ]
        factor_degrees.append(degrees)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SymPyDeprecationWarning)
        exceptional_factors = factor_list(
            exceptional,
            modulus=prime,
        )[1]
    exceptional_factor_degrees = [
        [factor.degree(), multiplicity]
        for factor, multiplicity in exceptional_factors
    ]
    linear_roots = sorted(
        {
            int(
                -factor.nth(0) * pow(int(factor.nth(1)), -1, prime)
            )
            % prime
            for factor, _multiplicity in exceptional_factors
            if factor.degree() == 1
        }
    )
    specialized_fibers = [
        specialized_relative_fiber(prime, value) for value in linear_roots
    ]
    for fiber in specialized_fibers:
        value = int(fiber["parameter_value"])
        components = []
        if int(unique[0].eval(value)) % prime == 0:
            components.append("common_degree_74")
        if int(extra_14.eval(value)) % prime == 0:
            components.append("extra_degree_14")
        if int(extra_20.eval(value)) % prime == 0:
            components.append("extra_degree_20")
        assert len(components) == 1
        fiber["denominator_component"] = components[0]
    classification_counts: dict[str, int] = {}
    component_classification_counts: dict[str, int] = {}
    for fiber in specialized_fibers:
        classification = str(fiber["classification"])
        classification_counts[classification] = (
            classification_counts.get(classification, 0) + 1
        )
        component_classification = (
            f"{fiber['denominator_component']}:{classification}"
        )
        component_classification_counts[component_classification] = (
            component_classification_counts.get(component_classification, 0)
            + 1
        )
    assert (
        component_classification_counts
        == EXPECTED_COMPONENT_CLASSIFICATIONS[prime]
    )

    return {
        "prime": prime,
        "saturation_exponent": 6,
        "relative_length": 18,
        "relations": relations,
        "unique_denominator_degrees": [74, 88, 94],
        "pairwise_denominator_gcd_degrees": gcd_degrees,
        "common_core_degree": 74,
        "extra_factor_degrees": [14, 20],
        "extra_factor_gcd_degree": 0,
        "exceptional_polynomial_degree": exceptional.degree(),
        "exceptional_polynomial_squarefree": (
            gcd(exceptional, exceptional.diff()).degree() == 0
        ),
        "denominator_factor_degree_multiplicities": factor_degrees,
        "exceptional_factor_degree_multiplicities": (
            exceptional_factor_degrees
        ),
        "distinct_exceptional_roots_in_base_field": (
            len(linear_roots)
        ),
        "base_field_exceptional_fibers": specialized_fibers,
        "base_field_exceptional_classification_counts": classification_counts,
        "base_field_component_classification_counts": (
            component_classification_counts
        ),
        "exceptional_polynomial_coefficients_low_to_high": [
            int(exceptional.nth(index)) % prime
            for index in range(exceptional.degree() + 1)
        ],
    }


def main() -> None:
    records = []
    for prime in PRIMES:
        output = pencil_border_basis(prime)
        records.append(parse_output(output, prime))
        print(
            f"PASS prime {prime}: length 18, six 19-term border "
            "relations, exceptional degree 108"
        )
        counts = records[-1]["base_field_exceptional_classification_counts"]
        print(
            f"PASS prime {prime}: classified "
            f"{records[-1]['distinct_exceptional_roots_in_base_field']} "
            f"base-field exceptional fibers as {counts}"
        )
    aggregate_profiles: dict[tuple[str, tuple[int, int, int]], int] = {}
    for record in records:
        for fiber in record["base_field_exceptional_fibers"]:
            key = (
                str(fiber["denominator_component"]),
                (
                    int(fiber["endpoint_t0_length"]),
                    int(fiber["endpoint_t1_length"]),
                    int(fiber["interior_length"]),
                ),
            )
            aggregate_profiles[key] = aggregate_profiles.get(key, 0) + 1
    assert aggregate_profiles == {
        ("common_degree_74", (2, 2, 14)): 4,
        ("extra_degree_14", (1, 2, 14)): 4,
        ("extra_degree_14", (2, 1, 14)): 4,
        ("extra_degree_14", (2, 2, 13)): 2,
        ("extra_degree_20", (2, 2, 13)): 1,
    }
    OUTPUT.write_text(
        json.dumps(
            {
                "format": "two-pair-sic-bidegree33-rank-two-pencil-border-v1",
                "status": (
                    "exact computations over three finite rational-function "
                    "fields on one rank-two pencil; not a characteristic-zero "
                    "universal parameter-space certificate"
                ),
                "pencil": "U(s)=U0+sU1, W(s)=W0+sW1",
                "records": records,
                "interpretation": {
                    "generic_pencil_chart": (
                        "the six monic logarithmic-Jacobian border relations "
                        "are regular away from a degree-108 polynomial"
                    ),
                    "exceptional_locus": (
                        "roots of the chart denominator require separate "
                        "specialization; every base-field root is classified "
                        "by coefficient rank, total relative length, and its "
                        "2+2+14 endpoint/interior profile"
                    ),
                },
            },
            indent=2,
        )
        + "\n"
    )
    print("PASS denominator degrees 74, 88, 94 have common gcd degree 74")
    print("PASS quotient degrees 14 and 20 are coprime")
    print("PASS border-basis chart denominator has degree 108")
    print("PASS every base-field exceptional root has an exact fiber classification")
    print("PASS result is modular one-pencil evidence, not a universal theorem")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

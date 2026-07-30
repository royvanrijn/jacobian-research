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
S = symbols("s")


def linear(left: int, right: int) -> str:
    return f"({left}+({right})*s)"


def q_expression() -> str:
    """Return Q=u^3 Phi(1,u,t,(1-t)/u) on U0+sU1,W0+sW1."""

    (u0, w0), (u1, w1) = POINTS
    summands = []
    for inner in range(2):
        dual = "+".join(
            f"{linear(u0[i][inner], u1[i][inner])}*u^{3-i}"
            for i in range(4)
        )
        coordinate = "+".join(
            (
                f"{linear(w0[inner][j], w1[inner][j])}"
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
    distinct_linear_factors = sum(
        1
        for factor, _multiplicity in exceptional_factors
        if factor.degree() == 1
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
            distinct_linear_factors
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
                        "specialization and may include coordinate-chart "
                        "artifacts"
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
    print("PASS result is modular one-pencil evidence, not a universal theorem")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

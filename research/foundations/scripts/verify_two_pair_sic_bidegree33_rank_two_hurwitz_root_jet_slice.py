#!/usr/bin/env python3
"""Classify a five-parameter exact-rank-two cubic Hurwitz slice.

On

    F=(1+x)B(y)+x^2*(lambda+x)D(y),

specialize b0=0 and b1=1.  Thus

    B=1+a1*y+a2*y^2+a3*y^3,
    D=y+b2*y^2+d3*y^3,
    d3=-1-a1/3-lambda*b2/3.

The coefficient matrix has rank exactly two everywhere.  Exact QQ msolve
calculations show that mu_2,...,mu_6 give a zero-dimensional quotient of
length 687 and that mu_2,...,mu_7 have support equal to the single rational
point

    (lambda,a1,a2,a3,b2)=(1,-4,5,-2,-2),

with local length 26.  A contraction-preserving flag change puts that point
in the strict chamber i>j.  This proves all-order pure vanishing, the
recurrence nu_(m+1)=0, a nonzero low-degree mixed value, and the sharp
eventual cutoff m>e for a degree-e multiplier.

This is a coefficient slice, not a diagonal-SL_2 orbit chart.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import sympy as sp

from research_two_pair_sic_bidegree33_rank_two_hurwitz import (
    exact_moment_polynomials,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_hurwitz_root_jet_slice.json"
)
ORDERS = tuple(range(2, 8))
VARIABLES = ("lambda", "a1", "a2", "a3", "b2")
POINT = (1, -4, 5, -2, -2)
SpecializedExponent = tuple[int, int, int, int, int]
SpecializedPolynomial = dict[SpecializedExponent, int]


def specialize(polynomial: dict[tuple[int, ...], int]) -> SpecializedPolynomial:
    """Set b0=0 and b1=1, retaining lambda,a1,a2,a3,b2."""

    answer: SpecializedPolynomial = {}
    for exponent, coefficient in polynomial.items():
        # Generic parameter order is lambda,a1,a2,a3,b0,b1,b2.
        if exponent[4]:
            continue
        specialized = (
            exponent[0],
            exponent[1],
            exponent[2],
            exponent[3],
            exponent[6],
        )
        answer[specialized] = answer.get(specialized, 0) + coefficient
    return {
        exponent: coefficient
        for exponent, coefficient in answer.items()
        if coefficient
    }


def primitive(polynomial: SpecializedPolynomial) -> tuple[SpecializedPolynomial, int]:
    content = 0
    for coefficient in polynomial.values():
        content = gcd(content, abs(coefficient))
    if content == 0:
        return {}, 0
    return (
        {
            exponent: coefficient // content
            for exponent, coefficient in polynomial.items()
        },
        content,
    )


def polynomial_string(polynomial: SpecializedPolynomial) -> str:
    polynomial, _ = primitive(polynomial)
    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponent]
        factors: list[str] = []
        if coefficient != 1 or not any(exponent):
            factors.append(str(coefficient))
        for name, power in zip(VARIABLES, exponent):
            if power == 1:
                factors.append(name)
            elif power:
                factors.append(f"{name}^{power}")
        terms.append("*".join(factors) or "1")
    return "+".join(terms).replace("+-", "-") or "0"


def evaluate(polynomial: SpecializedPolynomial) -> int:
    return sum(
        coefficient
        * int(sp.prod(value**power for value, power in zip(POINT, exponent)))
        for exponent, coefficient in polynomial.items()
    )


def msolve_source(
    moments: dict[int, SpecializedPolynomial], cutoff: int
) -> str:
    generators = [polynomial_string(moments[order]) for order in range(2, cutoff + 1)]
    return ",".join(VARIABLES) + "\n0\n" + ",\n".join(generators) + "\n"


def solve_exact(source: str) -> dict[str, object]:
    executable = shutil.which("msolve")
    if executable is None:
        raise RuntimeError("msolve is required")
    with tempfile.TemporaryDirectory(prefix="sic33-hurwitz-root-jet-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(source, encoding="utf-8")
        completed = subprocess.run(
            [
                executable,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-l",
                "2",
                "-v",
                "1",
                "--random-seed",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        result = output_path.read_text(encoding="utf-8") if output_path.exists() else ""
    if completed.returncode != 0:
        raise AssertionError(completed.stderr[-2000:])
    diagnostics = completed.stdout + "\n" + completed.stderr

    def diagnostic(pattern: str) -> int | None:
        match = re.search(pattern + r"\s+(\d+)", diagnostics)
        return int(match.group(1)) if match else None

    compact = " ".join(result.split())
    return {
        "returncode": completed.returncode,
        "result": compact,
        "result_sha256": sha256(compact.encode()).hexdigest(),
        "degree": diagnostic(r"degree of ideal"),
        "eliminant_degree": diagnostic(r"deg\. elim\. pol\."),
        "squarefree_eliminant_degree": diagnostic(
            r"deg\. sqfr\. elim\. pol\."
        ),
    }


def all_order_certificate() -> dict[str, object]:
    x, y = sp.symbols("x y")
    capital_w, capital_v, capital_z, capital_y = sp.symbols("W V Z Y")
    wp, vp, zp, yp = sp.symbols("Wp Vp Zp Yp")
    b_polynomial = 1 - 4 * y + 5 * y**2 - 2 * y**3
    d_polynomial = y - 2 * y**2 + y**3
    affine = sp.Poly(
        sp.expand((1 + x) * b_polynomial + x**2 * (1 + x) * d_polynomial),
        x,
        y,
    )
    homogeneous = sum(
        affine.coeff_monomial(x**i * y**j)
        * capital_w ** (3 - i)
        * capital_v**i
        * capital_z ** (3 - j)
        * capital_y**j
        for i in range(4)
        for j in range(4)
    )
    transformed = sp.expand(
        homogeneous.subs(
            {
                capital_w: wp,
                capital_v: vp - wp,
                capital_z: zp + yp,
                capital_y: yp,
            }
        )
    )
    expected = vp * zp**2 * (vp**2 * yp - 2 * vp * wp * yp + wp**2 * zp)
    if sp.expand(transformed - expected) != 0:
        raise AssertionError("fixed-flag factorization failed")

    entries = {(1, 0): 1, (2, 1): -2, (3, 1): 1}
    if min(i - j for i, j in entries) != 1:
        raise AssertionError("the transformed support is not strictly one-sided")
    u, t = sp.symbols("u t")
    period = sp.expand(
        sum(
            value
            * u ** (target_index - dual_index)
            * t**target_index
            * (1 - t) ** (3 - target_index)
            for (dual_index, target_index), value in entries.items()
        )
    )
    multiplier_numerator = u * t
    normalized_mixed = sp.integrate(
        sp.expand(multiplier_numerator * period).coeff(u, 0), (t, 0, 1)
    )
    raw_mixed = factorial(5) * normalized_mixed
    if normalized_mixed != Fraction(-1, 60) or raw_mixed != -2:
        raise AssertionError("unexpected degree-one mixed value")
    return {
        "point_form": str(sp.factor(affine.as_expr())),
        "flag_parameter": -1,
        "fixed_flag_factorization": str(expected),
        "fixed_flag_entries": {
            f"{i},{j}": value for (i, j), value in entries.items()
        },
        "period_laurent_polynomial": str(period),
        "creative_telescoping_recurrence": "nu_(m+1)=0 for m>=0",
        "initial_condition": "nu_1=0 by strict u-valuation",
        "degree_one_multiplier_period_numerator": "u*t",
        "degree_one_mixed_value_at_m1": int(raw_mixed),
        "mixed_tail": "every degree-e multiplier vanishes for m>e",
    }


def main() -> None:
    moments = {
        order: specialize(polynomial)
        for order, polynomial in exact_moment_polynomials(ORDERS).items()
    }
    profiles = {}
    for order, polynomial in moments.items():
        serialized = polynomial_string(polynomial)
        _, content = primitive(polynomial)
        profiles[str(order)] = {
            "terms": len(polynomial),
            "integer_content": str(content),
            "primitive_sha256": sha256(serialized.encode()).hexdigest(),
        }

    prefix_source = msolve_source(moments, 6)
    prefix = solve_exact(prefix_source)
    if (
        prefix["degree"] != 687
        or prefix["eliminant_degree"] != 668
        or prefix["squarefree_eliminant_degree"] != 658
    ):
        raise AssertionError(f"unexpected mu_2,...,mu_6 prefix: {prefix}")

    classification_source = msolve_source(moments, 7)
    classification = solve_exact(classification_source)
    expected_point_text = (
        "[[[1, 1], [-4, -4], [5, 5], [-2, -2], [-2, -2]]]"
    )
    if (
        classification["degree"] != 26
        or classification["squarefree_eliminant_degree"] != 1
        or expected_point_text not in classification["result"]
        or any(evaluate(moment) for moment in moments.values())
    ):
        raise AssertionError(
            f"unexpected mu_2,...,mu_7 classification: {classification}"
        )

    all_order = all_order_certificate()
    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-hurwitz-root-jet-slice-v2",
        "field": "characteristic zero",
        "chart": "F=(1+x)B(y)+x^2*(lambda+x)D(y)",
        "slice": {
            "B": "1+a1*y+a2*y^2+a3*y^3",
            "D": "y+b2*y^2+(-1-a1/3-lambda*b2/3)*y^3",
            "generic_coordinates": "b0=0 and b1=1",
            "coefficient_rank": 2,
        },
        "scaled_form": "3F",
        "moment_profiles": profiles,
        "prefix_mu2_through_mu6": {
            **{key: value for key, value in prefix.items() if key != "result"},
            "msolve_input_sha256": sha256(prefix_source.encode()).hexdigest(),
            "conclusion": "zero-dimensional quotient of length 687",
        },
        "classification_mu2_through_mu7": {
            **{
                key: value
                for key, value in classification.items()
                if key != "result"
            },
            "msolve_input_sha256": sha256(
                classification_source.encode()
            ).hexdigest(),
            "support": [list(POINT)],
            "local_length": 26,
        },
        "all_order_component_certificate": all_order,
        "conclusion": (
            "the only all-order point on the slice is fixed-flag one-sided "
            "and SIC-safe"
        ),
        "scope": (
            "exact five-parameter coefficient slice only; this is not a "
            "diagonal-SL_2 orbit chart or a generic Hurwitz classification"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS mu_2,...,mu_6 give the exact length-687 prefix")
    print("PASS mu_7 leaves one rational point of local length 26")
    print("PASS the point is fixed-flag one-sided with recurrence nu_(m+1)=0")
    print("PASS mixed value -2 at m=1 and degree-e cutoff m>e")


if __name__ == "__main__":
    main()

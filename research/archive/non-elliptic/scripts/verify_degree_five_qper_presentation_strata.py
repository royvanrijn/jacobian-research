#!/usr/bin/env python3
"""Audit the first coherent presentation charts for the hbar^5 period line.

The stable 16-monomial support gives 41 nonzero bounded conditions of generic
rank 15.  A fixed 15-row subsystem presents its generic dual-cocycle line.
After primitive row clearing, normalizing the first functional coordinate
produces the irreducible degree-34 determinant stored in the generated
artifact below.

This checker verifies the artifact and tests exact finite-field points on
that determinant.  At all of them the complete 41-condition fiber still has
rank 15, so the determinant is a normalization-chart boundary rather than
evidence for an extra fiberwise cocycle.  The characteristic-zero statement
that a second maximal minor removes this factor is checked fraction-freely by
``compute_degree_five_fifth_order_function_field.py --presentation-audit``.
"""

from __future__ import annotations

from functools import reduce
from math import gcd
from pathlib import Path
import shutil
import subprocess
import tempfile

import sympy as sp
from sympy.polys.domains import GF

from explore_degree_five_fifth_order_period_constraints import (
    constraint_pivots,
    constraint_system,
)


ROOT = Path(__file__).resolve().parents[1]
PIVOT_FACTOR = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_pivot_D34.sing"
)
ALTERNATE_FACTOR = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_pivot_E35.sing"
)
PRESENTATION = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_15x16_presentation.sing"
)


def main() -> None:
    a_symbol, tau_symbol = sp.symbols("a tau")
    expression = sp.sympify(
        PIVOT_FACTOR.read_text().strip().replace("^", "**")
    )
    polynomial = sp.Poly(expression, a_symbol, tau_symbol)
    assert len(polynomial.terms()) == 274
    assert polynomial.total_degree() == 34
    assert polynomial.degree(a_symbol) == 22
    assert polynomial.degree(tau_symbol) == 12
    assert reduce(gcd, (abs(int(value)) for value in polynomial.coeffs())) == 1
    alternate_expression = sp.sympify(
        ALTERNATE_FACTOR.read_text().strip().replace("^", "**")
    )
    alternate_polynomial = sp.Poly(
        alternate_expression,
        a_symbol,
        tau_symbol,
    )
    assert len(alternate_polynomial.terms()) == 283
    assert alternate_polynomial.total_degree() == 35
    assert alternate_polynomial.degree(a_symbol) == 23
    assert alternate_polynomial.degree(tau_symbol) == 12
    assert (
        reduce(
            gcd,
            (abs(int(value)) for value in alternate_polynomial.coeffs()),
        )
        == 1
    )

    terms = [
        (exponents, int(coefficient))
        for exponents, coefficient in polynomial.terms()
    ]

    def evaluate(a_value: int, tau_value: int, prime: int) -> int:
        return (
            sum(
                coefficient
                * pow(a_value, a_degree, prime)
                * pow(tau_value, tau_degree, prime)
                for (
                    a_degree,
                    tau_degree,
                ), coefficient in terms
            )
            % prime
        )

    samples = {
        101: ((2, 22), (7, 42), (8, 19), (9, 80), (10, 21)),
        103: ((3, 84), (6, 77), (7, 37)),
    }
    for prime, points in samples.items():
        field = GF(prime)
        for a_integer, tau_integer in points:
            assert evaluate(a_integer, tau_integer, prime) == 0
            a = field(a_integer)
            tau = field(tau_integer)
            assert a and a + field.one and tau

            # Auxiliary denominator of the chosen sparse lower-kernel basis.
            kernel_chart = (
                field(4) * a**3 * tau**2
                - field(24) * a**3 * tau
                - field(72) * a**3
                + field(8) * a**2 * tau**2
                - field(54) * a**2 * tau
                - field(216) * a**2
                + field(4) * a * tau**2
                - field(30) * a * tau
                - field(246) * a
                - field(105)
            )
            assert kernel_chart

            labels, constraints, lower_kernel_dimension = constraint_system(
                field,
                a,
                tau,
            )
            pivots = constraint_pivots(constraints, field)
            assert lower_kernel_dimension == 41
            assert len(labels) == 41
            assert len(pivots) == 15
            assert 16 - len(pivots) == 1

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for modular intersection audit")
    d_text = PIVOT_FACTOR.read_text().strip()
    e_text = ALTERNATE_FACTOR.read_text().strip()
    kernel_chart_text = (
        "4*a^3*tau^2-24*a^3*tau-72*a^3"
        "+8*a^2*tau^2-54*a^2*tau-216*a^2"
        "+4*a*tau^2-30*a*tau-246*a-105"
    )
    intersection_profiles = {}
    for prime in (31991, 32003, 65521):
        program = "\n".join(
            [
                'LIB "elim.lib";',
                f"ring r={prime},(a,tau),dp;",
                f"poly D={d_text};",
                f"poly E={e_text};",
                f"poly H={kernel_chart_text};",
                "ideal I=D,E;",
                "ideal H_boundary=H;",
                "list saturation_H=sat(I,H_boundary);",
                "ideal chart_interior=saturation_H[1];",
                "ideal G=std(chart_interior);",
                'print("DIM="+string(dim(G)));',
                'print("LENGTH="+string(vdim(G)));',
                "quit;",
            ]
        )
        with tempfile.TemporaryDirectory(
            prefix="degree-five-qper-chart-intersection-",
        ) as directory:
            path = Path(directory) / "intersection.sing"
            path.write_text(program)
            result = subprocess.run(
                [singular, "-q", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        values = dict(
            line.split("=", 1)
            for line in result.stdout.splitlines()
            if "=" in line
        )
        if "DIM" not in values:
            raise AssertionError(result.stdout + result.stderr)
        intersection_profiles[prime] = (
            int(values["DIM"]),
            int(values["LENGTH"]),
        )
    matrix_line = next(
        line
        for line in PRESENTATION.read_text().splitlines()
        if line.startswith("matrix M[15][16]=")
    )
    fitting_profiles = {}
    for prime in (31991, 32003, 65521):
        program = "\n".join(
            [
                'LIB "elim.lib";',
                f"ring r={prime},(a,tau),dp;",
                matrix_line,
                f"poly D={d_text};",
                f"poly E={e_text};",
                f"poly H={kernel_chart_text};",
                "ideal maximal_minors=minor(M,15);",
                "ideal a_boundary=a;",
                "ideal a1_boundary=a+1;",
                "list saturation_a=sat(maximal_minors,a_boundary);",
                "ideal away_a=saturation_a[1];",
                "list saturation_a1=sat(away_a,a1_boundary);",
                "ideal away_a1=saturation_a1[1];",
                "ideal H_boundary=H;",
                "list saturation_H=sat(away_a1,H_boundary);",
                "ideal interior=saturation_H[1];",
                "ideal G=std(interior);",
                'print("MINORS="+string(size(maximal_minors)));',
                'print("DIM="+string(dim(G)));',
                'print("LENGTH="+string(vdim(G)));',
                'print("GENERATORS="+string(size(G)));',
                'print("FIRST_DEGREE="+string(deg(G[1])));',
                'print("FIRST_TERMS="+string(size(G[1])));',
                'print("GCD_DEGREE="+string(deg(gcd(D,E))));',
                'print("D_REMAINDER_TERMS="+string(size(reduce(D,G))));',
                'print("E_REMAINDER_TERMS="+string(size(reduce(E,G))));',
                "quit;",
            ]
        )
        with tempfile.TemporaryDirectory(
            prefix="degree-five-qper-fitting-",
        ) as directory:
            path = Path(directory) / "fitting.sing"
            path.write_text(program)
            result = subprocess.run(
                [singular, "-q", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        values = dict(
            line.split("=", 1)
            for line in result.stdout.splitlines()
            if "=" in line
        )
        fitting_profiles[prime] = (
            int(values["MINORS"]),
            int(values["DIM"]),
            int(values["LENGTH"]),
            int(values["GENERATORS"]),
            int(values["FIRST_DEGREE"]),
            int(values["FIRST_TERMS"]),
            int(values["GCD_DEGREE"]),
            int(values["D_REMAINDER_TERMS"]),
            int(values["E_REMAINDER_TERMS"]),
        )
    assert fitting_profiles == {
        31991: (16, 0, 218, 21, 20, 219, 0, 0, 0),
        32003: (16, 0, 218, 21, 20, 219, 0, 0, 0),
        65521: (16, 0, 218, 21, 20, 219, 0, 0, 0),
    }

    print(
        "PASS: the primitive first-coordinate pivot factor has "
        "274 terms and total degree 34"
    )
    print(
        "PASS: eight exact points on D34 retain complete constraint "
        "rank 15 and one dual cocycle"
    )
    print(f"MODULAR_CHART_INTERSECTIONS={intersection_profiles}")
    print(f"MODULAR_FITTING_PROFILES={fitting_profiles}")
    print(
        "SCOPE: the modular fibers audit the alternate-chart phenomenon; "
        "the exact two-minor gcd is computed by --presentation-audit"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit the first biprojective graph calibration for the quartic seed.

The characteristic-zero part computes two independent rational linear
sections of the rational projective extension.  After removing its base
ideal, both sections are reduced of length 28.  The finite-field part
computes the full saturated biprojective graph over F_32003 and its four
projective multidegrees.

The modular graph calculation is deliberately reported separately: it is a
good-prime calibration for the equinormalization program, not by itself a
characteristic-zero flatness or normalization theorem.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from math import gcd, lcm
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    explicit_seed,
)


SINGULAR = shutil.which("Singular")
PRIME = 32003

X0, X1, X2, X3 = sp.symbols("X0 X1 X2 X3")
u, v = sp.symbols("u v")


@dataclass(frozen=True)
class SectionCard:
    name: str
    source_hyperplane: tuple[int, int, int, int]
    source_dehomogenizer: tuple[int, int, int, int]
    target_hyperplanes: tuple[
        tuple[int, int, int, int], tuple[int, int, int, int]
    ]


CARDS = (
    SectionCard(
        "A",
        (1, 2, 3, 5),
        (37, 41, 43, 47),
        ((1, 3, 7, 11), (13, 17, 19, 23)),
    ),
    SectionCard(
        "B",
        (2, 5, 7, 11),
        (31, 37, 41, 43),
        ((3, 11, 17, 29), (5, 13, 23, 31)),
    ),
)


def quartic_projective_coordinates() -> tuple[sp.Expr, ...]:
    """Return a common degree-twelve projective extension of [1:F]."""

    primitive = explicit_seed(4)
    assert sp.factor(primitive - w**2 * (w - 5) * (w - 1) / 4) == 0
    mapping = WeightedSeedModel(sp.diff(primitive, w)).mapping()
    degrees = tuple(sp.Poly(component, x, y, z).total_degree() for component in mapping)
    assert degrees == (12, 11, 4)

    coordinates: list[sp.Expr] = [X0**12]
    for component in mapping:
        polynomial = sp.Poly(component, x, y, z)
        homogenized = polynomial.homogenize(X0).as_expr().subs(
            {x: X1, y: X2, z: X3}
        )
        coordinates.append(
            sp.expand(X0 ** (12 - polynomial.total_degree()) * homogenized)
        )

    assert all(
        sp.Poly(coordinate, X0, X1, X2, X3).total_degree() == 12
        for coordinate in coordinates
    )
    assert tuple(
        len(sp.Poly(coordinate, X0, X1, X2, X3).terms())
        for coordinate in coordinates
    ) == (1, 16, 14, 3)
    return tuple(coordinates)


def verify_boundary_triangle(coordinates: tuple[sp.Expr, ...]) -> None:
    """Verify the reduced base triangle and its generic line orders."""

    assert sp.factor(
        coordinates[1].subs(X0, 0)
        - sp.Rational(3, 4) * X1**6 * X2**4 * X3**2
    ) == 0
    x0_valuations = tuple(
        min(monomial[0] for monomial, _ in sp.Poly(item, X0, X1, X2, X3).terms())
        for item in coordinates
    )
    assert x0_valuations == (12, 0, 1, 8)

    # Along L_i=(X0,Xi), the order is the minimum X0+Xi exponent at the
    # generic point of that line.  Together with the preceding restriction,
    # this proves that the reduced base support is the coordinate triangle
    # V(X0,X1*X2*X3) in the source hyperplane at infinity.
    expected_orders = {
        1: (12, 6, 7, 11),
        2: (12, 4, 4, 8),
        3: (12, 2, 3, 9),
    }
    variables = (X0, X1, X2, X3)
    for index, expected in expected_orders.items():
        actual = tuple(
            min(
                monomial[0] + monomial[index]
                for monomial, _ in sp.Poly(item, *variables).terms()
            )
            for item in coordinates
        )
        assert actual == expected


def primitive_integral(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    denominator = 1
    for coefficient in polynomial.coeffs():
        denominator = lcm(denominator, int(coefficient.q))
    integral_coefficients = [
        int(coefficient * denominator) for coefficient in polynomial.coeffs()
    ]
    content = 0
    for coefficient in integral_coefficients:
        content = gcd(content, abs(coefficient))
    if content == 0:
        raise AssertionError("zero polynomial has no primitive normalization")
    return sp.expand(polynomial.as_expr() * denominator / content)


def singular_polynomial(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def run_singular(program: str, timeout: int = 300) -> str:
    if SINGULAR is None:
        raise SystemExit("Singular is required for the quartic graph audit")
    completed = subprocess.run(
        [SINGULAR, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0 or "   ? " in completed.stdout:
        raise AssertionError(
            "Singular failed:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed.stdout


def linear_form(
    coefficients: tuple[int, int, int, int],
    variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
) -> sp.Expr:
    return sum(
        coefficient * variable
        for coefficient, variable in zip(coefficients, variables, strict=True)
    )


def exact_section_program(coordinates: tuple[sp.Expr, ...]) -> str:
    blocks: list[str] = [
        'LIB "primdec.lib";',
        "option(redSB);",
        "ring r=0,(u,v),dp;",
    ]
    for card in CARDS:
        source_hyperplane = linear_form(
            card.source_hyperplane, (X0, X1, X2, X3)
        )
        source_dehomogenizer = linear_form(
            card.source_dehomogenizer, (X0, X1, X2, X3)
        )
        solution = sp.solve(
            (source_hyperplane, source_dehomogenizer - 1),
            (X0, X3),
            dict=True,
        )
        assert len(solution) == 1
        substitution = {
            X0: solution[0][X0].subs({X1: u, X2: v}),
            X1: u,
            X2: v,
            X3: solution[0][X3].subs({X1: u, X2: v}),
        }
        restricted = tuple(sp.expand(item.subs(substitution)) for item in coordinates)
        target_sections = tuple(
            sp.expand(
                linear_form(target_hyperplane, restricted)  # type: ignore[arg-type]
            )
            for target_hyperplane in card.target_hyperplanes
        )

        base_generators = tuple(
            primitive_integral(item, (u, v)) for item in restricted
        )
        section_generators = tuple(
            primitive_integral(item, (u, v)) for item in target_sections
        )
        name = card.name
        blocks.extend(
            (
                f"ideal B{name}="
                + ",".join(singular_polynomial(item) for item in base_generators)
                + ";",
                f"ideal J{name}="
                + ",".join(singular_polynomial(item) for item in section_generators)
                + ";",
                f"list L{name}=sat(J{name},B{name});",
                f"ideal Q{name}=std(L{name}[1]);",
                f"ideal R{name}=std(radical(Q{name}));",
                f'print("CARD_{name} "+string(dim(Q{name}))+" "'
                f'+string(vdim(Q{name}))+" "+string(vdim(R{name})));',
            )
        )
    return "\n".join(blocks) + "\n"


def verify_exact_sections(coordinates: tuple[sp.Expr, ...]) -> None:
    output = run_singular(exact_section_program(coordinates))
    for card in CARDS:
        match = re.search(rf"^CARD_{card.name} (\d+) (\d+) (\d+)$", output, re.M)
        assert match is not None, output
        assert tuple(map(int, match.groups())) == (0, 28, 28)


def modular_graph_program(coordinates: tuple[sp.Expr, ...]) -> str:
    integral = tuple(
        primitive_integral(item, (X0, X1, X2, X3)) for item in coordinates
    )
    coordinate_strings = tuple(singular_polynomial(item) for item in integral)
    relations = (
        f"Y1*({coordinate_strings[0]})-Y0*({coordinate_strings[1]})",
        f"Y2*({coordinate_strings[0]})-Y0*({coordinate_strings[2]})",
        f"Y3*({coordinate_strings[0]})-Y0*({coordinate_strings[3]})",
    )

    source_hyperplanes = (
        "X0+2*X1+3*X2+5*X3",
        "7*X0+11*X1+13*X2+17*X3",
        "19*X0+23*X1+29*X2+31*X3",
    )
    target_hyperplanes = (
        "Y0+3*Y1+7*Y2+11*Y3",
        "13*Y0+17*Y1+19*Y2+23*Y3",
        "29*Y0+31*Y1+37*Y2+41*Y3",
    )
    source_dehomogenizer = "37*X0+41*X1+43*X2+47*X3-1"
    target_dehomogenizer = "43*Y0+47*Y1+53*Y2+59*Y3-1"

    blocks = [
        'LIB "elim.lib";',
        "option(redSB);",
        f"ring r={PRIME},(X0,X1,X2,X3,Y0,Y1,Y2,Y3),dp;",
        "ideal I=" + ",".join(relations) + ";",
        "ideal S=X0*Y0;",
        "list L=sat(I,S);",
        "ideal G=std(L[1]);",
        "ideal M=minbase(G);",
        'print("GRAPH "+string(dim(G))+" "+string(size(G))'
        '+" "+string(size(M)));',
    ]
    for target_count in range(4):
        forms = (
            source_hyperplanes[: 3 - target_count]
            + target_hyperplanes[:target_count]
            + (source_dehomogenizer, target_dehomogenizer)
        )
        blocks.extend(
            (
                f"ideal Q{target_count}=G," + ",".join(forms) + ";",
                f"Q{target_count}=std(Q{target_count});",
                f'print("DELTA_{target_count} "+string(dim(Q{target_count}))'
                f'+" "+string(vdim(Q{target_count})));',
            )
        )
    return "\n".join(blocks) + "\n"


def verify_modular_graph(coordinates: tuple[sp.Expr, ...]) -> None:
    output = run_singular(modular_graph_program(coordinates))
    graph = re.search(r"^GRAPH (\d+) (\d+) (\d+)$", output, re.M)
    assert graph is not None, output
    assert tuple(map(int, graph.groups())) == (5, 181, 38)

    expected = (1, 12, 28, 4)
    for target_count, length in enumerate(expected):
        match = re.search(rf"^DELTA_{target_count} (\d+) (\d+)$", output, re.M)
        assert match is not None, output
        assert tuple(map(int, match.groups())) == (0, length)


def main() -> None:
    coordinates = quartic_projective_coordinates()
    verify_boundary_triangle(coordinates)
    verify_exact_sections(coordinates)
    verify_modular_graph(coordinates)
    print(
        "PASS quartic biprojective graph: base triangle orders (6,4,2); "
        "two exact reduced length-28 sections; F_32003 graph dimension 3 "
        "and multidegrees (1,12,28,4)"
    )


if __name__ == "__main__":
    main()

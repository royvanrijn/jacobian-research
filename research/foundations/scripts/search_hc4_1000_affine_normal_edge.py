#!/usr/bin/env python3
"""Native-Singular search for the final chart-1000 affine-normal edge.

The exact boundary reductions leave only

    f=5*b^3/6-delta*b^2/2,
    g=gamma*a+g1*b+g2*b^2+g3*b^3,
    h=h1*b+h2*b^2+h3*b^3

up to Hessian-irrelevant affine terms.  This script adjoins every quartic
term containing at least two normal variables c,d.  Equivalently, it uses
arbitrary quadratic boundary data V_cc,V_cd,V_dd, arbitrary affine normal
third derivatives, and arbitrary constant normal fourth derivatives.

It constructs the graph determinant natively in Singular after restricting
to selected source hyperplanes, extracts every nonconstant coefficient,
and saturates by the common constant determinant.  A unit ideal over Q is
an exact obstruction; modular runs are screens only.
"""

from __future__ import annotations

import argparse
from itertools import product
import re
import runpy
import subprocess

import sympy as sp


native = runpy.run_path("scripts/verify_hc4_multislice_cubic_obstruction.py")
graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
graph_polynomial = native["graph_polynomial"]
hessian_entry = native["hessian_entry"]

h_variables = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
source_names = ("X", "Y", "W", "D")

a, b, c, d = sp.symbols("aa bb cc dd")
q_placeholders = sp.symbols("qq0:4")

delta, gamma = sp.symbols("delta gamma")
g_coefficients = sp.symbols("gb1:4")
h_coefficients = sp.symbols("hb1:4")

f = sp.Rational(5, 6) * b**3 - delta * b**2 / 2
g = gamma * a + sum(
    g_coefficients[degree - 1] * b**degree for degree in range(1, 4)
)
h = sum(
    h_coefficients[degree - 1] * b**degree for degree in range(1, 4)
)

quadratic_monomials = [
    a**a_degree * b**b_degree
    for a_degree, b_degree in product(range(3), repeat=2)
    if a_degree + b_degree <= 2
]
U_coefficients = sp.symbols("u0:6")
W_coefficients = sp.symbols("w0:6")
Z_coefficients = sp.symbols("zz0:6")
U = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        U_coefficients, quadratic_monomials, strict=True
    )
)
mixed_W = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        W_coefficients, quadratic_monomials, strict=True
    )
)
Z = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        Z_coefficients, quadratic_monomials, strict=True
    )
)

affine_monomials = (sp.Integer(1), a, b)
normal_third_coefficients = tuple(
    sp.symbols(f"n{index}_0:3") for index in range(4)
)
normal_thirds = tuple(
    sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            coefficients, affine_monomials, strict=True
        )
    )
    for coefficients in normal_third_coefficients
)
normal_fourths = sp.symbols("r0:5")

V = (
    f
    + c * g
    + d * h
    + c**2 * U / 2
    + c * d * mixed_W
    + d**2 * Z / 2
    + c**3 * normal_thirds[0] / 6
    + c**2 * d * normal_thirds[1] / 2
    + c * d**2 * normal_thirds[2] / 2
    + d**3 * normal_thirds[3] / 6
    + c**4 * normal_fourths[0] / 24
    + c**3 * d * normal_fourths[1] / 6
    + c**2 * d**2 * normal_fourths[2] / 4
    + c * d**3 * normal_fourths[3] / 6
    + d**4 * normal_fourths[4] / 24
)

coordinate_matrix = sp.Matrix(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 1],
    ]
)
hessian_adapted = sp.hessian(V, (a, b, c, d))
hessian_original = (
    coordinate_matrix.inv().T
    * hessian_adapted
    * coordinate_matrix.inv()
)
hessian_original = hessian_original.subs(
    {
        a: q_placeholders[0],
        b: q_placeholders[1],
        c: q_placeholders[2],
        d: q_placeholders[3] - q_placeholders[1],
    },
    simultaneous=True,
)

parameters = (
    (delta, gamma)
    + g_coefficients
    + h_coefficients
    + U_coefficients
    + W_coefficients
    + Z_coefficients
    + tuple(
        coefficient
        for coefficients in normal_third_coefficients
        for coefficient in coefficients
    )
    + normal_fourths
)


def singular_hessian_entry(expression: sp.Expr) -> str:
    """Print rational factors without ambiguous power/division parsing."""

    output = hessian_entry(sp.expand(expression))
    return re.sub(r"/([0-9]+)", r"*(1/\1)", output)


def singular_script(
    slices: tuple[str, ...], characteristic: int, algorithm: str
) -> str:
    mask = (1, 0, 0, 0)
    q = [M[index] if mask[index] else Q[index] for index in range(4)]
    m = [-Q[index] if mask[index] else M[index] for index in range(4)]
    saturation = "z_saturation"
    variables = source_names + tuple(map(str, parameters)) + (saturation,)
    lines = [f"ring ring_hc4={characteristic},({','.join(variables)}),dp;"]
    for index, expression in enumerate(q):
        lines.append(f"poly q{index}={graph_polynomial(expression)};")
    for index, expression in enumerate(m):
        lines.append(f"poly m{index}={graph_polynomial(expression)};")
    lines.append(
        "matrix H[4][4]="
        + ",".join(
            singular_hessian_entry(hessian_original[row, column])
            for row in range(4)
            for column in range(4)
        )
        + ";"
    )
    lines.append(
        "matrix Jq[4][4]="
        + ",".join(
            f"diff(q{row},{variable})"
            for row in range(4)
            for variable in source_names
        )
        + ";"
    )
    lines.append(
        "matrix Jm[4][4]="
        + ",".join(
            f"diff(m{row},{variable})"
            for row in range(4)
            for variable in source_names
        )
        + ";"
    )
    lines.extend(("matrix A=Jm+H*Jq;", "ideal I=0;", "poly constant_term=0;"))
    for index, zero_variable in enumerate(slices):
        coefficient_variables = "*".join(
            variable for variable in source_names if variable != zero_variable
        )
        lines.extend(
            (
                f"matrix A{index}=subst(A,{zero_variable},0);",
                f"poly determinant{index}=det(A{index});",
                f"matrix coefficients{index}=coef(determinant{index},{coefficient_variables});",
                f'print(\"SLICE {zero_variable} \"+string(ncols(coefficients{index})));',
                f"for (int i{index}=1; i{index}<=ncols(coefficients{index}); i{index}++) {{",
                f"  if (coefficients{index}[1,i{index}]==1) {{",
                f"    if ({index}==0) {{ constant_term=coefficients{index}[2,i{index}]; }}",
                "  } else {",
                f"    I=I,coefficients{index}[2,i{index}];",
                "  }",
                "}",
            )
        )
    lines.append(f"I=I,{saturation}*constant_term-1;")
    lines.extend(("option(redSB);", f"ideal G={algorithm}(I);", "reduce(1,G);", "size(G);"))
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slices",
        nargs="+",
        choices=source_names,
        default=("Y",),
    )
    parser.add_argument("--characteristic", type=int, default=0)
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="slimgb")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()

    script = singular_script(
        tuple(args.slices), args.characteristic, args.algorithm
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=args.timeout,
        check=False,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()

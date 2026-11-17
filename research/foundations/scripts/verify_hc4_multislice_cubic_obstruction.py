#!/usr/bin/env python3
"""Exact native-Singular slice ideals for bounded-degree PC(2) shears.

For a selected graph chart, take the unrestricted potential

    V = V_2 + ... + V_d

with every homogeneous coefficient free.  Degree three has ten quadratic
and twenty cubic coefficients; degree four adds 35 quartic coefficients.

Global constancy of

    det d(m + grad(V)(q))

implies constancy on every selected coordinate hyperplane.  This checker
substitutes each hyperplane into the Jacobian matrix *before* taking its
determinant, extracts all nonconstant source coefficients natively in
Singular, and saturates by the common constant term.

The final certificates use direct ``slimgb`` reduction over Q.  Prime-field
runs remain available as screens, but are explicitly not reported as
characteristic-zero proofs.

For the irreducible characteristic caustic branch in chart 1000, the checker
can either impose the five coefficients of

    det Hess(V_4(a,b,0,b)) = 0

or cover the resulting rational normal curve by either pure-power chart.
The latter eliminates five quartic coefficients before determinant
construction.
"""

from __future__ import annotations

import argparse
from itertools import product
import runpy
import subprocess

import sympy as sp


jet = runpy.run_path("scripts/verify_hc4_quadratic_cubic_jet_obstruction.py")
graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")

old_source = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
new_source = sp.symbols("X Y W D")
source_names = ("X", "Y", "W", "D")
source_rename = dict(zip(old_source, new_source, strict=True))

quadratic_parameters = tuple(jet["k"])
cubic_parameters = tuple(jet["c"])
z_saturation = "z_saturation"

q_placeholders = sp.symbols("qq0:4")
cubic_hessian = jet["cubic_hessian"].subs(
    dict(zip(jet["u"], q_placeholders, strict=True)), simultaneous=True
)


def potential_data(
    maximum_degree: int,
) -> tuple[tuple[sp.Symbol, ...], sp.Matrix]:
    """Return all potential parameters and the nonconstant Hessian part."""

    parameters = quadratic_parameters + cubic_parameters
    hessian = cubic_hessian
    u = jet["u"]
    for degree in range(4, maximum_degree + 1):
        exponents = [
            powers
            for powers in product(range(degree + 1), repeat=4)
            if sum(powers) == degree
        ]
        coefficients = sp.symbols(f"v{degree}_0:{len(exponents)}")
        homogeneous_part = sum(
            coefficient * sp.prod(u[index] ** powers[index] for index in range(4))
            for coefficient, powers in zip(coefficients, exponents, strict=True)
        )
        parameters += tuple(coefficients)
        hessian += sp.hessian(homogeneous_part, u).subs(
            dict(zip(u, q_placeholders, strict=True)),
            simultaneous=True,
        )
    return parameters, hessian


def chart_1000_binary_top_hessian_relations() -> tuple[sp.Expr, ...]:
    """Return Hess(f_4(a,b,0,b))=0 coefficient relations."""

    u = jet["u"]
    exponents = [
        powers
        for powers in product(range(5), repeat=4)
        if sum(powers) == 4
    ]
    coefficients = sp.symbols(f"v4_0:{len(exponents)}")
    quartic_part = sum(
        coefficient * sp.prod(u[index] ** powers[index] for index in range(4))
        for coefficient, powers in zip(coefficients, exponents, strict=True)
    )
    boundary_a, boundary_b = sp.symbols("boundary_a boundary_b")
    binary_top = sp.expand(
        quartic_part.subs(
            {
                u[0]: boundary_a,
                u[1]: boundary_b,
                u[2]: 0,
                u[3]: boundary_b,
            },
            simultaneous=True,
        )
    )
    binary_hessian = sp.det(sp.hessian(binary_top, (boundary_a, boundary_b)))
    polynomial = sp.Poly(sp.expand(binary_hessian), boundary_a, boundary_b)
    relations = (
        polynomial.coeff_monomial(
            boundary_a ** (4 - index) * boundary_b**index
        )
        for index in range(5)
    )
    return tuple(sp.Poly(relation).primitive()[1].as_expr() for relation in relations)


def chart_1000_pure_power_substitution(
    normalization: str,
) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    """Eliminate five quartic coefficients on one pure-power chart."""

    coefficients = sp.symbols("v4_0:35")
    top_lambda, top_slope = sp.symbols("top_lambda top_slope")
    pivots = (
        coefficients[34],
        coefficients[33],
        coefficients[30],
        coefficients[24],
        coefficients[14],
    )
    if normalization == "a":
        binary_coefficients = (
            top_lambda,
            4 * top_lambda * top_slope,
            6 * top_lambda * top_slope**2,
            4 * top_lambda * top_slope**3,
            top_lambda * top_slope**4,
        )
    elif normalization == "b":
        binary_coefficients = (
            top_lambda * top_slope**4,
            4 * top_lambda * top_slope**3,
            6 * top_lambda * top_slope**2,
            4 * top_lambda * top_slope,
            top_lambda,
        )
    else:
        raise ValueError("pure-power normalization must be a or b")

    substitution = {
        coefficients[34]: binary_coefficients[0],
        coefficients[33]: binary_coefficients[1] - coefficients[31],
        coefficients[30]: (
            binary_coefficients[2] - coefficients[25] - coefficients[28]
        ),
        coefficients[24]: (
            binary_coefficients[3]
            - coefficients[15]
            - coefficients[19]
            - coefficients[22]
        ),
        coefficients[14]: (
            binary_coefficients[4]
            - coefficients[0]
            - coefficients[5]
            - coefficients[9]
            - coefficients[12]
        ),
    }
    boundary_coefficients = (
        coefficients[34],
        coefficients[31] + coefficients[33],
        coefficients[25] + coefficients[28] + coefficients[30],
        (
            coefficients[15]
            + coefficients[19]
            + coefficients[22]
            + coefficients[24]
        ),
        (
            coefficients[0]
            + coefficients[5]
            + coefficients[9]
            + coefficients[12]
            + coefficients[14]
        ),
    )
    assert all(
        sp.factor(coefficient.subs(substitution) - expected) == 0
        for coefficient, expected in zip(
            boundary_coefficients, binary_coefficients, strict=True
        )
    )
    return substitution, pivots, (top_lambda, top_slope)


def graph_polynomial(expression: sp.Expr) -> str:
    """Format a rational graph polynomial safely for Singular."""

    polynomial = sp.Poly(sp.expand(expression.xreplace(source_rename)), *new_source)
    terms = []
    for exponents, coefficient in polynomial.terms():
        monomial = (
            "*".join(
                (
                    f"{source_names[index]}^{exponent}"
                    if exponent > 1
                    else source_names[index]
                )
                for index, exponent in enumerate(exponents)
                if exponent
            )
            or "1"
        )
        terms.append(f"({int(coefficient.p)}/{int(coefficient.q)})*{monomial}")
    return "+".join(terms).replace("+-", "-") or "0"


def hessian_entry(expression: sp.Expr) -> str:
    output = str(expression).replace("**", "^")
    for index in range(4):
        output = output.replace(f"qq{index}", f"(q{index})")
    return output


def singular_script(
    chart: str,
    slices: tuple[str, ...],
    characteristic: int,
    algorithm: str,
    potential_degree: int,
    chart_1000_zero_hessian_top: bool,
    chart_1000_pure_power_top: str | None,
) -> str:
    mask = tuple(int(bit) for bit in chart)
    q = [M[index] if mask[index] else Q[index] for index in range(4)]
    m = [-Q[index] if mask[index] else M[index] for index in range(4)]

    potential_parameters, potential_hessian = potential_data(potential_degree)
    if chart_1000_pure_power_top is not None:
        substitution, eliminated, introduced = (
            chart_1000_pure_power_substitution(chart_1000_pure_power_top)
        )
        potential_hessian = potential_hessian.subs(
            substitution, simultaneous=True
        )
        potential_parameters = tuple(
            parameter
            for parameter in potential_parameters
            if parameter not in eliminated
        ) + introduced
    variables = source_names + tuple(map(str, potential_parameters)) + (z_saturation,)
    lines = []
    lines.append(f"ring r={characteristic},({','.join(variables)}),dp;")
    for index, expression in enumerate(q):
        lines.append(f"poly q{index}={graph_polynomial(expression)};")
    for index, expression in enumerate(m):
        lines.append(f"poly m{index}={graph_polynomial(expression)};")

    K = jet["K"]
    lines.append(
        "matrix K[4][4]="
        + ",".join(str(K[row, column]) for row in range(4) for column in range(4))
        + ";"
    )
    lines.append(
        "matrix H[4][4]="
        + ",".join(
            hessian_entry(potential_hessian[row, column])
            for row in range(4)
            for column in range(4)
        )
        + ";"
    )
    lines.append(
        "matrix Jq[4][4]="
        + ",".join(
            f"diff(q{row},{variable})" for row in range(4) for variable in source_names
        )
        + ";"
    )
    lines.append(
        "matrix Jm[4][4]="
        + ",".join(
            f"diff(m{row},{variable})" for row in range(4) for variable in source_names
        )
        + ";"
    )
    lines.extend(
        (
            "matrix A=Jm+(K+H)*Jq;",
            "ideal I=0;",
            "poly constant_term=0;",
        )
    )
    if chart_1000_zero_hessian_top:
        for relation in chart_1000_binary_top_hessian_relations():
            lines.append(f"I=I,{str(relation).replace('**', '^')};")

    for index, zero_variable in enumerate(slices):
        coefficient_variables = "*".join(
            variable for variable in source_names if variable != zero_variable
        )
        lines.extend(
            (
                f"matrix A{index}=subst(A,{zero_variable},0);",
                f"poly F{index}=det(A{index});",
                f"matrix C{index}=coef(F{index},{coefficient_variables});",
                f'print("SLICE {zero_variable} "+string(ncols(C{index})));',
                f"for (int i{index}=1; i{index}<=ncols(C{index}); i{index}++) {{",
                f"  if (C{index}[1,i{index}]==1) {{",
                (f"    if ({index}==0) " f"{{ constant_term=C{index}[2,i{index}]; }}"),
                "  } else {",
                f"    I=I,C{index}[2,i{index}];",
                "  }",
                "}",
            )
        )

    lines.append(f"I=I,{z_saturation}*constant_term-1;")
    lines.extend(("option(redSB);", f"ideal G={algorithm}(I);"))
    lines.extend(("reduce(1,G);", "size(G);"))
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chart", required=True)
    parser.add_argument(
        "--slices",
        nargs="+",
        choices=source_names,
        required=True,
    )
    parser.add_argument("--characteristic", type=int, default=0)
    parser.add_argument(
        "--potential-degree",
        type=int,
        choices=(3, 4),
        default=3,
    )
    parser.add_argument(
        "--chart-1000-zero-hessian-top",
        action="store_true",
        help=(
            "impose Hess(V_4(a,b,0,b))=0 on the irreducible "
            "characteristic caustic branch"
        ),
    )
    parser.add_argument(
        "--chart-1000-pure-power-top",
        choices=("a", "b"),
        help=(
            "parameterize V_4(a,b,0,b) as lambda*(a+t*b)^4 "
            "or lambda*(t*a+b)^4 and eliminate five coefficients"
        ),
    )
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb"),
        default="slimgb",
    )
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    if len(args.chart) != 4 or any(bit not in "01" for bit in args.chart):
        raise ValueError("chart must be a four-bit mask")
    if len(set(args.slices)) != len(args.slices):
        raise ValueError("slice variables must be distinct")
    if args.chart_1000_zero_hessian_top and (
        args.chart != "1000" or args.potential_degree != 4
    ):
        raise ValueError(
            "--chart-1000-zero-hessian-top requires chart 1000 and degree 4"
        )
    if args.chart_1000_pure_power_top is not None and (
        args.chart != "1000" or args.potential_degree != 4
    ):
        raise ValueError(
            "--chart-1000-pure-power-top requires chart 1000 and degree 4"
        )
    if (
        args.chart_1000_zero_hessian_top
        and args.chart_1000_pure_power_top is not None
    ):
        raise ValueError("choose either zero-Hessian equations or a pure-power chart")
    script = singular_script(
        args.chart,
        tuple(args.slices),
        args.characteristic,
        args.algorithm,
        args.potential_degree,
        args.chart_1000_zero_hessian_top,
        args.chart_1000_pure_power_top,
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=args.timeout,
        check=True,
    )
    output = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if any("?" in line for line in output):
        raise RuntimeError("\n".join(output))
    slice_lines = output[: len(args.slices)]
    expected_prefixes = [f"SLICE {variable} " for variable in args.slices]
    assert all(
        line.startswith(prefix)
        for line, prefix in zip(slice_lines, expected_prefixes, strict=True)
    )
    assert output[-2:] == ["0", "1"], output[-10:]

    print(
        f"PASS: chart {args.chart}, potential degree <= "
        f"{args.potential_degree}, slice columns are "
        + ", ".join(line.removeprefix("SLICE ") for line in slice_lines)
    )
    if args.chart_1000_zero_hessian_top:
        print(
            "PASS: imposed the five coefficients of "
            "Hess(V_4(a,b,0,b))=0"
        )
    if args.chart_1000_pure_power_top is not None:
        print(
            "PASS: used pure-power top chart "
            f"{args.chart_1000_pure_power_top}"
        )
    if args.characteristic == 0:
        print(
            "PASS: saturation by the common determinant constant is the "
            "unit ideal over Q"
        )
        print(
            f"PASS: chart {args.chart} has no unrestricted degree-"
            f"{args.potential_degree} single-shear solution"
        )
    else:
        print(
            "SCREEN: the saturated slice ideal is unit in characteristic "
            f"{args.characteristic}"
        )
        print("SCOPE: a prime-field screen is not a characteristic-zero proof")


if __name__ == "__main__":
    main()

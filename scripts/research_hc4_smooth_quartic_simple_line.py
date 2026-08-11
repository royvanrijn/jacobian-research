#!/usr/bin/env python3
"""Research the squarefree-line simple-residual-line incidence in HC4NHM14.

This is an exploratory exact-arithmetic driver.  It constructs the full
Hessian-curl, residual-line divisibility, and scalar compatibility equations
for the boundary direction ``(x^2,y^2,0)`` and sends them to Singular.

The generic mode works over the rational function field in the residual-line
and surviving boundary parameters.  It therefore describes only a dense open
parameter stratum; a proof using it must separately record every denominator
appearing in a membership certificate.
"""

from __future__ import annotations

import argparse
import subprocess

import sympy as sp


def exact_quotient(numerator: sp.Expr, divisor: sp.Expr) -> sp.Expr:
    quotient = sp.cancel(numerator / divisor)
    assert sp.denom(quotient) == 1
    return sp.expand(quotient)


def coefficients(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> list[sp.Expr]:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    # Coefficient factorization is unnecessary for ideal generation and can
    # dominate runtime on the parameter-rich scalar-compatibility equations.
    return [value for _, value in polynomial.terms() if value != 0]


def singular(expression: sp.Expr) -> str:
    return sp.sstr(sp.expand(expression)).replace("**", "^")


def build_equations(
    case_name: str, all_b_unknown: bool
) -> tuple[list[sp.Expr], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    x, y, z = sp.symbols("x y z")
    tau, sigma = sp.symbols("tau sigma")
    p, q, r = sp.symbols("p q r")
    bs = sp.symbols("b0:18")
    u, v, w = sp.symbols("u v w")

    linear_entries = [
        bs[3 * index] * x + bs[3 * index + 1] * y + bs[3 * index + 2] * z
        for index in range(6)
    ]
    b00, b01, b02, b11, b12, b22 = linear_entries
    matrix_b = sp.Matrix(
        [[b00, b01, b02], [b01, b11, b12], [b02, b12, b22]]
    )
    if case_name == "squarefree-conic":
        cubic = (
            (x**3 + y**3) / 3
            + z * x * y
            + z**2 * (u * x + v * y)
            + w * z**3
        )
        matrix_0 = sp.Matrix(
            [
                [p * y**2, -q * x * y, y * (q * y - p * x)],
                [-q * x * y, r * x**2, x * (q * x - r * y)],
                [
                    y * (q * y - p * x),
                    x * (q * x - r * y),
                    p * x**2 - 2 * q * x * y + r * y**2,
                ],
            ]
        )
        lam = p * r - q**2
    elif case_name == "squarefree-line":
        cubic = (x**3 + y**3) / 3 + z**2 * (u * x + v * y) + w * z**3
        matrix_0 = sp.Matrix(
            [
                [0, 0, -y**2],
                [0, 0, x**2],
                [-y**2, x**2, p * x**2 + q * x * y + r * y**2],
            ]
        )
        lam = -1
    elif case_name == "double-conic":
        cubic = x**2 * y + z * y**2 + z**2 * (u * x + v * y) + w * z**3
        matrix_0 = sp.Matrix(
            [
                [
                    p * x**2 / 4 + q * x * y / 2 + r * y**2 / 4,
                    -y * (p * x + q * y) / 2,
                    -x * (q * x + r * y) / 2,
                ],
                [-y * (p * x + q * y) / 2, p * y**2, q * x * y],
                [-x * (q * x + r * y) / 2, q * x * y, r * x**2],
            ]
        )
        lam = (p * r - q**2) / 4
    elif case_name == "triple-line":
        cubic = x**3 / 3 + z * y**2 + z**2 * (u * x + v * y) + w * z**3
        matrix_0 = sp.Matrix(
            [
                [0, -y**2, 0],
                [-y**2, p * x**2 + q * x * y + r * y**2, x**2],
                [0, x**2, 0],
            ]
        )
        lam = -1
    else:
        raise ValueError(case_name)
    matrix_a = matrix_0 + z * matrix_b
    direction = sp.Matrix([sp.diff(cubic, variable) for variable in (x, y, z)])

    matrix_c = (matrix_a.adjugate() - lam * direction * direction.T).applyfunc(
        lambda entry: exact_quotient(entry, z)
    )
    vector_e = (matrix_a * direction).applyfunc(lambda entry: exact_quotient(entry, z))
    determinant_remainder = exact_quotient(matrix_a.det(), z)

    derivative = {variable: matrix_c.diff(variable) for variable in (x, y, z)}
    curls: list[sp.Expr] = []
    variables = (x, y, z)
    for row in range(3):
        for first in range(3):
            for second in range(first + 1, 3):
                curls.append(
                    derivative[variables[second]][row, first]
                    - derivative[variables[first]][row, second]
                )

    equations: list[sp.Expr] = []
    for curl in curls:
        equations.extend(coefficients(curl, variables))

    residual_line = y + tau * x + sigma * z
    line_remainder = determinant_remainder.subs(y, -tau * x - sigma * z)
    equations.extend(coefficients(line_remainder, (x, z)))

    scalar_compatibility = sp.expand(
        residual_line * (determinant_remainder - lam * vector_e.dot(direction))
        - z * determinant_remainder
    )
    equations.extend(coefficients(scalar_compatibility, variables))

    # Preserve first occurrence while removing exact duplicates.
    unique: list[sp.Expr] = []
    seen: set[str] = set()
    for equation in equations:
        normalized = sp.expand(equation)
        key = sp.srepr(normalized)
        if normalized != 0 and key not in seen:
            seen.add(key)
            unique.append(normalized)

    if all_b_unknown:
        unknowns = (*bs, u, v, w)
        parameters = (tau, sigma, p, q, r)
    else:
        unknowns = (*bs[:15], u, v, w)
        parameters = (tau, sigma, p, q, r, *bs[15:])
    return unique, unknowns, parameters


def singular_program(
    equations: list[sp.Expr],
    unknowns: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    algorithm: str,
    lift: bool,
    triangular: bool,
    reduced_linear_only: bool,
    known_uv_pivots: bool,
) -> str:
    field = ",".join(map(str, parameters))
    variables = ",".join(map(str, unknowns))
    generators = ",\n".join(singular(equation) for equation in equations)
    lines = [
        "option(redSB);",
        f"ring rr=(0,{field}),({variables}),dp;",
        f"ideal I={generators};",
        "timer=1;",
    ]
    if triangular:
        linear_generators = ",\n".join(
            singular(equation)
            for equation in equations
            if unknown_degree(equation, unknowns) <= 1
        )
        reduced_ideal = "jet(J,1)" if reduced_linear_only else "J"
        lines.extend(
            [
                f"ideal L={linear_generators};",
                "ideal GL=std(L);",
                "ideal J=reduce(I,GL);",
                f"ideal JR={reduced_ideal};",
                f"ideal H={algorithm}(JR);",
            ]
        )
        if known_uv_pivots:
            # In the squarefree-line chart H contains u and v.  Reducing GL
            # by these pivots before recombination avoids a costly but
            # mathematically redundant cross-reduction.
            lines.extend(
                [
                    "ideal K=subst(GL,u,0);",
                    "K=subst(K,v,0);",
                    "ideal G=std(K+H);",
                ]
            )
        else:
            lines.append("ideal G=std(GL+H);")
        lines.extend(
            [
                'print("linear_basis_size="+string(size(GL)));',
                'print("reduced_generator_count="+string(size(J)));',
                'print("active_reduced_generator_count="+string(size(JR)));',
                'print("nonlinear_basis_size="+string(size(H)));',
            ]
        )
    elif lift:
        lines.append(f'ideal G=liftstd(I,T,"{algorithm}");')
    else:
        lines.append(f"ideal G={algorithm}(I);")
    lines.extend(
        [
            "int elapsed=timer;",
            'print("RESULT_BEGIN");',
            f'print("equation_count={len(equations)}");',
            'print("basis_size="+string(size(G)));',
            'print("dimension="+string(dim(G)));',
            'print("elapsed_ticks="+string(elapsed));',
            'print("BASIS_BEGIN");',
            "print(G);",
            'print("BASIS_END");',
        ]
    )
    if lift:
        lines.extend(['print("LIFT_BEGIN");', "print(T);", 'print("LIFT_END");'])
    lines.extend(['print("RESULT_END");', "quit;"])
    return "\n".join(lines) + "\n"


def unknown_degree(
    equation: sp.Expr, unknowns: tuple[sp.Symbol, ...]
) -> int:
    return sp.Poly(equation, *unknowns).total_degree()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=("squarefree-conic", "squarefree-line", "double-conic", "triple-line"),
        default="squarefree-line",
    )
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="slimgb")
    parser.add_argument("--lift", action="store_true")
    parser.add_argument("--linear-only", action="store_true")
    parser.add_argument("--triangular", action="store_true")
    parser.add_argument("--print-stages", action="store_true")
    parser.add_argument("--reduced-linear-only", action="store_true")
    parser.add_argument("--all-b-unknown", action="store_true")
    parser.add_argument("--print-program", action="store_true")
    args = parser.parse_args()

    equations, unknowns, parameters = build_equations(args.case, args.all_b_unknown)
    degree_counts: dict[int, int] = {}
    for equation in equations:
        degree = unknown_degree(equation, unknowns)
        degree_counts[degree] = degree_counts.get(degree, 0) + 1
    if args.linear_only:
        equations = [
            equation for equation in equations if unknown_degree(equation, unknowns) <= 1
        ]
    program = singular_program(
        equations,
        unknowns,
        parameters,
        args.algorithm,
        args.lift,
        args.triangular,
        args.reduced_linear_only,
        args.case == "squarefree-line" and not args.reduced_linear_only,
    )
    if args.triangular and args.print_stages:
        program = program.replace(
            'print("nonlinear_basis_size="+string(size(H)));',
            'print("nonlinear_basis_size="+string(size(H)));\n'
            'print("LINEAR_BASIS_BEGIN"); print(GL); print("LINEAR_BASIS_END");\n'
            'print("NONLINEAR_BASIS_BEGIN"); print(H); print("NONLINEAR_BASIS_END");',
        )
    print(f"constructed {len(equations)} distinct nonzero coefficient equations")
    print(f"case: {args.case}")
    print(f"unknown-degree counts: {degree_counts}")
    print(f"unknowns: {','.join(map(str, unknowns))}")
    print(f"parameters: {','.join(map(str, parameters))}")
    if args.print_program:
        print(program)
        return

    result = subprocess.run(
        ["Singular", "--no-tty", "--quiet"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()

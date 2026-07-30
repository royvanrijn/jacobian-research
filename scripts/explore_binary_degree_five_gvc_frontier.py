#!/usr/bin/env python3
"""Explore the first full jet equations on the binary degree-five GVC frontier.

This is a research script, not a theorem checker.  It constructs the complete
operator jet that can occur in the first two pure equations for each of the
eight cubic-leading normal forms from Proposition 3.17 of
``extended-geometry/SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md``.  The first pure
equation is solved linearly, and the degree layers of the second pure equation
are serialized for subsequent triangular or finite-field analysis.

Only operator orders through ``10-r`` can occur in the second moment when the
lowest order is ``r`` and ``deg(P)=5``.  Higher jets enter later moments and
must be controlled by a separate face or recurrence argument before any
all-order conclusion is made.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "binary_degree_five_gvc_second_moments.json"
)

x, y = sp.symbols("x y")


@dataclass(frozen=True)
class NormalForm:
    name: str
    order: int
    operator: dict[tuple[int, int], sp.Expr]
    polynomial: sp.Expr


def normal_forms() -> tuple[NormalForm, ...]:
    cubic_triple = {(3, 0): sp.Integer(1)}
    cubic_double = {(2, 1): sp.Integer(1)}
    cubic_squarefree = {
        (2, 1): sp.Integer(1),
        (1, 2): sp.Integer(1),
    }
    return (
        NormalForm("r3_triple_x2y3", 3, cubic_triple, x**2 * y**3),
        NormalForm(
            "r3_triple_x_xminusy_y3",
            3,
            cubic_triple,
            x * (x - y) * y**3,
        ),
        NormalForm("r3_triple_xy4", 3, cubic_triple, x * y**4),
        NormalForm("r3_triple_y5", 3, cubic_triple, y**5),
        NormalForm("r3_double_x5", 3, cubic_double, x**5),
        NormalForm("r3_double_xy4", 3, cubic_double, x * y**4),
        NormalForm("r3_double_y5", 3, cubic_double, y**5),
        NormalForm("r3_squarefree_x5", 3, cubic_squarefree, x**5),
    )


def apply_operator(
    polynomial: sp.Expr,
    coefficients: dict[tuple[int, int], sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * sp.diff(polynomial, x, x_order, y, y_order)
            for (x_order, y_order), coefficient in coefficients.items()
        )
    )


def polynomial_coefficients(max_degree: int) -> dict[tuple[int, int], sp.Symbol]:
    return {
        (x_degree, total_degree - x_degree): sp.Symbol(
            f"p{x_degree}{total_degree - x_degree}"
        )
        for total_degree in range(max_degree + 1)
        for x_degree in range(total_degree + 1)
    }


def operator_coefficients(
    leading: dict[tuple[int, int], sp.Expr],
    lowest_order: int,
    max_order: int,
) -> dict[tuple[int, int], sp.Expr]:
    result = dict(leading)
    for total_order in range(lowest_order + 1, max_order + 1):
        for x_order in range(total_order + 1):
            y_order = total_order - x_order
            result[(x_order, y_order)] = sp.Symbol(
                f"l{x_order}{y_order}"
            )
    return result


def solve_first_equation(
    normal_form: NormalForm,
    operator: dict[tuple[int, int], sp.Expr],
    lower_polynomial_coefficients: dict[tuple[int, int], sp.Symbol],
) -> tuple[sp.Expr, dict[sp.Symbol, sp.Expr], list[sp.Expr]]:
    polynomial = sp.expand(
        normal_form.polynomial
        + sum(
            coefficient * x**x_degree * y**y_degree
            for (x_degree, y_degree), coefficient
            in lower_polynomial_coefficients.items()
        )
    )
    first = sp.Poly(apply_operator(polynomial, operator), x, y)
    equations = first.coeffs()
    # Only P_4 and P_3 can occur: the leading operator has order three and
    # every higher jet acting on lower terms is killed by total degree.
    candidates = [
        coefficient
        for (degree, coefficient) in sorted(
            (
                (sum(exponent), coefficient)
                for exponent, coefficient
                in lower_polynomial_coefficients.items()
            ),
            key=lambda item: (item[0], str(item[1])),
            reverse=True,
        )
        if degree >= normal_form.order
    ]
    solutions = sp.solve(
        equations,
        candidates,
        dict=True,
        simplify=False,
        manual=True,
    )
    if not solutions:
        raise RuntimeError(f"{normal_form.name}: first equation has no solution")
    solution = solutions[0]
    return sp.expand(polynomial.subs(solution)), solution, equations


def expression_record(expression: sp.Expr) -> dict[str, object]:
    expanded = sp.expand(expression)
    return {
        "expression": str(expanded),
        "factorization": str(sp.factor(expanded)),
        "symbols": sorted(str(symbol) for symbol in expanded.free_symbols),
        "total_degree_in_parameters": sp.Poly(
            expanded, *sorted(expanded.free_symbols, key=str)
        ).total_degree()
        if expanded.free_symbols
        else 0,
    }


def analyze(normal_form: NormalForm) -> dict[str, object]:
    max_operator_order = 10 - normal_form.order
    lower_coefficients = polynomial_coefficients(4)
    operator = operator_coefficients(
        normal_form.operator,
        normal_form.order,
        max_operator_order,
    )
    normalized_polynomial, first_solution, first_equations = (
        solve_first_equation(normal_form, operator, lower_coefficients)
    )
    second = normalized_polynomial**2
    for _ in range(2):
        second = apply_operator(second, operator)
    second_poly = sp.Poly(second, x, y)
    layers: dict[str, list[dict[str, object]]] = {}
    for monomial, coefficient in second_poly.terms():
        total_degree = sum(monomial)
        layers.setdefault(str(total_degree), []).append(
            {
                "monomial": [int(monomial[0]), int(monomial[1])],
                **expression_record(coefficient),
            }
        )
    return {
        "normal_form": normal_form.name,
        "lowest_operator_order": normal_form.order,
        "leading_operator": {
            f"{a},{b}": str(coefficient)
            for (a, b), coefficient in normal_form.operator.items()
        },
        "leading_polynomial": str(normal_form.polynomial),
        "complete_second_moment_operator_order": max_operator_order,
        "first_equation_count": len(first_equations),
        "first_equation_solution": {
            str(symbol): str(value)
            for symbol, value in sorted(
                first_solution.items(), key=lambda item: str(item[0])
            )
        },
        "normalized_polynomial": str(normalized_polynomial),
        "free_parameter_count_after_first_equation": len(
            (
                set(operator.values())
                | set(lower_coefficients.values())
            )
            - set(first_solution)
            - {sp.Integer(1)}
        ),
        "second_moment_equation_count": len(second_poly.terms()),
        "second_moment_layers": layers,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--branch",
        action="append",
        help="normal-form name; repeat to select several (default: all)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help="JSON output path",
    )
    args = parser.parse_args()

    forms = normal_forms()
    if args.branch:
        selected = tuple(form for form in forms if form.name in args.branch)
        missing = sorted(set(args.branch) - {form.name for form in selected})
        if missing:
            raise SystemExit(f"unknown branches: {', '.join(missing)}")
    else:
        selected = forms

    records = []
    for form in selected:
        print(f"analyzing {form.name}", flush=True)
        records.append(analyze(form))
    artifact = {
        "format": "binary-degree-five-gvc-second-moments-v1",
        "status": "experiment",
        "scope": (
            "complete first and second pure equations only; no all-order "
            "GVC conclusion"
        ),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")
    try:
        display_path = args.output.relative_to(ROOT)
    except ValueError:
        display_path = args.output
    print(f"wrote {display_path}")


if __name__ == "__main__":
    main()

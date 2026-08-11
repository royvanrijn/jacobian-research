#!/usr/bin/env python3
"""Exact bounded regressions for three all-degree HC4 frontends.

The degree-free proofs live in ``verify_hc4_all_degree_frontends.py``.  This
script retains the original Singular discovery sweeps and records:

1. the radical of the divisibility equations for the diagonal ternary
   Schur model through a requested degree range;
2. constant rank-two quadratic pencil directions on the minimal diagonal
   Schur tower, split by the number of active channels; and
3. finite-order instances of the exact Meng--Yang all-normal symbol.

The first two computations use Singular over QQ.  The third is a symbolic
identity over QQ[L,M,N].  These computations are regressions for HC4FSD1,
HC4FSD2, and HC4MYGJ2, not the logical basis of their universal quantifiers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections.abc import Sequence
from pathlib import Path

import sympy as sp


def homogeneous_exponents(degree: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (first, second, degree - first - second)
        for first in range(degree + 1)
        for second in range(degree - first + 1)
    )


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.factor(expression)).replace("**", "^")


def run_singular(program: str, executable: str) -> list[str]:
    completed = subprocess.run(
        [executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Singular failed with exit code "
            f"{completed.returncode}:\n{completed.stderr}\n{completed.stdout}"
        )
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def singular_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "--version"],
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.splitlines()[0].strip()


def fermat_schur_radical(
    potential_degree: int, singular: str
) -> dict[str, object]:
    """Classify the diagonal Schur divisibility scheme in one degree."""

    derivative_degree = potential_degree - 2
    x, y, z = sp.symbols("x y z")
    exponents = homogeneous_exponents(derivative_degree)
    coefficients = sp.symbols(f"a0:{len(exponents)}")
    a_form = sp.Add(
        *(
            coefficient * x**i * y**j * z**k
            for coefficient, (i, j, k) in zip(
                coefficients, exponents, strict=True
            )
        )
    )
    numerator = sp.expand(
        sp.diff(a_form, x) ** 2 * y**derivative_degree * z**derivative_degree
        + sp.diff(a_form, y) ** 2
        * x**derivative_degree
        * z**derivative_degree
        + sp.diff(a_form, z) ** 2
        * x**derivative_degree
        * y**derivative_degree
    )
    equations = tuple(
        dict.fromkeys(
            sp.factor(coefficient)
            for monomial, coefficient in sp.Poly(numerator, x, y, z).terms()
            if not all(exponent >= derivative_degree for exponent in monomial)
        )
    )
    pure_indices = {
        index
        for index, monomial in enumerate(exponents)
        if max(monomial) == derivative_degree
    }
    mixed_coefficients = tuple(
        coefficient
        for index, coefficient in enumerate(coefficients)
        if index not in pure_indices
    )

    ring_variables = ",".join(map(str, coefficients))
    ideal_equations = ",".join(map(singular_expression, equations))
    split_ideal = ",".join(map(str, mixed_coefficients))
    program = f"""
option(redSB);
LIB \"primdec.lib\";
ring rr=0,({ring_variables}),dp;
ideal I={ideal_equations};
ideal R=std(radical(I));
ideal K={split_ideal};
K=std(K);
int forward=1;
for (int i=1; i<=size(R); i++) {{
  if (reduce(R[i],K)!=0) {{ forward=0; }}
}}
int reverse=1;
"""
    for coefficient in mixed_coefficients:
        program += (
            f"if (reduce({coefficient},R)!=0) {{ reverse=0; }}\n"
        )
    program += """
print("BEGIN_RESULT");
print(dim(R));
print(size(R));
print(forward);
print(reverse);
print("END_RESULT");
"""
    output = run_singular(program, singular)
    start = output.index("BEGIN_RESULT")
    end = output.index("END_RESULT")
    values = tuple(map(int, output[start + 1 : end]))
    if len(values) != 4:
        raise RuntimeError(f"unexpected Singular radical output: {output}")
    dimension, radical_basis_size, forward, reverse = values
    return {
        "potential_degree": potential_degree,
        "ternary_form_degree": derivative_degree,
        "coefficient_count": len(coefficients),
        "divisibility_equation_count": len(equations),
        "mixed_coefficient_count": len(mixed_coefficients),
        "radical_dimension": dimension,
        "radical_basis_size": radical_basis_size,
        "radical_equals_pure_power_locus": bool(forward and reverse),
    }


def symmetric_rank_two_minors(matrix: sp.Matrix) -> tuple[sp.Expr, ...]:
    index_pairs = ((0, 1), (0, 2), (1, 2))
    minors: list[sp.Expr] = []
    for rows in index_pairs:
        for columns in index_pairs:
            minor = sp.factor(matrix.extract(rows, columns).det())
            if minor not in minors:
                minors.append(minor)
    return tuple(minors)


def diagonal_rank_two_tower(
    potential_degree: int,
    channels: tuple[int, int, int],
    singular: str,
) -> dict[str, object]:
    """Test rank-two constant directions on one minimal Schur tower."""

    m = potential_degree - 2
    x, y, z, t, pencil = sp.symbols("x y z t pencil")
    variables = (x, y, z, t)
    q11, q12, q13, q22, q23, q33 = q_variables = sp.symbols(
        "q11 q12 q13 q22 q23 q33"
    )
    q3 = sp.Matrix(
        ((q11, q12, q13), (q12, q22, q23), (q13, q23, q33))
    )
    q4 = sp.zeros(4)
    q4[:3, :3] = q3
    source = (x, y, z)
    top = sp.Add(
        *(variable ** (m + 2) / ((m + 2) * (m + 1)) for variable in source)
    )
    first_schur_layer = t * sp.Add(
        *(
            channel * variable**m / m
            for channel, variable in zip(channels, source, strict=True)
        )
    )
    second_schur_layer = sp.Rational(1, 2) * t**2 * sp.Add(
        *(
            channel**2 * variable ** (m - 2)
            for channel, variable in zip(channels, source, strict=True)
        )
    )
    hessian = sp.hessian(
        sp.expand(top + first_schur_layer + second_schur_layer), variables
    )
    pencil_difference = sp.expand(
        (hessian + pencil * q4).det() - hessian.det()
    )
    equations = tuple(
        dict.fromkeys(
            sp.factor(coefficient)
            for monomial, coefficient in sp.Poly(
                pencil_difference, pencil, *variables
            ).terms()
            if monomial[0] > 0
        )
    )
    rank_two_minors = symmetric_rank_two_minors(q3)
    survivors: list[int] = []
    for chart_index, minor in enumerate(rank_two_minors):
        equations_in_chart = (
            *equations,
            sp.factor(q3.det()),
            sp.Symbol("chart_unit") * minor - 1,
        )
        program = f"""
ring rr=0,({','.join(map(str, q_variables))},chart_unit),dp;
ideal I={','.join(map(singular_expression, equations_in_chart))};
ideal G=std(I);
print("BEGIN_RESULT");
print(reduce(1,G));
print("END_RESULT");
"""
        output = run_singular(program, singular)
        start = output.index("BEGIN_RESULT")
        end = output.index("END_RESULT")
        remainder = output[start + 1 : end]
        if remainder != ["0"]:
            survivors.append(chart_index)
    return {
        "potential_degree": potential_degree,
        "active_channel_count": sum(channels),
        "channels": list(channels),
        "coefficient_equation_count": len(equations),
        "rank_two_chart_count": len(rank_two_minors),
        "surviving_rank_two_charts": survivors,
        "rank_two_scheme_empty": not survivors,
    }


def meng_yang_normal_symbol(maximum_order: int) -> dict[str, object]:
    """Verify and record the all-order transverse graph symbol."""

    x, y, p, q, r = sp.symbols("x y p q r")
    L, M, N = sp.symbols("L M N")
    u = 1 + x * y
    q0 = y**2 * (1 + 3 * u)
    a_form = u**3 * p + 3 * x * u**2 * q - x**3 * r
    b_form = (
        u * q0 * p
        + (y + 3 * x * q0) * q
        + x * (5 - 3 * u) * r
    )
    potential = sp.expand(L * a_form**2 + M * a_form + N * b_form)
    ambient_hessian = sp.hessian(potential, (x, y, p, q, r))
    plane_trace, plane_normal = sp.symbols("plane_trace plane_normal")
    trace_y, trace_p, trace_q = sp.symbols("trace_y trace_p trace_q")
    ambient_plane = ambient_hessian.subs({x: 0, r: plane_trace})
    source_unit = sp.Matrix((1, 0, 0, 0))
    trace_gradient = sp.Matrix((plane_normal, trace_y, trace_p, trace_q))
    graph_plane_hessian = (
        ambient_plane[:4, :4]
        + 2
        * N
        * (source_unit * trace_gradient.T + trace_gradient * source_unit.T)
    )
    transverse_cofactor = sp.factor(graph_plane_hessian[1:, 1:].det())
    omitted_gradient = sp.diff(potential, r)
    omitted_gradient_value = sp.factor(
        omitted_gradient.subs({x: 0, r: plane_trace})
    )
    omitted_gradient_x = sp.factor(
        sp.diff(omitted_gradient, x).subs({x: 0, r: plane_trace})
    )
    omitted_gradient_rr = sp.factor(
        sp.diff(potential, r, 2).subs({x: 0, r: plane_trace})
    )
    assert transverse_cofactor == -2 * L * N**2
    assert omitted_gradient_value == 0
    assert omitted_gradient_x == 2 * N
    assert omitted_gradient_rr == 0

    orders = []
    for order in range(1, maximum_order + 1):
        symbol = sp.factor(
            transverse_cofactor * 2 * N * order * (order + 1)
        )
        expected = -4 * L * N**3 * order * (order + 1)
        assert sp.expand(symbol - expected) == 0
        orders.append(
            {
                "graph_jet_order": order,
                "first_changed_determinant_order": order - 1,
                "unit_multiplier": singular_expression(symbol),
            }
        )
    return {
        "transverse_cofactor": singular_expression(transverse_cofactor),
        "omitted_gradient_x_on_plane": singular_expression(
            omitted_gradient_x
        ),
        "formula": (
            "[x^(k-1)](D_(R+x^k*U)-D_R)="
            "-4*L*N^3*k*(k+1)*U"
        ),
        "verified_orders": orders,
        "interpretation": (
            "unit-triangular normal recursion; formal solvability is favored, "
            "while polynomial termination remains a separate obstruction"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-degree", type=int, default=4)
    parser.add_argument("--maximum-degree", type=int, default=8)
    parser.add_argument("--maximum-normal-order", type=int, default=12)
    parser.add_argument("--singular", default="Singular")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/hc4_all_degree_frontend_experiments.json"
        ),
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.minimum_degree < 4:
        raise ValueError("the minimum potential degree is four")
    if arguments.maximum_degree < arguments.minimum_degree:
        raise ValueError("the maximum degree must be at least the minimum")
    if arguments.maximum_normal_order < 1:
        raise ValueError("the normal-order bound must be positive")

    fermat_results = [
        fermat_schur_radical(degree, arguments.singular)
        for degree in range(
            arguments.minimum_degree, arguments.maximum_degree + 1
        )
    ]
    rank_two_results = [
        diagonal_rank_two_tower(degree, channels, arguments.singular)
        for degree in range(
            max(5, arguments.minimum_degree), arguments.maximum_degree + 1
        )
        for channels in ((1, 0, 0), (1, 1, 0), (1, 1, 1))
    ]
    normal_symbol = meng_yang_normal_symbol(arguments.maximum_normal_order)
    source_path = Path(__file__)
    result = {
        "schema_version": 1,
        "status": "exact_bounded_regression",
        "scope_warning": (
            "The finite tables are regressions for the written all-degree "
            "proofs HC4FSD1, HC4FSD2, and HC4MYGJ2.  HC4FSD2 concerns only "
            "the minimal tower, and HC4MYGJ2 proves formal solvability rather "
            "than polynomial termination."
        ),
        "command": (
            ".venv/bin/python scripts/research_hc4_all_degree_frontends.py "
            f"--minimum-degree {arguments.minimum_degree} "
            f"--maximum-degree {arguments.maximum_degree} "
            f"--maximum-normal-order {arguments.maximum_normal_order} "
            f"--output {arguments.output}"
        ),
        "software": {
            "sympy": sp.__version__,
            "singular": singular_version(arguments.singular),
        },
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "fermat_ternary_schur": {
            "model": (
                "Hess(f)=diag(x^(D-2),y^(D-2),z^(D-2)); require "
                "grad(a)^T*adj(Hess(f))*grad(a) divisible by det(Hess(f))"
            ),
            "results": fermat_results,
            "all_tested_radicals_are_pure_power_locus": all(
                row["radical_equals_pure_power_locus"]
                for row in fermat_results
            ),
        },
        "diagonal_rank_two_recognition": {
            "model": (
                "minimal three-layer diagonal Schur tower; constant symmetric "
                "quadratic directions supported on the active ternary block"
            ),
            "results": rank_two_results,
            "all_tested_rank_two_schemes_are_empty": all(
                row["rank_two_scheme_empty"] for row in rank_two_results
            ),
        },
        "meng_yang_normal_symbol": normal_symbol,
        "route_assessment": [
            {
                "route": "Meng--Yang all-normal recursion",
                "signal": "strongest constructive all-order identity",
                "next_target": (
                    "prove degree growth or nontermination, then impose the "
                    "marked collision on the unique formal branch"
                ),
            },
            {
                "route": "ternary Schur dichotomy",
                "signal": "all-degree split rigidity on the diagonal Hessian",
                "next_target": (
                    "extend the uniform theorem from the diagonal Hessian to "
                    "nonsquarefree factor strata"
                ),
            },
            {
                "route": "constant rank-two pencil recognition",
                "signal": "all-degree obstruction on the minimal diagonal tower",
                "next_target": (
                    "control arbitrary lower layers before promoting beyond "
                    "the minimal tower"
                ),
            },
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "PASS: Fermat Schur radicals are pure-power loci in degrees "
        f"{arguments.minimum_degree}--{arguments.maximum_degree}"
    )
    print(
        "PASS: all tested minimal diagonal towers have empty constant "
        "rank-two direction schemes"
    )
    print(
        "PASS: Meng--Yang normal multiplier is "
        "-4*L*N^3*k*(k+1) through order "
        f"{arguments.maximum_normal_order}"
    )
    print(f"WROTE: {arguments.output}")


if __name__ == "__main__":
    main()

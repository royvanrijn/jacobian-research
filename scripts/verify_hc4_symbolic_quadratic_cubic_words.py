#!/usr/bin/env python3
"""Exclude symbolic quadratic--cubic canonical words for HC(4).

For every shared-dual noncommuting support/sign pattern, put

    H1 = a*(q_i + epsilon1*p_j)^2,
    H2 = b*(q_k + epsilon2*p_j)^3

over Q[a,b].  The default order is the exact word T_H2 o T_H1.  The
``cubic-quadratic`` option reverses the two flows.  The signed search
HC4MCP4 found 54 support/sign patterns in each order and 216 signed
coefficient specializations with two coordinate affine pivots.

Parent constant-Hessian preservation is tested first.  Constancy would
force the Hessian determinants at every spatial point to equal the value at
the origin.  This checker takes exact differences at a deterministic list
of small axis points and adjoins

    z_saturation*a*b - 1.

A monomial difference is already a unit after localization at a*b.
Otherwise Singular computes an exact standard basis over Q and checks that
the sampled coefficient ideal is the unit ideal.  Since full polynomial
constancy implies all sampled equalities, a unit sampled ideal excludes
every a*b != 0 exactly; no interpolation claim is used.

The result is coefficient-uniform for this finite mixed-line alphabet.  It
does not classify oblique affine directions, zero coefficients, the reverse
degree order, longer words, or different Hamiltonian supports.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from functools import reduce
from itertools import product
from math import lcm
from pathlib import Path
from typing import Sequence

import sympy as sp

import search_hc4_mixed_canonical_pivots as base
import search_hc4_mixed_quadratic_words as words


a, b, z_saturation = sp.symbols("a b z_saturation")
parameter_variables = (a, b)
zero_point = (0,) * 6


@dataclass(frozen=True)
class Pattern:
    common_dual: int
    h1_source: int
    h2_source: int
    epsilon1: int
    epsilon2: int
    kappa: int
    incidence: str

    @property
    def pattern_id(self) -> str:
        return (
            f"p{self.common_dual}-q{self.h1_source}"
            f"-q{self.h2_source}-e{self.epsilon1:+d}"
            f"-e{self.epsilon2:+d}"
        )


def patterns() -> tuple[Pattern, ...]:
    result: list[Pattern] = []
    for common_dual, h1_source, h2_source, epsilon1, epsilon2 in product(
        range(3), range(3), range(3), (-1, 1), (-1, 1)
    ):
        kappa = (
            epsilon2 * int(h1_source == common_dual)
            - epsilon1 * int(common_dual == h2_source)
        )
        if kappa == 0:
            continue
        if h1_source == common_dual and h2_source == common_dual:
            incidence = "reciprocal"
        elif h1_source == common_dual:
            incidence = "h1_source_hits_h2_dual"
        else:
            assert h2_source == common_dual
            incidence = "h1_dual_hits_h2_source"
        result.append(
            Pattern(
                common_dual=common_dual,
                h1_source=h1_source,
                h2_source=h2_source,
                epsilon1=epsilon1,
                epsilon2=epsilon2,
                kappa=kappa,
                incidence=incidence,
            )
        )
    assert len(result) == 54
    return tuple(result)


def unit_axis(index: int, sign: int = 1) -> tuple[int, ...]:
    point = [0] * 6
    point[index] = sign
    return tuple(point)


def probe_points(pattern: Pattern) -> tuple[tuple[int, ...], ...]:
    """Put structurally relevant dual and source axes first."""

    candidates = [
        unit_axis(3 + pattern.h2_source),
        unit_axis(3 + pattern.h1_source),
        unit_axis(3 + pattern.common_dual),
        unit_axis(pattern.h1_source),
        unit_axis(pattern.h1_source, -1),
        unit_axis(pattern.h2_source),
        unit_axis(pattern.h2_source, -1),
        unit_axis(pattern.common_dual),
        unit_axis(pattern.common_dual, -1),
    ]
    candidates.extend(
        unit_axis(index, sign)
        for index in range(6)
        for sign in (1, -1)
    )
    unique: list[tuple[int, ...]] = []
    for point in candidates:
        if point not in unique:
            unique.append(point)
    return tuple(unique)


def two_axis_probe_points() -> tuple[tuple[int, ...], ...]:
    result: list[tuple[int, ...]] = []
    for left in range(6):
        for right in range(left + 1, 6):
            point = [0] * 6
            point[left] = 1
            point[right] = 1
            result.append(tuple(point))
    return tuple(result)


def transformed_potential(pattern: Pattern, order: str) -> sp.Expr:
    linear1 = (
        base.q[pattern.h1_source]
        + pattern.epsilon1 * base.p[pattern.common_dual]
    )
    linear2 = (
        base.q[pattern.h2_source]
        + pattern.epsilon2 * base.p[pattern.common_dual]
    )
    quadratic_flow = words.flow_of(a * linear1**2)
    cubic_flow = words.flow_of(b * linear2**3)
    if order == "quadratic-cubic":
        inner_flow = quadratic_flow
        outer_flow = cubic_flow
    elif order == "cubic-quadratic":
        inner_flow = cubic_flow
        outer_flow = quadratic_flow
    else:
        raise ValueError(f"unknown composition order: {order}")
    inner_substitution = dict(
        zip(base.variables, inner_flow, strict=True)
    )
    flow = tuple(
        sp.expand(
            expression.subs(inner_substitution, simultaneous=True)
        )
        for expression in outer_flow
    )
    return sp.expand(
        base.base_potential.subs(
            dict(zip(base.variables, flow, strict=True)),
            simultaneous=True,
        )
    )


def polynomial_hessian(
    expression: sp.Expr,
) -> tuple[tuple[sp.Poly, ...], ...]:
    polynomial = sp.Poly(
        expression,
        *(base.variables + parameter_variables),
        domain=sp.QQ,
    )
    return tuple(
        tuple(polynomial.diff(left).diff(right) for right in base.variables)
        for left in base.variables
    )


def determinant_at(
    hessian: Sequence[Sequence[sp.Poly]],
    point: Sequence[int],
) -> sp.Expr:
    matrix: list[list[sp.Expr]] = []
    for row in hessian:
        evaluated_row: list[sp.Expr] = []
        for entry in row:
            evaluated = entry
            for variable, value in zip(
                base.variables, point, strict=True
            ):
                evaluated = evaluated.eval(variable, value)
            evaluated_row.append(evaluated.as_expr())
        matrix.append(evaluated_row)
    return sp.expand(sp.Matrix(matrix).det(method="domain-ge"))


def primitive_integer_polynomial(expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(expression, *parameter_variables, domain=sp.QQ)
    _, cleared = polynomial.clear_denoms(convert=True)
    return sp.expand(cleared.primitive()[1].as_expr())


def is_monomial(expression: sp.Expr) -> bool:
    return len(sp.Poly(expression, *parameter_variables).terms()) == 1


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def singular_unit_ideal(
    equations: Sequence[sp.Expr],
    *,
    timeout: int,
) -> tuple[bool, int, str, tuple[str, ...]]:
    generators = [
        *(singular_expression(equation) for equation in equations),
        "z_saturation*a*b-1",
    ]
    program = "\n".join(
        [
            "ring R=0,(z_saturation,a,b),dp;",
            f"ideal I={','.join(generators)};",
            "ideal G=std(I);",
            'print("BASIS_SIZE="+string(size(G)));',
            'if (reduce(1,G)==0) { print("UNIT=1"); }'
            ' else { print("UNIT=0"); }',
            "ideal E=eliminate(G,z_saturation);",
            'print("ELIMINATION_BEGIN");',
            "E;",
            'print("ELIMINATION_END");',
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Singular failed:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    markers = dict(
        line.split("=", 1)
        for line in completed.stdout.splitlines()
        if "=" in line
    )
    elimination_text = completed.stdout.split(
        "ELIMINATION_BEGIN\n", 1
    )[1].split("\nELIMINATION_END", 1)[0]
    elimination = tuple(
        line.split("=", 1)[1].strip()
        for line in elimination_text.splitlines()
        if "=" in line
    )
    return (
        markers.get("UNIT") == "1",
        int(markers["BASIS_SIZE"]),
        hashlib.sha256(program.encode()).hexdigest(),
        elimination,
    )


def equation_record(point: Sequence[int], equation: sp.Expr) -> dict[str, object]:
    encoded = sp.sstr(equation)
    polynomial = sp.Poly(equation, *parameter_variables)
    return {
        "point": list(point),
        "terms": len(polynomial.terms()),
        "total_degree": polynomial.total_degree(),
        "factorization": sp.sstr(sp.factor(equation)),
        "sha256": hashlib.sha256(encoded.encode()).hexdigest(),
    }


def exceptional_a_value(elimination_basis: Sequence[str]) -> sp.Rational:
    if len(elimination_basis) != 1:
        raise AssertionError(
            f"expected one linear exceptional equation: "
            f"{elimination_basis}"
        )
    equation = sp.sympify(
        elimination_basis[0].replace("^", "**"),
        locals={"a": a, "b": b},
    )
    polynomial = sp.Poly(equation, a, domain=sp.QQ)
    if polynomial.degree() != 1:
        raise AssertionError(
            f"expected a linear exceptional equation: {equation}"
        )
    return sp.Rational(-polynomial.nth(0), polynomial.nth(1))


def singular_full_parent_identity(
    potential: sp.Expr,
    a_value: sp.Rational,
    reference: sp.Expr,
    *,
    timeout: int,
) -> dict[str, object]:
    specialized = sp.expand(potential.subs(a, a_value))
    hessian = sp.hessian(specialized, base.variables)
    denominators = [
        int(
            sp.Poly(
                hessian[row, column],
                *(base.variables + (b,)),
                domain=sp.QQ,
            ).clear_denoms(convert=True)[0]
        )
        for row in range(6)
        for column in range(6)
    ]
    scale = reduce(lcm, denominators, 1)
    entries = [
        singular_expression(sp.expand(scale * hessian[row, column]))
        for row in range(6)
        for column in range(6)
    ]
    target = sp.expand(
        scale**6 * reference.subs(a, a_value)
    )
    program = "\n".join(
        [
            "ring R=0,(x,y,z,u,v,w,b),dp;",
            f"matrix M[6][6]={','.join(entries)};",
            "poly D=det(M);",
            'print("DETERMINANT_TERMS="+string(size(D)));',
            (
                f"if (D-({singular_expression(target)})==0)"
                ' { print("IDENTITY=1"); }'
                ' else { print("IDENTITY=0"); }'
            ),
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Singular full determinant failed:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    markers = dict(
        line.split("=", 1)
        for line in completed.stdout.splitlines()
        if "=" in line
    )
    if markers.get("IDENTITY") != "1":
        raise AssertionError(
            f"exceptional line failed full determinant identity at "
            f"a={a_value}"
        )
    return {
        "a_value": str(a_value),
        "parent_hessian_determinant": sp.sstr(
            sp.expand(reference.subs(a, a_value))
        ),
        "hessian_denominator_scale": scale,
        "scaled_determinant_terms": int(
            markers["DETERMINANT_TERMS"]
        ),
        "singular_input_sha256": hashlib.sha256(
            program.encode()
        ).hexdigest(),
    }


def analyze_pattern(
    pattern: Pattern,
    *,
    order: str,
    singular_timeout: int,
) -> dict[str, object]:
    potential = transformed_potential(pattern, order)
    hessian = polynomial_hessian(potential)
    reference = determinant_at(hessian, zero_point)

    equations: list[sp.Expr] = []
    records: list[dict[str, object]] = []
    certificate_kind: str | None = None
    basis_size: int | None = None
    singular_input_hash: str | None = None
    elimination_basis: tuple[str, ...] = ()
    exact_parent_family: dict[str, object] | None = None

    def append_probe_equations(
        points: Sequence[Sequence[int]],
    ) -> None:
        nonlocal certificate_kind
        for point in points:
            difference = sp.expand(
                determinant_at(hessian, point) - reference
            )
            if difference == 0:
                continue
            equation = primitive_integer_polynomial(difference)
            equations.append(equation)
            records.append(equation_record(point, equation))
            if is_monomial(equation):
                certificate_kind = "localized_monomial"
                return

    append_probe_equations(probe_points(pattern))

    if certificate_kind is None and equations:
        # A partial axis set can have an unnecessarily expensive
        # positive-dimensional standard basis.  Collect every cheap
        # single-axis equation first.
        (
            unit,
            basis_size,
            singular_input_hash,
            elimination_basis,
        ) = singular_unit_ideal(equations, timeout=singular_timeout)
        if unit:
            certificate_kind = "singular_unit_ideal"

    if certificate_kind is None:
        append_probe_equations(two_axis_probe_points())
        if certificate_kind is None and equations:
            (
                unit,
                basis_size,
                singular_input_hash,
                elimination_basis,
            ) = singular_unit_ideal(
                equations,
                timeout=singular_timeout,
            )
            if unit:
                certificate_kind = "singular_unit_ideal"

    if certificate_kind is None and order == "quadratic-cubic":
        raise AssertionError(
            f"symbolic parent locus survived axis probes: "
            f"{pattern.pattern_id}"
        )
    if certificate_kind is None and order == "cubic-quadratic":
        a_value = exceptional_a_value(elimination_basis)
        exact_parent_family = singular_full_parent_identity(
            potential,
            a_value,
            reference,
            timeout=singular_timeout,
        )
        certificate_kind = "exact_parent_family"

    result: dict[str, object] = {
        "pattern_id": pattern.pattern_id,
        "composition_order": order,
        "common_dual": pattern.common_dual,
        "H1_source": pattern.h1_source,
        "H2_source": pattern.h2_source,
        "epsilon1": pattern.epsilon1,
        "epsilon2": pattern.epsilon2,
        "degree_pair": [2, 3],
        "poisson_incidence": pattern.incidence,
        "linear_pairing": pattern.kappa,
        "poisson_bracket_factor": (
            f"{6 * pattern.kappa}*a*b*L1*L2^2"
        ),
        "origin_parent_hessian_determinant": sp.sstr(reference),
        "certificate_kind": certificate_kind,
        "probe_equations": records,
        "parent_constant_locus_with_ab_nonzero": (
            "empty"
            if exact_parent_family is None
            else f"a={exact_parent_family['a_value']}, b!=0"
        ),
    }
    if basis_size is not None:
        result["singular_basis_size"] = basis_size
    if singular_input_hash is not None:
        result["singular_input_sha256"] = singular_input_hash
    if elimination_basis:
        result["sampled_parent_elimination_basis"] = elimination_basis
    if exact_parent_family is not None:
        result["exact_parent_family"] = exact_parent_family
    return result


def run(order: str, singular_timeout: int) -> dict[str, object]:
    if shutil.which("Singular") is None:
        raise RuntimeError("Singular is required for exact unit-ideal checks")

    rows: list[dict[str, object]] = []
    for index, pattern in enumerate(patterns(), start=1):
        print(
            f"progress={index}/54 pattern={pattern.pattern_id}",
            flush=True,
        )
        rows.append(
            analyze_pattern(
                pattern,
                order=order,
                singular_timeout=singular_timeout,
            )
        )

    incidence_census = Counter(
        str(row["poisson_incidence"]) for row in rows
    )
    certificate_census = Counter(
        str(row["certificate_kind"]) for row in rows
    )
    maximum_probes = max(len(row["probe_equations"]) for row in rows)
    parent_probe_survivors = sum(
        row["parent_constant_locus_with_ab_nonzero"]
        != "empty"
        for row in rows
    )
    assert incidence_census == {
        "h1_source_hits_h2_dual": 24,
        "h1_dual_hits_h2_source": 24,
        "reciprocal": 6,
    }
    if order == "quadratic-cubic":
        assert parent_probe_survivors == 0

    return {
        "status": "exact_symbolic_finite_support_theorem",
        "scope": {
            "base": (
                "collision-centred foundational cubic Keller doubling"
            ),
            "composition": (
                "T_H2 o T_H1"
                if order == "quadratic-cubic"
                else "T_H1 o T_H2"
            ),
            "degree_order": (
                [2, 3] if order == "quadratic-cubic" else [3, 2]
            ),
            "H1": "a*(q_i+epsilon1*p_j)^2",
            "H2": "b*(q_k+epsilon2*p_j)^3",
            "coefficient_ring": "Q[a,b]",
            "open_coefficient_locus": "a*b != 0",
            "support": "all 54 noncommuting shared-dual patterns",
            "parent_gate": (
                "exact Hessian-determinant differences at integer "
                "axis probes"
            ),
            "saturation": "z_saturation*a*b-1",
            "limitations": [
                "zero coefficients",
                "reverse degree order",
                "oblique affine directions",
                "longer words",
                "different Hamiltonian supports",
            ],
        },
        "summary": {
            "patterns": len(rows),
            "parent_constant_locus_empty": len(rows)
            - parent_probe_survivors,
            "parent_probe_survivors": parent_probe_survivors,
            "parent_survivors": parent_probe_survivors,
            "maximum_probe_equations": maximum_probes,
            "incidence_census": dict(sorted(incidence_census.items())),
            "certificate_census": dict(
                sorted(certificate_census.items())
            ),
            "reduced_rank_audits_needed_after_parent_gate": (
                parent_probe_survivors
            ),
        },
        "patterns": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--order",
        choices=("quadratic-cubic", "cubic-quadratic"),
        default="quadratic-cubic",
        help="inner-to-outer degree order",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete symbolic census as JSON",
    )
    parser.add_argument(
        "--singular-timeout",
        type=int,
        default=30,
        help="seconds allowed for each incremental exact standard basis",
    )
    args = parser.parse_args()

    result = run(args.order, args.singular_timeout)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(f"artifact={args.output}")
        print(
            "artifact_sha256="
            f"{hashlib.sha256(encoded.encode()).hexdigest()}"
        )
    print(
        "HC4_SYMBOLIC_QUADRATIC_CUBIC_WORD_SUMMARY "
        f"order={args.order}"
    )
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

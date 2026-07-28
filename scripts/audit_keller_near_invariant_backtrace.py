#!/usr/bin/env python3
"""Backtrace the quadratic near-invariant through the frozen BCW circuit.

The quotient coordinates satisfy

    Q = X_18*X_20 - X_6*X_8.

This audit reconstructs the 17-step degree-lowering trace and proves that in
the rank-compressed 24D homogenization,

    Q = c_4*s - v_3*v_5,
    c_4 = -x*(v_3*y + v_5*z).

On the stable source section v_3=-x*z and v_5=-x*y, so Q=x^2*y*z=M.
The stored 21D map satisfies Q(V)-Q=-M*s^2.  Thus Q is a shared-factor gate
residual: it vanishes after the map on the stable source section.  This
identifies its circuit provenance without claiming it is a boundary invariant.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from rank_compressed_bcw_homogenization import (  # noqa: E402
    extract_quadratic_cubic,
    factor_cubic_output,
)
from search_rank_aware_bcw import State, initial_state, support  # noqa: E402
from verify_shared_bcw_33_route import (  # noqa: E402
    apply_shared_step,
    dense_factor,
)


SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "essential_bcw_21_counterexample.json"
)
ROUTE = SCRIPTS / "verify_essential_bcw_21_route.py"


def frozen_plan() -> list[tuple[object, ...]]:
    tree = ast.parse(ROUTE.read_text())
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(getattr(target, "id", None) == "PLAN" for target in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("frozen PLAN not found")


def replay_trace() -> State:
    state = initial_state()
    for component, first_support, second_support in frozen_plan():
        dimension = len(state.variables)
        first = dense_factor(first_support, dimension)
        second = dense_factor(second_support, dimension)
        removed = tuple(a + b for a, b in zip(first, second))
        polynomial = sp.Poly(
            state.expressions[component], *state.variables, domain=sp.QQ
        )
        coefficient = polynomial.coeff_monomial(removed)
        assert coefficient
        selected = (component, removed, coefficient, sum(removed))
        result = apply_shared_step(
            state.expressions,
            state.variables,
            state.registry,
            selected,
            (first, second),
            state.introduced,
        )
        assert result is not None
        state = State(
            result[0],
            result[1],
            result[2],
            result[3],
            state.plan + ((component, support(first), support(second)),),
        )
    assert len(state.variables) == 17
    return state


def dense_row(stored: list[list[object]], columns: int) -> list[sp.Rational]:
    row = [sp.Rational(0)] * columns
    for index, value in stored:
        row[index] = sp.Rational(value)
    return row


def decode_component(
    stored: list[dict[str, object]], variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    answer = sp.Integer(0)
    for term in stored:
        monomial = sp.Integer(1)
        for variable, power in term["monomial"]:
            monomial *= variables[variable] ** power
        answer += sp.Rational(term["coefficient"]) * monomial
    return sp.expand(answer)


def main() -> None:
    artifact = json.loads(SOURCE.read_text())
    state = replay_trace()
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    assert factorization.basis_components == (0, 1, 2, 3, 4, 5)

    x, y, z = state.variables[:3]
    v3 = state.variables[6]
    v5 = state.variables[8]
    c4 = sp.factor(factorization.c[4].as_expr())
    assert c4 == -x * (v3 * y + v5 * z)

    # The quotient matrix identifies q_6,q_8,q_18,q_20 with ambient
    # coordinates X_6,X_8,c_4,s.
    b_rows = artifact["quotient_factorization"]["B_rows"]
    assert dense_row(b_rows[6], 24) == [
        sp.Rational(int(index == 6)) for index in range(24)
    ]
    assert dense_row(b_rows[8], 24) == [
        sp.Rational(int(index == 8)) for index in range(24)
    ]
    assert dense_row(b_rows[18], 24) == [
        sp.Rational(int(index == 21)) for index in range(24)
    ]
    assert dense_row(b_rows[20], 24) == [
        sp.Rational(int(index == 23)) for index in range(24)
    ]

    q_trace = sp.factor(c4 - v3 * v5)
    assert q_trace == -v3 * v5 - v3 * x * y - v5 * x * z
    source_section = {v3: -x * z, v5: -x * y}
    m = x**2 * y * z
    assert sp.expand(q_trace.subs(source_section)) == m

    variables = sp.symbols("X_0:21")
    h = [
        decode_component(component, variables)
        for component in artifact["H"]
    ]
    image = [variable + correction for variable, correction in zip(variables, h)]
    s = variables[20]
    q = variables[18] * s - variables[6] * variables[8]
    q_image = image[18] * image[20] - image[6] * image[8]
    stored_m = variables[0] ** 2 * variables[1] * variables[2]
    assert sp.factor(q_image - q) == -stored_m * s**2

    points = [
        [sp.Rational(value) for value in point]
        for point in artifact["collision_points"]
    ]
    q_values = [
        point[18] * point[20] - point[6] * point[8]
        for point in points
    ]
    m_values = [point[0] ** 2 * point[1] * point[2] for point in points]
    assert q_values == m_values == [
        sp.Rational(0),
        -sp.Rational(39, 16),
        sp.Rational(39, 16),
    ]

    print("PASS near-invariant backtrace: q_18*s-q_6*q_8 lifts to c_4*s-v_3*v_5")
    print("PASS near-invariant backtrace: c_4=-x*(v_3*y+v_5*z)")
    print("PASS near-invariant backtrace: stable source section gives Q=M=x^2*y*z")
    print("PASS near-invariant backtrace: Q(V)-Q=-M*s^2")
    print("PASS near-invariant backtrace: Q is a shared-factor gate residual")


if __name__ == "__main__":
    main()

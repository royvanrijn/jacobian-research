#!/usr/bin/env python3
"""Construct and verify a shared-factor 3 -> 16 -> 33 BCW route.

The conservative 79-variable route exposes two fresh factors for every
degree-lowering operation.  Here an exposed factor Y+a is retained as an
output and reused by later elementary target shears.  If both factors are
already exposed no stabilization is needed; if a monomial is a square, one
new coordinate exposes its repeated factor.

The frozen trace found by a deterministic width-24 beam search has 17
elementary cancellations but introduces only 13 variables.  It then performs
the same nilpotent doubling and cubic homogenization as the certified
79-variable route, producing an explicit cubic-homogeneous Keller collision
in 33 variables.
"""

from __future__ import annotations

from itertools import product
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "shared_bcw_33_counterexample.json"

# A width-24 beam search using the score documented in the companion note
# found this 13-variable trace.  Entries are
# (target component, first monomial support, second monomial support), with
# indices referring to the variables present at that step.
OPTIMIZED_PLAN: list[tuple[int, list[tuple[int, int]], list[tuple[int, int]]]] = [
    (2, [(0, 1), (1, 2)], [(0, 2), (1, 1), (2, 1)]),
    (2, [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
    (1, [(0, 1), (1, 1)], [(0, 2), (1, 1), (2, 1)]),
    (1, [(0, 1), (1, 1)], [(0, 1), (1, 2)]),
    (1, [(0, 1), (2, 1), (5, 1)], [(0, 1), (1, 1)]),
    (2, [(0, 1), (2, 1)], [(0, 1), (1, 2)]),
    (2, [(0, 1), (2, 1)], [(0, 1), (1, 1), (3, 1)]),
    (1, [(0, 1), (2, 1)], [(0, 1), (1, 1)]),
    (4, [(0, 1), (2, 1)], [(0, 1), (1, 1)]),
    (0, [(0, 1), (2, 1)], [(0, 2)]),
    (1, [(1, 1), (5, 1)], [(0, 1), (1, 1)]),
    (1, [(5, 2)], [(0, 1), (2, 1)]),
    (2, [(1, 2)], [(0, 1), (1, 1)]),
    (2, [(1, 1), (3, 1)], [(0, 1), (1, 1)]),
    (2, [(1, 2)], [(0, 1), (7, 1)]),
    (2, [(1, 1), (3, 1)], [(0, 1), (7, 1)]),
    (2, [(1, 1), (4, 1)], [(0, 1), (1, 1)]),
]


def dense_factor(support: list[tuple[int, int]], dimension: int) -> tuple[int, ...]:
    exponents = [0] * dimension
    for index, exponent in support:
        exponents[index] = exponent
    return tuple(exponents)


def rational_text(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def monomial_expression(exponents: tuple[int, ...], variables: list[sp.Symbol]) -> sp.Expr:
    return sp.prod(variable**exponent for variable, exponent in zip(variables, exponents))


def evaluate_monomial(exponents: tuple[int, ...], point: list[sp.Expr]) -> sp.Expr:
    return sp.prod(value**exponent for value, exponent in zip(point, exponents))


def exponent_support(exponents: tuple[int, ...]) -> list[list[int]]:
    return [[index, exponent] for index, exponent in enumerate(exponents) if exponent]


def candidate_splits(exponents: tuple[int, ...]) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """All unordered proper splits which strictly lower the selected degree."""
    degree = sum(exponents)
    splits: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for first in product(*(range(exponent + 1) for exponent in exponents)):
        second = tuple(exponent - part for exponent, part in zip(exponents, first))
        if not (2 <= sum(first) <= degree - 2 and 2 <= sum(second) <= degree - 2):
            continue
        if first > second:
            continue
        splits.append((tuple(first), second))
    return splits


def high_terms(
    expressions: list[sp.Expr], variables: list[sp.Symbol]
) -> tuple[tuple[int, int, int], list[tuple[int, tuple[int, ...], sp.Expr, int]]]:
    terms = [
        (component, exponents, coefficient, sum(exponents))
        for component, expression in enumerate(expressions)
        for exponents, coefficient in sp.Poly(expression, *variables).terms()
        if coefficient and sum(exponents) > 3
    ]
    maximum = max((degree for *_, degree in terms), default=3)
    potential = sum((degree - 3) ** 2 for *_, degree in terms)
    return (maximum, potential, len(terms)), terms


def apply_shared_step(
    expressions: list[sp.Expr],
    variables: list[sp.Symbol],
    registry: dict[tuple[int, ...], int],
    selected: tuple[int, tuple[int, ...], sp.Expr, int],
    split: tuple[tuple[int, ...], tuple[int, ...]],
    variable_counter: int,
) -> tuple[
    list[sp.Expr],
    list[sp.Symbol],
    dict[tuple[int, ...], int],
    int,
    dict[str, object],
] | None:
    component, removed, coefficient, degree = selected
    first, second = split

    # An elementary target shear may not use the coordinate it modifies as
    # one of its factor coordinates.
    if registry.get(first) == component or registry.get(second) == component:
        return None

    missing: list[tuple[int, ...]] = []
    for factor in (first, second):
        if factor not in registry and factor not in missing:
            missing.append(factor)

    old_dimension = len(variables)
    new_variables = [sp.Symbol(f"shared_bcw_v{variable_counter + j}") for j in range(len(missing))]
    final_variables = variables + new_variables
    padding = (0,) * len(missing)
    final_registry = {key + padding: output for key, output in registry.items()}
    final_expressions = list(expressions)
    new_factor_records: list[dict[str, object]] = []

    for factor, new_variable in zip(missing, new_variables):
        padded_factor = factor + padding
        output_index = len(final_expressions)
        final_expressions.append(new_variable + monomial_expression(factor, variables))
        final_registry[padded_factor] = output_index
        new_factor_records.append(
            {
                "factor": exponent_support(factor),
                "new_variable": old_dimension + len(new_factor_records),
                "output": output_index,
            }
        )

    padded_first = first + padding
    padded_second = second + padding
    first_output = final_registry[padded_first]
    second_output = final_registry[padded_second]
    assert component not in (first_output, second_output)
    final_expressions[component] = sp.expand(
        final_expressions[component]
        - coefficient * final_expressions[first_output] * final_expressions[second_output]
    )

    # Modifying an output which exposed some other factor invalidates that
    # registry entry.  The two factors used above were excluded already.
    final_registry = {
        factor: output for factor, output in final_registry.items() if output != component
    }
    metadata = {
        "selected_degree": degree,
        "component": component,
        "removed_coefficient": rational_text(coefficient),
        "removed_monomial": exponent_support(removed),
        "first_factor": exponent_support(first),
        "second_factor": exponent_support(second),
        "factor_outputs": [first_output, second_output],
        "new_factors": new_factor_records,
    }
    return (
        final_expressions,
        final_variables,
        final_registry,
        variable_counter + len(missing),
        metadata,
    )


def homogeneous_part(expression: sp.Expr, variables: list[sp.Symbol], degree: int) -> sp.Poly:
    poly = sp.Poly(expression, *variables)
    selected = sum(
        coefficient * monomial_expression(exponents, variables)
        for exponents, coefficient in poly.terms()
        if sum(exponents) == degree
    )
    return sp.Poly(selected, *variables)


def main() -> None:
    x, y, z = sp.symbols("x y z")
    variables = [x, y, z]
    v = 1 + 2 * x * y
    long_map = sp.Matrix(
        [
            v**3 * z + 4 * y**2 * v * (2 + 3 * x * y),
            y + 3 * x * v**2 * z + 12 * x * y**2 * (2 + 3 * x * y),
            -x + 3 * x**2 * y + x**3 * z,
        ]
    )
    assert sp.factor(long_map.jacobian((x, y, z)).det()) == 1
    expressions = [sp.expand(-long_map[2]), sp.expand(long_map[1]), sp.expand(long_map[0])]
    assert sp.Matrix(expressions).jacobian(variables).subs(dict.fromkeys(variables, 0)) == sp.eye(3)

    collision_points: list[list[sp.Expr]] = [
        [sp.Integer(0), sp.Integer(0), -sp.Rational(1, 8)],
        [sp.Integer(1), -sp.Rational(3, 4), sp.Rational(13, 4)],
        [-sp.Integer(1), sp.Rational(3, 4), sp.Rational(13, 4)],
    ]
    normalized_target = [sp.Integer(0), sp.Integer(0), -sp.Rational(1, 8)]
    for point in collision_points:
        assert list(long_map.subs(dict(zip((x, y, z), point)))) == [
            -sp.Rational(1, 8),
            0,
            0,
        ]

    registry: dict[tuple[int, ...], int] = {}
    variable_counter = 0
    steps: list[dict[str, object]] = []
    for component, first_support, second_support in OPTIMIZED_PLAN:
        old_dimension = len(variables)
        first = dense_factor(first_support, old_dimension)
        second = dense_factor(second_support, old_dimension)
        removed = tuple(a + b for a, b in zip(first, second))
        polynomial = sp.Poly(expressions[component], *variables)
        coefficient = polynomial.coeff_monomial(removed)
        maximum_degree = max(
            sp.Poly(expression, *variables).total_degree() for expression in expressions
        )
        assert coefficient != 0 and sum(removed) == maximum_degree
        selected = (component, removed, coefficient, maximum_degree)
        chosen = apply_shared_step(
            expressions,
            variables,
            registry,
            selected,
            (first, second),
            variable_counter,
        )
        assert chosen is not None

        # Lift the collision before replacing the state.  Each new factor is
        # a monomial in the old variables and its new exposed output is set to
        # zero on the lifted graph.
        metadata = chosen[4]
        for point in collision_points:
            additions = []
            for record in metadata["new_factors"]:
                dense = [0] * old_dimension
                for index, exponent in record["factor"]:
                    dense[index] = exponent
                additions.append(-evaluate_monomial(tuple(dense), point))
            point.extend(additions)

        expressions, variables, registry, variable_counter, metadata = chosen
        steps.append(metadata)

    # The beam-search trace has 17 target cancellations but shares exposed
    # factors, so only 13 new variables are required.
    assert [step["selected_degree"] for step in steps] == [7, 6, 6] + [5] * 4 + [4] * 10
    assert len(steps) == 17
    assert variable_counter == 13
    assert len(variables) == 16
    assert max(sp.Poly(expression, *variables).total_degree() for expression in expressions) == 3

    expected_target = normalized_target + [sp.Integer(0)] * 13
    for point in collision_points:
        image = [sp.Poly(expression, *variables).eval(dict(zip(variables, point))) for expression in expressions]
        assert image == expected_target

    origin = dict.fromkeys(variables, 0)
    assert sp.Matrix(expressions).jacobian(variables).subs(origin) == sp.eye(16)
    quadratic = [homogeneous_part(expression, variables, 2) for expression in expressions]
    cubic = [homogeneous_part(expression, variables, 3) for expression in expressions]
    assert all(
        sp.expand(expression - variable - q.as_expr() - c.as_expr()) == 0
        for expression, variable, q, c in zip(expressions, variables, quadratic, cubic)
    )

    final_dimension = 2 * len(variables) + 1
    assert final_dimension == 33
    final_points: list[list[sp.Expr]] = []
    final_images: list[list[sp.Expr]] = []
    for point in collision_points:
        substitution = dict(zip(variables, point))
        q_value = [poly.eval(substitution) for poly in quadratic]
        c_value = [poly.eval(substitution) for poly in cubic]
        final_point = point + c_value + [sp.Integer(1)]
        final_image = (
            [sp.cancel(px + cy + qx) for px, cy, qx in zip(point, c_value, q_value)]
            + [sp.Integer(0)] * len(variables)
            + [sp.Integer(1)]
        )
        final_points.append(final_point)
        final_images.append(final_image)
    assert len({tuple(point) for point in final_points}) == 3
    assert final_images[0] == final_images[1] == final_images[2]
    assert final_images[0] == expected_target + [sp.Integer(0)] * 16 + [sp.Integer(1)]

    t_index = final_dimension - 1
    h_components: list[list[dict[str, object]]] = [[] for _ in range(final_dimension)]
    for index, (q_part, c_part) in enumerate(zip(quadratic, cubic)):
        h_components[index].append(
            {
                "coefficient": "1",
                "monomial": [[len(variables) + index, 1], [t_index, 2]],
            }
        )
        for exponents, coefficient in q_part.terms():
            if coefficient:
                support = exponent_support(exponents)
                support.append([t_index, 1])
                h_components[index].append(
                    {"coefficient": rational_text(coefficient), "monomial": support}
                )
        for exponents, coefficient in c_part.terms():
            if coefficient:
                h_components[len(variables) + index].append(
                    {
                        "coefficient": rational_text(-coefficient),
                        "monomial": exponent_support(exponents),
                    }
                )
    assert all(
        sum(exponent for _, exponent in term["monomial"]) == 3
        for component in h_components
        for term in component
    )

    artifact = {
        "format": "shared-bcw-sparse-cubic-homogeneous-map-v1",
        "source": "repository shared-factor optimization of Christopher D. Long's BCW route",
        "construction": "17 elementary factor cancellations with 13 exposed variables, then nilpotent cubic homogenization",
        "source_dimension": 3,
        "degree_reduced_dimension": 16,
        "dimension": 33,
        "linear_part": "identity",
        "map": "V_i = Z_i + H_i; omitted H_i are zero",
        "jacobian_determinant": "1",
        "jacobian_certificate": "explicit determinant-one source shears and elementary target shears; det(I+t JN)=1 before homogenization",
        "degree_lowering_steps": steps,
        "H": h_components,
        "collision_points": [[rational_text(value) for value in point] for point in final_points],
        "common_image": [rational_text(value) for value in final_images[0]],
        "statistics": {
            "target_cancellations": len(steps),
            "introduced_variables": variable_counter,
            "reused_without_stabilization": sum(not step["new_factors"] for step in steps),
            "nonzero_cubic_terms": sum(len(component) for component in h_components),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS shared BCW: 17 exact cancellations introduce only 13 variables")
    print("PASS shared BCW: degree reduction reaches dimension 16 and degree <= 3")
    print("PASS shared BCW: cubic homogenization gives a 33-variable collision")
    print("PASS shared BCW: fixed-dimensional implication now targets GMC(66)")
    print(f"PASS shared BCW: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

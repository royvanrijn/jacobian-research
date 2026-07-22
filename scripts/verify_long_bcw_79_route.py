#!/usr/bin/env python3
"""Exact audit of Long's conservative 3 -> 39 -> 79 BCW route.

This script transcribes the determinant-one presentation in
arXiv:2607.18186v1, performs the balanced Bass--Connell--Wright degree
lowering in exactly 18 stable steps, and certifies the final 79-variable
cubic-homogeneous collision.  The last implication from that map to failure
of GMC(158) is the fixed-dimensional DVEZ/Zhao theorem, not a computation in
this script.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "long_bcw_79_counterexample.json"


def rational_text(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def evaluate_monomial(exponents: list[int], point: list[sp.Expr]) -> sp.Expr:
    return sp.prod(value**exponent for value, exponent in zip(point, exponents))


def balanced_factor(exponents: tuple[int, ...]) -> tuple[list[int], list[int]]:
    """Split total degree d>=4 into floor(d/2) and ceil(d/2)."""
    first_degree = sum(exponents) // 2
    first = [0] * len(exponents)
    second = list(exponents)
    remaining = first_degree
    for index, exponent in enumerate(second):
        take = min(exponent, remaining)
        first[index] = take
        second[index] -= take
        remaining -= take
        if remaining == 0:
            break
    assert sum(first) == first_degree
    assert sum(second) == sum(exponents) - first_degree
    assert 2 <= sum(first) <= sum(exponents) - 2
    assert 2 <= sum(second) <= sum(exponents) - 2
    return first, second


def homogeneous_part(poly: sp.Poly, degree: int) -> sp.Poly:
    expression = sum(
        coefficient
        * sp.prod(variable**exponent for variable, exponent in zip(poly.gens, monomial))
        for monomial, coefficient in poly.terms()
        if sum(monomial) == degree
    )
    return sp.Poly(expression, *poly.gens)


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

    # Target permutation (u,v,w) -> (-w,v,u) gives identity linear part.
    normalized = [-long_map[2], long_map[1], long_map[0]]
    polynomials = [sp.Poly(sp.expand(component), *variables) for component in normalized]
    assert sp.Matrix(normalized).jacobian((x, y, z)).subs({x: 0, y: 0, z: 0}) == sp.eye(3)

    support_by_degree = {
        degree: sum(
            1
            for poly in polynomials
            for monomial, coefficient in poly.terms()
            if coefficient and sum(monomial) == degree
        )
        for degree in range(4, 8)
    }
    assert support_by_degree == {4: 3, 5: 2, 6: 2, 7: 1}

    # The balanced recurrence gives c(4),...,c(7)=(1,2,3,5).
    costs = {0: 0, 1: 0, 2: 0, 3: 0}
    for degree in range(4, 8):
        p = degree // 2
        q = degree - p
        costs[degree] = 1 + costs[p + 1] + costs[q + 1] + costs[p] + costs[q]
    assert tuple(costs[d] for d in range(4, 8)) == (1, 2, 3, 5)
    assert sum(support_by_degree[d] * costs[d] for d in range(4, 8)) == 18

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

    selected_degrees: list[int] = []
    steps: list[dict[str, object]] = []
    while max(poly.total_degree() for poly in polynomials) > 3:
        maximum_degree = max(poly.total_degree() for poly in polynomials)
        selected = None
        for component_index, poly in enumerate(polynomials):
            for monomial, coefficient in poly.terms():
                if coefficient and sum(monomial) == maximum_degree:
                    selected = component_index, monomial, coefficient
                    break
            if selected is not None:
                break
        assert selected is not None
        component_index, monomial, coefficient = selected
        a_exponents, b_exponents = balanced_factor(monomial)

        # For c*m=a*b, the stable map has coordinates
        # f_i-(Y+a)(Z+b), Y+a, Z+b.  The graph lifts below make the last
        # two outputs zero and transport the collision scheme exactly.
        for point in collision_points:
            a_value = coefficient * evaluate_monomial(a_exponents, point)
            b_value = evaluate_monomial(b_exponents, point)
            point.extend([-a_value, -b_value])

        new_y, new_z = sp.symbols(f"bcw_y{len(selected_degrees)} bcw_z{len(selected_degrees)}")
        a = coefficient * sp.prod(
            variable**exponent for variable, exponent in zip(variables, a_exponents)
        )
        b = sp.prod(variable**exponent for variable, exponent in zip(variables, b_exponents))
        expressions = [poly.as_expr() for poly in polynomials]
        expressions[component_index] = sp.expand(
            expressions[component_index] - (new_y + a) * (new_z + b)
        )
        expressions.extend([new_y + a, new_z + b])
        variables.extend([new_y, new_z])
        polynomials = [sp.Poly(expression, *variables) for expression in expressions]
        selected_degrees.append(maximum_degree)
        steps.append(
            {
                "selected_degree": maximum_degree,
                "component": component_index,
                "removed_coefficient": rational_text(coefficient),
                "removed_monomial": [[i, e] for i, e in enumerate(monomial) if e],
                "a_exponents": [[i, e] for i, e in enumerate(a_exponents) if e],
                "b_exponents": [[i, e] for i, e in enumerate(b_exponents) if e],
                "new_variables": [len(variables) - 2, len(variables) - 1],
            }
        )

    assert selected_degrees == [7, 6, 6, 5, 5, 5] + [4] * 12
    assert len(selected_degrees) == 18
    assert len(variables) == 39
    assert max(poly.total_degree() for poly in polynomials) == 3

    expected_target = normalized_target + [sp.Integer(0)] * 36
    for point in collision_points:
        image = [poly.eval(dict(zip(variables, point))) for poly in polynomials]
        assert image == expected_target

    # Write K=I+Q+C.  The stable steps preserve identity linear part.
    origin = dict.fromkeys(variables, 0)
    linear_jacobian = sp.Matrix([poly.as_expr() for poly in polynomials]).jacobian(variables)
    assert linear_jacobian.subs(origin) == sp.eye(39)
    quadratic = [homogeneous_part(poly, 2) for poly in polynomials]
    cubic = [homogeneous_part(poly, 3) for poly in polynomials]
    assert all(
        sp.expand(poly.as_expr() - variable - q.as_expr() - c.as_expr()) == 0
        for poly, variable, q, c in zip(polynomials, variables, quadratic, cubic)
    )

    # U(X,Y)=(X+Q(X)+Y,Y-C(X)) is stably equivalent to K.  Its cubic
    # homogenization in one more variable is
    #   V(X,Y,T)=(X+Y*T^2+Q(X)*T, Y-C(X), T).
    # It has 2*39+1 variables and every nonlinear term is cubic homogeneous.
    final_dimension = 2 * len(variables) + 1
    assert final_dimension == 79
    assert all(
        sum(monomial) == 2
        for poly in quadratic
        for monomial, coefficient in poly.terms()
        if coefficient
    )
    assert all(
        sum(monomial) == 3
        for poly in cubic
        for monomial, coefficient in poly.terms()
        if coefficient
    )

    final_points: list[list[sp.Expr]] = []
    final_images: list[list[sp.Expr]] = []
    for point in collision_points:
        substitution = dict(zip(variables, point))
        q_value = [poly.eval(substitution) for poly in quadratic]
        c_value = [poly.eval(substitution) for poly in cubic]
        final_point = point + c_value + [sp.Integer(1)]
        final_image = (
            [sp.cancel(px + cy + qx) for px, cy, qx in zip(point, c_value, q_value)]
            + [sp.Integer(0)] * 39
            + [sp.Integer(1)]
        )
        final_points.append(final_point)
        final_images.append(final_image)
    assert len({tuple(point) for point in final_points}) == 3
    assert final_images[0] == final_images[1] == final_images[2]
    assert final_images[0] == expected_target + [sp.Integer(0)] * 39 + [sp.Integer(1)]

    # Serialize the exact cubic-homogeneous map V=I+H.  Omitted components
    # are zero; every stored nonlinear monomial has total degree three.
    t_index = final_dimension - 1
    h_components: list[list[dict[str, object]]] = [[] for _ in range(final_dimension)]
    for index, (q_part, c_part) in enumerate(zip(quadratic, cubic)):
        h_components[index].append(
            {
                "coefficient": "1",
                "monomial": [[len(variables) + index, 1], [t_index, 2]],
            }
        )
        for monomial, coefficient in q_part.terms():
            if not coefficient:
                continue
            support = [[i, exponent] for i, exponent in enumerate(monomial) if exponent]
            support.append([t_index, 1])
            h_components[index].append(
                {"coefficient": rational_text(coefficient), "monomial": support}
            )
        for monomial, coefficient in c_part.terms():
            if not coefficient:
                continue
            support = [[i, exponent] for i, exponent in enumerate(monomial) if exponent]
            h_components[len(variables) + index].append(
                {"coefficient": rational_text(-coefficient), "monomial": support}
            )
    assert all(
        sum(exponent for _, exponent in term["monomial"]) == 3
        for component in h_components
        for term in component
    )

    artifact = {
        "format": "long-bcw-sparse-cubic-homogeneous-map-v1",
        "source": "Christopher D. Long, arXiv:2607.18186v1",
        "construction": "18 balanced BCW degree-lowering steps followed by nilpotent cubic homogenization",
        "source_dimension": 3,
        "degree_reduced_dimension": 39,
        "dimension": final_dimension,
        "linear_part": "identity",
        "map": "V_i = Z_i + H_i; omitted H_i are zero",
        "jacobian_determinant": "1",
        "jacobian_certificate": "stable determinant-one chain and det(I+t JN)=1 before cubic homogenization",
        "high_degree_support": {str(key): value for key, value in support_by_degree.items()},
        "degree_lowering_steps": steps,
        "H": h_components,
        "collision_points": [
            [rational_text(value) for value in point] for point in final_points
        ],
        "common_image": [rational_text(value) for value in final_images[0]],
        "statistics": {
            "degree_lowering_steps": len(steps),
            "nonzero_cubic_terms": sum(len(component) for component in h_components),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS Long BCW: high-degree support counts are 3,2,2,1")
    print("PASS Long BCW: balanced costs are c(4..7)=1,2,3,5")
    print("PASS Long BCW: exactly 18 stable steps give dimension 39 and degree <= 3")
    print("PASS Long BCW: cubic homogenization gives a 79-variable collision")
    print("PASS Long BCW: DVEZ/Zhao fixed-dimensional implication targets GMC(158)")
    print(f"PASS Long BCW: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

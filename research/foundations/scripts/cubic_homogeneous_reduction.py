#!/usr/bin/env python3
"""Construct an explicit BCW/Yagzhev cubic-homogeneous counterexample.

The construction follows the four-step reduction recorded in Campbell's
implementation of the Bass--Connell--Wright/Yagzhev reduction:

1. Remove every term of degree > 3 by polynomial stable equivalences.
2. Keep the identity linear part (the input map is normalized in advance).
3. Apply the Segre homogenization with a new variable t.
4. Remove the quadratic part with a second stable equivalence.

The generated map has the form I + H on A^95, with H cubic homogeneous.
The output is sparse JSON so that the result is explicit without printing 95
large coordinate polynomials in a manuscript.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "cubic_homogeneous_counterexample.json"


def rational_text(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def eval_monomial(exponents: tuple[int, ...] | list[int], point: list[sp.Expr]) -> sp.Expr:
    result = sp.Integer(1)
    for exponent, coordinate in zip(exponents, point):
        result *= coordinate**exponent
    return sp.cancel(result)


def sparse_terms(poly: sp.Poly) -> list[dict[str, object]]:
    terms = []
    for exponents, coefficient in poly.terms():
        if coefficient == 0:
            continue
        support = [[index, exponent] for index, exponent in enumerate(exponents) if exponent]
        terms.append({"coefficient": rational_text(coefficient), "monomial": support})
    return terms


def choose_degree_two_factor(exponents: tuple[int, ...]) -> tuple[list[int], list[int]]:
    """Split a monomial of degree >= 4 into factors of degrees 2 and d-2."""
    first = [0] * len(exponents)
    second = list(exponents)
    remaining = 2
    for index, exponent in enumerate(second):
        take = min(exponent, remaining)
        first[index] = take
        second[index] -= take
        remaining -= take
        if remaining == 0:
            break
    assert sum(first) == 2 and sum(second) >= 2
    return first, second


def homogeneous_part(poly: sp.Poly, degree: int) -> sp.Poly:
    expression = sum(
        coefficient * sp.prod(variable**exponent for variable, exponent in zip(poly.gens, monomial))
        for monomial, coefficient in poly.terms()
        if sum(monomial) == degree
    )
    return sp.Poly(expression, *poly.gens)


def main() -> None:
    x, y, z = sp.symbols("x y z")
    variables = [x, y, z]
    u = 1 + x * y
    announced = [
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    ]

    # Reorder and scale the target so that the linear part is the identity and
    # the Jacobian determinant is 1.
    normalized = [announced[2] / 2, announced[1], announced[0]]
    polynomials = [sp.Poly(sp.expand(component), *variables) for component in normalized]
    assert sp.Matrix(normalized).jacobian((x, y, z)).subs({x: 0, y: 0, z: 0}) == sp.eye(3)
    assert sp.factor(sp.Matrix(normalized).jacobian((x, y, z)).det()) == 1

    collision_points = [
        [sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)],
        [sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)],
        [-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)],
    ]
    normalized_target = [sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)]
    steps: list[dict[str, object]] = []

    while True:
        selected = None
        maximum_degree = max(poly.total_degree() for poly in polynomials)
        if maximum_degree <= 3:
            break
        for component_index, poly in enumerate(polynomials):
            for monomial, coefficient in poly.terms():
                if sum(monomial) == maximum_degree:
                    selected = (component_index, monomial, coefficient)
                    break
            if selected is not None:
                break
        assert selected is not None
        component_index, monomial, coefficient = selected
        a_exponents, b_exponents = choose_degree_two_factor(monomial)

        # Lift every certified source point to the stable extension.  On this
        # graph slice, the two new output coordinates are zero.
        for point in collision_points:
            a_value = coefficient * eval_monomial(a_exponents, point)
            b_value = eval_monomial(b_exponents, point)
            point.extend([-a_value, -b_value])

        new_y, new_z = sp.symbols(f"v{len(steps)} w{len(steps)}")
        a = coefficient * sp.prod(
            variable**exponent for variable, exponent in zip(variables, a_exponents)
        )
        b = sp.prod(variable**exponent for variable, exponent in zip(variables, b_exponents))
        expressions = [poly.as_expr() for poly in polynomials]
        expressions[component_index] = sp.expand(
            expressions[component_index] - (new_y + a) * (new_z + b)
        )
        expressions.extend([new_y + a, new_z + b])

        steps.append(
            {
                "component": component_index,
                "removed_coefficient": rational_text(coefficient),
                "removed_monomial": [[i, e] for i, e in enumerate(monomial) if e],
                "a_exponents": [[i, e] for i, e in enumerate(a_exponents) if e],
                "b_exponents": [[i, e] for i, e in enumerate(b_exponents) if e],
                "new_variables": [len(variables), len(variables) + 1],
            }
        )
        variables.extend([new_y, new_z])
        polynomials = [sp.Poly(expression, *variables) for expression in expressions]

    reduced_dimension = len(variables)
    assert reduced_dimension == 47
    assert len(steps) == 22
    assert max(poly.total_degree() for poly in polynomials) == 3
    assert all(poly.eval(dict(zip(variables, [0] * reduced_dimension))) == 0 for poly in polynomials)
    linear_jacobian = sp.Matrix([poly.as_expr() for poly in polynomials]).jacobian(variables)
    assert linear_jacobian.subs(dict.fromkeys(variables, 0)) == sp.eye(reduced_dimension)

    quadratic = [homogeneous_part(poly, 2) for poly in polynomials]
    cubic = [homogeneous_part(poly, 3) for poly in polynomials]
    reconstructed = [
        sp.Poly(variables[index] + quadratic[index].as_expr() + cubic[index].as_expr(), *variables)
        for index in range(reduced_dimension)
    ]
    assert reconstructed == polynomials

    # Verify the transported collision before homogenization.
    reduced_images = []
    for point in collision_points:
        image = [poly.eval(dict(zip(variables, point))) for poly in polynomials]
        reduced_images.append(image)
    expected_reduced_target = normalized_target + [sp.Integer(0)] * (reduced_dimension - 3)
    assert all(image == expected_reduced_target for image in reduced_images)

    # Final cubic-homogeneous map K=I+H in variables (X,Y,t):
    #   K_X = X - t^2 Y + t Q(X)
    #   K_Y = Y + C(X)
    #   K_t = t.
    final_dimension = 2 * reduced_dimension + 1
    h_components: list[list[dict[str, object]]] = [[] for _ in range(final_dimension)]
    for index in range(reduced_dimension):
        h_components[index].append(
            {"coefficient": "-1", "monomial": [[reduced_dimension + index, 1], [final_dimension - 1, 2]]}
        )
        for term in sparse_terms(quadratic[index]):
            term = {"coefficient": term["coefficient"], "monomial": list(term["monomial"])}
            term["monomial"].append([final_dimension - 1, 1])
            h_components[index].append(term)
        h_components[reduced_dimension + index] = sparse_terms(cubic[index])

    # Lift the collision through the Segre and final stable-equivalence steps.
    final_collision_points: list[list[sp.Expr]] = []
    final_images: list[list[sp.Expr]] = []
    for point in collision_points:
        c_value = [poly.eval(dict(zip(variables, point))) for poly in cubic]
        lifted = point + [-value for value in c_value] + [sp.Integer(1)]
        final_collision_points.append(lifted)

        image = list(lifted)
        for component_index, terms in enumerate(h_components):
            correction = sp.Integer(0)
            for term in terms:
                monomial_value = sp.Integer(1)
                for variable_index, exponent in term["monomial"]:
                    monomial_value *= lifted[variable_index] ** exponent
                correction += sp.Rational(term["coefficient"]) * monomial_value
            image[component_index] = sp.cancel(image[component_index] + correction)
        final_images.append(image)
    assert final_images[0] == final_images[1] == final_images[2]
    assert len({tuple(point) for point in final_collision_points}) == 3

    # Every nonlinear term is cubic homogeneous.  The determinant-one
    # certificate is the recorded chain of stable equivalences plus
    # det D(X+tQ+t^2C)=det D(reduced_map)(tX)=1.
    assert all(
        sum(exponent for _, exponent in term["monomial"]) == 3
        for component in h_components
        for term in component
    )

    artifact = {
        "format": "sparse-cubic-homogeneous-map-v1",
        "construction": "Bass--Connell--Wright/Yagzhev via Campbell four-step reduction",
        "dimension": final_dimension,
        "linear_part": "identity",
        "map": "K_i = X_i + H_i; omitted H_i are zero",
        "jacobian_determinant": "1",
        "nilpotency_certificate": "H is cubic homogeneous and det(I+JH)=1, hence JH is nilpotent",
        "source_dimension": 3,
        "degree_reduced_dimension": reduced_dimension,
        "degree_lowering_steps": steps,
        "H": h_components,
        "collision_points": [
            [rational_text(value) for value in point] for point in final_collision_points
        ],
        "common_image": [rational_text(value) for value in final_images[0]],
        "statistics": {
            "nonzero_cubic_terms": sum(len(component) for component in h_components),
            "degree_lowering_steps": len(steps),
            "quadratic_terms_before_homogenization": sum(len(sparse_terms(poly)) for poly in quadratic),
            "cubic_terms_before_homogenization": sum(len(sparse_terms(poly)) for poly in cubic),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"PASS: reduced degree 7 to degree 3 in dimension {reduced_dimension}")
    print(f"PASS: constructed K=I+H with cubic homogeneous H in dimension {final_dimension}")
    print("PASS: transported three distinct rational points to one common image")
    print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

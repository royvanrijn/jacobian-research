#!/usr/bin/env python3
"""Dependency-free replay of the 20/40-dimensional identity-slice witnesses."""

from __future__ import annotations

from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
ARTIFACT = ROOT / "artifacts" / "generated-results" / "image_vanishing_counterexamples_20_40.json"

Complex = tuple[Q, Q]
Poly = dict[tuple[int, ...], Complex]
ZERO: Complex = (Q(0), Q(0))
ONE: Complex = (Q(1), Q(0))
I: Complex = (Q(0), Q(1))


def cadd(a: Complex, b: Complex) -> Complex:
    return (a[0] + b[0], a[1] + b[1])


def cmul(a: Complex, b: Complex) -> Complex:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cscale(a: Complex, scalar: Q) -> Complex:
    return (a[0] * scalar, a[1] * scalar)


def cpow(a: Complex, exponent: int) -> Complex:
    answer = ONE
    for _ in range(exponent):
        answer = cmul(answer, a)
    return answer


def add_term(poly: Poly, exponent: tuple[int, ...], coefficient: Complex) -> None:
    value = cadd(poly.get(exponent, ZERO), coefficient)
    if value != ZERO:
        poly[exponent] = value
    elif exponent in poly:
        del poly[exponent]


def multiply(left: Poly, right: Poly) -> Poly:
    answer: Poly = {}
    for a, ca in left.items():
        for b, cb in right.items():
            add_term(answer, tuple(x + y for x, y in zip(a, b)), cmul(ca, cb))
    return answer


def dense(support: list[list[int]], dimension: int) -> tuple[int, ...]:
    exponent = [0] * dimension
    for variable, power in support:
        exponent[variable] = power
    return tuple(exponent)


def restrict_components(source: dict[str, object]) -> list[Poly]:
    n = source["dimension"] - 1
    components: list[Poly] = []
    for stored_component in source["H"][:n]:
        component: Poly = {}
        for term in stored_component:
            exponent = [0] * n
            for variable, power in term["monomial"]:
                if variable < n:
                    exponent[variable] = power
                else:
                    assert variable == n
            add_term(component, tuple(exponent), (Q(term["coefficient"]), Q(0)))
        components.append(component)
    return components


def decode_rational_terms(terms: list[dict[str, object]], dimension: int) -> Poly:
    return {
        dense(term["monomial"], dimension): (Q(term["coefficient"]), Q(0))
        for term in terms
    }


def decode_complex_terms(terms: list[dict[str, object]], dimension: int) -> Poly:
    return {
        dense(term["monomial"], dimension): (
            Q(term["coefficient"]["real"]), Q(term["coefficient"]["imaginary"])
        )
        for term in terms
    }


def evaluate_real(poly: Poly, point: list[Q]) -> Q:
    answer = Q(0)
    for exponent, coefficient in poly.items():
        assert coefficient[1] == 0
        monomial = Q(1)
        for power, value in zip(exponent, point):
            monomial *= value**power
        answer += coefficient[0] * monomial
    return answer


def contraction_polynomial(components: list[Poly]) -> Poly:
    n = len(components)
    answer: Poly = {}
    for output, component in enumerate(components):
        for old, coefficient in component.items():
            exponent = [0] * (2 * n)
            exponent[output] = 1
            exponent[n:] = old
            add_term(answer, tuple(exponent), cscale(coefficient, Q(-1)))
    return answer


def transformed_polynomial(components: list[Poly], symmetric: bool) -> Poly:
    n = len(components)
    answer: Poly = {}
    for output, component in enumerate(components):
        for old, coefficient in component.items():
            partial: Poly = {
                (0,) * (2 * n): coefficient if symmetric else cscale(coefficient, Q(-1))
            }
            output_factor: Poly = {}
            for variable, factor_coefficient in (
                (output, ONE), (n + output, cscale(I, Q(-1)) if symmetric else I)
            ):
                exponent = [0] * (2 * n)
                exponent[variable] = 1
                output_factor[tuple(exponent)] = factor_coefficient
            partial = multiply(partial, output_factor)

            input_i = I if symmetric else cscale(I, Q(-1))
            for variable, power in enumerate(old):
                if not power:
                    continue
                factor: Poly = {}
                for v_power in range(power + 1):
                    exponent = [0] * (2 * n)
                    exponent[variable] = power - v_power
                    exponent[n + variable] = v_power
                    factor[tuple(exponent)] = cscale(cpow(input_i, v_power), Q(comb(power, v_power)))
                partial = multiply(partial, factor)
            for exponent, partial_coefficient in partial.items():
                add_term(answer, exponent, partial_coefficient)
    return answer


def reflected_negative(poly: Poly, n: int) -> Poly:
    return {
        exponent: cscale(coefficient, Q(-((-1) ** sum(exponent[n:]))))
        for exponent, coefficient in poly.items()
    }


def rational_rank(matrix: list[list[Q]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for index in range(len(work)):
            if index != row and work[index][column]:
                factor = work[index][column]
                work[index] = [a - factor * b for a, b in zip(work[index], work[row])]
        row += 1
        if row == len(work):
            break
    return row


def main() -> None:
    source = json.loads(SOURCE.read_text())
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "identity-slice-image-vanishing-counterexamples-20-40-v1"
    assert source["dimension"] == 21 and source["H"][20] == []

    components = restrict_components(source)
    n = len(components)
    assert n == stored["source_dimension"] == 20
    decoded_components = [
        decode_rational_terms(component, n)
        for component in stored["restricted_nonlinear_part"]["components"]
    ]
    assert components == decoded_components
    monomials = sorted({exponent for component in components for exponent in component})
    coefficient_matrix = [
        [component.get(exponent, ZERO)[0] for exponent in monomials]
        for component in components
    ]
    assert rational_rank(coefficient_matrix) == n
    assert stored["restricted_nonlinear_part"]["output_span_rank"] == n
    assert stored["restricted_nonlinear_part"]["linear_identity_output_functionals"] == 0

    points = [[Q(value) for value in point[:n]] for point in source["collision_points"]]
    assert points == [[Q(value) for value in point] for point in stored["collision_points"]]
    images = [
        [point[i] + evaluate_real(components[i], point) for i in range(n)]
        for point in points
    ]
    assert images[0] == images[1] == images[2]
    assert images[0] == [Q(value) for value in stored["common_image"]]
    assert [point[0] for point in points] == [Q(0), Q(1), Q(-1)]

    p = contraction_polynomial(components)
    stored_p = decode_rational_terms(stored["special_image_counterexample"]["p_terms"], 2 * n)
    assert p == stored_p and len(p) == stored["statistics"]["contraction_terms"]

    laplacian = transformed_polynomial(components, symmetric=False)
    stored_laplacian = decode_complex_terms(stored["laplacian_counterexample"]["R_terms"], 2 * n)
    assert laplacian == stored_laplacian
    assert len(laplacian) == stored["statistics"]["expanded_laplacian_terms"]
    assert sorted({sum(exponent) for exponent in laplacian}) == [2, 3, 4]

    symmetric = transformed_polynomial(components, symmetric=True)
    assert symmetric == reflected_negative(laplacian, n)

    print("PASS (stdlib) identity slice: reconstructed all 20 restricted components")
    print("PASS (stdlib) identity slice: excluded further linear identity outputs")
    print("PASS (stdlib) identity slice: replayed the exact collision and coordinate separation")
    print("PASS (stdlib) identity slice: rebuilt the 40-variable contraction polynomial")
    print("PASS (stdlib) identity slice: independently expanded all 628 Laplacian terms")
    print("PASS (stdlib) identity slice: verified P_20(u,v)=-R_20(u,-v)")


if __name__ == "__main__":
    main()

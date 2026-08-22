#!/usr/bin/env python3
"""Exact generic relation between the two diameter-235 affine sections.

The exact rational component is checked separately.  Here we use its compact
Mestre recursion over Q(p,T), triangular cubic ordinates, and the binary
quartic covariant map to test the sparse relation discovered in the p=-294,
T=2 fibre.  No eliminated residual is formed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from screen_mestre_fermigier_two_section_escape import square_root
from verify_mestre_transverse_conic_component_generic_relations import (
    add,
    covariant_point,
    negative,
)


Q = Fraction
P, T, X = sp.symbols("p T X")
SEED_P = -294
SEED_T = 2
# Coefficients refer to the canonical ``(-,+)`` displayed-root order.
SPECIAL_RELATION = (-1, 0, -1, -1, -1, -1, 0, -1, -1, 0, -1, 0)


def polynomial(value: object, coefficients: tuple[int, ...]) -> object:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def setup() -> tuple[object, object, list[object], object, list[object], object]:
    """Return Q(p), its raw Mestre recursion, roots, and selected lines."""

    ground = sp.QQ.frac_field(P)
    b = (P - 66) * (P + 54) * (P**2 + 18 * P + 456)
    k = 3 * P**2 + 4 * P + 1068
    c1 = -polynomial(P, (-7245936, -246096, 4704, -876, 79)) / (2 * b)
    c2 = 5 * polynomial(
        P,
        (14277841000704, 888999349248, -30000699648, 2783206656,
         -278543520, -14802624, 293232, -40272, 1849),
    ) / (16 * b**2)
    c3 = -polynomial(
        P,
        (-32402742765539106816, -2677382748374243328,
         191800563007140864, -2244816893422080, 657002839322880,
         128514275004672, -3626525084544, 276573028032,
         -9211864320, -761052280, 16418244, -1903778, 59329),
    ) / (16 * b**3)
    c4 = (
        25 * (P - 26) * (P - 6) * (P + 6) * (P + 14)
        * (P**2 - 12 * P + 276) * k * (7 * P**2 - 204 * P - 2628)
        * (29 * P**3 + 378 * P**2 + 3132 * P + 177336)
        * (37 * P**3 - 126 * P**2 + 8316 * P - 269352) / (64 * b**4)
    )
    roots = [
        ground.zero,
        ground.one,
        ground.from_sympy(8 - (P + 294) * (3 * P**3 - 314 * P**2 - 7356 * P - 161208) / (4 * b)),
        ground.from_sympy(k * (7 * P**2 - 204 * P - 2628) / (2 * b)),
        ground.from_sympy((P + 14) * (37 * P**3 - 126 * P**2 + 8316 * P - 269352) / (4 * b)),
        ground.from_sympy(sp.Rational(235, 17) - 45 * (P + 294) * (P**3 - 118 * P**2 - 2292 * P - 57416) / (34 * b)),
    ]
    lines = [
        (
            ground.from_sympy(-(29 * P**3 + 378 * P**2 + 3132 * P + 177336) / (3 * (P - 66) * k)),
            ground.from_sympy(-(13 * P**2 + 204 * P + 1908) / (3 * k)),
        ),
        (
            ground.from_sympy(polynomial(P, (-2491329312, -167485968, -6454080, -770760, 9870, -1053, 53)) / (3 * b * k)),
            ground.from_sympy(-2 * (P - 6) * (P + 54) / (3 * k)),
        ),
    ]
    coefficient_ring = ground.poly_ring(T)
    q = sp.Poly(
        X * (X - 1) * (X**4 + c1 * X**3 + c2 * X**2 + c3 * X + c4),
        X,
        domain=coefficient_ring,
    )
    product = (
        sp.Poly(q.as_expr().subs(X, X - T), X, domain=coefficient_ring)
        * sp.Poly(q.as_expr().subs(X, X + T), X, domain=coefficient_ring)
    )
    approximant = sp.Poly(X**6, X, domain=coefficient_ring)
    for lower_degree in range(5, -1, -1):
        target = product.coeff_monomial(X ** (6 + lower_degree))
        square = (approximant * approximant).coeff_monomial(X ** (6 + lower_degree))
        approximant += sp.Poly((target - square) * X**lower_degree / 2, X, domain=coefficient_ring)
    remainder = approximant * approximant - product
    if any(remainder.coeff_monomial(X**degree) for degree in range(5, 13)):
        raise AssertionError("the generic Mestre recursion left a degree-five remainder")
    quartic = [
        sp.cancel(remainder.coeff_monomial(X**degree).as_expr() / T**2)
        for degree in range(5)
    ]
    square_root_d = ground.from_sympy(
        15 * (P - 26) * (P - 6) * (P + 6) * (P + 14) * k
        / ((P - 66) * (P + 54) * (P**2 + 18 * P + 456) ** 2)
    )
    return ground, approximant, quartic, remainder, roots, lines, square_root_d


def triangular_ordinate(
    ground: object, remainder: object, intercept: object, slope: object, root_d: object
) -> object:
    """Selected rational cubic ordinate, constructed from the triangular recursion."""

    line = sp.Poly(
        remainder.as_expr().subs(X, intercept.as_expr() + slope.as_expr() * T) / T**2,
        T,
        domain=ground,
    )
    f = [ground.from_sympy(line.coeff_monomial(T**degree)) for degree in range(7)]
    y3 = (ground.one - slope**2) * root_d / ground.convert(2)
    y2 = f[5] / (ground.convert(2) * y3)
    y1 = (f[4] - y2**2) / (ground.convert(2) * y3)
    y0 = (f[3] - ground.convert(2) * y1 * y2) / (ground.convert(2) * y3)
    ordinate = y0.as_expr() + y1.as_expr() * T + y2.as_expr() * T**2 + y3.as_expr() * T**3
    values = (y0, y1, y2, y3)
    if any(
        f[degree] != sum(values[left] * values[degree - left] for left in range(max(0, degree - 3), min(3, degree) + 1))
        for degree in range(7)
    ):
        raise AssertionError("a selected triangular ordinate failed")
    return ordinate


def specialize(value: object) -> Q:
    return Q(sp.cancel(value.as_expr()).subs({P: SEED_P, T: SEED_T}))


def orientation(ordinate: object) -> int:
    value = specialize(ordinate)
    if value == 0 or square_root(value * value) is None:
        raise AssertionError("the orientation specialization is not a nonzero rational ordinate")
    return 1 if value > 0 else -1


def signed_sum(field: object, coefficient_a: object, terms: list[tuple[tuple[object, object], int]]) -> tuple[object, object] | None:
    """Use a balanced exact addition tree to limit rational-function swell."""

    level = [point if sign > 0 else negative(point) for point, sign in terms]
    while len(level) > 1:
        next_level = []
        for index in range(0, len(level) - 1, 2):
            next_level.append(add(field, coefficient_a, level[index], level[index + 1]))
        if len(level) % 2:
            next_level.append(level[-1])
        level = next_level
    return level[0] if level else None


def replay() -> dict[str, object]:
    ground, approximant, quartic, remainder, roots, lines, root_d = setup()
    field = sp.QQ.frac_field(P, T)
    raw_quartic = [field.from_sympy(value) for value in quartic]
    visible = []
    visible_orientations = []
    for root in roots:
        for sign in (-1, 1):
            x_value = field.convert(root) + field.convert(sign) * field.convert(T)
            y_value = field.from_sympy(approximant.as_expr().subs(X, x_value.as_expr()) / T)
            visible.append(covariant_point(field, raw_quartic, x_value, y_value))
            visible_orientations.append(orientation(y_value))
    affine = []
    affine_orientations = []
    for intercept, slope in lines:
        ordinate = field.from_sympy(triangular_ordinate(ground, remainder, intercept, slope, root_d))
        x_value = field.convert(intercept) + field.convert(slope) * field.convert(T)
        point = covariant_point(
            field, raw_quartic, x_value, ordinate
        )
        affine.append(point)
        affine_orientations.append(orientation(ordinate))
    invariant_i = 12 * raw_quartic[4] * raw_quartic[0] - 3 * raw_quartic[3] * raw_quartic[1] + raw_quartic[2] ** 2
    coefficient_a = -27 * invariant_i
    left = signed_sum(
        field,
        coefficient_a,
        [(affine[0], affine_orientations[0]), (affine[1], affine_orientations[1])],
    )
    right = signed_sum(
        field,
        coefficient_a,
        [(point, coefficient * orientation_sign) for point, coefficient, orientation_sign in zip(visible, SPECIAL_RELATION, visible_orientations) if coefficient],
    )
    if left != right:
        raise AssertionError("the candidate generic P1+P2 visible relation failed")
    return {
        "status": "exact generic relation between the diameter-235 affine sections verified",
        "base_field": "Q(p,T)",
        "method": "recursive Mestre square root, triangular affine ordinates, covariant map, and short-Weierstrass addition",
        "specialization": {"p": str(SEED_P), "T": str(SEED_T)},
        "specialized_visible_coefficients_for_P1_plus_P2": list(SPECIAL_RELATION),
        "visible_orientation_signs": visible_orientations,
        "affine_orientation_signs": affine_orientations,
        "conclusion": "P1 and P2 differ by a generic visible-subgroup combination, so the displayed pair supplies at most one new generic direction",
        "not_established": [
            "whether P1 is generically independent of the visible subgroup",
            "generic rank at least 14, saturation, pair intersections, or a Shioda Gram matrix",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

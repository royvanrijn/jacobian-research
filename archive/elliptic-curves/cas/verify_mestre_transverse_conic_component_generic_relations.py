#!/usr/bin/env python3
"""Exact generic visible-subgroup relation on the second two-section component.

The computation is over ``Q(s,T)``.  It reconstructs the Mestre quartic by
the recursive monic-square procedure and obtains the selected cubic ordinates
from the same triangular recursion.  It then applies the binary-quartic
covariant map and short-Weierstrass addition directly in the fraction field.
No universal two-section residual is expanded or materialized.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


S, T, X = sp.symbols("s T X")


def setup() -> tuple[object, object, object, list[object], list[object], object]:
    """Return the raw quartic recursion data over Q(s)[T,X]."""

    ground = sp.QQ.frac_field(S)
    u = ground.from_sympy(-12 * (47 * S + 357) / (S**2 - 49))
    c1 = -ground.convert(sp.Rational(425, 7)) + u
    c2 = ground.convert(sp.Rational(66335, 49)) - ground.convert(sp.Rational(929, 21)) * u + ground.convert(sp.Rational(13, 36)) * u**2
    c3 = -ground.convert(sp.Rational(4501495, 343)) + ground.convert(sp.Rational(93718, 147)) * u - ground.convert(sp.Rational(2599, 252)) * u**2 + u**3 / ground.convert(18)
    c4 = ground.convert(sp.Rational(112212864, 2401)) - ground.convert(sp.Rational(1029264, 343)) * u + ground.convert(sp.Rational(63671, 882)) * u**2 - ground.convert(sp.Rational(583, 756)) * u**3 + u**4 / ground.convert(324)
    roots = [
        ground.zero,
        ground.one,
        (ground.convert(411) - 7 * u) / ground.convert(21),
        (ground.convert(474) - 7 * u) / ground.convert(42),
        (ground.convert(1254) - 21 * u - (ground.convert(282) + ground.convert(S) * u)) / ground.convert(84),
        (ground.convert(1254) - 21 * u + (ground.convert(282) + ground.convert(S) * u)) / ground.convert(84),
    ]
    coefficient_ring = ground.poly_ring(T)
    q = sp.Poly(
        X * (X - 1) * (X**4 + c1.as_expr() * X**3 + c2.as_expr() * X**2 + c3.as_expr() * X + c4.as_expr()),
        X, domain=coefficient_ring,
    )
    shift_minus = sp.Poly(q.as_expr().subs(X, X - T), X, domain=coefficient_ring)
    shift_plus = sp.Poly(q.as_expr().subs(X, X + T), X, domain=coefficient_ring)
    product = shift_minus * shift_plus
    approximant = sp.Poly(X**6, X, domain=coefficient_ring)
    for lower_degree in range(5, -1, -1):
        target = product.coeff_monomial(X ** (6 + lower_degree))
        square = (approximant * approximant).coeff_monomial(X ** (6 + lower_degree))
        approximant += sp.Poly((target - square) * X**lower_degree / 2, X, domain=coefficient_ring)
    remainder = approximant * approximant - product
    quartic = []
    for degree in range(5):
        coefficient = remainder.coeff_monomial(X**degree)
        quartic.append(sp.cancel(coefficient.as_expr() / T**2))
    if any(remainder.coeff_monomial(X**degree) for degree in range(5, 13)):
        raise AssertionError("the generic Mestre recursion left a degree-five remainder")
    return ground, u, approximant, roots, quartic, remainder


def section_ordinate(
    ground: object,
    u: object,
    remainder: object,
    intercept: object,
    slope: object,
) -> object:
    """Return the triangular-recursion cubic ordinate on a chosen line."""

    line = sp.Poly(
        remainder.as_expr().subs(X, intercept.as_expr() + slope.as_expr() * T) / T**2,
        T,
        domain=ground,
    )
    f = [
        ground.from_sympy(line.coeff_monomial(T**degree))
        for degree in range(7)
    ]
    root_d = 4 * (7 * u - ground.convert(432)) / ground.convert(21)
    y3 = (ground.one - slope**2) * root_d / ground.convert(2)
    y2 = f[5] / (ground.convert(2) * y3)
    y1 = (f[4] - y2**2) / (ground.convert(2) * y3)
    y0 = (f[3] - ground.convert(2) * y1 * y2) / (ground.convert(2) * y3)
    ordinate = y0.as_expr() + y1.as_expr() * T + y2.as_expr() * T**2 + y3.as_expr() * T**3
    failures = [
        degree for degree in range(7)
        if f[degree] != sum(
            (y0, y1, y2, y3)[left] * (y0, y1, y2, y3)[degree - left]
            for left in range(max(0, degree - 3), min(3, degree) + 1)
        )
    ]
    if failures:
        raise AssertionError(f"the triangular selected-section square failed at {failures}")
    return ordinate


def covariant_point(field: object, quartic: list[object], x_value: object, y_value: object) -> tuple[object, object]:
    """Map an exact raw quartic point to its raw short Jacobian."""

    e, d, c, b, a = quartic
    g0 = b**2 / field.convert(16) - a * c / field.convert(6)
    g1 = b * c / field.convert(12) - a * d / field.convert(2)
    g2 = c**2 / field.convert(12) - b * d / field.convert(8) - a * e
    g3 = c * d / field.convert(12) - b * e / field.convert(2)
    g4 = d**2 / field.convert(16) - c * e / field.convert(6)
    g_value = g0 * x_value**4 + g1 * x_value**3 + g2 * x_value**2 + g3 * x_value + g4
    g_x = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
    g_y = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
    u_x = 4 * a * x_value**3 + 3 * b * x_value**2 + 2 * c * x_value + d
    u_y = b * x_value**3 + 2 * c * x_value**2 + 3 * d * x_value + 4 * e
    h_value = (u_x * g_y - u_y * g_x) / field.convert(8)
    return 36 * g_value / y_value**2, 108 * h_value / y_value**3


def add(field: object, coefficient_a: object, left: tuple[object, object] | None, right: tuple[object, object] | None) -> tuple[object, object] | None:
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    if x_left == x_right:
        slope = (3 * x_left**2 + coefficient_a) / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_sum = slope**2 - x_left - x_right
    return x_sum, -y_left + slope * (x_left - x_sum)


def negative(point: tuple[object, object]) -> tuple[object, object]:
    return point[0], -point[1]


def replay(*, include_second: bool = False) -> dict[str, object]:
    ground, u, approximant, roots, quartic, remainder = setup()
    field = sp.QQ.frac_field(S, T)
    raw_quartic = [field.from_sympy(value) for value in quartic]
    root_values = [field.convert(value) for value in roots]
    g = approximant
    relations = [((0, -1), (1, 1), (3, -1), (4, -1), (5, 1))]
    signs = [1]
    if include_second:
        relations.append(
            ((0, -1), (1, 1), (3, -1))
        )
        signs.append(-1)
    required_visible = set(term for relation in relations for term in relation)
    visible: dict[tuple[int, int], tuple[object, object]] = {}
    for index, root in enumerate(root_values):
        for sign in (-1, 1):
            if (index, sign) not in required_visible:
                continue
            x_value = root + field.convert(sign) * field.convert(T)
            y_value = field.from_sympy(g.as_expr().subs(X, x_value.as_expr()) / T)
            visible[index, sign] = covariant_point(field, raw_quartic, x_value, y_value)
    invariant_i = 12 * raw_quartic[4] * raw_quartic[0] - 3 * raw_quartic[3] * raw_quartic[1] + raw_quartic[2] ** 2
    coefficient_a = -27 * invariant_i
    denominator = field.from_sympy(65 * S**2 + 658 * S + 1813)
    x01 = field.from_sympy((137 * S**2 + 1316 * S + 3283) / denominator.as_expr())
    x11 = field.from_sympy(-(47 * S**2 + 714 * S + 2303) / denominator.as_expr())
    x02 = field.from_sympy((4 * S + 21) * (9 * S + 35) * (137 * S**2 + 1316 * S + 3283) / (7 * (S - 7) * (S + 7) * denominator.as_expr()))
    y01 = field.from_sympy(section_ordinate(ground, u, remainder, ground.from_sympy(x01.as_expr()), ground.from_sympy(x11.as_expr())))
    affine = [covariant_point(field, raw_quartic, x01 + x11 * field.convert(T), y01)]
    if include_second:
        y02 = field.from_sympy(section_ordinate(ground, u, remainder, ground.from_sympy(x02.as_expr()), ground.zero))
        affine.append(covariant_point(field, raw_quartic, x02, y02))
    for relation_index, (target, terms, coefficient) in enumerate(
        zip(affine, relations, signs), start=1
    ):
        total = None
        for term in terms:
            summand = visible[term] if coefficient > 0 else negative(visible[term])
            total = add(field, coefficient_a, total, summand)
        if total != target:
            raise AssertionError(f"the declared generic visible-subgroup relation {relation_index} failed")
    result = {
        "status": "exact generic visible-subgroup relations verified",
        "base_field": "Q(s,T)",
        "method": "recursive Mestre square root, triangular selected ordinates, covariant map, and exact short-Weierstrass addition",
        "relation": "P1=V(0,-)+V(1,+)+V(r3,-)+V(r4,-)+V(r5,+)",
        "conclusion": "the first selected affine section lies in the generic visible subgroup on this component",
        "not_established": [
            "relations for every other two-section component",
            "the visible subgroup's generic rank or saturation",
            "pair intersections or a Shioda Gram matrix",
        ],
    }
    if include_second:
        result["second_relation"] = (
            "P2=-V(0,-)-V(1,+)-V(r3,-)"
        )
        result["conclusion"] = "both selected affine sections lie in the generic visible subgroup on this component"
    else:
        result["not_established"].insert(0, "the visible-subgroup status of the second selected affine section")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-second", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(include_second=args.include_second), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

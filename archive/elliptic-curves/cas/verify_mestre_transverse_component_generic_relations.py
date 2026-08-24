#!/usr/bin/env python3
"""Exact generic visible relations on the component through (0,1,7,8,9,11)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from verify_mestre_transverse_conic_component_generic_relations import (
    add,
    covariant_point,
    negative,
)


R, T, X = sp.symbols("r T X")


def setup() -> tuple[object, object, object, list[object], list[object], object]:
    ground = sp.QQ.frac_field(R)
    ring = ground.poly_ring(T)
    z = ground.from_sympy(12 * (R + 3) / (1 - R**2))
    c1 = -ground.convert(35) + z
    c2 = ground.convert(455) - ground.convert(sp.Rational(77, 3)) * z + ground.convert(sp.Rational(13, 36)) * z**2
    c3 = -ground.convert(2605) + ground.convert(sp.Rational(652, 3)) * z - ground.convert(sp.Rational(217, 36)) * z**2 + z**3 / ground.convert(18)
    c4 = ground.convert(5544) - ground.convert(608) * z + ground.convert(sp.Rational(449, 18)) * z**2 - ground.convert(sp.Rational(49, 108)) * z**3 + z**4 / ground.convert(324)
    roots = [
        ground.zero, ground.one, (ground.convert(33) - z) / ground.convert(3),
        (ground.convert(42) - z) / ground.convert(6),
        (ground.convert(102) - 3 * z - (ground.convert(6) + ground.convert(R) * z)) / ground.convert(12),
        (ground.convert(102) - 3 * z + (ground.convert(6) + ground.convert(R) * z)) / ground.convert(12),
    ]
    q = sp.Poly(X * (X - 1) * (X**4 + c1.as_expr() * X**3 + c2.as_expr() * X**2 + c3.as_expr() * X + c4.as_expr()), X, domain=ring)
    product = sp.Poly(q.as_expr().subs(X, X - T), X, domain=ring) * sp.Poly(q.as_expr().subs(X, X + T), X, domain=ring)
    approximant = sp.Poly(X**6, X, domain=ring)
    for lower in range(5, -1, -1):
        correction = (product.coeff_monomial(X ** (6 + lower)) - (approximant * approximant).coeff_monomial(X ** (6 + lower))) / 2
        approximant += sp.Poly(correction * X**lower, X, domain=ring)
    remainder = approximant * approximant - product
    if any(remainder.coeff_monomial(X**degree) for degree in range(5, 13)):
        raise AssertionError("the generic Mestre recursion left a degree-five remainder")
    return ground, z, approximant, roots, [sp.cancel(remainder.coeff_monomial(X**i).as_expr() / T**2) for i in range(5)], remainder


def ordinate(ground: object, z: object, remainder: object, intercept: object, slope: object) -> object:
    line = sp.Poly(remainder.as_expr().subs(X, intercept.as_expr() + slope.as_expr() * T) / T**2, T, domain=ground)
    f = [ground.from_sympy(line.coeff_monomial(T**degree)) for degree in range(7)]
    y3 = (ground.one - slope**2) * 2 * (ground.convert(36) - z) / 3
    y2 = f[5] / (2 * y3)
    y1 = (f[4] - y2**2) / (2 * y3)
    y0 = (f[3] - 2 * y1 * y2) / (2 * y3)
    coefficients = (y0, y1, y2, y3)
    if any(f[d] != sum(coefficients[i] * coefficients[d-i] for i in range(max(0, d-3), min(3, d)+1)) for d in range(7)):
        raise AssertionError("the triangular selected-section square failed")
    return sum(value.as_expr() * T**degree for degree, value in enumerate(coefficients))


def replay() -> dict[str, object]:
    ground, z, approximant, roots, quartic, remainder = setup()
    field = sp.QQ.frac_field(R, T)
    raw_quartic = [field.from_sympy(value) for value in quartic]
    root_values = [field.convert(value) for value in roots]
    required = {(0, -1), (0, 1), (1, 1), (3, -1)}
    visible = {}
    for index, root in enumerate(root_values):
        for sign in (-1, 1):
            if (index, sign) not in required:
                continue
            x_value = root + field.convert(sign) * field.convert(T)
            y_value = field.from_sympy(approximant.as_expr().subs(X, x_value.as_expr()) / T)
            visible[index, sign] = covariant_point(field, raw_quartic, x_value, y_value)
    invariant_i = 12 * raw_quartic[4] * raw_quartic[0] - 3 * raw_quartic[3] * raw_quartic[1] + raw_quartic[2]**2
    coefficient_a = -27 * invariant_i
    denominator = field.from_sympy(30 - z.as_expr())
    x01 = field.from_sympy((z.as_expr()**2 - 66*z.as_expr() + 1098) / (3 * denominator.as_expr()))
    x11 = field.from_sympy((42 - z.as_expr()) / denominator.as_expr())
    x02 = field.from_sympy((z.as_expr()**2 - 69*z.as_expr() + 1188) / (6 * denominator.as_expr()))
    p1 = covariant_point(field, raw_quartic, x01 + x11 * field.convert(T), field.from_sympy(ordinate(ground, z, remainder, ground.from_sympy(x01.as_expr()), ground.from_sympy(x11.as_expr()))))
    p2 = covariant_point(field, raw_quartic, x02, field.from_sympy(ordinate(ground, z, remainder, ground.from_sympy(x02.as_expr()), ground.zero)))
    relations = (
        ((visible[0, -1], -1), (visible[0, 1], 1), (visible[1, 1], -1)),
        ((visible[0, -1], 1), (visible[1, 1], 1), (visible[3, -1], 1)),
    )
    for index, (target, terms) in enumerate(zip((p1, p2), relations), start=1):
        total = None
        for point, sign in terms:
            total = add(field, coefficient_a, total, point if sign > 0 else negative(point))
        if total != target:
            if total == negative(target):
                raise AssertionError(f"relation {index} has the opposite triangular-ordinate sign")
            raise AssertionError(f"the declared generic visible-subgroup relation {index} failed")
    return {
        "status": "exact generic visible-subgroup relations verified",
        "base_field": "Q(r,T)",
        "relations": ["P1=-V(0,-)+V(0,+)-V(1,+)", "P2=V(0,-)+V(1,+)+V(r3,-)"],
        "conclusion": "both selected affine sections lie in the generic visible subgroup on this component",
        "not_established": ["relations for other components", "the visible subgroup's rank or saturation", "pair intersections or a Shioda Gram matrix"],
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

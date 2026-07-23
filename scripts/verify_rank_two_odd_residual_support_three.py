#!/usr/bin/env python3
"""Eliminate exact-support-three directions inside the residual L2 space."""

from __future__ import annotations

from collections import Counter
from itertools import combinations

import sympy as sp
from sympy.polys.domains import QQ

from explore_rank_two_odd_mixed_quantization import (
    all_scale_fourth_order_span_audit,
    essential_problem,
    extend_sparse_poly,
    linear_combination,
)
from explore_rank_two_odd_quantization import third_order_axis_test
from explore_rank_two_odd_residual_five_space import (
    RESIDUAL_VECTORS,
    exact_scale_locus,
)


def residual_support_three_candidates():
    """Intersect V(bQ) with every exact-support-three residual chart."""

    linear, quadric = exact_scale_locus()
    scale_polynomial = sp.expand(linear * quadric)
    z = sp.symbols("z0:5")
    parameter = sp.Symbol("t")
    x0_coefficients = [
        vector[0] for vector in RESIDUAL_VECTORS
    ]
    contained_lines = []
    candidates = []

    def intersect(kind, indices, coordinates):
        restriction = sp.Poly(
            sp.expand(
                scale_polynomial.subs(dict(zip(z, coordinates)))
            ),
            parameter,
            domain=sp.QQ,
        )
        if restriction.is_zero:
            contained_lines.append((kind, indices, coordinates))
            return
        selected = [coordinates[index] for index in indices]
        x0_value = sum(
            sp.Rational(
                x0_coefficients[index].numerator,
                x0_coefficients[index].denominator,
            )
            * coordinates[index]
            for index in range(5)
        )
        for factor, _multiplicity in sp.factor_list(restriction)[1]:
            # Coordinate or x0 zeros are lower-support boundary points.
            if any(
                sp.rem(
                    sp.Poly(expression, parameter, domain=sp.QQ),
                    factor,
                ).is_zero
                for expression in selected + [x0_value]
            ):
                continue
            candidates.append(
                {
                    "kind": kind,
                    "indices": indices,
                    "factor": factor.monic(),
                    "coordinates": coordinates,
                }
            )

    # Two nonzero residual coordinates and nonzero x0 give exact support 3.
    for left, right in combinations(range(5), 2):
        coordinates = [sp.S.Zero] * 5
        coordinates[left] = sp.S.One
        coordinates[right] = parameter
        intersect("pair", (left, right), coordinates)

    # Three nonzero residual coordinates give exact support 3 only when x0
    # cancels.  Normalize the first coordinate and solve that linear equation.
    for first, second, third in combinations(range(5), 3):
        coordinates = [sp.S.Zero] * 5
        coordinates[first] = sp.S.One
        coordinates[second] = parameter
        first_x0 = sp.Rational(
            x0_coefficients[first].numerator,
            x0_coefficients[first].denominator,
        )
        second_x0 = sp.Rational(
            x0_coefficients[second].numerator,
            x0_coefficients[second].denominator,
        )
        third_x0 = sp.Rational(
            x0_coefficients[third].numerator,
            x0_coefficients[third].denominator,
        )
        coordinates[third] = -(
            first_x0 + second_x0 * parameter
        ) / third_x0
        intersect(
            "x0-cancelled triple",
            (first, second, third),
            coordinates,
        )
    return contained_lines, candidates


def extend_problem(S, T, pairs, field):
    if field == QQ:
        return S, T, pairs
    return (
        extend_sparse_poly(S, field),
        extend_sparse_poly(T, field),
        [
            (
                extend_sparse_poly(s_part, field),
                extend_sparse_poly(t_part, field),
            )
            for s_part, t_part in pairs
        ],
    )


def main() -> None:
    contained_lines, candidates = residual_support_three_candidates()
    assert len(contained_lines) == 1
    assert contained_lines[0][0:2] == ("pair", (0, 1))
    assert len(candidates) == 12
    assert Counter(
        candidate["factor"].degree()
        for candidate in candidates
    ) == Counter({1: 3, 2: 9})
    assert Counter(
        candidate["kind"] for candidate in candidates
    ) == Counter({"pair": 12})

    parameter = sp.Symbol("t")
    S, T, _, pairs, _, _ = essential_problem()
    rank_profile = Counter()
    fourth_profile = Counter()
    for candidate in candidates:
        factor = candidate["factor"]
        if factor.degree() == 1:
            root_expression = -factor.nth(0) / factor.nth(1)
            field = QQ
        else:
            root_expression = sp.RootOf(factor.as_expr(), 0)
            field = QQ.algebraic_field(root_expression)
        residual_coordinates = [
            field.from_sympy(
                sp.sympify(coordinate).subs(
                    parameter,
                    root_expression,
                )
            )
            for coordinate in candidate["coordinates"]
        ]
        global_coefficients = {}
        for vector, residual_coefficient in zip(
            RESIDUAL_VECTORS,
            residual_coordinates,
        ):
            for index, coefficient in vector.items():
                global_coefficients[index] = (
                    global_coefficients.get(index, field.zero)
                    + residual_coefficient * field.convert(coefficient)
                )
        global_coefficients = {
            index: coefficient
            for index, coefficient in global_coefficients.items()
            if coefficient
        }
        assert len(global_coefficients) == 3

        field_S, field_T, field_pairs = extend_problem(
            S,
            T,
            pairs,
            field,
        )
        direction = linear_combination(
            field_pairs,
            global_coefficients,
            field,
        )
        third_result = third_order_axis_test(
            field_S,
            field_T,
            direction,
            field,
        )
        rank_profile[
            (
                factor.degree(),
                third_result[0],
                third_result[1],
                third_result[2],
                third_result[3],
            )
        ] += 1
        fourth_result = all_scale_fourth_order_span_audit(
            field_S,
            field_T,
            direction,
            field,
        )
        normalized_fourth = (
            fourth_result[0] == field.zero,
            *fourth_result[1:],
        )
        fourth_profile[
            (factor.degree(), *normalized_fourth)
        ] += 1
        print(
            f"residual pair {candidate['indices']} "
            f"degree={factor.degree()}: third={third_result}; "
            f"fourth={normalized_fourth}"
        )

    assert rank_profile == Counter(
        {
            (1, True, "u is free", 2572, 63): 3,
            (2, True, "u is free", 2572, 63): 9,
        }
    )
    assert fourth_profile == Counter(
        {
            (1, True, 63, 2079, 626, 627): 3,
            (2, True, 63, 2079, 626, 627): 9,
        }
    )
    print(
        "PASS: the only contained residual support-three line is the "
        "already-eliminated (e0,e1) line"
    )
    print(
        "PASS: the other 12 closed-point classes all have exact "
        "next-obstruction rank 626->627"
    )
    print(
        "PASS: every exact-support-three branch inside the residual L2 "
        "space is eliminated"
    )


if __name__ == "__main__":
    main()

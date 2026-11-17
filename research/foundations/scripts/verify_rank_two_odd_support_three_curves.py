#!/usr/bin/env python3
"""Continue every closed point left by the support-three curve relaxation."""

from __future__ import annotations

from collections import Counter

import sympy as sp
from sympy.polys.domains import QQ

from explore_rank_two_odd_mixed_quantization import (
    all_scale_fourth_order_span_audit,
    essential_problem,
    extend_sparse_poly,
    linear_combination,
)
from explore_rank_two_odd_quantization import third_order_axis_test
from explore_rank_two_odd_support_three_curves import (
    support_three_curve_candidates,
)


def candidate_direction(candidate):
    """Realize a candidate over its exact rational or quadratic residue field."""

    parameter = sp.Symbol("t")
    kind = candidate["kind"]
    factor = sp.Poly(candidate["factor"], parameter, domain=sp.QQ)
    if kind == "finite":
        if factor.degree() == 1:
            root_expression = -factor.nth(0) / factor.nth(1)
            field = QQ
        else:
            assert factor.degree() == 2
            root_expression = sp.RootOf(factor.as_expr(), 0)
            field = QQ.algebraic_field(root_expression)
        coordinates = tuple(
            field.from_sympy(
                sp.sympify(coordinate).subs(parameter, root_expression)
            )
            for coordinate in candidate["coordinates"]
        )
    else:
        field = QQ
        coordinates = tuple(
            field.convert(coordinate)
            for coordinate in candidate["coordinates"]
        )
    assert all(coordinates)
    coefficients = {
        global_index: coordinate
        for global_index, coordinate in zip(
            candidate["triple"],
            coordinates,
        )
    }
    return field, coefficients


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
    dimension_profile, candidates, contained = (
        support_three_curve_candidates()
    )
    assert dimension_profile == Counter({0: 18, 1: 110, 2: 21})
    assert not contained
    assert len(candidates) == 23

    S, T, _, pairs, _, _ = essential_problem()
    residue_degree_profile = Counter()
    rank_profile = Counter()
    survivors = []
    for candidate in candidates:
        field, coefficients = candidate_direction(candidate)
        degree = field.mod.degree() if field != QQ else 1
        residue_degree_profile[degree] += 1
        field_S, field_T, field_pairs = extend_problem(
            S,
            T,
            pairs,
            field,
        )
        direction = linear_combination(
            field_pairs,
            coefficients,
            field,
        )
        survives, reason, rank, nullity = third_order_axis_test(
            field_S,
            field_T,
            direction,
            field,
        )
        rank_profile[(degree, rank, nullity, reason)] += 1
        print(
            f"{candidate['triple']} degree={degree}: "
            f"rank={rank}, nullity={nullity}; {reason}"
        )
        if survives:
            survivors.append(
                (
                    candidate,
                    field,
                    coefficients,
                    direction,
                    rank,
                    nullity,
                    reason,
                )
            )

    assert residue_degree_profile == Counter({1: 4, 2: 19})
    assert rank_profile == Counter(
        {
            (1, 2561, 74, "u is free"): 4,
            (2, 2571, 64, "u is forced to zero"): 19,
        }
    )
    expected_survivors = {
        (
            (1, 11, 12),
            (
                (1, QQ(331) / QQ(126)),
                (11, QQ.one),
                (12, -QQ.one / QQ(2)),
            ),
        ),
        (
            (2, 13, 14),
            (
                (2, -QQ(2948) / QQ(189)),
                (13, -QQ(8) / QQ(3)),
                (14, QQ.one),
            ),
        ),
        (
            (3, 15, 16),
            (
                (3, -QQ(305) / QQ(9)),
                (15, -QQ(10) / QQ(3)),
                (16, QQ.one),
            ),
        ),
        (
            (8, 21, 22),
            (
                (8, QQ(2) / QQ(3)),
                (21, QQ.one),
                (22, -QQ(3) / QQ(4)),
            ),
        ),
    }
    assert {
        (
            candidate["triple"],
            tuple(sorted(coefficients.items())),
        )
        for candidate, _, coefficients, _, _, _, _ in survivors
    } == expected_survivors

    next_order_profile = Counter()
    for (
        _candidate,
        field,
        _coefficients,
        direction,
        _rank,
        _nullity,
        _reason,
    ) in survivors:
        assert field == QQ
        next_order_profile[
            all_scale_fourth_order_span_audit(
                S,
                T,
                direction,
                QQ,
            )
        ] += 1
    assert next_order_profile == Counter(
        {(QQ.zero, 74, 2849, 646, 647): 4}
    )

    print(f"rank profile={sorted(rank_profile.items(), key=str)}")
    print(
        "surviving closed-point candidates="
        f"{[(candidate['triple'], reason) for candidate, _, _, _, _, _, reason in survivors]}"
    )
    print(
        "PASS: all four rational survivors have 74 lower-lift parameters "
        "and exact next-obstruction rank 646->647"
    )
    print(
        "PASS: every positive-dimensional exact-support-three chart is "
        "eliminated"
    )


if __name__ == "__main__":
    main()

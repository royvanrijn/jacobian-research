#!/usr/bin/env python3
"""Eliminate the complete residual odd line, including its exceptional point."""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import QQ

from explore_rank_two_odd_mixed_quantization import (
    all_scale_fourth_order_span_audit,
    all_scale_fourth_order_span_data,
    essential_problem,
    extend_sparse_poly,
    linear_combination,
)
from explore_rank_two_odd_quantization import third_order_axis_test


RESIDUE_WEIGHTS = {
    (18, 0, 3): sp.Rational(1, 211680),
    (14, 0, 1): sp.Rational(19, 1270080),
    (16, 0, 2): sp.Rational(1, 105840),
}


def residue(poly, field):
    return sum(
        (
            field.convert(weight) * poly.get(monomial, field.zero)
            for monomial, weight in RESIDUE_WEIGHTS.items()
        ),
        field.zero,
    )


def denominator_factors(polys, kernel, parameter):
    """Return the irreducible parameter poles of the chosen generic basis."""

    factors = {}

    def include(value):
        if not value:
            return
        denominator = sp.Poly(
            value.denom.as_expr(),
            parameter,
            domain=sp.QQ,
        ).monic()
        for factor, exponent in sp.factor_list(denominator)[1]:
            normalized = sp.Poly(
                factor,
                parameter,
                domain=sp.QQ,
            ).monic()
            factors[normalized.as_expr()] = max(
                factors.get(normalized.as_expr(), 0),
                exponent,
            )

    for second_pair, third_pair, u_value in kernel:
        for poly in second_pair + third_pair:
            for value in poly.values():
                include(value)
        include(u_value)
    for poly in polys:
        for value in poly.values():
            include(value)
    return factors


def main() -> None:
    parameter = sp.Symbol("r")
    field = QQ.frac_field(parameter)
    parameter_value = field.from_sympy(parameter)

    S, T, _, pairs, _, _ = essential_problem()
    field_S = extend_sparse_poly(S, field)
    field_T = extend_sparse_poly(T, field)
    field_pairs = [
        (
            extend_sparse_poly(s_part, field),
            extend_sparse_poly(t_part, field),
        )
        for s_part, t_part in pairs
    ]

    # d_1 + r*d_2, where
    # d_1=x_9+(9/5)x_0 and d_2=x_10+(26/15)x_0.
    coefficients = {
        0: field(9) / field(5)
        + parameter_value * field(26) / field(15),
        9: field.one,
        10: parameter_value,
    }
    direction = linear_combination(field_pairs, coefficients, field)
    (
        base_u,
        lower_kernel,
        nonconstant,
        correction_columns,
        constant,
    ) = all_scale_fourth_order_span_data(
        field_S,
        field_T,
        direction,
        field,
    )
    assert base_u == field.zero
    assert len(lower_kernel) == 63
    assert len(nonconstant) == 2079
    assert len(correction_columns) == 614
    assert all(
        residue(column, field) == field.zero
        for column in correction_columns
    )
    assert all(
        residue(column, field) == field.zero
        for column in nonconstant
    )
    assert residue(constant, field) == field.one

    factors = denominator_factors(
        nonconstant + [constant],
        lower_kernel,
        parameter,
    )
    assert factors == {parameter + sp.Rational(3, 4): 2}

    # The unique pole of the generic lower-lift basis is r=-3/4.
    # There the direction is (1/2)x_0+x_9-(3/4)x_10 and the lower
    # system gains eleven parameters.  Audit it without specialization
    # of the generic basis.
    exceptional_coefficients = {
        0: QQ(1) / QQ(2),
        9: QQ.one,
        10: -QQ(3) / QQ(4),
    }
    rational_S, rational_T, _, rational_pairs, _, _ = essential_problem()
    exceptional_direction = linear_combination(
        rational_pairs,
        exceptional_coefficients,
        QQ,
    )
    survives, reason, rank, nullity = third_order_axis_test(
        rational_S,
        rational_T,
        exceptional_direction,
        QQ,
    )
    assert survives
    assert reason == "u is free"
    assert (rank, nullity) == (2561, 74)
    exceptional_audit = all_scale_fourth_order_span_audit(
        rational_S,
        rational_T,
        exceptional_direction,
        QQ,
    )
    assert exceptional_audit == (
        QQ.zero,
        74,
        2849,
        646,
        647,
    )

    print(
        "PASS: a fixed three-monomial residue annihilates all 2693 "
        "generic enlarged-span columns over QQ(r)"
    )
    print(
        "PASS: its obstruction value is 1 and its only lower-basis pole is "
        "r=-3/4 with multiplicity 2"
    )
    print(
        "PASS: at r=-3/4 the exact lower-lift nullity is 74 and the "
        "next-obstruction rank jumps 646->647"
    )
    print(
        "PASS: every projective parameter on the residual line is eliminated"
    )


if __name__ == "__main__":
    main()

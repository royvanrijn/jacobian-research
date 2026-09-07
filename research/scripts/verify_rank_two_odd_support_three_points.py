#!/usr/bin/env python3
"""Eliminate the isolated exact-support-three odd quantization branches."""

from __future__ import annotations

from collections import Counter

import sympy as sp
from sympy.polys.domains import QQ

from explore_degree_five_a2_subprincipal import (
    add,
    filtered_monomials,
    pi_power,
    poisson,
    scale,
)
from explore_degree_five_quantum_residue import column_rank, solve_affine
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    essential_problem,
    extend_sparse_poly,
    linear_combination,
    multiprime_all_scale_audit,
    split_correction,
)
from explore_rank_two_odd_quantization import third_order_axis_test
from explore_rank_two_odd_support_three import exact_support_three_audit


def lower_lift_data(S, T):
    """Compute the rational u=0 affine base, kernel, and third differential."""

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    second_columns = [
        poisson({monomial: QQ.one}, T) for monomial in s2_monomials
    ]
    second_columns += [
        poisson(S, {monomial: QQ.one}) for monomial in t2_monomials
    ]
    base_rhs = scale(pi_power(S, T, 3), -QQ.one / QQ(24))
    particular, kernel, rank = solve_affine(
        second_columns,
        base_rhs,
        QQ,
    )
    assert rank == 1527
    parity_pair = split_correction(
        particular,
        s2_monomials,
        t2_monomials,
    )
    kernel_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in kernel
    ]
    assert len(kernel_pairs) == 42

    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    third_columns = [
        poisson({monomial: QQ.one}, T) for monomial in s3_monomials
    ]
    third_columns += [
        poisson(S, {monomial: QQ.one}) for monomial in t3_monomials
    ]
    assert len(third_columns) == 1065
    return parity_pair, kernel_pairs, third_columns


def closed_point_directions(points):
    """Yield one exact direction over each residue field of a closed point."""

    r_symbol, s_symbol = sp.symbols("r s")
    for triple, equations in points.items():
        univariate = max(
            (
                sp.Poly(equation, r_symbol, domain=QQ)
                for equation in equations
                if equation.free_symbols <= {r_symbol}
            ),
            key=lambda polynomial: polynomial.degree(),
        )
        for factor, _multiplicity in sp.factor_list(univariate)[1]:
            if factor.degree() == 1:
                field = QQ
                r_expression = -factor.nth(0) / factor.nth(1)
                r_value = field.convert(r_expression)
            else:
                assert factor.degree() == 2
                r_expression = sp.RootOf(factor.as_expr(), 0)
                field = QQ.algebraic_field(r_expression)
                r_value = field.from_sympy(r_expression)

            specialized = [
                sp.Poly(
                    equation.subs(r_symbol, r_expression),
                    s_symbol,
                    domain=field,
                )
                for equation in equations
                if equation.subs(r_symbol, r_expression) != 0
            ]
            common = specialized[0]
            for polynomial in specialized[1:]:
                common = sp.gcd(common, polynomial)
            common = common.sqf_part().monic()
            assert common.degree() == 1
            s_value = field.convert(
                -common.nth(0) / common.nth(1)
            )
            yield (
                triple,
                factor.degree(),
                field,
                {
                    triple[0]: field.one,
                    triple[1]: r_value,
                    triple[2]: s_value,
                },
            )


def extend_pair(pair, field):
    return (
        extend_sparse_poly(pair[0], field),
        extend_sparse_poly(pair[1], field),
    )


def zero_scale_rank_test(
    S,
    T,
    pairs,
    parity_pair,
    kernel_pairs,
    third_columns,
    coefficients,
    field,
):
    """Test the third-order equation at u=0 over an exact residue field."""

    direction = linear_combination(pairs, coefficients, field)
    available = list(third_columns)
    available.extend(
        coupling(direction, kernel_pair)
        for kernel_pair in kernel_pairs
    )
    rhs = scale(
        add(
            pi_power(direction[0], T, 3),
            pi_power(S, direction[1], 3),
        ),
        -field.one / field(24),
    )
    rhs = add(rhs, coupling(direction, parity_pair), -field.one)
    rank = column_rank(available)
    augmented_rank = column_rank(available + [rhs])
    return direction, rank, augmented_rank


def main() -> None:
    S, T, _, pairs, _, _ = essential_problem()
    parity_pair, kernel_pairs, third_columns = lower_lift_data(S, T)
    (
        contained,
        monomial_empty,
        saturated_empty,
        curves,
        points,
    ) = exact_support_three_audit()
    assert (
        len(contained),
        len(monomial_empty),
        len(saturated_empty),
        len(curves),
        len(points),
    ) == (66, 7924, 268, 149, 29)

    directions = list(closed_point_directions(points))
    degree_profile = Counter(degree for _, degree, _, _ in directions)
    assert degree_profile == Counter({1: 38, 2: 28})

    rank_profile = Counter()
    compatible = []
    for triple, degree, field, coefficients in directions:
        if field == QQ:
            field_S = S
            field_T = T
            field_pairs = pairs
            field_parity = parity_pair
            field_kernel = kernel_pairs
            field_third = third_columns
        else:
            field_S = extend_sparse_poly(S, field)
            field_T = extend_sparse_poly(T, field)
            field_pairs = [extend_pair(pair, field) for pair in pairs]
            field_parity = extend_pair(parity_pair, field)
            field_kernel = [
                extend_pair(pair, field) for pair in kernel_pairs
            ]
            field_third = [
                extend_sparse_poly(column, field)
                for column in third_columns
            ]
        direction, rank, augmented_rank = zero_scale_rank_test(
            field_S,
            field_T,
            field_pairs,
            field_parity,
            field_kernel,
            field_third,
            coefficients,
            field,
        )
        rank_profile[(degree, rank, augmented_rank)] += 1
        if rank == augmented_rank:
            compatible.append((triple, field, coefficients, direction))

    assert rank_profile == Counter(
        {
            (1, 1045, 1046): 28,
            (1, 1046, 1047): 9,
            (1, 1034, 1034): 1,
            (2, 1045, 1046): 19,
            (2, 1046, 1047): 9,
        }
    )
    assert len(compatible) == 1
    triple, field, coefficients, direction = compatible[0]
    assert triple == (4, 17, 18)
    assert field == QQ
    assert coefficients == {
        4: QQ.one,
        17: QQ(126) / QQ(1963),
        18: -QQ(63) / QQ(3926),
    }

    survives, reason, rank, nullity = third_order_axis_test(
        S,
        T,
        direction,
        QQ,
    )
    assert survives
    assert reason == "u is free"
    assert (rank, nullity) == (2561, 74)

    (
        lower_nullity,
        coefficient_count,
        prime_results,
    ) = multiprime_all_scale_audit(
        direction,
        coefficients,
        {},
    )
    assert (lower_nullity, coefficient_count) == (74, 2849)
    assert prime_results == [
        (31991, 646, 647),
        (32003, 646, 647),
        (65521, 646, 647),
    ]

    print(
        "PASS: the isolated support-three locus has 38 rational and "
        "28 quadratic closed-point classes"
    )
    print(
        "PASS: 65 classes are inconsistent already at u=0; rank profile "
        f"{sorted(rank_profile.items())}"
    )
    print(
        "PASS: the sole survivor "
        "x4+(126/1963)x17-(63/3926)x18 has u free and lower-lift "
        "nullity 74"
    )
    print(
        "PASS: its all-scale hbar^4 span has 2849 coefficients and rank "
        f"646->647 at {[prime for prime, _, _ in prime_results]}"
    )


if __name__ == "__main__":
    main()

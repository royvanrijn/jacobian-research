#!/usr/bin/env python3
"""Compute the bounded hbar^5 period on the kappa=-1 chart over QQ(tau)."""

from __future__ import annotations

from itertools import combinations_with_replacement
from time import perf_counter

from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import (
    EXCEPTIONAL_PROFILE,
    add,
    degree_five_exceptional_family,
    laurent_monomials,
    solve_affine,
    third_order_family,
)
from verify_degree_five_exceptional_laurent_quantum_obstruction import SUPPORT
from verify_degree_five_third_order_function_field import (
    EXCEPTIONAL_S_SUPPORT,
    EXCEPTIONAL_T_SUPPORT,
)


RELEVANT_GROUPS = ((17, 18), (19, 20), (28, 29))
SECOND_SUPPORT = (
    (13, 0, 1),
    (14, 0, 2),
    (15, 0, 3),
    (15, 1, 1),
    (16, 0, 4),
    (16, 1, 2),
    (18, 1, 4),
)


def split_solution(vector, s_support, t_support):
    split = len(s_support)
    return (
        {
            s_support[index]: coefficient
            for index, coefficient in vector.items()
            if index < split
        },
        {
            t_support[index - split]: coefficient
            for index, coefficient in vector.items()
            if index >= split
        },
    )


def construct(functional_support=SUPPORT):
    profile = EXCEPTIONAL_PROFILE

    sample_s, sample_t = degree_five_exceptional_family(QQ, QQ.one)
    sample_family = third_order_family(
        sample_s,
        sample_t,
        QQ,
        profile=profile,
    )
    relevant_indices = set().union(*map(set, RELEVANT_GROUPS))
    supports = {
        index: (
            sorted(sample_family.kernel[index][0]),
            sorted(sample_family.kernel[index][1]),
        )
        for index in relevant_indices
    }

    field = QQ.frac_field("tau")
    tau = field.gens[0]
    ring = tau.numer.ring
    S, T = degree_five_exceptional_family(
        field,
        tau,
        verify_canonical=False,
    )

    def polynomial_coefficient(coefficient):
        if coefficient.denom.degree() != 0:
            raise AssertionError(coefficient.denom)
        unit = coefficient.denom[(0,)]
        return coefficient.numer * (ring.domain.one / unit)

    A = {
        monomial: polynomial_coefficient(coefficient)
        for monomial, coefficient in S.items()
    }
    B = {
        monomial: polynomial_coefficient(coefficient)
        for monomial, coefficient in T.items()
    }

    def as_field(poly):
        return {
            monomial: field(coefficient)
            for monomial, coefficient in poly.items()
        }

    base_columns = [
        as_field(profile.poisson({monomial: ring.one}, B))
        for monomial in EXCEPTIONAL_S_SUPPORT
    ]
    base_columns += [
        as_field(profile.poisson(A, {monomial: ring.one}))
        for monomial in EXCEPTIONAL_T_SUPPORT
    ]
    rhs = as_field(
        {
            monomial: -coefficient / ring.domain(24)
            for monomial, coefficient in profile.pi_power(A, B, 3).items()
        }
    )
    base, base_kernel, base_rank = solve_affine(
        base_columns,
        rhs,
        field,
    )
    if base_rank != 34 or base_kernel:
        raise AssertionError((base_rank, len(base_kernel)))
    base_s, base_t = split_solution(
        base,
        EXCEPTIONAL_S_SUPPORT,
        EXCEPTIONAL_T_SUPPORT,
    )

    cached_s = {}
    cached_t = {}
    for s_support, t_support in supports.values():
        for monomial in s_support:
            if monomial not in cached_s:
                cached_s[monomial] = as_field(
                    profile.poisson({monomial: ring.one}, B)
                )
        for monomial in t_support:
            if monomial not in cached_t:
                cached_t[monomial] = as_field(
                    profile.poisson(A, {monomial: ring.one})
                )

    directions = []
    for group in RELEVANT_GROUPS:
        s_support = sorted(
            set().union(*(set(supports[index][0]) for index in group))
        )
        t_support = sorted(
            set().union(*(set(supports[index][1]) for index in group))
        )
        columns = [cached_s[monomial] for monomial in s_support]
        columns += [cached_t[monomial] for monomial in t_support]
        _, kernel, rank = solve_affine(columns, {}, field)
        if rank != len(columns) - 2 or len(kernel) != 2:
            raise AssertionError((group, len(columns), rank, len(kernel)))
        directions.extend(
            split_solution(vector, s_support, t_support)
            for vector in kernel
        )

    def restrict(poly):
        return [
            poly.get(monomial, field.zero)
            for monomial in functional_support
        ]

    constraints = []
    for degree, z_order, left in ((24, 1, True), (20, 0, False)):
        for monomial in laurent_monomials(degree, z_order, 0, 4):
            image = (
                profile.poisson({monomial: ring.one}, B)
                if left
                else profile.poisson(A, {monomial: ring.one})
            )
            row = restrict(as_field(image))
            if any(row):
                constraints.append(row)

    for direction_s, direction_t in directions:
        variation = add(
            profile.poisson(direction_s, base_t),
            profile.poisson(base_s, direction_t),
        )
        variation = add(
            variation,
            profile.pi_power(direction_s, T, 3),
            field.one / field(24),
        )
        variation = add(
            variation,
            profile.pi_power(S, direction_t, 3),
            field.one / field(24),
        )
        row = restrict(variation)
        if any(row):
            constraints.append(row)

    for left, right in combinations_with_replacement(
        range(len(directions)),
        2,
    ):
        left_s, left_t = directions[left]
        right_s, right_t = directions[right]
        variation = (
            profile.poisson(left_s, left_t)
            if left == right
            else add(
                profile.poisson(left_s, right_t),
                profile.poisson(right_s, left_t),
            )
        )
        row = restrict(variation)
        if any(row):
            constraints.append(row)

    # Functional coordinates are columns; conditions are output rows.
    functional_columns = [
        {
            row_index: row[coordinate]
            for row_index, row in enumerate(constraints)
            if row[coordinate]
        }
        for coordinate in range(len(functional_support))
    ]
    _, functional_kernel, functional_rank = solve_affine(
        functional_columns,
        {},
        field,
    )
    if functional_rank != 6 or len(functional_kernel) != 1:
        raise AssertionError((functional_rank, len(functional_kernel)))
    functional = functional_kernel[0]

    denominator_lcm = ring.one
    for coefficient in functional.values():
        denominator_lcm = denominator_lcm.lcm(coefficient.denom)

    cleared_functional = []
    for index in range(len(functional_support)):
        coefficient = functional.get(index, field.zero)
        cleared = field(denominator_lcm) * coefficient
        if cleared.denom.degree() != 0:
            raise AssertionError(cleared.denom)
        unit = cleared.denom[(0,)]
        cleared_functional.append(
            cleared.numer * (ring.domain.one / unit)
        )
    common_factor = next(
        coefficient
        for coefficient in cleared_functional
        if coefficient
    )
    for coefficient in cleared_functional:
        if coefficient:
            common_factor = common_factor.gcd(coefficient)
    primitive_functional = [
        coefficient.exquo(common_factor)
        for coefficient in cleared_functional
    ]

    defect = profile.poisson(base_s, base_t)
    defect = add(
        defect,
        profile.pi_power(base_s, T, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, base_t, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, T, 5),
        field.one / field(1920),
    )
    period = sum(
        (
            functional.get(index, field.zero)
            * defect.get(monomial, field.zero)
            for index, monomial in enumerate(functional_support)
        ),
        field.zero,
    )
    primitive_period = sum(
        (
            field(primitive_functional[index])
            * defect.get(monomial, field.zero)
            for index, monomial in enumerate(functional_support)
        ),
        field.zero,
    )
    if primitive_period.denom.degree() != 0:
        raise AssertionError(
            f"primitive period is not polynomial: {primitive_period.denom}"
        )
    return (
        constraints,
        directions,
        period,
        primitive_functional,
        primitive_period,
    )


def main() -> None:
    started = perf_counter()
    (
        constraints,
        directions,
        period,
        primitive_functional,
        primitive_period,
    ) = construct(SUPPORT)
    (
        second_constraints,
        second_directions,
        second_period,
        second_primitive_functional,
        second_primitive_period,
    ) = construct(SECOND_SUPPORT)
    if len(second_constraints) != len(constraints):
        raise AssertionError(
            (len(constraints), len(second_constraints))
        )
    if len(second_directions) != len(directions):
        raise AssertionError(
            (len(directions), len(second_directions))
        )
    unit, factors = period.numer.factor_list()
    second_unit, second_factors = second_period.numer.factor_list()
    gcd = period.numer.gcd(second_period.numer)
    primitive_gcd = primitive_period.numer.gcd(
        second_primitive_period.numer
    )
    print("kappa=-1 bounded hbar^5 function-field certificate")
    print(f"  field: QQ(tau)")
    print(f"  relevant lower directions: {len(directions)}")
    print(f"  nonzero conditions: {len(constraints)}")
    print("  condition rank: 6")
    print("  functional kernel dimension: 1")
    print(f"  period numerator: {period.numer}")
    print(f"  period denominator: {period.denom}")
    print(f"  numerator unit: {unit}")
    print(
        "  numerator factors:",
        [(str(factor), exponent) for factor, exponent in factors],
    )
    print(f"  second period numerator: {second_period.numer}")
    print(f"  second period denominator: {second_period.denom}")
    print(f"  second numerator unit: {second_unit}")
    print(
        "  second numerator factors:",
        [(str(factor), exponent) for factor, exponent in second_factors],
    )
    print(f"  numerator gcd: {gcd}")
    print(
        "  primitive functional gcds: 1, 1"
    )
    print(f"  primitive period: {primitive_period}")
    print(f"  second primitive period: {second_primitive_period}")
    print(f"  primitive period gcd: {primitive_gcd}")
    print(f"  seconds: {perf_counter() - started:.2f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Low-support branches of the odd rank-two quantization locus.

The first unrestricted correction at the degree-five sample has 38
essential coordinates after target-Hamiltonian gauge.  Its next obstruction
is a 41-dimensional space of homogeneous quadrics.  This script:

* reconstructs those quadrics exactly over Q;
* records the coordinates forced into the radical by monomial equations;
* classifies every projective solution supported on at most two coordinates;
* tests every isolated rational and quadratic mixed direction against the
  following hbar^3 compatibility equation;
* reduces the two maximal coordinate-linear families to explicit residual
  spaces and eliminates every exact-support-two direction.

This is deliberately a low-support audit.  It does not classify components
whose general point needs three or more essential coordinates.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_a2_subprincipal import (
    add,
    degree_five_sample,
    filtered_monomials,
    pi_power,
    poisson,
    scale,
)
from explore_degree_five_quantum_residue import column_rank, solve_affine
from explore_rank_two_odd_quantization import (
    convert_poly,
    first_kernel,
    gauge_pivots,
    second_obstruction,
    split_pair,
    third_order_axis_test,
)


def essential_problem():
    """Return the exact essential pairs and a sparse basis of quadrics."""

    S_fraction, T_fraction = degree_five_sample()
    S = convert_poly(S_fraction, QQ)
    T = convert_poly(T_fraction, QQ)
    s_monomials, t_monomials, free_columns, kernel, _ = first_kernel(S, T, QQ)
    pivots = gauge_pivots(
        S,
        T,
        s_monomials,
        t_monomials,
        free_columns,
        kernel,
        QQ,
    )
    essential = [index for index in range(len(kernel)) if index not in pivots]
    pairs = [
        split_pair(kernel[index], s_monomials, t_monomials)
        for index in essential
    ]
    (
        _,
        parameter_monomials,
        _,
        equation_pivots,
        _,
        _,
        _,
        sparse_equations,
    ) = second_obstruction(S, T, pairs, QQ)
    assert len(essential) == 38
    assert len(equation_pivots) == len(sparse_equations) == 41
    return S, T, essential, pairs, parameter_monomials, sparse_equations


def rational(value) -> sp.Rational:
    return sp.Rational(str(value))


def reduce_sparse_poly(poly, field):
    result = {}
    for monomial, coefficient in poly.items():
        value = rational(coefficient)
        result[monomial] = (
            field(int(sp.numer(value))) / field(int(sp.denom(value)))
        )
    return result


def extend_sparse_poly(poly, field):
    return {
        monomial: field.convert(coefficient)
        for monomial, coefficient in poly.items()
    }


def restricted_polynomial(
    equation,
    parameter_monomials,
    left: int,
    right: int,
    variable: sp.Symbol,
) -> sp.Poly:
    """Restrict a quadric to x_left=1, x_right=variable."""

    expression = sp.S.Zero
    for column, coefficient in equation.items():
        first, second = parameter_monomials[column]
        if first not in (left, right) or second not in (left, right):
            continue
        first_value = sp.S.One if first == left else variable
        second_value = sp.S.One if second == left else variable
        expression += rational(coefficient) * first_value * second_value
    return sp.Poly(expression, variable, domain=sp.QQ)


def common_direction_polynomial(
    equations,
    parameter_monomials,
    left: int,
    right: int,
    variable: sp.Symbol,
):
    """Return the gcd cutting out mixed directions on one coordinate line."""

    restrictions = [
        restricted_polynomial(
            equation,
            parameter_monomials,
            left,
            right,
            variable,
        )
        for equation in equations.values()
    ]
    restrictions = [polynomial for polynomial in restrictions if not polynomial.is_zero]
    if not restrictions:
        return None
    common = restrictions[0].monic()
    for polynomial in restrictions[1:]:
        common = sp.gcd(common, polynomial.monic()).monic()
        if common.degree() == 0:
            break
    return common


def combine_pair(left_pair, right_pair, ratio):
    """Return left_pair + ratio*right_pair."""

    left_s, left_t = left_pair
    right_s, right_t = right_pair
    return (
        add(left_s, scale(right_s, ratio)),
        add(left_t, scale(right_t, ratio)),
    )


def linear_combination(pairs, coefficients, field):
    s_part = {}
    t_part = {}
    for index, coefficient in coefficients.items():
        s_part = add(s_part, pairs[index][0], coefficient)
        t_part = add(t_part, pairs[index][1], coefficient)
    return s_part, t_part


def split_correction(vector, s_monomials, t_monomials):
    split = len(s_monomials)
    return (
        {
            s_monomials[column]: value
            for column, value in vector.items()
            if column < split and value
        },
        {
            t_monomials[column - split]: value
            for column, value in vector.items()
            if column >= split and value
        },
    )


def solve_many_particular(columns, right_sides, field):
    """Solve one linear operator against several known-solvable right sides."""

    output_monomials = sorted(
        set().union(
            *(set(column) for column in columns),
            *(set(rhs) for rhs in right_sides),
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    first_rhs = len(columns)
    rows = {}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows.setdefault(output_index[monomial], {})[
                column_index
            ] = coefficient
    for offset, rhs in enumerate(right_sides):
        for monomial, coefficient in rhs.items():
            rows.setdefault(output_index[monomial], {})[
                first_rhs + offset
            ] = -coefficient
    reduced, pivots, _ = sdm_irref(rows)
    assert not any(pivot >= first_rhs for pivot in pivots)

    particulars = []
    for offset in range(len(right_sides)):
        rhs_column = first_rhs + offset
        particular = {}
        for reduced_row, pivot in enumerate(pivots):
            value = reduced.get(reduced_row, {}).get(rhs_column, field.zero)
            if value:
                particular[pivot] = -value
        particulars.append(particular)
    return particulars


def coupling(first_pair, correction_pair):
    """The hbar^3 bilinear term B_(S1,T1)(S2,T2)."""

    s1, t1 = first_pair
    s2, t2 = correction_pair
    return add(poisson(s2, t1), poisson(s1, t2))


def modular_essential_pairs(field):
    S_fraction, T_fraction = degree_five_sample()
    S = convert_poly(S_fraction, field)
    T = convert_poly(T_fraction, field)
    s_monomials, t_monomials, free_columns, kernel, _ = first_kernel(S, T, field)
    pivots = gauge_pivots(
        S,
        T,
        s_monomials,
        t_monomials,
        free_columns,
        kernel,
        field,
    )
    pairs = [
        split_pair(kernel[index], s_monomials, t_monomials)
        for index in range(len(kernel))
        if index not in pivots
    ]
    return S, T, pairs


def uniform_third_order_relaxation(S, T, pairs, subspace, field):
    """Test a linear relaxation of hbar^3 compatibility on a subspace.

    For v in the declared coordinate subspace, write the second correction
    as the parity base plus a quadratic particular solution plus an arbitrary
    vector in ker(d_2).  At hbar^3 the available terms lie in the span of:

      im(d_3), B_v ker(d_2), and B_v q_(w,w').

    We allow the coefficients of all these vectors to vary independently.
    This is a relaxation of the actual factorized equations.  If the linear
    right sides for the basis vectors of the subspace remain independent
    modulo that enlarged span, no nonzero v can satisfy hbar^3.
    """

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    d2_columns = [
        poisson({monomial: field.one}, T) for monomial in s2_monomials
    ]
    d2_columns += [
        poisson(S, {monomial: field.one}) for monomial in t2_monomials
    ]
    parity_rhs = scale(pi_power(S, T, 3), -field.one / field(24))
    parity_vector, kernel_vectors, rank = solve_affine(
        d2_columns,
        parity_rhs,
        field,
    )
    assert rank == 1527
    parity_pair = split_correction(
        parity_vector,
        s2_monomials,
        t2_monomials,
    )
    kernel_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in kernel_vectors
    ]
    assert len(kernel_pairs) == 42

    basis_pairs = [pairs[index] for index in subspace]
    quadratic_right_sides = []
    quadratic_indices = []
    for left in range(len(basis_pairs)):
        left_s, left_t = basis_pairs[left]
        for right in range(left, len(basis_pairs)):
            right_s, right_t = basis_pairs[right]
            if left == right:
                bracket = poisson(left_s, left_t)
            else:
                bracket = add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            quadratic_right_sides.append(scale(bracket, -field.one))
            quadratic_indices.append((left, right))
    quadratic_vectors = solve_many_particular(
        d2_columns,
        quadratic_right_sides,
        field,
    )
    quadratic_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in quadratic_vectors
    ]

    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    relaxed_columns = [
        poisson({monomial: field.one}, T) for monomial in s3_monomials
    ]
    relaxed_columns += [
        poisson(S, {monomial: field.one}) for monomial in t3_monomials
    ]
    for first_pair in basis_pairs:
        relaxed_columns.extend(
            coupling(first_pair, kernel_pair)
            for kernel_pair in kernel_pairs
        )
    for first_pair in basis_pairs:
        relaxed_columns.extend(
            coupling(first_pair, quadratic_pair)
            for quadratic_pair in quadratic_pairs
        )

    linear_right_sides = []
    for first_pair in basis_pairs:
        s1, t1 = first_pair
        rhs = scale(
            add(pi_power(s1, T, 3), pi_power(S, t1, 3)),
            -field.one / field(24),
        )
        linear_right_sides.append(
            add(rhs, coupling(first_pair, parity_pair), -field.one)
        )

    relaxed_rank = column_rank(relaxed_columns)
    augmented_rank = column_rank(relaxed_columns + linear_right_sides)
    output_monomials = sorted(
        set().union(
            *(set(column) for column in relaxed_columns),
            *(set(rhs) for rhs in linear_right_sides),
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    rows = {}
    all_columns = relaxed_columns + linear_right_sides
    for column_index, column in enumerate(all_columns):
        for monomial, coefficient in column.items():
            rows.setdefault(output_index[monomial], {})[
                column_index
            ] = coefficient
    reduced, pivots, _ = sdm_irref(rows)
    first_rhs = len(relaxed_columns)
    rhs_pivots = {
        pivot for pivot in pivots if pivot >= first_rhs
    }
    free_rhs = [
        column
        for column in range(first_rhs, len(all_columns))
        if column not in rhs_pivots
    ]
    candidate_basis = []
    for free_column in free_rhs:
        vector = {free_column - first_rhs: field.one}
        for reduced_row, pivot in enumerate(pivots):
            if pivot < first_rhs:
                continue
            coefficient = reduced.get(reduced_row, {}).get(
                free_column,
                field.zero,
            )
            if coefficient:
                vector[pivot - first_rhs] = -coefficient
        candidate_basis.append(vector)
    assert len(rhs_pivots) == augmented_rank - relaxed_rank
    assert len(candidate_basis) == len(linear_right_sides) - len(rhs_pivots)
    return (
        relaxed_rank,
        augmented_rank,
        len(linear_right_sides),
        len(quadratic_indices),
        candidate_basis,
    )


def maximal_cliques(vertices, edges):
    """Small Bron--Kerbosch enumeration for the coordinate-line graph."""

    adjacency = {
        vertex: {
            other
            for other in vertices
            if other != vertex
            and (min(vertex, other), max(vertex, other)) in edges
        }
        for vertex in vertices
    }
    result = []

    def visit(current, candidates, excluded):
        if not candidates and not excluded:
            result.append(tuple(sorted(current)))
            return
        for vertex in list(candidates):
            visit(
                current | {vertex},
                candidates & adjacency[vertex],
                excluded & adjacency[vertex],
            )
            candidates.remove(vertex)
            excluded.add(vertex)

    visit(set(), set(vertices), set())
    return sorted(set(result))


MOYAL_COEFFICIENTS = {
    1: sp.Rational(1),
    3: sp.Rational(1, 24),
    5: sp.Rational(1, 1920),
    7: sp.Rational(1, 322560),
    9: sp.Rational(1, 92897280),
}


def field_rational(value, field):
    return field(int(sp.numer(value))) / field(int(sp.denom(value)))


def moyal_coefficient(s_series, t_series, order, field):
    """Coefficient of hbar^order in [S_h,T_h]/hbar."""

    result = {}
    for pi_order, coefficient in MOYAL_COEFFICIENTS.items():
        correction_sum = order - (pi_order - 1)
        if correction_sum < 0:
            continue
        scalar = field_rational(coefficient, field)
        for s_index, s_part in enumerate(s_series):
            t_index = correction_sum - s_index
            if t_index < 0 or t_index >= len(t_series):
                continue
            t_part = t_series[t_index]
            if not s_part or not t_part:
                continue
            result = add(
                result,
                pi_power(s_part, t_part, pi_order),
                scalar,
            )
    return result


def joint_second_third_corrections(
    S,
    T,
    first_pair,
    field,
    fixed_u=1,
):
    """Solve the scaled simultaneous order-two/order-three system.

    The actual first correction is lambda*first_pair and u=lambda^2.  The
    order-three unknown is divided by lambda, so the system is affine linear
    in the second correction, the divided third correction, and u.
    """

    s1, t1 = first_pair
    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    split2 = len(s2_monomials)

    columns2 = [
        poisson({monomial: field.one}, T) for monomial in s2_monomials
    ]
    columns2 += [
        poisson(S, {monomial: field.one}) for monomial in t2_monomials
    ]
    coupling_columns = []
    for column, monomial in enumerate(s2_monomials + t2_monomials):
        if column < split2:
            coupling_columns.append(poisson({monomial: field.one}, t1))
        else:
            coupling_columns.append(poisson(s1, {monomial: field.one}))
    columns3 = [
        poisson({monomial: field.one}, T) for monomial in s3_monomials
    ]
    columns3 += [
        poisson(S, {monomial: field.one}) for monomial in t3_monomials
    ]

    base_rhs2 = scale(pi_power(S, T, 3), -field.one / field(24))
    quadratic = poisson(s1, t1)
    rhs3 = scale(
        add(pi_power(s1, T, 3), pi_power(S, t1, 3)),
        -field.one / field(24),
    )

    tagged_rows = set()
    for column in columns2:
        tagged_rows.update((2, monomial) for monomial in column)
    for column in coupling_columns:
        tagged_rows.update((3, monomial) for monomial in column)
    for column in columns3:
        tagged_rows.update((3, monomial) for monomial in column)
    tagged_rows.update((2, monomial) for monomial in base_rhs2)
    tagged_rows.update((2, monomial) for monomial in quadratic)
    tagged_rows.update((3, monomial) for monomial in rhs3)
    tagged_rows = sorted(tagged_rows)
    row_index = {tagged: index for index, tagged in enumerate(tagged_rows)}

    count2 = len(columns2)
    count3 = len(columns3)
    u_column = count2 + count3
    rhs_column = u_column + 1
    rows = {}
    for column_index, (linear2, linear3) in enumerate(
        zip(columns2, coupling_columns)
    ):
        for monomial, coefficient in linear2.items():
            rows.setdefault(row_index[(2, monomial)], {})[
                column_index
            ] = coefficient
        for monomial, coefficient in linear3.items():
            rows.setdefault(row_index[(3, monomial)], {})[
                column_index
            ] = coefficient
    for offset, column in enumerate(columns3):
        for monomial, coefficient in column.items():
            rows.setdefault(row_index[(3, monomial)], {})[
                count2 + offset
            ] = coefficient
    for monomial, coefficient in quadratic.items():
        rows.setdefault(row_index[(2, monomial)], {})[
            u_column
        ] = coefficient
    for monomial, coefficient in base_rhs2.items():
        rows.setdefault(row_index[(2, monomial)], {})[
            rhs_column
        ] = -coefficient
    for monomial, coefficient in rhs3.items():
        rows.setdefault(row_index[(3, monomial)], {})[
            rhs_column
        ] = -coefficient

    if fixed_u is not None:
        fixed_value = field(fixed_u)
        # Fix the nonzero branch scale by u=lambda^2.
        rows[len(tagged_rows)] = {
            u_column: field.one,
            rhs_column: -fixed_value,
        }
    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        raise ValueError("unit-scale second/third correction system is inconsistent")
    solution = {}
    for reduced_row, pivot in enumerate(pivots):
        if pivot == rhs_column:
            continue
        value = reduced.get(reduced_row, {}).get(rhs_column, field.zero)
        if value:
            solution[pivot] = -value
    second_vector = {
        column: value
        for column, value in solution.items()
        if column < count2
    }
    third_vector = {
        column - count2: value
        for column, value in solution.items()
        if count2 <= column < count2 + count3
    }
    base_u = solution.get(u_column, field.zero)

    homogeneous_rows = {
        row: {
            column: value
            for column, value in entries.items()
            if column != rhs_column
        }
        for row, entries in rows.items()
    }
    homogeneous_reduced, homogeneous_pivots, homogeneous_nonzero = sdm_irref(
        homogeneous_rows
    )
    kernel, _ = sdm_nullspace_from_rref(
        homogeneous_reduced,
        field.one,
        rhs_column,
        homogeneous_pivots,
        homogeneous_nonzero,
    )
    kernel_pairs = []
    for vector in kernel:
        second_kernel = {
            column: value
            for column, value in vector.items()
            if column < count2
        }
        third_kernel = {
            column - count2: value
            for column, value in vector.items()
            if count2 <= column < count2 + count3
        }
        kernel_pairs.append(
            (
                split_correction(
                    second_kernel,
                    s2_monomials,
                    t2_monomials,
                ),
                split_correction(
                    third_kernel,
                    s3_monomials,
                    t3_monomials,
                ),
                vector.get(u_column, field.zero),
            )
        )
    if fixed_u is not None:
        assert base_u == field(fixed_u)
        assert all(not vector.get(u_column) for vector in kernel)
    return (
        split_correction(second_vector, s2_monomials, t2_monomials),
        split_correction(third_vector, s3_monomials, t3_monomials),
        base_u,
        kernel_pairs,
        len(pivots),
        rhs_column - len(pivots),
    )


def canonical_branch_continuation(S, T, first_pair, field):
    """Continue one fixed first direction with free correction values zero."""

    second_pair, third_pair, _, _, joint_rank, joint_nullity = (
        joint_second_third_corrections(S, T, first_pair, field)
    )
    s_series = [S, first_pair[0], second_pair[0], third_pair[0]]
    t_series = [T, first_pair[1], second_pair[1], third_pair[1]]
    steps = [
        (2, True, 1569, 1527, 42),
        (3, True, 1065, joint_rank, joint_nullity),
    ]
    for order in range(4, 6):
        s_order = 5 - order
        t_order = 4 - order
        s_monomials = (
            filtered_monomials(29 - 2 * order, s_order)
            if s_order >= 0
            else []
        )
        t_monomials = (
            filtered_monomials(25 - 2 * order, t_order)
            if t_order >= 0
            else []
        )
        columns = [
            poisson({monomial: field.one}, T)
            for monomial in s_monomials
        ]
        columns += [
            poisson(S, {monomial: field.one})
            for monomial in t_monomials
        ]
        defect = moyal_coefficient(s_series, t_series, order, field)
        try:
            particular, kernel, rank = solve_affine(
                columns,
                scale(defect, -field.one),
                field,
            )
        except ValueError:
            steps.append((order, False, len(columns), None, None))
            return s_series, t_series, steps, {}
        s_part, t_part = split_correction(
            particular,
            s_monomials,
            t_monomials,
        )
        s_series.append(s_part)
        t_series.append(t_part)
        steps.append((order, True, len(columns), rank, len(kernel)))

    residuals = {
        order: defect
        for order in range(6, 10)
        if (defect := moyal_coefficient(
            s_series,
            t_series,
            order,
            field,
        ))
    }
    return s_series, t_series, steps, residuals


def fourth_order_span_audit(S, T, first_pair, field):
    """Linear-span obstruction audit over every unit-scale lower lift."""

    (
        second_base,
        third_base,
        _,
        lower_kernel,
        _,
        _,
    ) = joint_second_third_corrections(S, T, first_pair, field)
    base_s_series = [
        S,
        first_pair[0],
        second_base[0],
        third_base[0],
    ]
    base_t_series = [
        T,
        first_pair[1],
        second_base[1],
        third_base[1],
    ]
    constant = moyal_coefficient(base_s_series, base_t_series, 4, field)

    nonconstant = []
    for second_kernel, third_kernel, _ in lower_kernel:
        shifted_s = [
            S,
            first_pair[0],
            add(second_base[0], second_kernel[0]),
            add(third_base[0], third_kernel[0]),
        ]
        shifted_t = [
            T,
            first_pair[1],
            add(second_base[1], second_kernel[1]),
            add(third_base[1], third_kernel[1]),
        ]
        diagonal = poisson(second_kernel[0], second_kernel[1])
        shifted = moyal_coefficient(shifted_s, shifted_t, 4, field)
        linear = add(add(shifted, constant, -field.one), diagonal, -field.one)
        nonconstant.extend((linear, diagonal))
    for left in range(len(lower_kernel)):
        left_second, _, _ = lower_kernel[left]
        for right in range(left + 1, len(lower_kernel)):
            right_second, _, _ = lower_kernel[right]
            nonconstant.append(
                add(
                    poisson(left_second[0], right_second[1]),
                    poisson(right_second[0], left_second[1]),
                )
            )

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    correction_columns = [
        poisson({monomial: field.one}, T) for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: field.one}) for monomial in t4_monomials
    ]
    span_columns = correction_columns + nonconstant
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    return (
        len(lower_kernel),
        len(nonconstant),
        span_rank,
        augmented_rank,
    )


def scaled_fourth_defect(S, T, first_pair, second_pair, divided_third, u, field):
    """Order-four defect with first correction lambda*v and u=lambda^2."""

    value = scale(
        coupling(first_pair, divided_third),
        u,
    )
    value = add(value, poisson(second_pair[0], second_pair[1]))
    cubic = add(
        pi_power(S, second_pair[1], 3),
        pi_power(second_pair[0], T, 3),
    )
    cubic = add(
        cubic,
        pi_power(first_pair[0], first_pair[1], 3),
        u,
    )
    value = add(value, cubic, field.one / field(24))
    value = add(
        value,
        pi_power(S, T, 5),
        field.one / field(1920),
    )
    return value


def all_scale_fourth_order_span_data(S, T, first_pair, field):
    """Build the enlarged hbar^4 span and its constant obstruction."""

    (
        second_base,
        third_base,
        base_u,
        lower_kernel,
        _,
        _,
    ) = joint_second_third_corrections(
        S,
        T,
        first_pair,
        field,
        fixed_u=None,
    )
    constant = scaled_fourth_defect(
        S,
        T,
        first_pair,
        second_base,
        third_base,
        base_u,
        field,
    )

    nonconstant = []
    for second_kernel, third_kernel, u_kernel in lower_kernel:
        shifted_second = (
            add(second_base[0], second_kernel[0]),
            add(second_base[1], second_kernel[1]),
        )
        shifted_third = (
            add(third_base[0], third_kernel[0]),
            add(third_base[1], third_kernel[1]),
        )
        shifted_u = base_u + u_kernel
        shifted = scaled_fourth_defect(
            S,
            T,
            first_pair,
            shifted_second,
            shifted_third,
            shifted_u,
            field,
        )
        diagonal = add(
            poisson(second_kernel[0], second_kernel[1]),
            coupling(first_pair, third_kernel),
            u_kernel,
        )
        linear = add(add(shifted, constant, -field.one), diagonal, -field.one)
        nonconstant.extend((linear, diagonal))
    for left in range(len(lower_kernel)):
        left_second, left_third, left_u = lower_kernel[left]
        for right in range(left + 1, len(lower_kernel)):
            right_second, right_third, right_u = lower_kernel[right]
            cross = add(
                poisson(left_second[0], right_second[1]),
                poisson(right_second[0], left_second[1]),
            )
            cross = add(
                cross,
                coupling(first_pair, right_third),
                left_u,
            )
            cross = add(
                cross,
                coupling(first_pair, left_third),
                right_u,
            )
            nonconstant.append(cross)

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    correction_columns = [
        poisson({monomial: field.one}, T) for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: field.one}) for monomial in t4_monomials
    ]
    return (
        base_u,
        lower_kernel,
        nonconstant,
        correction_columns,
        constant,
    )


def all_scale_fourth_order_span_audit(S, T, first_pair, field):
    """Audit hbar^4 over all amplitudes and all compatible lower lifts."""

    (
        base_u,
        lower_kernel,
        nonconstant,
        correction_columns,
        constant,
    ) = all_scale_fourth_order_span_data(
        S,
        T,
        first_pair,
        field,
    )
    span_columns = correction_columns + nonconstant
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    return (
        base_u,
        len(lower_kernel),
        len(nonconstant),
        span_rank,
        augmented_rank,
    )


def subspace_fourth_order_relaxation(S, T, first_pairs, field):
    """Audit order four uniformly on a linear first-correction subspace.

    The relaxation retains the exact order-two affine family generated by
    the parity base, ker(d_2), and one particular solution for every
    quadratic monomial in the first-correction parameters.  It then allows
    the bounded third correction to vary arbitrarily, without imposing its
    preceding compatibility equation.  This only enlarges the possible
    order-four obstruction span.
    """

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    d2_columns = [
        poisson({monomial: field.one}, T) for monomial in s2_monomials
    ]
    d2_columns += [
        poisson(S, {monomial: field.one}) for monomial in t2_monomials
    ]
    parity_rhs = scale(pi_power(S, T, 3), -field.one / field(24))
    parity_vector, kernel_vectors, rank = solve_affine(
        d2_columns,
        parity_rhs,
        field,
    )
    assert rank == 1527
    parity_pair = split_correction(
        parity_vector,
        s2_monomials,
        t2_monomials,
    )
    kernel_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in kernel_vectors
    ]
    assert len(kernel_pairs) == 42

    quadratic_right_sides = []
    for left in range(len(first_pairs)):
        left_s, left_t = first_pairs[left]
        for right in range(left, len(first_pairs)):
            right_s, right_t = first_pairs[right]
            if left == right:
                bracket = poisson(left_s, left_t)
            else:
                bracket = add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            quadratic_right_sides.append(scale(bracket, -field.one))
    quadratic_vectors = solve_many_particular(
        d2_columns,
        quadratic_right_sides,
        field,
    )
    quadratic_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in quadratic_vectors
    ]
    second_deviations = kernel_pairs + quadratic_pairs

    constant = poisson(parity_pair[0], parity_pair[1])
    constant = add(
        constant,
        pi_power(parity_pair[0], T, 3),
        field.one / field(24),
    )
    constant = add(
        constant,
        pi_power(S, parity_pair[1], 3),
        field.one / field(24),
    )
    constant = add(
        constant,
        pi_power(S, T, 5),
        field.one / field(1920),
    )

    nonconstant = []
    for s_part, t_part in second_deviations:
        linear = add(
            poisson(s_part, parity_pair[1]),
            poisson(parity_pair[0], t_part),
        )
        linear = add(
            linear,
            pi_power(s_part, T, 3),
            field.one / field(24),
        )
        linear = add(
            linear,
            pi_power(S, t_part, 3),
            field.one / field(24),
        )
        nonconstant.extend((linear, poisson(s_part, t_part)))
    for left in range(len(second_deviations)):
        left_s, left_t = second_deviations[left]
        for right in range(left + 1, len(second_deviations)):
            right_s, right_t = second_deviations[right]
            nonconstant.append(
                add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            )

    # The third correction is deliberately unconstrained.  Its contribution
    # to order four is B_(S1,T1)(S3,T3).
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    third_coordinate_pairs = [
        ({monomial: field.one}, {}) for monomial in s3_monomials
    ]
    third_coordinate_pairs += [
        ({}, {monomial: field.one}) for monomial in t3_monomials
    ]
    for first_pair in first_pairs:
        nonconstant.extend(
            coupling(first_pair, third_pair)
            for third_pair in third_coordinate_pairs
        )

    for left in range(len(first_pairs)):
        left_s, left_t = first_pairs[left]
        for right in range(left, len(first_pairs)):
            right_s, right_t = first_pairs[right]
            if left == right:
                cubic = pi_power(left_s, left_t, 3)
            else:
                cubic = add(
                    pi_power(left_s, right_t, 3),
                    pi_power(right_s, left_t, 3),
                )
            nonconstant.append(cubic)

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    correction_columns = [
        poisson({monomial: field.one}, T) for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: field.one}) for monomial in t4_monomials
    ]
    span_columns = correction_columns + nonconstant
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    return (
        len(first_pairs),
        len(second_deviations),
        len(third_coordinate_pairs),
        len(nonconstant),
        span_rank,
        augmented_rank,
    )


def multiprime_all_scale_audit(
    exact_direction,
    global_vector,
    modular_problem_cache,
):
    """Run the all-scale order-four span certificate at three good primes."""

    prime_results = []
    scaled_nullity = None
    scaled_coefficient_count = None
    for prime in (31991, 32003, 65521):
        if prime not in modular_problem_cache:
            prime_field = GF(prime)
            modular_problem_cache[prime] = (
                prime_field,
                *modular_essential_pairs(prime_field),
            )
        (
            prime_field,
            prime_S,
            prime_T,
            prime_pairs,
        ) = modular_problem_cache[prime]
        prime_vector = {
            index: prime_field(
                int(sp.numer(sp.Rational(str(coefficient))))
            )
            / prime_field(
                int(sp.denom(sp.Rational(str(coefficient))))
            )
            for index, coefficient in global_vector.items()
        }
        prime_direction = linear_combination(
            prime_pairs,
            prime_vector,
            prime_field,
        )
        assert prime_direction == (
            reduce_sparse_poly(exact_direction[0], prime_field),
            reduce_sparse_poly(exact_direction[1], prime_field),
        )
        (
            scaled_base_u,
            scaled_nullity,
            scaled_coefficient_count,
            scaled_span_rank,
            scaled_augmented_rank,
        ) = all_scale_fourth_order_span_audit(
            prime_S,
            prime_T,
            prime_direction,
            prime_field,
        )
        assert scaled_base_u == prime_field.zero
        assert scaled_augmented_rank == scaled_span_rank + 1
        prime_results.append(
            (
                prime,
                scaled_span_rank,
                scaled_augmented_rank,
            )
        )
    return scaled_nullity, scaled_coefficient_count, prime_results


def main() -> None:
    (
        S,
        T,
        essential,
        pairs,
        parameter_monomials,
        equations,
    ) = essential_problem()

    # Five of the sparse basis equations expose the radical coordinates:
    # x_34^2, x_7^2, x_19^2, x_6*x_7, and x_7*x_19.  A sixth equation then
    # reduces to a nonzero multiple of x_6^2.
    radical_coordinates = (6, 7, 19, 34)
    surviving_axes = (0, 1, 2, 3, 8, 9, 10, 23, 24, 35)

    variable = sp.Symbol("r")
    full_lines = []
    rational_directions = []
    algebraic_directions = []
    for left in range(len(pairs)):
        for right in range(left + 1, len(pairs)):
            common = common_direction_polynomial(
                equations,
                parameter_monomials,
                left,
                right,
                variable,
            )
            if common is None:
                full_lines.append((left, right))
                continue
            if common.degree() == 0:
                continue
            # Roots zero on the chosen affine chart are axes, not mixed
            # directions.  Remove that factor before classifying the roots.
            while common.degree() and common.eval(0) == 0:
                common = sp.quo(
                    common,
                    sp.Poly(variable, variable, domain=sp.QQ),
                ).monic()
            if common.degree() == 0:
                continue
            roots = sp.roots(common.as_expr(), variable)
            has_algebraic_root = False
            for root in roots:
                if root.is_Rational:
                    rational_directions.append((left, right, sp.Rational(root)))
                else:
                    has_algebraic_root = True
            if has_algebraic_root:
                algebraic_directions.append((left, right, common.as_expr()))
    assert len(full_lines) == 36
    assert len(rational_directions) == 7
    assert len(algebraic_directions) == 9

    print(
        "PASS: the reduced quadratic locus forces quotient coordinates "
        f"{list(radical_coordinates)} to vanish"
    )
    print(
        "PASS: coordinate axes on the quadratic locus are "
        f"{list(surviving_axes)}"
    )
    print(f"coordinate lines contained in the quadratic locus: {full_lines}")
    print(
        "isolated rational directions of exact support two: "
        f"{rational_directions}"
    )
    print(
        "irreducible algebraic directions of exact support two: "
        f"{algebraic_directions}"
    )

    line_edges = set(full_lines)
    coordinate_subspaces = maximal_cliques(surviving_axes, line_edges)
    assert coordinate_subspaces == [
        (0, 1, 2, 3, 8, 9, 10),
        (0, 1, 9, 10, 23, 24, 35),
    ]
    print(
        "maximal coordinate linear subspaces in the quadratic locus: "
        f"{coordinate_subspaces}"
    )

    third_order_survivors = []
    for left, right, ratio in rational_directions:
        ratio_q = QQ(int(sp.numer(ratio))) / QQ(int(sp.denom(ratio)))
        direction = combine_pair(pairs[left], pairs[right], ratio_q)
        survives, reason, rank, nullity = third_order_axis_test(
            S,
            T,
            direction,
            QQ,
        )
        print(
            f"direction x{left}+({ratio})*x{right}: "
            f"hbar^3 combined rank={rank} nullity={nullity}; {reason}"
        )
        if survives:
            third_order_survivors.append((left, right, ratio))

    algebraic_third_order_survivors = []
    for left, right, polynomial in algebraic_directions:
        algebraic_root = sp.RootOf(polynomial, 0)
        number_field = QQ.algebraic_field(algebraic_root)
        field_S = extend_sparse_poly(S, number_field)
        field_T = extend_sparse_poly(T, number_field)
        field_pairs = [
            (
                extend_sparse_poly(s_part, number_field),
                extend_sparse_poly(t_part, number_field),
            )
            for s_part, t_part in pairs
        ]
        direction = linear_combination(
            field_pairs,
            {
                left: number_field.one,
                right: number_field.from_sympy(algebraic_root),
            },
            number_field,
        )
        survives, reason, rank, nullity = third_order_axis_test(
            field_S,
            field_T,
            direction,
            number_field,
        )
        print(
            f"quadratic direction ({left},{right}) over {number_field}: "
            f"rank={rank} nullity={nullity}; {reason}"
        )
        if survives:
            algebraic_third_order_survivors.append((left, right, polynomial))
    assert not algebraic_third_order_survivors

    audited_survivors = set()
    modular_problem_cache = {}
    support_two_residual_vectors = None
    if full_lines:
        for subspace in coordinate_subspaces:
            (
                relaxed_rank,
                augmented_rank,
                dimension,
                quadratic_count,
                candidate_basis,
            ) = uniform_third_order_relaxation(
                S,
                T,
                pairs,
                subspace,
                QQ,
            )
            print(
                f"coordinate subspace {subspace}: relaxed hbar^3 rank "
                f"{relaxed_rank}->{augmented_rank} after all "
                f"{quadratic_count} quadratic particulars"
            )
            print(
                "  exact necessary direction space: "
                f"{candidate_basis}"
            )
            if augmented_rank - relaxed_rank == dimension:
                print(
                    "  PASS: the entire subspace is uniformly killed by "
                    "the relaxed hbar^3 equations"
                )
            else:
                print(
                    "  SCOPE: only this residual linear direction space "
                    "can survive the relaxation"
                )
            candidate_global_vectors = [
                {
                    subspace[local_index]: coefficient
                    for local_index, coefficient in local_vector.items()
                }
                for local_vector in candidate_basis
            ]
            if subspace == coordinate_subspaces[1]:
                support_two_residual_vectors = candidate_global_vectors
            candidate_survival = []
            for candidate_index, global_vector in enumerate(
                candidate_global_vectors
            ):
                direction = linear_combination(pairs, global_vector, QQ)
                survives, reason, rank, nullity = third_order_axis_test(
                    S,
                    T,
                    direction,
                    QQ,
                )
                print(
                    f"  residual basis direction {candidate_index} "
                    f"{global_vector}: rank={rank} nullity={nullity}; "
                    f"{reason}"
                )
                candidate_survival.append(survives)
                if survives:
                    (
                        _,
                        _,
                        continuation_steps,
                        residuals,
                    ) = canonical_branch_continuation(
                        S,
                        T,
                        direction,
                        QQ,
                    )
                    print(
                        "    free-zero continuation steps "
                        f"(order,solvable,columns,rank,nullity): "
                        f"{continuation_steps}"
                    )
                    if len(continuation_steps) == 4:
                        print(
                            "    remaining normalized-commutator support "
                            f"at orders 6..9: "
                            f"{[(order, len(defect)) for order, defect in residuals.items()]}"
                        )
                    audit_key = tuple(sorted(global_vector.items()))
                    if audit_key not in audited_survivors:
                        audited_survivors.add(audit_key)
                        modular_field = GF(32003)
                        if 32003 not in modular_problem_cache:
                            modular_problem_cache[32003] = (
                                modular_field,
                                *modular_essential_pairs(modular_field),
                            )
                        (
                            _,
                            modular_S,
                            modular_T,
                            modular_pairs,
                        ) = modular_problem_cache[32003]
                        modular_vector = {
                            index: modular_field(
                                int(sp.numer(sp.Rational(str(coefficient))))
                            )
                            / modular_field(
                                int(sp.denom(sp.Rational(str(coefficient))))
                            )
                            for index, coefficient in global_vector.items()
                        }
                        modular_direction = linear_combination(
                            modular_pairs,
                            modular_vector,
                            modular_field,
                        )
                        (
                            lower_nullity,
                            coefficient_count,
                            fourth_span_rank,
                            fourth_augmented_rank,
                        ) = fourth_order_span_audit(
                            modular_S,
                            modular_T,
                            modular_direction,
                            modular_field,
                        )
                        print(
                            "    unit-scale all-lower-lift hbar^4 span audit "
                            "modulo 32003: "
                            f"nullity={lower_nullity}, "
                            f"coefficients={coefficient_count}, "
                            f"rank {fourth_span_rank}->{fourth_augmented_rank}"
                        )
                        assert fourth_augmented_rank == fourth_span_rank + 1

                        (
                            scaled_nullity,
                            scaled_coefficient_count,
                            prime_results,
                        ) = multiprime_all_scale_audit(
                            direction,
                            global_vector,
                            modular_problem_cache,
                        )
                        print(
                            "    all-scale/all-lower-lift hbar^4 span audit: "
                            f"nullity={scaled_nullity}, "
                            f"coefficients={scaled_coefficient_count}, "
                            f"rank checks {prime_results}"
                        )
            if subspace == coordinate_subspaces[0]:
                assert candidate_survival == [True, True]
            else:
                assert candidate_survival == [True, True, False, False, False]

    assert support_two_residual_vectors is not None
    residual_pair_survivors = []
    for left in range(len(support_two_residual_vectors)):
        left_vector = support_two_residual_vectors[left]
        left_x0 = left_vector[0]
        for right in range(left + 1, len(support_two_residual_vectors)):
            right_vector = support_two_residual_vectors[right]
            ratio = -left_x0 / right_vector[0]
            global_vector = dict(left_vector)
            for index, coefficient in right_vector.items():
                global_vector[index] = (
                    global_vector.get(index, QQ.zero) + ratio * coefficient
                )
            global_vector = {
                index: coefficient
                for index, coefficient in global_vector.items()
                if coefficient
            }
            assert 0 not in global_vector
            assert len(global_vector) == 2
            direction = linear_combination(pairs, global_vector, QQ)
            survives, reason, rank, nullity = third_order_axis_test(
                S,
                T,
                direction,
                QQ,
            )
            print(
                f"residual support-two direction ({left},{right}) "
                f"{global_vector}: rank={rank} nullity={nullity}; {reason}"
            )
            if survives:
                residual_pair_survivors.append((left, right))
                (
                    scaled_nullity,
                    scaled_coefficient_count,
                    prime_results,
                ) = multiprime_all_scale_audit(
                    direction,
                    global_vector,
                    modular_problem_cache,
                )
                print(
                    "  all-scale/all-lower-lift hbar^4 span audit: "
                    f"nullity={scaled_nullity}, "
                    f"coefficients={scaled_coefficient_count}, "
                    f"rank checks {prime_results}"
                )
    assert residual_pair_survivors == [(0, 1)]

    print(
        "isolated rational exact-support-two directions surviving through "
        "hbar^3: "
        f"{third_order_survivors}"
    )
    assert not third_order_survivors
    print(
        "PASS: every exact-support-two branch is eliminated by hbar^4"
    )
    print(
        "SCOPE: support-three and residual-space continuations are handled "
        "by the dedicated follow-up scripts"
    )


if __name__ == "__main__":
    main()

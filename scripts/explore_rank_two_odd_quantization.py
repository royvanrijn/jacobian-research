#!/usr/bin/env python3
"""Odd-hbar continuation of the localized rank-two Weyl-symbol lift.

This calculation uses the same degree-five sample and filtered spaces as
``explore_degree_five_a2_subprincipal.py``.  It computes the 57-dimensional
first-order kernel, removes the 19-dimensional target-Hamiltonian gauge
subspace, and projects the second-order Maurer--Cartan obstruction to the
cokernel of the next linearized commutator operator.

The kernel, gauge quotient, quadratic obstruction, and coordinate-axis
continuation are computed exactly over Q.  The large quadratic projection is
also repeated over the good prime 32003 as an independent rank check.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_a2_subprincipal import (
    add,
    degree_five_sample,
    filtered_monomials,
    multiply,
    pi_power,
    poisson,
    scale,
)


PRIME = 32003


def convert_poly(poly, field):
    result = {}
    for monomial, value in poly.items():
        value = Fraction(value)
        result[monomial] = field(value.numerator) / field(value.denominator)
    return result


def operator_rref(S, T, s_monomials, t_monomials, field):
    columns = [poisson({monomial: field.one}, T) for monomial in s_monomials]
    columns += [poisson(S, {monomial: field.one}) for monomial in t_monomials]
    output_monomials = sorted(set().union(*(set(column) for column in columns)))
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    rows = {}
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows.setdefault(output_index[monomial], {})[column_index] = coefficient
    reduced, pivots, nonzero = sdm_irref(rows)
    return columns, output_monomials, reduced, pivots, nonzero


def first_kernel(S, T, field):
    s_monomials = filtered_monomials(27, 4)
    t_monomials = filtered_monomials(23, 3)
    columns, outputs, reduced, pivots, nonzero = operator_rref(
        S, T, s_monomials, t_monomials, field
    )
    kernel, _ = sdm_nullspace_from_rref(
        reduced, field.one, len(columns), pivots, nonzero
    )
    assert len(columns) == 2132
    assert len(pivots) == 2075
    assert len(kernel) == 57
    free_columns = [column for column in range(len(columns)) if column not in pivots]
    assert all(vector[free] == field.one for vector, free in zip(kernel, free_columns))
    return s_monomials, t_monomials, free_columns, kernel, len(outputs)


def poly_power(poly, exponent, field):
    result = {(0, 0, 0): field.one}
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def gauge_pivots(S, T, s_monomials, t_monomials, free_columns, kernel, field):
    """Return pivots for the 19 target-Hamiltonian gauge directions.

    The Hamiltonians are R^a*T (0<=a<=9), T^2/2, and R^a*S
    (0<=a<=7).  Their induced first corrections are respectively an
    R-dependent translation of S, the shear S -> S+hbar*T, and an
    R-dependent translation of T.  Signs do not affect their span.
    """

    R = {(1, 0, 0): field(2), (2, 1, 0): field(-3)}
    s_index = {monomial: index for index, monomial in enumerate(s_monomials)}
    t_index = {
        monomial: len(s_monomials) + index
        for index, monomial in enumerate(t_monomials)
    }
    ambient_gauge = []
    for exponent in range(10):
        ambient_gauge.append(
            {s_index[monomial]: value for monomial, value in poly_power(R, exponent, field).items()}
        )
    ambient_gauge.append({s_index[monomial]: value for monomial, value in T.items()})
    for exponent in range(8):
        ambient_gauge.append(
            {t_index[monomial]: value for monomial, value in poly_power(R, exponent, field).items()}
        )

    gauge_coordinates = []
    for gauge in ambient_gauge:
        coordinates = {
            index: gauge[free]
            for index, free in enumerate(free_columns)
            if gauge.get(free)
        }
        reconstruction = defaultdict(lambda: field.zero)
        for index, coefficient in coordinates.items():
            for column, value in kernel[index].items():
                reconstruction[column] += coefficient * value
        assert {column: value for column, value in reconstruction.items() if value} == gauge
        gauge_coordinates.append(coordinates)

    _, pivots, _ = sdm_irref(
        {index: row for index, row in enumerate(gauge_coordinates)}
    )
    assert len(pivots) == 19
    return pivots


def split_pair(vector, s_monomials, t_monomials):
    split = len(s_monomials)
    s_part = {
        s_monomials[column]: value
        for column, value in vector.items()
        if column < split and value
    }
    t_part = {
        t_monomials[column - split]: value
        for column, value in vector.items()
        if column >= split and value
    }
    return s_part, t_part


def second_obstruction(S, T, quotient_pairs, field):
    """Project {S1,T1} to coker(d_2) as homogeneous quadrics."""

    quadratic_vectors = []
    parameter_monomials = []
    for index, (s_part, t_part) in enumerate(quotient_pairs):
        quadratic_vectors.append(poisson(s_part, t_part))
        parameter_monomials.append((index, index))
    for left in range(len(quotient_pairs)):
        for right in range(left + 1, len(quotient_pairs)):
            s_left, t_left = quotient_pairs[left]
            s_right, t_right = quotient_pairs[right]
            quadratic_vectors.append(
                add(poisson(s_left, t_right), poisson(s_right, t_left))
            )
            parameter_monomials.append((left, right))

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    columns = [poisson({monomial: field.one}, T) for monomial in s2_monomials]
    columns += [poisson(S, {monomial: field.one}) for monomial in t2_monomials]
    output_monomials = sorted(
        set().union(
            *(set(column) for column in columns),
            *(set(vector) for vector in quadratic_vectors),
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose_rows)
    left_kernel, _ = sdm_nullspace_from_rref(
        reduced, field.one, len(output_monomials), pivots, nonzero
    )
    incidence = defaultdict(list)
    for functional_index, functional in enumerate(left_kernel):
        for output_row, coefficient in functional.items():
            incidence[output_row].append((functional_index, coefficient))

    projected = [defaultdict(lambda: field.zero) for _ in left_kernel]
    for coefficient_index, vector in enumerate(quadratic_vectors):
        for monomial, value in vector.items():
            for functional_index, functional_value in incidence[output_index[monomial]]:
                projected[functional_index][coefficient_index] += functional_value * value
    projected_rows = {
        row: {column: value for column, value in equation.items() if value}
        for row, equation in enumerate(projected)
        if any(equation.values())
    }
    obstruction_rref, obstruction_pivots, _ = sdm_irref(projected_rows)
    # A basis selected from the original projected equations is dramatically
    # sparser than the reduced basis and generates the same quadratic ideal.
    ordered_rows = sorted(projected_rows.values(), key=len)
    transposed_projected = {}
    for equation_index, equation in enumerate(ordered_rows):
        for monomial_index, coefficient in equation.items():
            transposed_projected.setdefault(monomial_index, {})[
                equation_index
            ] = coefficient
    _, sparse_pivots, _ = sdm_irref(transposed_projected)
    sparse_basis = {
        index: ordered_rows[pivot] for index, pivot in enumerate(sparse_pivots)
    }
    assert len(sparse_basis) == len(obstruction_pivots)
    return (
        quadratic_vectors,
        parameter_monomials,
        obstruction_rref,
        obstruction_pivots,
        len(pivots),
        len(left_kernel),
        len(projected_rows),
        sparse_basis,
    )


def evaluate_equations(equations, parameter_monomials, values, field):
    for equation in equations.values():
        total = field.zero
        for column, coefficient in equation.items():
            left, right = parameter_monomials[column]
            total += coefficient * values.get(left, field.zero) * values.get(
                right, field.zero
            )
        if total:
            return False
    return True


def third_order_axis_test(S, T, first_pair, field):
    """Test a whole nonzero coordinate axis through the hbar^3 equation.

    If the first correction is lambda*(s1,t1), put u=lambda^2 and divide
    the odd third-order equation by lambda.  The hbar^2 and hbar^3 equations
    are then one affine *linear* system in (S2,T2), (S3/lambda,T3/lambda),
    and u.  We report whether the system permits u != 0.
    """

    s1, t1 = first_pair
    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    split2 = len(s2_monomials)

    columns2 = [poisson({monomial: field.one}, T) for monomial in s2_monomials]
    columns2 += [poisson(S, {monomial: field.one}) for monomial in t2_monomials]
    coupling = []
    for column, monomial in enumerate(s2_monomials + t2_monomials):
        if column < split2:
            coupling.append(poisson({monomial: field.one}, t1))
        else:
            coupling.append(poisson(s1, {monomial: field.one}))
    columns3 = [poisson({monomial: field.one}, T) for monomial in s3_monomials]
    columns3 += [poisson(S, {monomial: field.one}) for monomial in t3_monomials]

    base_rhs2 = scale(pi_power(S, T, 3), -field.one / field(24))
    quadratic = poisson(s1, t1)
    rhs3 = scale(
        add(pi_power(s1, T, 3), pi_power(S, t1, 3)),
        -field.one / field(24),
    )

    tagged_rows = set()
    for column in columns2:
        tagged_rows.update((2, monomial) for monomial in column)
    for column in coupling:
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
    for column_index, (linear2, linear3) in enumerate(zip(columns2, coupling)):
        for monomial, coefficient in linear2.items():
            rows.setdefault(row_index[(2, monomial)], {})[column_index] = coefficient
        for monomial, coefficient in linear3.items():
            rows.setdefault(row_index[(3, monomial)], {})[column_index] = coefficient
    for offset, column in enumerate(columns3):
        column_index = count2 + offset
        for monomial, coefficient in column.items():
            rows.setdefault(row_index[(3, monomial)], {})[column_index] = coefficient
    for monomial, coefficient in quadratic.items():
        rows.setdefault(row_index[(2, monomial)], {})[u_column] = coefficient
    for monomial, coefficient in base_rhs2.items():
        rows.setdefault(row_index[(2, monomial)], {})[rhs_column] = -coefficient
    for monomial, coefficient in rhs3.items():
        rows.setdefault(row_index[(3, monomial)], {})[rhs_column] = -coefficient

    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        return False, "inconsistent even at u=0", len(pivots), 0
    nullity = rhs_column - len(pivots)
    if u_column not in pivots:
        return True, "u is free", len(pivots), nullity
    pivot_row = pivots.index(u_column)
    equation = reduced[pivot_row]
    other_unknowns = {
        column: value
        for column, value in equation.items()
        if column not in (u_column, rhs_column) and value
    }
    if other_unknowns:
        return True, "u depends on free correction parameters", len(pivots), nullity
    forced_u = -equation.get(rhs_column, field.zero)
    if forced_u:
        return True, f"u is forced to {int(forced_u)}", len(pivots), nullity
    return False, "u is forced to zero", len(pivots), nullity


def main():
    # Compute the complete kernel, gauge quotient, and next obstruction over Q.
    S_fraction, T_fraction = degree_five_sample()
    S_q = convert_poly(S_fraction, QQ)
    T_q = convert_poly(T_fraction, QQ)
    s_q_monomials, t_q_monomials, free_q, kernel_q, output_count_q = first_kernel(
        S_q, T_q, QQ
    )
    assert len(kernel_q) == 57
    exact_gauge_pivots = gauge_pivots(
        S_q,
        T_q,
        s_q_monomials,
        t_q_monomials,
        free_q,
        kernel_q,
        QQ,
    )
    exact_essential = [
        index for index in range(57) if index not in exact_gauge_pivots
    ]
    exact_pairs = [
        split_pair(kernel_q[index], s_q_monomials, t_q_monomials)
        for index in exact_essential
    ]
    (
        _,
        exact_parameter_monomials,
        exact_equations,
        exact_equation_pivots,
        _,
        _,
        _,
        _,
    ) = second_obstruction(S_q, T_q, exact_pairs, QQ)
    exact_surviving_axes = [
        index
        for index in range(len(exact_essential))
        if evaluate_equations(
            exact_equations,
            exact_parameter_monomials,
            {index: QQ.one},
            QQ,
        )
    ]
    print(
        "PASS: exact first operator has 2132 columns, rank 2075, "
        f"nullity 57, and {output_count_q} occupied output monomials"
    )
    print(
        "PASS: exact target-Hamiltonian gauge rank=19 and exact quadratic "
        f"obstruction rank={len(exact_equation_pivots)}"
    )

    field = GF(PRIME)
    S = convert_poly(S_fraction, field)
    T = convert_poly(T_fraction, field)
    s_monomials, t_monomials, free_columns, kernel, _ = first_kernel(
        S, T, field
    )
    pivots = gauge_pivots(
        S, T, s_monomials, t_monomials, free_columns, kernel, field
    )
    assert pivots == exact_gauge_pivots
    essential = [index for index in range(57) if index not in pivots]
    assert len(essential) == 38
    quotient_pairs = [
        split_pair(kernel[index], s_monomials, t_monomials) for index in essential
    ]
    print(f"PASS: target-Hamiltonian gauge rank=19; quotient dimension={len(essential)}")
    print(f"gauge pivots={pivots}")
    print(f"quotient kernel indices={essential}")

    (
        _,
        parameter_monomials,
        equations,
        equation_pivots,
        second_rank,
        second_cokernel,
        raw_equations,
        sparse_equations,
    ) = second_obstruction(S, T, quotient_pairs, field)
    print(
        f"PASS: next linear operator rank={second_rank}; "
        f"ambient left-cokernel dimension={second_cokernel}"
    )
    print(
        f"PASS: quadratic obstruction map has {raw_equations} nonzero "
        f"cokernel coordinates and rank {len(equation_pivots)} modulo {PRIME}"
    )
    term_counts = sorted(len(equation) for equation in equations.values())
    print(f"reduced quadratic term-count range={term_counts[0]}..{term_counts[-1]}")
    sparse_counts = sorted(len(equation) for equation in sparse_equations.values())
    print(f"sparse quadratic term-count range={sparse_counts[0]}..{sparse_counts[-1]}")

    surviving_axes = []
    for index in range(len(essential)):
        if evaluate_equations(equations, parameter_monomials, {index: field.one}, field):
            surviving_axes.append(index)
    assert surviving_axes == exact_surviving_axes
    print(
        f"coordinate-axis branches surviving the next obstruction: "
        f"{surviving_axes} ({len(surviving_axes)} of {len(essential)})"
    )
    third_order_survivors = []
    for axis in surviving_axes:
        survives, reason, rank, nullity = third_order_axis_test(
            S_q, T_q, exact_pairs[axis], QQ
        )
        print(
            f"axis {axis}: hbar^3 combined rank={rank} nullity={nullity}; "
            f"{reason}"
        )
        if survives:
            third_order_survivors.append(axis)
    print(f"coordinate-axis branches surviving through hbar^3: {third_order_survivors}")


if __name__ == "__main__":
    main()

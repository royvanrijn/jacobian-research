#!/usr/bin/env python3
"""Certify the first restricted quantization stages of the (4,3) symbols.

This is a bounded, symbol-specific calculation.  It does not construct an
endomorphism of A_2 and it does not obstruct quantizations outside the stated
Bernstein filtration or the normal ordering encoded by ``pi_power``.

For the specialization (a,tau)=(-4/3,0), the classical symbols have

    (deg_B(S), ord_Z(S)) = (22,4),
    (deg_B(T), ord_Z(T)) = (18,3).

At hbar^n the inherited filtered ansatz lowers Bernstein degree by 2*n and
Z-order by n.  The script:

* solves the complete parity-preserving hbar^3 affine equation;
* constructs an exact dual cocycle which annihilates every allowed hbar^5
  correction and every constant/linear/quadratic variation over that affine
  family, but evaluates to one on the hbar^5 defect;
* rebuilds the unrestricted hbar^1 kernel and its complete admissible
  target-Hamiltonian gauge subspace; and
* projects the hbar^2 Maurer--Cartan quadrics to the next cokernel.

The exact calculation is independently repeated over GF(32003).  A JSON
certificate containing the rational dual cocycle can optionally be written
under artifacts/generated-results/.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_a2_subprincipal import (
    add,
    filtered_monomials,
    pi_power,
    poisson,
    scale,
)
from explore_degree_five_quantum_residue import (
    column_rank,
    degree_five_family,
    poly_power,
    solve_affine,
    split_pair,
)
from explore_rank_two_odd_quantization import operator_rref
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    moyal_coefficient,
    solve_many_particular,
    split_correction,
)


PRIME = 32003
SPECIAL_A = Fraction(-4, 3)
SPECIAL_TAU = Fraction(0)


@dataclass(frozen=True)
class SymbolBounds:
    s_degree: int
    s_order: int
    t_degree: int
    t_order: int

    def correction(self, order: int) -> tuple[list[tuple[int, int, int]], ...]:
        return (
            filtered_monomials(
                self.s_degree - 2 * order,
                self.s_order - order,
            )
            if self.s_order >= order
            else [],
            filtered_monomials(
                self.t_degree - 2 * order,
                self.t_order - order,
            )
            if self.t_order >= order
            else [],
        )


BOUNDS = SymbolBounds(22, 4, 18, 3)


def field_fraction(field, value: Fraction):
    return field(value.numerator) / field(value.denominator)


def pairing(functional, vector, field):
    return sum(
        (coefficient * vector.get(monomial, field.zero)
         for monomial, coefficient in functional.items()),
        field.zero,
    )


def dual_witness(span_columns, constant, field):
    """Return a functional killing ``span_columns`` and pairing 1 with constant."""

    output_monomials = sorted(
        set(constant).union(*(set(column) for column in span_columns))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(span_columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose_rows)
    left_kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(output_monomials),
        pivots,
        nonzero,
    )
    for vector in left_kernel:
        functional = {
            output_monomials[index]: coefficient
            for index, coefficient in vector.items()
            if coefficient
        }
        value = pairing(functional, constant, field)
        if value:
            functional = {
                monomial: coefficient / value
                for monomial, coefficient in functional.items()
            }
            assert all(
                pairing(functional, column, field) == field.zero
                for column in span_columns
            )
            assert pairing(functional, constant, field) == field.one
            return functional
    raise AssertionError("rank jump has no dual witness")


def parity_audit(S, T, field):
    s2_monomials, t2_monomials = BOUNDS.correction(2)
    columns = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    columns += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    rhs = scale(pi_power(S, T, 3), -field.one / field(24))
    particular, kernel, rank = solve_affine(columns, rhs, field)
    base = split_pair(particular, s2_monomials, t2_monomials)
    kernel_pairs = [
        split_pair(vector, s2_monomials, t2_monomials)
        for vector in kernel
    ]

    def fifth_defect(pair):
        s2, t2 = pair
        value = poisson(s2, t2)
        value = add(value, pi_power(s2, T, 3), field.one / field(24))
        value = add(value, pi_power(S, t2, 3), field.one / field(24))
        return add(
            value,
            pi_power(S, T, 5),
            field.one / field(1920),
        )

    constant = fifth_defect(base)
    nonconstant = []
    for basis_s, basis_t in kernel_pairs:
        diagonal = poisson(basis_s, basis_t)
        shifted = fifth_defect(
            (add(base[0], basis_s), add(base[1], basis_t))
        )
        linear = add(add(shifted, constant, -1), diagonal, -1)
        nonconstant.extend((linear, diagonal))
    for left, right in combinations(range(len(kernel_pairs)), 2):
        left_s, left_t = kernel_pairs[left]
        right_s, right_t = kernel_pairs[right]
        nonconstant.append(
            add(
                poisson(left_s, right_t),
                poisson(right_s, left_t),
            )
        )

    s4_monomials, t4_monomials = BOUNDS.correction(4)
    correction_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t4_monomials
    ]
    span_columns = correction_columns + nonconstant
    correction_rank = column_rank(correction_columns)
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    if augmented_rank != span_rank + 1:
        raise AssertionError("the expected exact hbar^5 rank jump failed")
    witness = dual_witness(span_columns, constant, field)
    return {
        "h3_s_columns": len(s2_monomials),
        "h3_t_columns": len(t2_monomials),
        "h3_rank": rank,
        "h3_nullity": len(kernel),
        "h5_s_columns": len(s4_monomials),
        "h5_t_columns": len(t4_monomials),
        "h5_correction_rank": correction_rank,
        "h5_parameter_coefficients": len(nonconstant),
        "h5_span_rank": span_rank,
        "h5_augmented_rank": augmented_rank,
        "witness": witness,
    }


def target_gauge_pivots(
    S,
    T,
    s_monomials,
    t_monomials,
    free_columns,
    kernel,
    field,
):
    """Compute all target-Hamiltonian directions fitting the hbar^1 bounds."""

    R = {(1, 0, 0): field(2), (2, 1, 0): field(-3)}
    s_index = {monomial: index for index, monomial in enumerate(s_monomials)}
    t_index = {
        monomial: len(s_monomials) + index
        for index, monomial in enumerate(t_monomials)
    }
    ambient_gauge = []
    s1_degree = BOUNDS.s_degree - 2
    t1_degree = BOUNDS.t_degree - 2
    for exponent in range(s1_degree // 3 + 1):
        power = poly_power(R, exponent, field.one)
        if all(monomial in s_index for monomial in power):
            ambient_gauge.append(
                {s_index[monomial]: value for monomial, value in power.items()}
            )
    if all(monomial in s_index for monomial in T):
        ambient_gauge.append(
            {s_index[monomial]: value for monomial, value in T.items()}
        )
    for exponent in range(t1_degree // 3 + 1):
        power = poly_power(R, exponent, field.one)
        if all(monomial in t_index for monomial in power):
            ambient_gauge.append(
                {t_index[monomial]: value for monomial, value in power.items()}
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
        if {
            column: value
            for column, value in reconstruction.items()
            if value
        } != gauge:
            raise AssertionError("target gauge did not reconstruct in kernel")
        gauge_coordinates.append(coordinates)
    _, pivots, _ = sdm_irref(
        {index: row for index, row in enumerate(gauge_coordinates)}
    )
    return pivots, len(ambient_gauge)


def project_second_obstruction(S, T, quotient_pairs, field):
    quadratic_vectors = []
    parameter_monomials = []
    for index, (s_part, t_part) in enumerate(quotient_pairs):
        quadratic_vectors.append(poisson(s_part, t_part))
        parameter_monomials.append((index, index))
    for left, right in combinations(range(len(quotient_pairs)), 2):
        s_left, t_left = quotient_pairs[left]
        s_right, t_right = quotient_pairs[right]
        quadratic_vectors.append(
            add(
                poisson(s_left, t_right),
                poisson(s_right, t_left),
            )
        )
        parameter_monomials.append((left, right))

    s2_monomials, t2_monomials = BOUNDS.correction(2)
    columns = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    columns += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
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
        reduced,
        field.one,
        len(output_monomials),
        pivots,
        nonzero,
    )
    incidence = defaultdict(list)
    for functional_index, functional in enumerate(left_kernel):
        for output_row, coefficient in functional.items():
            incidence[output_row].append((functional_index, coefficient))
    projected = [defaultdict(lambda: field.zero) for _ in left_kernel]
    for coefficient_index, vector in enumerate(quadratic_vectors):
        for monomial, value in vector.items():
            for functional_index, functional_value in incidence[
                output_index[monomial]
            ]:
                projected[functional_index][coefficient_index] += (
                    functional_value * value
                )
    projected_rows = {
        row: {
            column: value
            for column, value in equation.items()
            if value
        }
        for row, equation in enumerate(projected)
        if any(equation.values())
    }
    obstruction_rref, obstruction_pivots, _ = sdm_irref(projected_rows)
    surviving_axes = []
    for parameter in range(len(quotient_pairs)):
        diagonal = parameter_monomials.index((parameter, parameter))
        if all(
            equation.get(diagonal, field.zero) == field.zero
            for equation in obstruction_rref.values()
        ):
            surviving_axes.append(parameter)
    summary = {
        "next_operator_columns": len(columns),
        "next_operator_rank": len(pivots),
        "ambient_cokernel": len(left_kernel),
        "quadratic_coefficients": len(quadratic_vectors),
        "raw_nonzero_equations": len(projected_rows),
        "quadratic_rank": len(obstruction_pivots),
        "surviving_axes": surviving_axes,
    }
    return summary, parameter_monomials, obstruction_rref


def third_order_axis_test(S, T, first_pair, field):
    """Test one nonzero hbar^1 coordinate axis through hbar^3."""

    s1, t1 = first_pair
    s2_monomials, t2_monomials = BOUNDS.correction(2)
    s3_monomials, t3_monomials = BOUNDS.correction(3)
    split2 = len(s2_monomials)

    columns2 = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    columns2 += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    coupling = []
    for column, monomial in enumerate(s2_monomials + t2_monomials):
        if column < split2:
            coupling.append(poisson({monomial: field.one}, t1))
        else:
            coupling.append(poisson(s1, {monomial: field.one}))
    columns3 = [
        poisson({monomial: field.one}, T)
        for monomial in s3_monomials
    ]
    columns3 += [
        poisson(S, {monomial: field.one})
        for monomial in t3_monomials
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
    for column in coupling:
        tagged_rows.update((3, monomial) for monomial in column)
    for column in columns3:
        tagged_rows.update((3, monomial) for monomial in column)
    tagged_rows.update((2, monomial) for monomial in base_rhs2)
    tagged_rows.update((2, monomial) for monomial in quadratic)
    tagged_rows.update((3, monomial) for monomial in rhs3)
    tagged_rows = sorted(tagged_rows)
    row_index = {
        tagged: index for index, tagged in enumerate(tagged_rows)
    }

    count2 = len(columns2)
    count3 = len(columns3)
    u_column = count2 + count3
    rhs_column = u_column + 1
    rows = {}
    for column_index, (linear2, linear3) in enumerate(
        zip(columns2, coupling)
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
        column_index = count2 + offset
        for monomial, coefficient in column.items():
            rows.setdefault(row_index[(3, monomial)], {})[
                column_index
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

    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        return {
            "allows_nonzero_u": False,
            "u_status": "combined_system_inconsistent_for_all_u",
            "rank": len(pivots),
            "nullity": 0,
        }
    nullity = rhs_column - len(pivots)
    if u_column not in pivots:
        status = "free"
        allows_nonzero = True
    else:
        pivot_row = pivots.index(u_column)
        equation = reduced[pivot_row]
        other_unknowns = {
            column: value
            for column, value in equation.items()
            if column not in (u_column, rhs_column) and value
        }
        if other_unknowns:
            status = "depends_on_free_corrections"
            allows_nonzero = True
        elif equation.get(rhs_column, field.zero):
            status = "forced_nonzero"
            allows_nonzero = True
        else:
            status = "forced_zero"
            allows_nonzero = False
    return {
        "allows_nonzero_u": allows_nonzero,
        "u_status": status,
        "rank": len(pivots),
        "nullity": nullity,
    }


def sympy_rational(value):
    return sp.Rational(str(value))


def restricted_quadratic(
    equation,
    parameter_monomials,
    left,
    right,
    variable,
):
    expression = sp.S.Zero
    for column, coefficient in equation.items():
        first, second = parameter_monomials[column]
        if first not in (left, right) or second not in (left, right):
            continue
        first_value = sp.S.One if first == left else variable
        second_value = sp.S.One if second == left else variable
        expression += (
            sympy_rational(coefficient) * first_value * second_value
        )
    return sp.Poly(expression, variable, domain=sp.QQ)


def common_direction_polynomial(
    equations,
    parameter_monomials,
    left,
    right,
    variable,
):
    restrictions = [
        restricted_quadratic(
            equation,
            parameter_monomials,
            left,
            right,
            variable,
        )
        for equation in equations.values()
    ]
    restrictions = [
        polynomial
        for polynomial in restrictions
        if not polynomial.is_zero
    ]
    if not restrictions:
        return None
    common = restrictions[0].monic()
    for polynomial in restrictions[1:]:
        common = sp.gcd(common, polynomial.monic()).monic()
        if common.degree() == 0:
            break
    while common.degree() and common.eval(0) == 0:
        common = sp.quo(
            common,
            sp.Poly(variable, variable, domain=sp.QQ),
        ).monic()
    return common


def linear_combination(pairs, coefficients, field):
    s_part = {}
    t_part = {}
    for index, coefficient in coefficients.items():
        s_part = add(s_part, pairs[index][0], coefficient)
        t_part = add(t_part, pairs[index][1], coefficient)
    return s_part, t_part


def uniform_third_order_relaxation(S, T, basis_pairs, field):
    """Return the necessary linear direction space after relaxed hbar^3."""

    s2_monomials, t2_monomials = BOUNDS.correction(2)
    d2_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    d2_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    parity_rhs = scale(pi_power(S, T, 3), -field.one / field(24))
    parity_vector, kernel_vectors, _ = solve_affine(
        d2_columns,
        parity_rhs,
        field,
    )
    parity_pair = split_correction(
        parity_vector,
        s2_monomials,
        t2_monomials,
    )
    kernel_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in kernel_vectors
    ]

    quadratic_right_sides = []
    for left in range(len(basis_pairs)):
        left_s, left_t = basis_pairs[left]
        for right in range(left, len(basis_pairs)):
            right_s, right_t = basis_pairs[right]
            bracket = (
                poisson(left_s, left_t)
                if left == right
                else add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            )
            quadratic_right_sides.append(scale(bracket, -field.one))
    quadratic_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in solve_many_particular(
            d2_columns,
            quadratic_right_sides,
            field,
        )
    ]

    s3_monomials, t3_monomials = BOUNDS.correction(3)
    relaxed_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s3_monomials
    ]
    relaxed_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t3_monomials
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
    for s1, t1 in basis_pairs:
        rhs = scale(
            add(pi_power(s1, T, 3), pi_power(S, t1, 3)),
            -field.one / field(24),
        )
        linear_right_sides.append(
            add(rhs, coupling((s1, t1), parity_pair), -field.one)
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
    return {
        "quadratic_particulars": len(quadratic_pairs),
        "relaxed_columns": len(relaxed_columns),
        "relaxed_rank": relaxed_rank,
        "augmented_rank": augmented_rank,
        "directions_killed": augmented_rank - relaxed_rank,
    }, candidate_basis


def residual_plane_compatibility_data(S, T, residual_pairs):
    """Project the genuine hbar^3 equation for the residual plane."""

    field = QQ
    s2_monomials, t2_monomials = BOUNDS.correction(2)
    d2_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    d2_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    parity_rhs = scale(pi_power(S, T, 3), -field.one / field(24))
    parity_vector, kernel_vectors, _ = solve_affine(
        d2_columns,
        parity_rhs,
        field,
    )
    parity_pair = split_correction(
        parity_vector,
        s2_monomials,
        t2_monomials,
    )
    kernel_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in kernel_vectors
    ]

    quadratic_indices = list(combinations_with_replacement(range(3), 2))
    quadratic_right_sides = []
    for left, right in quadratic_indices:
        left_s, left_t = residual_pairs[left]
        right_s, right_t = residual_pairs[right]
        bracket = (
            poisson(left_s, left_t)
            if left == right
            else add(
                poisson(left_s, right_t),
                poisson(right_s, left_t),
            )
        )
        quadratic_right_sides.append(scale(bracket, -field.one))
    quadratic_pairs = [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in solve_many_particular(
            d2_columns,
            quadratic_right_sides,
            field,
        )
    ]

    s3_monomials, t3_monomials = BOUNDS.correction(3)
    d3_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s3_monomials
    ]
    d3_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t3_monomials
    ]
    kernel_couplings = [
        [
            coupling(first_pair, kernel_pair)
            for kernel_pair in kernel_pairs
        ]
        for first_pair in residual_pairs
    ]
    linear_right_sides = []
    for s1, t1 in residual_pairs:
        rhs = scale(
            add(pi_power(s1, T, 3), pi_power(S, t1, 3)),
            -field.one / field(24),
        )
        linear_right_sides.append(
            add(rhs, coupling((s1, t1), parity_pair), -field.one)
        )
    cubic_couplings = [
        [
            coupling(first_pair, quadratic_pair)
            for quadratic_pair in quadratic_pairs
        ]
        for first_pair in residual_pairs
    ]

    all_vectors = (
        d3_columns
        + [
            vector
            for row in kernel_couplings
            for vector in row
        ]
        + linear_right_sides
        + [
            vector
            for row in cubic_couplings
            for vector in row
        ]
    )
    output_monomials = sorted(
        set().union(*(set(vector) for vector in all_vectors))
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }
    transpose_rows = {
        column_index: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(d3_columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose_rows)
    cokernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(output_monomials),
        pivots,
        nonzero,
    )
    incidence = defaultdict(list)
    for functional_index, functional in enumerate(cokernel):
        for output_row, coefficient in functional.items():
            incidence[output_row].append(
                (functional_index, coefficient)
            )

    def project(vector):
        projected = defaultdict(lambda: field.zero)
        for monomial, coefficient in vector.items():
            for functional_index, functional_value in incidence[
                output_index[monomial]
            ]:
                projected[functional_index] += (
                    functional_value * coefficient
                )
        return {
            index: coefficient
            for index, coefficient in projected.items()
            if coefficient
        }

    return {
        "d3_columns": len(d3_columns),
        "d3_rank": len(pivots),
        "d3_cokernel": len(cokernel),
        "kernel_columns": [
            [project(vector) for vector in row]
            for row in kernel_couplings
        ],
        "linear_rhs": [project(vector) for vector in linear_right_sides],
        "cubic_columns": [
            [project(vector) for vector in row]
            for row in cubic_couplings
        ],
        "quadratic_indices": quadratic_indices,
    }


def residual_plane_resonance_audit(data):
    """Prove the exact linear resonance cutting out hbar^3 compatibility."""

    a, b, c = sp.symbols("a b c")
    variables = (a, b, c)
    resonance = 21 * a + 28 * b + 64 * c
    kernel_columns = data["kernel_columns"]
    linear_rhs = data["linear_rhs"]
    cubic_columns = data["cubic_columns"]
    quadratic_indices = data["quadratic_indices"]
    fixed_span = [
        vector for row in kernel_columns for vector in row
    ]
    if column_rank(fixed_span) != 6:
        raise AssertionError("unexpected residual coupling span")
    if column_rank(fixed_span + linear_rhs) != 6:
        raise AssertionError("linear hbar^3 right side left the fixed span")

    def rational(value):
        return sp.Rational(str(value))

    for kernel_index in range(len(kernel_columns[0])):
        rows = set().union(
            *(
                set(kernel_columns[coordinate][kernel_index])
                for coordinate in range(3)
            )
        )
        for row in rows:
            expression = sum(
                rational(
                    kernel_columns[coordinate][kernel_index].get(
                        row,
                        QQ.zero,
                    )
                )
                * variables[coordinate]
                for coordinate in range(3)
            )
            _, remainder = sp.div(
                sp.Poly(expression, a, b, c),
                sp.Poly(resonance, a, b, c),
            )
            if not remainder.is_zero:
                raise AssertionError("kernel coupling missed the resonance")

    linear_rows = set().union(*(set(vector) for vector in linear_rhs))
    for row in linear_rows:
        expression = sum(
            rational(linear_rhs[coordinate].get(row, QQ.zero))
            * variables[coordinate]
            for coordinate in range(3)
        )
        _, remainder = sp.div(
            sp.Poly(expression, a, b, c),
            sp.Poly(resonance, a, b, c),
        )
        if not remainder.is_zero:
            raise AssertionError("linear right side missed the resonance")

    cubic_by_row = defaultdict(lambda: sp.S.Zero)
    for coordinate in range(3):
        for quadratic_index, (left, right) in enumerate(
            quadratic_indices
        ):
            monomial = (
                variables[coordinate]
                * variables[left]
                * variables[right]
            )
            for row, coefficient in cubic_columns[coordinate][
                quadratic_index
            ].items():
                cubic_by_row[row] += rational(coefficient) * monomial
    for expression in cubic_by_row.values():
        _, remainder = sp.div(
            sp.Poly(expression, a, b, c),
            sp.Poly(resonance, a, b, c),
        )
        if not remainder.is_zero:
            raise AssertionError("cubic right side missed the resonance")

    witness = dual_witness(
        fixed_span,
        cubic_columns[0][0],
        QQ,
    )
    obstruction = sp.S.Zero
    for row, expression in cubic_by_row.items():
        obstruction += rational(witness.get(row, QQ.zero)) * expression
    obstruction = sp.factor(obstruction)
    expected = sp.expand(resonance**3 / 21**3)
    if sp.expand(obstruction - expected) != 0:
        raise AssertionError(
            f"unexpected residual obstruction: {obstruction}"
        )
    return {
        "fixed_coupling_rank": 6,
        "resonance_linear_form": "21*a + 28*b + 64*c",
        "obstruction": "(21*a + 28*b + 64*c)^3 / 21^3",
        "coefficientwise_divisibility": {
            "kernel_couplings": True,
            "linear_right_side": True,
            "cubic_right_side": True,
        },
        "hbar3_locus": "21*a + 28*b + 64*c = 0",
        "rational_line_basis": [
            [4, -3, 0],
            [0, 16, -7],
        ],
    }


def fixed_unit_fourth_order_audit(
    S,
    T,
    first_pair,
    field=QQ,
    collect_exceptional_denominator=False,
    direction_label=None,
):
    """Audit hbar^4 over every lower lift above one unit-scale direction."""

    s1, t1 = first_pair
    s2_monomials, t2_monomials = BOUNDS.correction(2)
    s3_monomials, t3_monomials = BOUNDS.correction(3)
    split2 = len(s2_monomials)
    columns2 = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    columns2 += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    coupling_columns = []
    for column, monomial in enumerate(s2_monomials + t2_monomials):
        if column < split2:
            coupling_columns.append(
                poisson({monomial: field.one}, t1)
            )
        else:
            coupling_columns.append(
                poisson(s1, {monomial: field.one})
            )
    columns3 = [
        poisson({monomial: field.one}, T)
        for monomial in s3_monomials
    ]
    columns3 += [
        poisson(S, {monomial: field.one})
        for monomial in t3_monomials
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
    row_index = {
        tagged: index for index, tagged in enumerate(tagged_rows)
    }

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
    rows[len(tagged_rows)] = {
        u_column: field.one,
        rhs_column: -field.one,
    }

    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        raise AssertionError("the certified resonant direction failed hbar^3")
    solution = {}
    for reduced_row, pivot in enumerate(pivots):
        if pivot == rhs_column:
            continue
        value = reduced.get(reduced_row, {}).get(
            rhs_column,
            field.zero,
        )
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
    second_base = split_correction(
        second_vector,
        s2_monomials,
        t2_monomials,
    )
    third_base = split_correction(
        third_vector,
        s3_monomials,
        t3_monomials,
    )

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
    lower_kernel = []
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
        lower_kernel.append(
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
            )
        )

    constant = moyal_coefficient(
        [S, s1, second_base[0], third_base[0]],
        [T, t1, second_base[1], third_base[1]],
        4,
        field,
    )
    nonconstant = []
    for second_kernel, third_kernel in lower_kernel:
        shifted = moyal_coefficient(
            [
                S,
                s1,
                add(second_base[0], second_kernel[0]),
                add(third_base[0], third_kernel[0]),
            ],
            [
                T,
                t1,
                add(second_base[1], second_kernel[1]),
                add(third_base[1], third_kernel[1]),
            ],
            4,
            field,
        )
        diagonal = poisson(second_kernel[0], second_kernel[1])
        linear = add(
            add(shifted, constant, -field.one),
            diagonal,
            -field.one,
        )
        nonconstant.extend((linear, diagonal))
    for left, right in combinations(range(len(lower_kernel)), 2):
        left_second, _ = lower_kernel[left]
        right_second, _ = lower_kernel[right]
        nonconstant.append(
            add(
                poisson(left_second[0], right_second[1]),
                poisson(right_second[0], left_second[1]),
            )
        )

    s4_monomials, t4_monomials = BOUNDS.correction(4)
    correction_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t4_monomials
    ]
    span_columns = correction_columns + nonconstant
    span_rank = column_rank(span_columns)
    augmented_rank = column_rank(span_columns + [constant])
    if augmented_rank != span_rank + 1:
        raise AssertionError("the fixed resonant hbar^4 rank jump failed")
    witness = dual_witness(span_columns, constant, field)
    result = {
        "direction": (
            [4, -3, 0] if direction_label is None else direction_label
        ),
        "unit_scale": True,
        "joint_hbar2_hbar3_rank": len(pivots),
        "lower_lift_dimension": len(lower_kernel),
        "hbar4_correction_columns": len(correction_columns),
        "hbar4_parameter_coefficients": len(nonconstant),
        "hbar4_span_rank": span_rank,
        "hbar4_augmented_rank": augmented_rank,
        "dual_cocycle_support": len(witness),
    }
    if collect_exceptional_denominator:
        denominator = field.gens[0].denom.ring.one

        def collect_poly(poly):
            nonlocal denominator
            for coefficient in poly.values():
                denominator = denominator.lcm(coefficient.denom)

        for matrix in (reduced, homogeneous_reduced):
            for entries in matrix.values():
                collect_poly(entries)
        for poly in (
            second_base
            + third_base
            + tuple(witness for _ in range(1))
        ):
            collect_poly(poly)
        for second_kernel, third_kernel in lower_kernel:
            for poly in second_kernel + third_kernel:
                collect_poly(poly)
        _, factors = denominator.factor_list()
        result["exceptional_denominator"] = str(denominator)
        result["exceptional_factors"] = [
            {
                "factor": str(factor),
                "multiplicity": multiplicity,
                "degree": factor.degree(),
            }
            for factor, multiplicity in factors
        ]
    return result


def resonance_line_fourth_order_audit(
    S,
    T,
    residual_pairs,
    zero_audit,
):
    """Certify hbar^4 obstruction on the complete resonant projective line."""

    zero_direction = linear_combination(
        residual_pairs,
        {0: QQ(4), 1: QQ(-3)},
        QQ,
    )
    infinity_direction = linear_combination(
        residual_pairs,
        {1: QQ(16), 2: QQ(-7)},
        QQ,
    )
    # Guard the direction passed into ``zero_audit`` against bookkeeping
    # changes in the residual-plane basis.
    if zero_audit["direction"] != [4, -3, 0]:
        raise AssertionError("unexpected zero-section direction")
    infinity_audit = fixed_unit_fourth_order_audit(
        S,
        T,
        infinity_direction,
        direction_label=[0, 16, -7],
    )

    parameter = sp.Symbol("t")
    function_field = QQ.frac_field(parameter)
    t = function_field.gens[0]

    def extend(poly):
        return {
            monomial: function_field.convert(coefficient)
            for monomial, coefficient in poly.items()
        }

    generic_s = extend(zero_direction[0])
    generic_t = extend(zero_direction[1])
    for monomial, coefficient in infinity_direction[0].items():
        generic_s[monomial] = (
            generic_s.get(monomial, function_field.zero)
            + t * function_field.convert(coefficient)
        )
    for monomial, coefficient in infinity_direction[1].items():
        generic_t[monomial] = (
            generic_t.get(monomial, function_field.zero)
            + t * function_field.convert(coefficient)
        )
    generic_direction = (
        {
            monomial: coefficient
            for monomial, coefficient in generic_s.items()
            if coefficient
        },
        {
            monomial: coefficient
            for monomial, coefficient in generic_t.items()
            if coefficient
        },
    )
    generic_audit = fixed_unit_fourth_order_audit(
        extend(S),
        extend(T),
        generic_direction,
        function_field,
        collect_exceptional_denominator=True,
        direction_label="e0 + t*e1",
    )
    if generic_audit["exceptional_factors"] != [
        {"factor": "t", "multiplicity": 1, "degree": 1}
    ]:
        raise AssertionError(
            "unexpected exceptional parameter set on the resonance line: "
            f"{generic_audit['exceptional_factors']}"
        )
    for audit in (zero_audit, infinity_audit, generic_audit):
        if audit["hbar4_augmented_rank"] != audit["hbar4_span_rank"] + 1:
            raise AssertionError("resonance-line hbar^4 obstruction failed")
    return {
        "parameterization": {
            "e0": [4, -3, 0],
            "e1": [0, 16, -7],
            "affine_chart": "e0 + t*e1",
        },
        "generic_function_field_audit": generic_audit,
        "exceptional_t_zero_audit": zero_audit,
        "projective_infinity_audit": infinity_audit,
        "coverage": (
            "The Q(t) cocycle has sole denominator factor t; exact audits "
            "at t=0 and t=infinity cover the complete projective line."
        ),
    }


def low_support_audit(
    S,
    T,
    quotient_pairs,
    parameter_monomials,
    equations,
):
    """Classify exact support at most two and relax the remaining P4."""

    variable = sp.Symbol("r")
    full_lines = []
    rational_directions = []
    algebraic_directions = []
    for left, right in combinations(range(len(quotient_pairs)), 2):
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
        for root in sp.roots(common.as_expr(), variable):
            if root.is_Rational:
                rational_directions.append(
                    (left, right, sp.Rational(root))
                )
            else:
                algebraic_directions.append(
                    (left, right, common.as_expr())
                )

    expected_vertices = (0, 1, 7, 8, 17)
    expected_lines = list(combinations(expected_vertices, 2))
    if full_lines != expected_lines:
        raise AssertionError(f"unexpected full coordinate lines: {full_lines}")
    if algebraic_directions:
        raise AssertionError(
            f"unexpected algebraic support-two points: {algebraic_directions}"
        )

    isolated_tests = []
    for left, right, ratio in rational_directions:
        ratio_q = QQ(int(sp.numer(ratio))) / QQ(int(sp.denom(ratio)))
        direction = linear_combination(
            quotient_pairs,
            {left: QQ.one, right: ratio_q},
            QQ,
        )
        result = third_order_axis_test(S, T, direction, QQ)
        isolated_tests.append(
            {
                "coordinates": [left, right],
                "ratio": str(ratio),
                **result,
            }
        )
    if any(test["allows_nonzero_u"] for test in isolated_tests):
        raise AssertionError("an isolated support-two branch reached hbar^3")

    p4_pairs = [quotient_pairs[index] for index in expected_vertices]
    p4_relaxation, p4_candidates = uniform_third_order_relaxation(
        S,
        T,
        p4_pairs,
        QQ,
    )
    expected_candidates = [
        {2: QQ.one, 0: QQ(2)},
        {3: QQ.one, 0: QQ(28) / QQ(9)},
        {4: QQ.one, 0: QQ(824) / QQ(81)},
    ]
    if p4_candidates != expected_candidates:
        raise AssertionError(f"unexpected residual plane: {p4_candidates}")

    residual_global = []
    residual_pairs = []
    residual_basis_tests = []
    for local_vector in p4_candidates:
        global_vector = {
            expected_vertices[index]: coefficient
            for index, coefficient in local_vector.items()
        }
        residual_global.append(global_vector)
        direction = linear_combination(quotient_pairs, global_vector, QQ)
        residual_pairs.append(direction)
        residual_basis_tests.append(
            third_order_axis_test(S, T, direction, QQ)
        )
    if any(
        test["allows_nonzero_u"] for test in residual_basis_tests
    ):
        raise AssertionError("a residual-plane basis direction reached hbar^3")

    p2_relaxation, p2_candidates = uniform_third_order_relaxation(
        S,
        T,
        residual_pairs,
        QQ,
    )
    if p2_candidates != [
        {0: QQ.one},
        {1: QQ.one},
        {2: QQ.one},
    ]:
        raise AssertionError(
            f"the residual plane unexpectedly shrank: {p2_candidates}"
        )
    compatibility_data = residual_plane_compatibility_data(
        S,
        T,
        residual_pairs,
    )
    resonance = residual_plane_resonance_audit(compatibility_data)
    first_resonant_direction = None
    for coordinates in (
        (4, -3, 0),
        (0, 16, -7),
        (4, 13, -7),
    ):
        direction = linear_combination(
            residual_pairs,
            {
                index: QQ(value)
                for index, value in enumerate(coordinates)
                if value
            },
            QQ,
        )
        result = third_order_axis_test(S, T, direction, QQ)
        if not result["allows_nonzero_u"]:
            raise AssertionError(
                f"resonant direction failed hbar^3: {coordinates}"
            )
        if coordinates == (4, -3, 0):
            first_resonant_direction = direction
    fourth_order = fixed_unit_fourth_order_audit(
        S,
        T,
        first_resonant_direction,
    )
    resonance_line_fourth = resonance_line_fourth_order_audit(
        S,
        T,
        residual_pairs,
        fourth_order,
    )

    def encode_vector(vector):
        return {
            str(index): encode_rational(coefficient)
            for index, coefficient in sorted(vector.items())
        }

    return {
        "full_coordinate_subspace": list(expected_vertices),
        "full_coordinate_lines": [list(line) for line in full_lines],
        "isolated_rational_directions": isolated_tests,
        "isolated_algebraic_directions": [],
        "p4_relaxation": p4_relaxation,
        "necessary_residual_plane_basis": [
            encode_vector(vector) for vector in residual_global
        ],
        "residual_basis_tests": residual_basis_tests,
        "p2_relaxation": p2_relaxation,
        "resonance_audit": resonance,
        "fixed_resonant_fourth_order_audit": fourth_order,
        "resonance_line_fourth_order_audit": resonance_line_fourth,
        "scope": (
            "The displayed resonance line is the exact hbar^3 locus. "
            "A Q(t) dual cocycle together with exact exceptional audits "
            "obstructs the complete projective resonance line at hbar^4."
        ),
    }


def unrestricted_audit(S, T, field):
    s1_monomials, t1_monomials = BOUNDS.correction(1)
    columns, outputs, reduced, pivots, nonzero = operator_rref(
        S,
        T,
        s1_monomials,
        t1_monomials,
        field,
    )
    kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(columns),
        pivots,
        nonzero,
    )
    free_columns = [
        column for column in range(len(columns)) if column not in pivots
    ]
    gauge_pivots, gauge_directions = target_gauge_pivots(
        S,
        T,
        s1_monomials,
        t1_monomials,
        free_columns,
        kernel,
        field,
    )
    essential = [
        index for index in range(len(kernel)) if index not in gauge_pivots
    ]
    quotient_pairs = [
        split_pair(kernel[index], s1_monomials, t1_monomials)
        for index in essential
    ]
    second, parameter_monomials, equations = project_second_obstruction(
        S,
        T,
        quotient_pairs,
        field,
    )
    axis_tests = {
        str(axis): third_order_axis_test(
            S,
            T,
            quotient_pairs[axis],
            field,
        )
        for axis in second["surviving_axes"]
    }
    audit = {
        "h1_s_columns": len(s1_monomials),
        "h1_t_columns": len(t1_monomials),
        "h1_rank": len(pivots),
        "h1_nullity": len(kernel),
        "h1_output_monomials": len(outputs),
        "gauge_directions": gauge_directions,
        "gauge_rank": len(gauge_pivots),
        "quotient_dimension": len(essential),
        **second,
        "third_order_axis_tests": axis_tests,
        "third_order_surviving_axes": [
            int(axis)
            for axis, result in axis_tests.items()
            if result["allows_nonzero_u"]
        ],
    }
    problem = {
        "quotient_pairs": quotient_pairs,
        "parameter_monomials": parameter_monomials,
        "equations": equations,
    }
    return audit, problem


def poly_stats(poly):
    return {
        "terms": len(poly),
        "bernstein_degree": max(
            x_degree + q_degree + 3 * z_degree
            for x_degree, q_degree, z_degree in poly
        ),
        "z_order": max(z_degree for _, _, z_degree in poly),
    }


def encode_rational(value):
    numerator, denominator = value.numerator, value.denominator
    return str(numerator) if denominator == 1 else f"{numerator}/{denominator}"


def certificate_payload(symbols, parity, unrestricted, low_support):
    return {
        "claim": (
            "The parity-preserving normal-ordered lift in the inherited "
            "(4,3) Bernstein filtration is obstructed at hbar^5."
        ),
        "specialization": {"a": "-4/3", "tau": "0"},
        "symbol_stats": symbols,
        "bounds": {
            "S": {"bernstein_degree": 22, "z_order": 4},
            "T": {"bernstein_degree": 18, "z_order": 3},
            "rule": "hbar^n lowers Bernstein degree by 2n and Z-order by n",
        },
        "parity_audit": {
            key: value for key, value in parity.items() if key != "witness"
        },
        "dual_cocycle": [
            {
                "monomial": list(monomial),
                "coefficient": encode_rational(coefficient),
            }
            for monomial, coefficient in sorted(parity["witness"].items())
        ],
        "dual_pairing": {
            "allowed_hbar5_corrections": "0",
            "all_affine_hbar3_parameter_coefficients": "0",
            "constant_hbar5_defect": "1",
        },
        "unrestricted_audit": unrestricted,
        "low_support_audit": low_support,
        "software": {"python_dependencies": "requirements.txt"},
    }


def run(field):
    a = field_fraction(field, SPECIAL_A)
    tau = field_fraction(field, SPECIAL_TAU)
    S, T = degree_five_family(field, a, tau)
    symbols = {"S": poly_stats(S), "T": poly_stats(T)}
    if symbols != {
        "S": {"terms": 33, "bernstein_degree": 22, "z_order": 4},
        "T": {"terms": 23, "bernstein_degree": 18, "z_order": 3},
    }:
        raise AssertionError(f"unexpected specialized supports: {symbols}")
    unrestricted, problem = unrestricted_audit(S, T, field)
    return symbols, parity_audit(S, T, field), unrestricted, problem


def comparable(audit):
    return {
        key: value
        for key, value in audit.items()
        if key != "witness"
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--certificate",
        type=Path,
        help="write the exact rational dual-cocycle certificate as JSON",
    )
    args = parser.parse_args()

    (
        symbols,
        exact_parity,
        exact_unrestricted,
        exact_problem,
    ) = run(QQ)
    (
        finite_symbols,
        finite_parity,
        finite_unrestricted,
        _,
    ) = run(GF(PRIME))
    assert symbols == finite_symbols
    assert comparable(exact_parity) == comparable(finite_parity)
    assert exact_unrestricted == finite_unrestricted
    exact_S, exact_T = degree_five_family(
        QQ,
        -QQ(4) / QQ(3),
        QQ.zero,
    )
    low_support = low_support_audit(
        exact_S,
        exact_T,
        exact_problem["quotient_pairs"],
        exact_problem["parameter_monomials"],
        exact_problem["equations"],
    )

    print(
        "PASS: specialized symbols "
        f"S=(degree {symbols['S']['bernstein_degree']}, "
        f"order {symbols['S']['z_order']}), "
        f"T=(degree {symbols['T']['bernstein_degree']}, "
        f"order {symbols['T']['z_order']})"
    )
    print(
        "PASS: exact hbar^3 affine space has "
        f"{exact_parity['h3_s_columns']}+"
        f"{exact_parity['h3_t_columns']} columns, "
        f"rank {exact_parity['h3_rank']}, "
        f"nullity {exact_parity['h3_nullity']}"
    )
    print(
        "PASS: exact hbar^5 restricted obstruction has correction rank "
        f"{exact_parity['h5_correction_rank']} and span rank jump "
        f"{exact_parity['h5_span_rank']}->"
        f"{exact_parity['h5_augmented_rank']}; "
        f"dual cocycle support={len(exact_parity['witness'])}"
    )
    print(
        "PASS: unrestricted hbar^1 operator has "
        f"{exact_unrestricted['h1_s_columns']}+"
        f"{exact_unrestricted['h1_t_columns']} columns, "
        f"rank {exact_unrestricted['h1_rank']}, "
        f"nullity {exact_unrestricted['h1_nullity']}; "
        f"gauge rank {exact_unrestricted['gauge_rank']}; "
        f"quotient dimension {exact_unrestricted['quotient_dimension']}"
    )
    print(
        "PASS: projected hbar^2 quadratic obstruction has rank "
        f"{exact_unrestricted['quadratic_rank']} on "
        f"{exact_unrestricted['quadratic_coefficients']} coefficients; "
        f"surviving coordinate axes="
        f"{exact_unrestricted['surviving_axes']}"
    )
    print(
        "PASS: coordinate axes surviving the coupled hbar^2/hbar^3 test="
        f"{exact_unrestricted['third_order_surviving_axes']}"
    )
    print(
        "PASS: exact support <=2 is one coordinate P4 plus "
        f"{len(low_support['isolated_rational_directions'])} isolated "
        "rational directions; every isolated direction fails by hbar^3"
    )
    print(
        "PASS: uniform hbar^3 relaxation reduces the P4 to the explicit "
        "residual P2; exact compatibility on that plane is the resonance "
        "line 21*a+28*b+64*c=0"
    )
    fixed_fourth = low_support["fixed_resonant_fourth_order_audit"]
    print(
        "PASS: the complete projective resonance line is obstructed at "
        "hbar^4 over every compatible unit-scale lower lift; the generic "
        "Q(t) rank jump is "
        f"{fixed_fourth['hbar4_span_rank']}->"
        f"{fixed_fourth['hbar4_augmented_rank']}, with sole affine "
        "exception t=0 checked exactly together with projective infinity"
    )
    print(f"PASS: all discrete ranks agree over Q and GF({PRIME})")

    if args.certificate:
        args.certificate.parent.mkdir(parents=True, exist_ok=True)
        args.certificate.write_text(
            json.dumps(
                certificate_payload(
                    symbols,
                    exact_parity,
                    exact_unrestricted,
                    low_support,
                ),
                indent=2,
            )
            + "\n"
        )
        print(f"WROTE: {args.certificate}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Analyze the genuine bounded hbar^5 equations on the cubic period locus.

Four bounded fifth-order periods leave the reduced characteristic-zero
component

    94*a^3 + 335*a^2 + 400*a + 160 = 0,
    8*tau + 658*a^2 + 1593*a + 976 = 0.

At its finite-field points the constant fifth-order defect enters the
*linear span* of the correction and lower-lift coefficient vectors.  That
only makes the span obstruction inconclusive: the linear, square, and
cross-term coefficients must still be evaluated at one common lower-lift
parameter vector.

This script constructs those actual quadratic equations.  It projects the
defect modulo the bounded hbar^5 correction image using a basis of the left
kernel, records the exact parameter-monomial structure, and extracts every
linear consequence obtained by cancelling all quadratic terms.  An optional
Singular calculation tests the resulting quadratic ideal.  Over the exact
cubic field it also verifies a sparse K-rational hbar^5 solution and can
test its affine-linear hbar^7 extension problem.

The default point is one root of the cubic component over GF(32003).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from itertools import combinations
from pathlib import Path

from sympy import Poly, symbols
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    GENERIC_PROFILE,
    add,
    degree_five_family,
    fifth_order_coefficients,
    laurent_monomials,
    scale,
    solve_affine,
    third_order_family,
)


DEFAULT_PRIME = 32003
DEFAULT_A = -7418
DEFAULT_TAU = 25070


def sparse_nullspace(rows, column_count, field):
    """Return the right nullspace of a sparse row matrix."""

    reduced, pivots, nonzero = sdm_irref(rows)
    kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        column_count,
        pivots,
        nonzero,
    )
    return kernel, len(pivots)


def affine_particular(columns, rhs, field):
    """Solve columns*x=rhs without constructing its homogeneous kernel."""

    coordinates = sorted(
        set(rhs).union(*(set(column) for column in columns))
    )
    coordinate_index = {
        coordinate: index for index, coordinate in enumerate(coordinates)
    }
    rhs_column = len(columns)
    rows = {}
    for column_index, column in enumerate(columns):
        for coordinate, coefficient in column.items():
            rows.setdefault(coordinate_index[coordinate], {})[
                column_index
            ] = coefficient
    for coordinate, coefficient in rhs.items():
        rows.setdefault(coordinate_index[coordinate], {})[
            rhs_column
        ] = -coefficient
    reduced, pivots, _ = sdm_irref(rows)
    if rhs_column in pivots:
        raise ValueError("affine system is inconsistent")
    particular = {}
    for reduced_row, pivot in enumerate(pivots):
        rhs_value = reduced.get(reduced_row, {}).get(
            rhs_column,
            field.zero,
        )
        if rhs_value:
            particular[pivot] = -rhs_value
    return particular, len(pivots)


def dot(functional, polynomial, output_index, field):
    value = field.zero
    for monomial, coefficient in polynomial.items():
        coordinate = output_index.get(monomial)
        if coordinate is not None:
            value += functional.get(coordinate, field.zero) * coefficient
    return value


def modular_int(value, prime: int) -> int:
    """Return the least nonnegative integer representing a field element."""

    return int(value) % prime


def singular_polynomial(
    coefficients,
    variable_count: int,
    coefficient_text,
) -> str:
    constant, linear, diagonal, cross = coefficients
    terms = []
    if constant:
        terms.append(coefficient_text(constant))
    for index, value in enumerate(linear, start=1):
        if value:
            terms.append(f"({coefficient_text(value)})*u{index}")
    for index, value in enumerate(diagonal, start=1):
        if value:
            terms.append(f"({coefficient_text(value)})*u{index}^2")
    for (left, right), value in zip(
        combinations(range(1, variable_count + 1), 2),
        cross,
    ):
        if value:
            terms.append(
                f"({coefficient_text(value)})*u{left}*u{right}"
            )
    return "+".join(terms) if terms else "0"


def cubic_number_field():
    indeterminate = symbols("x")
    polynomial = Poly(
        94 * indeterminate**3
        + 335 * indeterminate**2
        + 400 * indeterminate
        + 160,
        indeterminate,
    )
    field = QQ.algebraic_field((polynomial, "a"))
    a = field.unit
    tau = -(
        field(658) * a**2
        + field(1593) * a
        + field(976)
    ) / field(8)
    return field, a, tau


def construct_projected_system(
    prime: int | None,
    a_value: int,
    tau_value: int,
):
    if prime is None:
        field, a, tau = cubic_number_field()
        field_label = "Q(a)"
    else:
        field = GF(prime)
        a = field(a_value)
        tau = field(tau_value)
        field_label = f"GF({prime})"
    cubic = 94 * a**3 + 335 * a**2 + 400 * a + 160
    tau_relation = 8 * tau + 658 * a**2 + 1593 * a + 976
    if cubic or tau_relation:
        raise ValueError(
            "the selected finite-field point is not on the rational "
            f"cubic locus: cubic={cubic}, tau relation={tau_relation}"
        )
    print(f"coefficient field: {field_label}", flush=True)

    S, T = degree_five_family(field, a, tau)
    third = third_order_family(S, T, field)
    print(
        "h3:",
        f"columns={third.column_count}",
        f"rank={third.operator_rank}",
        f"kernel={len(third.kernel)}",
        flush=True,
    )
    constant, nonconstant = fifth_order_coefficients(
        S,
        T,
        third,
        field,
    )
    variable_count = len(third.kernel)
    expected_coefficients = (
        2 * variable_count
        + variable_count * (variable_count - 1) // 2
    )
    if len(nonconstant) != expected_coefficients:
        raise AssertionError((len(nonconstant), expected_coefficients))

    s4_monomials = laurent_monomials(
        GENERIC_PROFILE.s4_degree,
        1,
        0,
        GENERIC_PROFILE.z_weight,
    )
    t4_monomials = laurent_monomials(
        GENERIC_PROFILE.t4_degree,
        0,
        0,
        GENERIC_PROFILE.z_weight,
    )
    corrections = [
        GENERIC_PROFILE.poisson({monomial: field.one}, T)
        for monomial in s4_monomials
    ]
    corrections += [
        GENERIC_PROFILE.poisson(S, {monomial: field.one})
        for monomial in t4_monomials
    ]
    output_monomials = sorted(
        set(constant).union(
            *(set(column) for column in corrections + nonconstant)
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }

    # Rows of C^T.  Its right kernel is the left kernel of the correction
    # operator C, hence a basis of linear functionals on coker(C).
    transpose_rows = {
        row: {
            output_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for row, column in enumerate(corrections)
        if column
    }
    functionals, correction_rank = sparse_nullspace(
        transpose_rows,
        len(output_monomials),
        field,
    )
    print(
        "h5 correction:",
        f"columns={len(corrections)}",
        f"outputs={len(output_monomials)}",
        f"rank={correction_rank}",
        f"cokernel={len(functionals)}",
        flush=True,
    )

    linear_columns = nonconstant[0 : 2 * variable_count : 2]
    diagonal_columns = nonconstant[1 : 2 * variable_count : 2]
    cross_columns = nonconstant[2 * variable_count :]
    equations = []
    for functional in functionals:
        coefficients = (
            dot(functional, constant, output_index, field),
            [
                dot(functional, column, output_index, field)
                for column in linear_columns
            ],
            [
                dot(functional, column, output_index, field)
                for column in diagonal_columns
            ],
            [
                dot(functional, column, output_index, field)
                for column in cross_columns
            ],
        )
        if any(
            value
            for part in coefficients
            for value in (part if isinstance(part, list) else [part])
        ):
            equations.append(coefficients)

    quadratic_rows = {}
    quadratic_width = (
        variable_count
        + variable_count * (variable_count - 1) // 2
    )
    for quadratic_coordinate in range(quadratic_width):
        entries = {}
        for equation_index, (_, _, diagonal, cross) in enumerate(equations):
            values = diagonal + cross
            value = values[quadratic_coordinate]
            if value:
                entries[equation_index] = value
        if entries:
            quadratic_rows[quadratic_coordinate] = entries
    linear_consequence_weights, quadratic_rank = sparse_nullspace(
        quadratic_rows,
        len(equations),
        field,
    )

    linear_consequences = []
    for weights in linear_consequence_weights:
        constant_value = field.zero
        linear_values = [field.zero] * variable_count
        for equation_index, weight in weights.items():
            equation_constant, equation_linear, _, _ = equations[
                equation_index
            ]
            constant_value += weight * equation_constant
            for variable_index, value in enumerate(equation_linear):
                linear_values[variable_index] += weight * value
        if constant_value or any(linear_values):
            linear_consequences.append((constant_value, linear_values))

    linear_rows = {
        row: {
            coordinate: value
            for coordinate, value in enumerate(values)
            if value
        }
        for row, (_, values) in enumerate(linear_consequences)
    }
    _, linear_rank = sparse_nullspace(
        linear_rows,
        variable_count,
        field,
    )
    augmented_linear_rows = {
        row: {
            **{
                coordinate: value
                for coordinate, value in enumerate(values)
                if value
            },
            **(
                {variable_count: -constant_value}
                if constant_value
                else {}
            ),
        }
        for row, (constant_value, values) in enumerate(linear_consequences)
    }
    _, augmented_linear_rank = sparse_nullspace(
        augmented_linear_rows,
        variable_count + 1,
        field,
    )
    active_variables = {
        index
        for _, linear, diagonal, cross in equations
        for index, value in enumerate(linear)
        if value
    }
    active_variables.update(
        index
        for _, _, diagonal, _ in equations
        for index, value in enumerate(diagonal)
        if value
    )
    for _, _, _, cross in equations:
        active_variables.update(
            index
            for (index, other), value in zip(
                combinations(range(variable_count), 2),
                cross,
            )
            if value
        )
        active_variables.update(
            other
            for (index, other), value in zip(
                combinations(range(variable_count), 2),
                cross,
            )
            if value
        )
    print(
        "projected equations:",
        f"nonzero={len(equations)}",
        f"quadratic-rank={quadratic_rank}",
        f"linear-consequences={len(linear_consequences)}",
        f"linear-rank={linear_rank}",
        f"augmented-linear-rank={augmented_linear_rank}",
        f"active-variables={len(active_variables)}/{variable_count}",
        flush=True,
    )
    if augmented_linear_rank > linear_rank:
        print("VERDICT=OBSTRUCTED_BY_LINEAR_CONSEQUENCE", flush=True)
    else:
        print("VERDICT=NONLINEAR_SYSTEM_REMAINS", flush=True)
    return (
        equations,
        variable_count,
        field,
        S,
        T,
        third,
        constant,
        nonconstant,
        corrections,
        s4_monomials,
        t4_monomials,
    )


def verify_explicit_cubic_branch(
    equations,
    variable_count: int,
    field,
    S,
    T,
    third,
    constant,
    nonconstant,
    corrections,
    s4_monomials,
    t4_monomials,
    seventh_order: bool = False,
    branch_u36=None,
) -> None:
    """Verify a sparse K-rational solution of all projected quadrics."""

    if variable_count != 41:
        raise AssertionError(variable_count)
    a = field.unit
    parameters = [field.zero] * variable_count
    parameters[33] = -(
        field(101668771215) * a**2 / field(2097152)
        + field(487549466415) * a / field(4194304)
        + field(18474132105) / field(262144)
    )
    parameters[18] = (
        field(243)
        * (
            field(122116896574) * a**2
            + field(292049112895) * a
            + field(176583624080)
        )
        / field(10485760)
    )
    if branch_u36 is not None and branch_u36:
        parameters[35] = branch_u36
        parameters[34] = -field(4) * branch_u36 / field(3)
        parameters[24] = -(
            field(32) * a + field(24)
        ) * branch_u36 / field(9)
    nonzero_residuals = 0
    for equation_constant, linear, diagonal, cross in equations:
        value = equation_constant
        value += sum(
            coefficient * parameter
            for coefficient, parameter in zip(linear, parameters)
        )
        value += sum(
            coefficient * parameter**2
            for coefficient, parameter in zip(diagonal, parameters)
        )
        value += sum(
            coefficient * parameters[left] * parameters[right]
            for coefficient, (left, right) in zip(
                cross,
                combinations(range(variable_count), 2),
            )
        )
        nonzero_residuals += bool(value)
    if nonzero_residuals:
        raise AssertionError(
            f"explicit cubic branch leaves {nonzero_residuals} residuals"
        )
    print(
        "EXPLICIT_K_POINT_NONZERO_PARAMETERS="
        + ",".join(
            f"u{index + 1}"
            for index, value in enumerate(parameters)
            if value
        ),
        flush=True,
    )
    print(
        "EXPLICIT_K_POINT_RESIDUALS=0/"
        + str(len(equations)),
        flush=True,
    )

    # Reassemble the actual fifth-order defect at the same parameter point
    # and solve d5(S4,T4)=-O5.  This independently turns the cokernel
    # vanishing statement into an explicit correction-vector certificate.
    defect = dict(constant)
    for index, parameter in enumerate(parameters):
        defect = add(
            defect,
            nonconstant[2 * index],
            parameter,
        )
        defect = add(
            defect,
            nonconstant[2 * index + 1],
            parameter**2,
        )
    cross_offset = 2 * variable_count
    for cross_index, (left, right) in enumerate(
        combinations(range(variable_count), 2)
    ):
        defect = add(
            defect,
            nonconstant[cross_offset + cross_index],
            parameters[left] * parameters[right],
        )
    correction, correction_kernel, correction_rank = solve_affine(
        corrections,
        scale(defect, -field.one),
        field,
    )
    reconstructed = dict(defect)
    for column, coefficient in correction.items():
        reconstructed = add(
            reconstructed,
            corrections[column],
            coefficient,
        )
    if reconstructed:
        raise AssertionError(
            f"explicit hbar^5 correction leaves {len(reconstructed)} terms"
        )
    print(
        "EXPLICIT_H5_CORRECTION="
        f"rank{correction_rank},kernel{len(correction_kernel)},"
        f"nonzero{len(correction)}",
        flush=True,
    )
    if not seventh_order:
        return

    def split_correction(vector):
        split = len(s4_monomials)
        return (
            {
                s4_monomials[index]: coefficient
                for index, coefficient in vector.items()
                if index < split
            },
            {
                t4_monomials[index - split]: coefficient
                for index, coefficient in vector.items()
                if index >= split
            },
        )

    base_s2, base_t2 = third.base
    s2 = dict(base_s2)
    t2 = dict(base_t2)
    for index, parameter in enumerate(parameters):
        direction_s, direction_t = third.kernel[index]
        s2 = add(s2, direction_s, parameter)
        t2 = add(t2, direction_t, parameter)
    s4, t4 = split_correction(correction)

    def seventh_defect(left4, right4):
        value = GENERIC_PROFILE.poisson(s2, right4)
        value = add(
            value,
            GENERIC_PROFILE.poisson(left4, t2),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(left4, T, 3),
            field.one / field(24),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(s2, t2, 3),
            field.one / field(24),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(S, right4, 3),
            field.one / field(24),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(s2, T, 5),
            field.one / field(1920),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(S, t2, 5),
            field.one / field(1920),
        )
        value = add(
            value,
            GENERIC_PROFILE.pi_power(S, T, 7),
            field.one / field(322560),
        )
        return value

    seventh_constant = seventh_defect(s4, t4)
    seventh_variations = []
    for kernel_vector in correction_kernel:
        direction_s4, direction_t4 = split_correction(kernel_vector)
        variation = GENERIC_PROFILE.poisson(s2, direction_t4)
        variation = add(
            variation,
            GENERIC_PROFILE.poisson(direction_s4, t2),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(direction_s4, T, 3),
            field.one / field(24),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(S, direction_t4, 3),
            field.one / field(24),
        )
        seventh_variations.append(variation)
    try:
        seventh_parameters, seventh_kernel, seventh_rank = solve_affine(
            seventh_variations,
            scale(seventh_constant, -field.one),
            field,
        )
    except ValueError:
        augmented_rank = None
        # Recover the two ranks for a precise obstruction report.
        from explore_degree_five_quantum_residue import column_rank

        seventh_rank = column_rank(seventh_variations)
        augmented_rank = column_rank(
            seventh_variations + [seventh_constant]
        )
        print(
            "EXPLICIT_H7_EXTENSION=OBSTRUCTED,"
            f"variables={len(correction_kernel)},"
            f"ranks={seventh_rank}->{augmented_rank},"
            f"constant_terms={len(seventh_constant)}",
            flush=True,
        )
        output_monomials = sorted(
            set(seventh_constant).union(
                *(set(column) for column in seventh_variations)
            )
        )
        constant_equation = len(seventh_variations)
        dual_columns = []
        for monomial in output_monomials:
            column = {
                index: variation.get(monomial, field.zero)
                for index, variation in enumerate(seventh_variations)
                if variation.get(monomial, field.zero)
            }
            constant_value = seventh_constant.get(
                monomial,
                field.zero,
            )
            if constant_value:
                column[constant_equation] = constant_value
            dual_columns.append(column)
        dual_vector, dual_rank = affine_particular(
            dual_columns,
            {constant_equation: field.one},
            field,
        )
        dual_functional = {
            output_monomials[index]: coefficient
            for index, coefficient in dual_vector.items()
            if coefficient
        }
        print(
            "EXPLICIT_H7_PERIOD="
            f"support{len(dual_functional)},rank{dual_rank},value1",
            flush=True,
        )
        for monomial, coefficient in dual_functional.items():
            print(
                f"  H7_PERIOD_TERM={monomial}:{field.to_sympy(coefficient)}",
                flush=True,
            )
        return

    reconstructed_seventh = dict(seventh_constant)
    for index, coefficient in seventh_parameters.items():
        reconstructed_seventh = add(
            reconstructed_seventh,
            seventh_variations[index],
            coefficient,
        )
    if reconstructed_seventh:
        raise AssertionError(
            f"hbar^7 extension leaves {len(reconstructed_seventh)} terms"
        )
    print(
        "EXPLICIT_H7_EXTENSION=SOLVABLE,"
        f"variables={len(correction_kernel)},"
        f"rank={seventh_rank},"
        f"kernel={len(seventh_kernel)},"
        f"nonzero={len(seventh_parameters)}",
        flush=True,
    )


def analyze_seventh_order_solution_line(
    field,
    S,
    T,
    third,
    constant,
    nonconstant,
    corrections,
    s4_monomials,
    t4_monomials,
) -> None:
    """Prove or isolate exceptions to hbar^7 obstruction on the u36-line."""

    variable_count = len(third.kernel)
    if variable_count != 41:
        raise AssertionError(variable_count)
    a = field.unit

    def parameters_at(value):
        parameters = [field.zero] * variable_count
        parameters[33] = -(
            field(101668771215) * a**2 / field(2097152)
            + field(487549466415) * a / field(4194304)
            + field(18474132105) / field(262144)
        )
        parameters[18] = (
            field(243)
            * (
                field(122116896574) * a**2
                + field(292049112895) * a
                + field(176583624080)
            )
            / field(10485760)
        )
        parameters[35] = value
        parameters[34] = -field(4) * value / field(3)
        parameters[24] = -(
            field(32) * a + field(24)
        ) * value / field(9)
        return parameters

    def fifth_defect(parameters):
        defect = dict(constant)
        for index, parameter in enumerate(parameters):
            defect = add(defect, nonconstant[2 * index], parameter)
            defect = add(
                defect,
                nonconstant[2 * index + 1],
                parameter**2,
            )
        cross_offset = 2 * variable_count
        for cross_index, (left, right) in enumerate(
            combinations(range(variable_count), 2)
        ):
            defect = add(
                defect,
                nonconstant[cross_offset + cross_index],
                parameters[left] * parameters[right],
            )
        return defect

    zero_parameters = parameters_at(field.zero)
    plus_parameters = parameters_at(field.one)
    minus_parameters = parameters_at(-field.one)
    defect_zero = fifth_defect(zero_parameters)
    defect_plus = fifth_defect(plus_parameters)
    defect_minus = fifth_defect(minus_parameters)
    defect_linear = scale(
        add(defect_plus, defect_minus, -field.one),
        field.one / field(2),
    )
    defect_quadratic = add(
        scale(add(defect_plus, defect_minus), field.one / field(2)),
        defect_zero,
        -field.one,
    )

    correction_coefficients = []
    correction_kernel = None
    for coefficient in (
        defect_zero,
        defect_linear,
        defect_quadratic,
    ):
        solution, kernel, rank = solve_affine(
            corrections,
            scale(coefficient, -field.one),
            field,
        )
        if rank != 594 or len(kernel) != 20:
            raise AssertionError((rank, len(kernel)))
        correction_coefficients.append(solution)
        if correction_kernel is None:
            correction_kernel = kernel
    print(
        "H7_LINE_H5_CORRECTION_DEGREES=0,1,2;"
        "rank594,kernel20",
        flush=True,
    )

    def split_correction(vector):
        split = len(s4_monomials)
        return (
            {
                s4_monomials[index]: coefficient
                for index, coefficient in vector.items()
                if index < split
            },
            {
                t4_monomials[index - split]: coefficient
                for index, coefficient in vector.items()
                if index >= split
            },
        )

    correction_pairs = [
        split_correction(vector)
        for vector in correction_coefficients
    ]
    kernel_pairs = [
        split_correction(vector)
        for vector in correction_kernel
    ]
    base_s2, base_t2 = third.base
    parameter_zero = zero_parameters
    parameter_direction = [
        plus - zero
        for plus, zero in zip(plus_parameters, zero_parameters)
    ]

    def lower_pair(parameters):
        s2 = dict(base_s2)
        t2 = dict(base_t2)
        for index, parameter in enumerate(parameters):
            direction_s, direction_t = third.kernel[index]
            s2 = add(s2, direction_s, parameter)
            t2 = add(t2, direction_t, parameter)
        return s2, t2

    s2_zero, t2_zero = lower_pair(parameter_zero)
    s2_direction, t2_direction = lower_pair(parameter_direction)
    # ``lower_pair`` includes the affine base.  Remove it for the direction.
    s2_direction = add(s2_direction, base_s2, -field.one)
    t2_direction = add(t2_direction, base_t2, -field.one)

    def evaluate_correction(value):
        s4 = {}
        t4 = {}
        power = field.one
        for coefficient_s4, coefficient_t4 in correction_pairs:
            s4 = add(s4, coefficient_s4, power)
            t4 = add(t4, coefficient_t4, power)
            power *= value
        return s4, t4

    def evaluate_lower(value):
        return (
            add(s2_zero, s2_direction, value),
            add(t2_zero, t2_direction, value),
        )

    def seventh_defect_at(value):
        s2, t2 = evaluate_lower(value)
        s4, t4 = evaluate_correction(value)
        result = GENERIC_PROFILE.poisson(s2, t4)
        result = add(result, GENERIC_PROFILE.poisson(s4, t2))
        result = add(
            result,
            GENERIC_PROFILE.pi_power(s4, T, 3),
            field.one / field(24),
        )
        result = add(
            result,
            GENERIC_PROFILE.pi_power(s2, t2, 3),
            field.one / field(24),
        )
        result = add(
            result,
            GENERIC_PROFILE.pi_power(S, t4, 3),
            field.one / field(24),
        )
        result = add(
            result,
            GENERIC_PROFILE.pi_power(s2, T, 5),
            field.one / field(1920),
        )
        result = add(
            result,
            GENERIC_PROFILE.pi_power(S, t2, 5),
            field.one / field(1920),
        )
        result = add(
            result,
            GENERIC_PROFILE.pi_power(S, T, 7),
            field.one / field(322560),
        )
        return result

    values = {
        integer: seventh_defect_at(field(integer))
        for integer in (0, 1, -1, 2)
    }
    constant_coefficients = [None] * 4
    constant_coefficients[0] = values[0]
    odd_part = scale(
        add(values[1], values[-1], -field.one),
        field.one / field(2),
    )
    constant_coefficients[2] = add(
        scale(add(values[1], values[-1]), field.one / field(2)),
        values[0],
        -field.one,
    )
    constant_coefficients[3] = scale(
        add(
            add(
                values[2],
                values[0],
                -field.one,
            ),
            constant_coefficients[2],
            -field(4),
        ),
        field.one / field(6),
    )
    constant_coefficients[3] = add(
        constant_coefficients[3],
        odd_part,
        -field.one / field(3),
    )
    # The previous two lines equal
    # (f(2)-f(0)-4*c2-2*(c1+c3))/6.
    constant_coefficients[1] = add(
        odd_part,
        constant_coefficients[3],
        -field.one,
    )
    for integer, expected in values.items():
        reconstructed = {}
        power = field.one
        for coefficient in constant_coefficients:
            reconstructed = add(reconstructed, coefficient, power)
            power *= field(integer)
        if reconstructed != expected:
            raise AssertionError(
                f"cubic interpolation failed at r={integer}"
            )

    def seventh_variation(pair, value):
        direction_s4, direction_t4 = pair
        s2, t2 = evaluate_lower(value)
        variation = GENERIC_PROFILE.poisson(s2, direction_t4)
        variation = add(
            variation,
            GENERIC_PROFILE.poisson(direction_s4, t2),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(direction_s4, T, 3),
            field.one / field(24),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(S, direction_t4, 3),
            field.one / field(24),
        )
        return variation

    variation_coefficients = []
    for pair in kernel_pairs:
        value_zero = seventh_variation(pair, field.zero)
        value_one = seventh_variation(pair, field.one)
        variation_coefficients.append(
            (value_zero, add(value_one, value_zero, -field.one))
        )

    fraction_field = field.frac_field("r")
    r = fraction_field.gens[0]
    polynomial_ring = fraction_field.get_ring()

    def promote(coefficients):
        result = {}
        monomials = set().union(*(set(value) for value in coefficients))
        for monomial in monomials:
            coefficient = fraction_field.zero
            power = fraction_field.one
            for value in coefficients:
                coefficient += fraction_field(
                    value.get(monomial, field.zero)
                ) * power
                power *= r
            if coefficient:
                result[monomial] = coefficient
        return result

    variation_columns = [
        promote(coefficients)
        for coefficients in variation_coefficients
    ]
    line_constant = promote(constant_coefficients)
    output_monomials = sorted(
        set(line_constant).union(
            *(set(column) for column in variation_columns)
        )
    )
    output_index = {
        monomial: index for index, monomial in enumerate(output_monomials)
    }

    def generic_column_pivots(column_order):
        rows = {}
        for reordered_index, original_index in enumerate(column_order):
            for monomial, coefficient in variation_columns[
                original_index
            ].items():
                rows.setdefault(output_index[monomial], {})[
                    reordered_index
                ] = coefficient
        _, pivots, _ = sdm_irref(rows)
        return [column_order[index] for index in pivots]

    def determinant(matrix):
        matrix = [list(row) for row in matrix]
        value = fraction_field.one
        sign = 1
        for column in range(len(matrix)):
            pivot = next(
                (
                    row
                    for row in range(column, len(matrix))
                    if matrix[row][column]
                ),
                None,
            )
            if pivot is None:
                return fraction_field.zero
            if pivot != column:
                matrix[column], matrix[pivot] = (
                    matrix[pivot],
                    matrix[column],
                )
                sign *= -1
            pivot_value = matrix[column][column]
            value *= pivot_value
            for row in range(column + 1, len(matrix)):
                ratio = matrix[row][column] / pivot_value
                for index in range(column + 1, len(matrix)):
                    matrix[row][index] -= (
                        ratio * matrix[column][index]
                    )
        return value if sign == 1 else -value

    def chart(column_order, monomial_order):
        basis_indices = generic_column_pivots(column_order)
        if len(basis_indices) != 6:
            raise AssertionError(len(basis_indices))
        monomial_position = {
            monomial: index
            for index, monomial in enumerate(monomial_order)
        }
        transpose_rows = {
            row: {
                monomial_position[monomial]: coefficient
                for monomial, coefficient in variation_columns[
                    basis_index
                ].items()
            }
            for row, basis_index in enumerate(basis_indices)
        }
        _, pivot_rows, _ = sdm_irref(transpose_rows)
        selected_monomials = [
            monomial_order[index] for index in pivot_rows
        ]
        restricted_columns = [
            {
                row: coefficient
                for row, monomial in enumerate(selected_monomials)
                if (
                    coefficient := variation_columns[basis_index].get(
                        monomial,
                        fraction_field.zero,
                    )
                )
            }
            for basis_index in basis_indices
        ]
        restricted_rhs = {
            row: coefficient
            for row, monomial in enumerate(selected_monomials)
            if (
                coefficient := line_constant.get(
                    monomial,
                    fraction_field.zero,
                )
            )
        }
        solution, _, rank = solve_affine(
            restricted_columns,
            restricted_rhs,
            fraction_field,
        )
        if rank != 6:
            raise AssertionError(rank)
        residual = dict(line_constant)
        for index, coefficient in solution.items():
            residual = add(
                residual,
                variation_columns[basis_indices[index]],
                -coefficient,
            )
        nonzero_residuals = list(residual.values())
        if not nonzero_residuals:
            raise AssertionError("generic line unexpectedly extends")
        residual_gcd = None
        for value in nonzero_residuals:
            numerator = polynomial_ring.convert(value.numer)
            residual_gcd = (
                numerator
                if residual_gcd is None
                else polynomial_ring.gcd(residual_gcd, numerator)
            )
        matrix = [
            [
                variation_columns[basis_index].get(
                    monomial,
                    fraction_field.zero,
                )
                for basis_index in basis_indices
            ]
            for monomial in selected_monomials
        ]
        chart_determinant = determinant(matrix)
        determinant_numerator = polynomial_ring.convert(
            chart_determinant.numer
        )
        return (
            residual_gcd,
            determinant_numerator,
            basis_indices,
            selected_monomials,
        )

    column_orders = (
        list(range(len(variation_columns))),
        list(reversed(range(len(variation_columns)))),
        list(range(5, len(variation_columns)))
        + list(range(5)),
    )
    monomial_orders = (
        output_monomials,
        list(reversed(output_monomials)),
        output_monomials[::2] + output_monomials[1::2],
    )
    determinant_gcd = None
    all_charts_unit = True
    for chart_index, (column_order, monomial_order) in enumerate(
        zip(column_orders, monomial_orders),
        start=1,
    ):
        (
            residual_gcd,
            determinant_numerator,
            basis_indices,
            selected_monomials,
        ) = chart(column_order, monomial_order)
        determinant_gcd = (
            determinant_numerator
            if determinant_gcd is None
            else polynomial_ring.gcd(
                determinant_gcd,
                determinant_numerator,
            )
        )
        residual_degree = residual_gcd.degree()
        determinant_degree = determinant_numerator.degree()
        all_charts_unit &= residual_degree == 0
        print(
            f"H7_LINE_CHART={chart_index},"
            f"basis={basis_indices},"
            f"rows={selected_monomials},"
            f"residual_gcd_degree={residual_degree},"
            f"pivot_degree={determinant_degree}",
            flush=True,
        )
        if residual_degree:
            print(
                "H7_LINE_RESIDUAL_GCD="
                + str(polynomial_ring.to_sympy(residual_gcd)),
                flush=True,
            )
    determinant_gcd_degree = determinant_gcd.degree()
    print(
        "H7_LINE_PIVOT_GCD_DEGREE="
        + str(determinant_gcd_degree),
        flush=True,
    )
    if all_charts_unit and determinant_gcd_degree == 0:
        print(
            "H7_LINE_VERDICT=OBSTRUCTED_FOR_ALL_R",
            flush=True,
        )
    else:
        print(
            "H7_LINE_VERDICT=EXCEPTIONAL_PARAMETERS_REMAIN",
            flush=True,
        )


def reduced_component_parameters(field, a):
    """Return an affine parametrization of the reduced hbar^5 lift scheme.

    The exact radical computation leaves 27 free coordinates.  This helper
    records the six linear generators of the 16-variable nonlinear core
    together with the eight linear equations at the head of the original
    standard basis.
    """

    free_indices = (
        *range(10),
        11,
        13,
        15,
        17,
        19,
        21,
        23,
        26,
        28,
        30,
        32,
        35,
        36,
        37,
        38,
        39,
        40,
    )
    if len(free_indices) != 27:
        raise AssertionError(len(free_indices))

    def parameters(free_values):
        result = [field.zero] * 41
        for index, value in zip(free_indices, free_values, strict=True):
            result[index] = value

        # Six generators of the reduced 16-variable core.
        result[34] = -field(4) * result[35] / field(3)
        result[33] = -(
            field(203337542430) * a**2
            + field(487549466415) * a
            + field(295586113680)
        ) / field(4194304)
        result[31] = -field(10) * result[32] / field(3)
        result[29] = -(
            field(24) * result[30]
            + (field(872) * a + field(1856)) * result[40]
        ) / field(9)
        result[27] = -(
            field(18) * result[28]
            + (field(392) * a + field(752)) * result[38]
            + (field(848) * a + field(1496)) * result[39]
        ) / field(9)
        result[25] = -(
            field(36) * result[26]
            + (field(456) * a + field(672)) * result[36]
            + (field(896) * a + field(1184)) * result[37]
        ) / field(27)

        # Eight linear generators at the head of the full standard basis.
        result[24] = -(
            (field(108) * a + field(108)) * result[34]
            + (field(176) * a + field(168)) * result[35]
        ) / field(9)
        result[22] = -field(14) * result[23] / field(3)
        result[20] = -field(4) * result[21]
        result[18] = -(
            field(314572800) * result[19]
            + (
                field(8766095360) * a
                + field(19587399680)
            )
            * result[33]
            + field(186954838846362) * a**2
            + field(447513347427885) * a
            + field(270836273654640)
        ) / field(94371840)
        result[16] = -(
            field(72) * result[17]
            + (field(1308) * a + field(2784)) * result[31]
            + (field(3200) * a + field(6500)) * result[32]
        ) / field(27)
        result[14] = -(
            field(162) * result[15]
            + (field(1764) * a + field(3384)) * result[29]
            + (field(3816) * a + field(6732)) * result[30]
            + (
                field(121712) * a**2
                + field(481504) * a
                + field(469520)
            )
            * result[40]
        ) / field(81)
        result[12] = -(
            field(324) * result[13]
            + (field(2052) * a + field(3024)) * result[27]
            + (field(4032) * a + field(5328)) * result[28]
            + (
                field(74112) * a**2
                + field(250320) * a
                + field(200832)
            )
            * result[38]
            + (
                field(152768) * a**2
                + field(493024) * a
                + field(384320)
            )
            * result[39]
        ) / field(243)
        result[10] = -(
            field(54) * result[11]
            + (field(324) * a + field(324)) * result[25]
            + (field(504) * a + field(504)) * result[26]
            + (
                field(4512) * a**2
                + field(12048) * a
                + field(7536)
            )
            * result[36]
            + (
                field(8896) * a**2
                + field(21824) * a
                + field(12928)
            )
            * result[37]
        ) / field(81)
        return result

    zero = [field.zero] * len(free_indices)
    base = parameters(zero)
    directions = []
    for free_index in range(len(free_indices)):
        values = list(zero)
        values[free_index] = field.one
        point = parameters(values)
        directions.append(
            [
                value - base_value
                for value, base_value in zip(point, base, strict=True)
            ]
        )
    return free_indices, base, directions


def profile_seventh_order_reduced_component(
    field,
    a,
    S,
    T,
    third,
    constant,
    nonconstant,
    corrections,
    s4_monomials,
    t4_monomials,
) -> None:
    """Find the effective parameter count in the full reduced hbar^7 matrix."""

    from explore_degree_five_quantum_residue import column_rank

    free_indices, base_parameters, parameter_directions = (
        reduced_component_parameters(field, a)
    )
    variable_count = len(third.kernel)
    if variable_count != 41:
        raise AssertionError(variable_count)

    def fifth_defect(parameters):
        defect = dict(constant)
        for index, parameter in enumerate(parameters):
            defect = add(defect, nonconstant[2 * index], parameter)
            defect = add(
                defect,
                nonconstant[2 * index + 1],
                parameter**2,
            )
        cross_offset = 2 * variable_count
        for cross_index, (left, right) in enumerate(
            combinations(range(variable_count), 2)
        ):
            defect = add(
                defect,
                nonconstant[cross_offset + cross_index],
                parameters[left] * parameters[right],
            )
        return defect

    base_defect = fifth_defect(base_parameters)
    correction, correction_kernel, correction_rank = solve_affine(
        corrections,
        scale(base_defect, -field.one),
        field,
    )
    if correction_rank != 594 or len(correction_kernel) != 20:
        raise AssertionError((correction_rank, len(correction_kernel)))

    def split_correction(vector):
        split = len(s4_monomials)
        return (
            {
                s4_monomials[index]: coefficient
                for index, coefficient in vector.items()
                if index < split
            },
            {
                t4_monomials[index - split]: coefficient
                for index, coefficient in vector.items()
                if index >= split
            },
        )

    correction_kernel_pairs = [
        split_correction(vector) for vector in correction_kernel
    ]
    base_s2, base_t2 = third.base

    def lower_pair(parameters):
        s2 = dict(base_s2)
        t2 = dict(base_t2)
        for index, parameter in enumerate(parameters):
            direction_s, direction_t = third.kernel[index]
            s2 = add(s2, direction_s, parameter)
            t2 = add(t2, direction_t, parameter)
        return s2, t2

    lower_base = lower_pair(base_parameters)
    lower_directions = []
    for direction in parameter_directions:
        point_s, point_t = lower_pair(direction)
        lower_directions.append(
            (
                add(point_s, base_s2, -field.one),
                add(point_t, base_t2, -field.one),
            )
        )

    def variation_columns(lower):
        s2, t2 = lower
        result = []
        for direction_s4, direction_t4 in correction_kernel_pairs:
            variation = GENERIC_PROFILE.poisson(s2, direction_t4)
            variation = add(
                variation,
                GENERIC_PROFILE.poisson(direction_s4, t2),
            )
            variation = add(
                variation,
                GENERIC_PROFILE.pi_power(direction_s4, T, 3),
                field.one / field(24),
            )
            variation = add(
                variation,
                GENERIC_PROFILE.pi_power(S, direction_t4, 3),
                field.one / field(24),
            )
            result.append(variation)
        return result

    base_variations = variation_columns(lower_base)
    active_parameters = []
    direction_rank_changes = []
    for parameter_index, lower_direction in enumerate(lower_directions):
        direction_columns = variation_columns(
            (
                add(lower_base[0], lower_direction[0]),
                add(lower_base[1], lower_direction[1]),
            )
        )
        direction_columns = [
            add(column, base_column, -field.one)
            for column, base_column in zip(
                direction_columns,
                base_variations,
                strict=True,
            )
        ]
        if any(direction_columns):
            active_parameters.append(parameter_index)
            direction_rank_changes.append(column_rank(direction_columns))
    print(
        "H7_REDUCED_COMPONENT="
        "dimension27,nonlinear_core16,radical_generators6,"
        f"variation_active_parameters={len(active_parameters)},"
        f"base_variation_rank={column_rank(base_variations)}",
        flush=True,
    )
    for parameter_index, direction_rank in zip(
        active_parameters,
        direction_rank_changes,
        strict=True,
    ):
        print(
            "H7_ACTIVE_PARAMETER="
            f"v{parameter_index + 1},"
            f"source=u{free_indices[parameter_index] + 1},"
            f"matrix_direction_rank={direction_rank}",
            flush=True,
        )


def sample_seventh_order_reduced_component(
    sample_count,
    field,
    a,
    S,
    T,
    third,
    constant,
    nonconstant,
    corrections,
    s4_monomials,
    t4_monomials,
) -> None:
    """Sample the exact reduced component over a finite coefficient field."""

    import random

    from explore_degree_five_quantum_residue import column_rank

    if not hasattr(field, "mod"):
        raise ValueError("component sampling is intended for a finite field")
    free_indices, base_parameters, parameter_directions = (
        reduced_component_parameters(field, a)
    )
    variable_count = len(third.kernel)

    def parameters_at(values):
        result = list(base_parameters)
        for value, direction in zip(
            values,
            parameter_directions,
            strict=True,
        ):
            for index, coefficient in enumerate(direction):
                result[index] += value * coefficient
        return result

    def fifth_defect(parameters):
        defect = dict(constant)
        for index, parameter in enumerate(parameters):
            defect = add(defect, nonconstant[2 * index], parameter)
            defect = add(
                defect,
                nonconstant[2 * index + 1],
                parameter**2,
            )
        cross_offset = 2 * variable_count
        for cross_index, (left, right) in enumerate(
            combinations(range(variable_count), 2)
        ):
            defect = add(
                defect,
                nonconstant[cross_offset + cross_index],
                parameters[left] * parameters[right],
            )
        return defect

    def split_correction(vector):
        split = len(s4_monomials)
        return (
            {
                s4_monomials[index]: coefficient
                for index, coefficient in vector.items()
                if index < split
            },
            {
                t4_monomials[index - split]: coefficient
                for index, coefficient in vector.items()
                if index >= split
            },
        )

    random_source = random.Random(20260723)
    samples = [[field.zero] * 27]
    samples.extend(
        [
            field(random_source.randrange(0, field.mod))
            for _ in range(27)
        ]
        for _ in range(sample_count - 1)
    )
    rank_counts = {}
    solvable = []
    correction_kernel_pairs = None
    for sample_index, values in enumerate(samples):
        parameters = parameters_at(values)
        defect = fifth_defect(parameters)
        correction, correction_kernel, correction_rank = solve_affine(
            corrections,
            scale(defect, -field.one),
            field,
        )
        if correction_rank != 594 or len(correction_kernel) != 20:
            raise AssertionError((correction_rank, len(correction_kernel)))
        if correction_kernel_pairs is None:
            correction_kernel_pairs = [
                split_correction(vector) for vector in correction_kernel
            ]
        s2, t2 = third.base
        for index, parameter in enumerate(parameters):
            direction_s, direction_t = third.kernel[index]
            s2 = add(s2, direction_s, parameter)
            t2 = add(t2, direction_t, parameter)
        s4, t4 = split_correction(correction)
        seventh_constant = GENERIC_PROFILE.poisson(s2, t4)
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.poisson(s4, t2),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(s4, T, 3),
            field.one / field(24),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(s2, t2, 3),
            field.one / field(24),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(S, t4, 3),
            field.one / field(24),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(s2, T, 5),
            field.one / field(1920),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(S, t2, 5),
            field.one / field(1920),
        )
        seventh_constant = add(
            seventh_constant,
            GENERIC_PROFILE.pi_power(S, T, 7),
            field.one / field(322560),
        )
        seventh_variations = []
        for direction_s4, direction_t4 in correction_kernel_pairs:
            variation = GENERIC_PROFILE.poisson(s2, direction_t4)
            variation = add(
                variation,
                GENERIC_PROFILE.poisson(direction_s4, t2),
            )
            variation = add(
                variation,
                GENERIC_PROFILE.pi_power(direction_s4, T, 3),
                field.one / field(24),
            )
            variation = add(
                variation,
                GENERIC_PROFILE.pi_power(S, direction_t4, 3),
                field.one / field(24),
            )
            seventh_variations.append(variation)
        rank = column_rank(seventh_variations)
        augmented_rank = column_rank(
            seventh_variations + [seventh_constant]
        )
        rank_counts[(rank, augmented_rank)] = (
            rank_counts.get((rank, augmented_rank), 0) + 1
        )
        if rank == augmented_rank:
            solvable.append(sample_index)
    print(
        "H7_REDUCED_COMPONENT_SAMPLES="
        f"{sample_count},rank_counts={rank_counts},"
        f"solvable={solvable},free_sources="
        + ",".join(f"u{index + 1}" for index in free_indices),
        flush=True,
    )


def eliminate_seventh_order_reduced_component(
    field,
    a,
    S,
    T,
    third,
    constant,
    nonconstant,
    corrections,
    s4_monomials,
    t4_monomials,
) -> None:
    """Test the complete reduced hbar^7 consistency ideal."""

    free_indices, base_parameters, parameter_directions = (
        reduced_component_parameters(field, a)
    )
    variable_count = len(third.kernel)
    parameter_names = tuple(f"v{index}" for index in range(1, 28))
    parameter_field = field.frac_field(*parameter_names)
    parameter_ring = parameter_field.get_ring()
    ground_polynomial = parameter_ring.ring.ground_new
    variables = parameter_field.gens

    parameters = []
    for coordinate in range(variable_count):
        value = parameter_field(base_parameters[coordinate])
        for variable, direction in zip(
            variables,
            parameter_directions,
            strict=True,
        ):
            coefficient = direction[coordinate]
            if coefficient:
                value += parameter_field(coefficient) * variable
        parameters.append(value)

    def promote(poly):
        return {
            monomial: parameter_field(coefficient)
            for monomial, coefficient in poly.items()
            if coefficient
        }

    defect = promote(constant)
    for index, parameter in enumerate(parameters):
        defect = add(
            defect,
            promote(nonconstant[2 * index]),
            parameter,
        )
        defect = add(
            defect,
            promote(nonconstant[2 * index + 1]),
            parameter**2,
        )
    cross_offset = 2 * variable_count
    for cross_index, (left, right) in enumerate(
        combinations(range(variable_count), 2)
    ):
        defect = add(
            defect,
            promote(nonconstant[cross_offset + cross_index]),
            parameters[left] * parameters[right],
        )

    promoted_corrections = [promote(column) for column in corrections]
    correction, correction_rank = affine_particular(
        promoted_corrections,
        scale(defect, -parameter_field.one),
        parameter_field,
    )
    if correction_rank != 594:
        raise AssertionError(correction_rank)

    def as_polynomial(value):
        value = parameter_field.convert(value)
        numerator = parameter_ring.convert(value.numer)
        denominator = parameter_ring.convert(value.denom)
        denominator_terms = denominator.to_dict()
        zero_exponent = (0,) * len(parameter_names)
        if set(denominator_terms) != {zero_exponent}:
            raise AssertionError("unexpected parameter denominator")
        inverse = field.one / denominator_terms[zero_exponent]
        return numerator * ground_polynomial(inverse)

    polynomial_parameters = [
        as_polynomial(parameter) for parameter in parameters
    ]

    def promote_polynomial(poly):
        return {
            monomial: ground_polynomial(coefficient)
            for monomial, coefficient in poly.items()
            if coefficient
        }

    correction = {
        index: as_polynomial(coefficient)
        for index, coefficient in correction.items()
    }

    base_defect = dict(constant)
    for index, parameter in enumerate(base_parameters):
        base_defect = add(
            base_defect,
            nonconstant[2 * index],
            parameter,
        )
        base_defect = add(
            base_defect,
            nonconstant[2 * index + 1],
            parameter**2,
        )
    for cross_index, (left, right) in enumerate(
        combinations(range(variable_count), 2)
    ):
        base_defect = add(
            base_defect,
            nonconstant[cross_offset + cross_index],
            base_parameters[left] * base_parameters[right],
        )
    _, correction_kernel, base_rank = solve_affine(
        corrections,
        scale(base_defect, -field.one),
        field,
    )
    if base_rank != 594 or len(correction_kernel) != 20:
        raise AssertionError((base_rank, len(correction_kernel)))

    def split_correction(vector, converter):
        split = len(s4_monomials)
        return (
            {
                s4_monomials[index]: converter(coefficient)
                for index, coefficient in vector.items()
                if index < split
            },
            {
                t4_monomials[index - split]: converter(coefficient)
                for index, coefficient in vector.items()
                if index >= split
            },
        )

    s4, t4 = split_correction(correction, parameter_ring.convert)
    s2 = promote_polynomial(third.base[0])
    t2 = promote_polynomial(third.base[1])
    for index, parameter in enumerate(polynomial_parameters):
        direction_s, direction_t = third.kernel[index]
        s2 = add(s2, promote_polynomial(direction_s), parameter)
        t2 = add(t2, promote_polynomial(direction_t), parameter)

    seventh_constant = GENERIC_PROFILE.poisson(s2, t4)
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.poisson(s4, t2),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(s4, promote_polynomial(T), 3),
        ground_polynomial(field.one / field(24)),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(s2, t2, 3),
        ground_polynomial(field.one / field(24)),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(promote_polynomial(S), t4, 3),
        ground_polynomial(field.one / field(24)),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(s2, promote_polynomial(T), 5),
        ground_polynomial(field.one / field(1920)),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(promote_polynomial(S), t2, 5),
        ground_polynomial(field.one / field(1920)),
    )
    seventh_constant = add(
        seventh_constant,
        GENERIC_PROFILE.pi_power(
            promote_polynomial(S),
            promote_polynomial(T),
            7,
        ),
        ground_polynomial(field.one / field(322560)),
    )

    seventh_variations = []
    for kernel_vector in correction_kernel:
        direction_s4, direction_t4 = split_correction(
            kernel_vector,
            ground_polynomial,
        )
        variation = GENERIC_PROFILE.poisson(s2, direction_t4)
        variation = add(
            variation,
            GENERIC_PROFILE.poisson(direction_s4, t2),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(
                direction_s4,
                promote_polynomial(T),
                3,
            ),
            ground_polynomial(field.one / field(24)),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(
                promote_polynomial(S),
                direction_t4,
                3,
            ),
            ground_polynomial(field.one / field(24)),
        )
        seventh_variations.append(variation)

    pivot_columns = (9, 11, 13, 15, 17, 19)
    pivot_rows = tuple((degree, 0, 0) for degree in range(8, 14))
    restricted_columns = [
        {
            row: parameter_field(
                seventh_variations[column].get(
                    monomial,
                    parameter_ring.zero,
                )
            )
            for row, monomial in enumerate(pivot_rows)
            if seventh_variations[column].get(
                monomial,
                parameter_ring.zero,
            )
        }
        for column in pivot_columns
    ]
    restricted_rhs = {
        row: -parameter_field(
            seventh_constant.get(monomial, parameter_ring.zero)
        )
        for row, monomial in enumerate(pivot_rows)
        if seventh_constant.get(monomial, parameter_ring.zero)
    }
    pivot_solution, pivot_kernel, pivot_rank = solve_affine(
        restricted_columns,
        restricted_rhs,
        parameter_field,
    )
    if pivot_rank != 6 or pivot_kernel:
        raise AssertionError((pivot_rank, len(pivot_kernel)))
    pivot_solution = {
        index: as_polynomial(coefficient)
        for index, coefficient in pivot_solution.items()
    }
    pivot_matrix = [
        [
            restricted_columns[column].get(
                row,
                parameter_field.zero,
            )
            for column in range(6)
        ]
        for row in range(6)
    ]
    pivot_determinant = parameter_field.one
    determinant_sign = 1
    for column in range(6):
        pivot = next(
            (
                row
                for row in range(column, 6)
                if pivot_matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            pivot_determinant = parameter_field.zero
            break
        if pivot != column:
            pivot_matrix[column], pivot_matrix[pivot] = (
                pivot_matrix[pivot],
                pivot_matrix[column],
            )
            determinant_sign *= -1
        pivot_value = pivot_matrix[column][column]
        pivot_determinant *= pivot_value
        for row in range(column + 1, 6):
            ratio = pivot_matrix[row][column] / pivot_value
            for index in range(column + 1, 6):
                pivot_matrix[row][index] -= (
                    ratio * pivot_matrix[column][index]
                )
    if determinant_sign == -1:
        pivot_determinant = -pivot_determinant
    def total_degree(polynomial):
        return max(
            (sum(exponents) for exponents in polynomial.to_dict()),
            default=-1,
        )

    constant_pivot = bool(pivot_determinant) and (
        total_degree(parameter_ring.convert(pivot_determinant.numer)) == 0
        and total_degree(parameter_ring.convert(pivot_determinant.denom))
        == 0
    )
    if not constant_pivot:
        raise AssertionError("the selected hbar^7 pivot is not constant")

    for column_index, column in enumerate(seventh_variations):
        if column_index in pivot_columns:
            continue
        column_rhs = {
            row: parameter_field(
                column.get(monomial, parameter_ring.zero)
            )
            for row, monomial in enumerate(pivot_rows)
            if column.get(monomial, parameter_ring.zero)
        }
        span_solution, span_kernel, span_rank = solve_affine(
            restricted_columns,
            column_rhs,
            parameter_field,
        )
        if span_rank != 6 or span_kernel:
            raise AssertionError((column_index, span_rank, len(span_kernel)))
        span_solution = {
            index: as_polynomial(coefficient)
            for index, coefficient in span_solution.items()
        }
        column_residual = dict(column)
        for local_index, coefficient in span_solution.items():
            column_residual = add(
                column_residual,
                seventh_variations[pivot_columns[local_index]],
                -coefficient,
            )
        if column_residual:
            raise AssertionError(
                f"hbar^7 column {column_index} leaves the global rank-six span"
            )

    residual = dict(seventh_constant)
    for local_index, coefficient in pivot_solution.items():
        residual = add(
            residual,
            seventh_variations[pivot_columns[local_index]],
            coefficient,
        )
    residual_numerators = [
        parameter_ring.convert(value)
        for value in residual.values()
        if value
    ]
    if not residual_numerators:
        raise AssertionError("generic reduced component unexpectedly extends")

    active_variables = sorted(
        {
            str(generator)
            for polynomial in residual_numerators
            for generator in parameter_ring.to_sympy(
                polynomial
            ).free_symbols
        }
    )
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for component elimination")
    if hasattr(field, "mod"):
        ring_declaration = (
            f"ring r={int(field.mod)},"
            f"({','.join(active_variables)}),dp;"
        )

        def coefficient_text(value):
            return str(int(value) % int(field.mod))

        field_label = f"GF({int(field.mod)})"
    else:
        ring_declaration = (
            f"ring r=(0,a),({','.join(active_variables)}),dp;\n"
            "minpoly=94*a^3+335*a^2+400*a+160;"
        )

        def coefficient_text(value):
            terms = []
            coefficients = value.to_list()
            for index, coefficient in enumerate(coefficients):
                if not coefficient:
                    continue
                exponent = len(coefficients) - index - 1
                scalar = (
                    f"({coefficient.numerator}/"
                    f"{coefficient.denominator})"
                )
                if exponent == 0:
                    term = scalar
                elif exponent == 1:
                    term = f"{scalar}*a"
                else:
                    term = f"{scalar}*a^{exponent}"
                terms.append(term)
            return "+".join(terms) if terms else "0"

        field_label = "Q(a)"

    def serialize_polynomial(polynomial):
        terms = []
        for exponents, coefficient in polynomial.to_dict().items():
            factors = [f"({coefficient_text(coefficient)})"]
            for name, exponent in zip(
                parameter_names,
                exponents,
                strict=True,
            ):
                if exponent == 1:
                    factors.append(name)
                elif exponent:
                    factors.append(f"{name}^{exponent}")
            terms.append("*".join(factors))
        return "+".join(terms) if terms else "0"

    polynomial_text = [
        serialize_polynomial(polynomial)
        for polynomial in residual_numerators
    ]
    program = f"""\
{ring_declaration}
option(redSB);
ideal I={','.join(polynomial_text)};
ideal G=std(I);
print("H7_CONSISTENCY_FIELD={field_label}");
print("H7_CONSISTENCY_GENERATORS="+string(size(I)));
print("H7_CONSISTENCY_ACTIVE_PARAMETERS={len(active_variables)}/27");
print("H7_CONSISTENCY_GB_SIZE="+string(size(G)));
print("H7_CONSISTENCY_UNIT="+string(reduce(1,G)==0));
if(reduce(1,G)!=0)
{{
  print("H7_CONSISTENCY_DIMENSION="+string(dim(G)));
}}
quit;
"""
    with tempfile.TemporaryDirectory(
        prefix="degree-five-cubic-seventh-component-",
    ) as directory:
        script = Path(directory) / "eliminate.sing"
        script.write_text(program)
        result = subprocess.run(
            [singular, "-q", str(script)],
            check=True,
            capture_output=True,
            text=True,
        )
    print(
        "H7_GLOBAL_PIVOT="
        f"columns={pivot_columns},rows={pivot_rows},constant=1",
        flush=True,
    )
    print("H7_GLOBAL_VARIATION_RANK=6", flush=True)
    print(
        "H7_CONSISTENCY_ACTIVE_SOURCES="
        + ",".join(
            f"u{free_indices[int(name[1:]) - 1] + 1}"
            for name in active_variables
        ),
        flush=True,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")


def run_singular(
    equations,
    variable_count: int,
    prime: int | None,
    field,
    print_basis: bool = False,
    radical_analysis: bool = False,
    print_radical_basis: bool = False,
) -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for --groebner")
    active_variables = set()
    variable_pairs = tuple(combinations(range(variable_count), 2))
    for _, linear, diagonal, cross in equations:
        active_variables.update(
            index for index, value in enumerate(linear) if value
        )
        active_variables.update(
            index for index, value in enumerate(diagonal) if value
        )
        for (left, right), value in zip(variable_pairs, cross):
            if value:
                active_variables.update((left, right))
    inactive_count = variable_count - len(active_variables)
    ring_variables = (
        sorted(active_variables)
        if radical_analysis
        else range(variable_count)
    )
    variables = ",".join(f"u{index + 1}" for index in ring_variables)
    if prime is None:
        ring_declaration = (
            f"ring r=(0,a),({variables}),dp;\n"
            "minpoly=94*a^3+335*a^2+400*a+160;"
        )

        def coefficient_text(value):
            # Singular can parse ``(large_integer/denominator)*a^k`` but
            # may misclassify ``large_integer*a^k/denominator`` as a
            # number-to-number power.  Emit the algebraic-field vector
            # directly, with every rational scalar parenthesized.
            terms = []
            coefficients = value.to_list()
            for index, coefficient in enumerate(coefficients):
                if not coefficient:
                    continue
                exponent = len(coefficients) - index - 1
                scalar = (
                    f"({coefficient.numerator}/"
                    f"{coefficient.denominator})"
                )
                if exponent == 0:
                    term = scalar
                elif exponent == 1:
                    term = f"{scalar}*a"
                else:
                    term = f"{scalar}*a^{exponent}"
                terms.append(term)
            return "+".join(terms) if terms else "0"
    else:
        ring_declaration = f"ring r={prime},({variables}),dp;"

        def coefficient_text(value):
            return str(modular_int(value, prime))

    generators = ",\n".join(
        singular_polynomial(equation, variable_count, coefficient_text)
        for equation in equations
    )
    radical_free_dimension = inactive_count + 7
    radical_variables = ",".join(
        f"u{index}" for index in range(26, 42)
    )
    if prime is None:
        radical_ring = (
            f"ring q=(0,a),({radical_variables}),dp;\n"
            "minpoly=94*a^3+335*a^2+400*a+160;"
        )
    else:
        radical_ring = f"ring q={prime},({radical_variables}),dp;"
    radical_generators = ",".join(
        f"G[{index}]" for index in range(9, 20)
    )
    radical_program = """
LIB "primdec.lib";
if(size(G)!=19)
{
  ERROR("unexpected standard basis size in compact radical reduction");
}
ideal J=%s;
%s
ideal Jq=imap(r,J);
ideal R=radical(Jq);
ideal GR=std(R);
list MP=minAssGTZ(R);
print("RADICAL_GB_SIZE="+string(size(GR)));
print("RADICAL_DIMENSION="+string(dim(GR)));
print("RADICAL_AMBIENT_DIMENSION="+string(dim(GR)+%d));
print("MINIMAL_PRIME_COUNT="+string(size(MP)));
%s
int component_index;
for(component_index=1;component_index<=size(MP);component_index++)
{
  ideal PC=std(MP[component_index]);
  print(
    "MINIMAL_PRIME_"+string(component_index)
    +"_SIZE="+string(size(PC))
    +",DIMENSION="+string(dim(PC))
    +",AMBIENT_DIMENSION="+string(dim(PC)+%d)
  );
}
""" % (
        radical_generators,
        radical_ring,
        radical_free_dimension,
        (
            'print("RADICAL_STANDARD_BASIS="); print(GR);'
            if print_radical_basis
            else ""
        ),
        radical_free_dimension,
    ) if radical_analysis else ""
    program = f"""\
{ring_declaration}
option(redSB);
ideal I=
{generators};
ideal G=std(I);
print("ACTIVE_VARIABLES={len(active_variables)}/{variable_count}");
print("IDEAL_GENERATORS="+string(size(I)));
print("GB_SIZE="+string(size(G)));
print("UNIT_REMAINDER="+string(reduce(1,G)));
print("DIMENSION="+string(dim(G)));
if(dim(G)==0)
{{
  print("VECTOR_DIMENSION="+string(vdim(G)));
}}
{'''
print("STANDARD_BASIS=");
print(G);
''' if print_basis else ''}
{radical_program}
quit;
"""
    with tempfile.TemporaryDirectory(
        prefix="degree-five-cubic-fifth-order-",
    ) as directory:
        script = Path(directory) / "solve.sing"
        script.write_text(program)
        result = subprocess.run(
            [singular, "-q", str(script)],
            check=True,
            capture_output=True,
            text=True,
        )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument("--a", type=int, default=DEFAULT_A)
    parser.add_argument("--tau", type=int, default=DEFAULT_TAU)
    parser.add_argument(
        "--exact-cubic",
        action="store_true",
        help="work over Q(a) on the reconstructed cubic locus",
    )
    parser.add_argument(
        "--groebner",
        action="store_true",
        help="compute a Singular standard basis of the projected quadrics",
    )
    parser.add_argument(
        "--print-basis",
        action="store_true",
        help="print the standard basis (implies --groebner)",
    )
    parser.add_argument(
        "--seventh-order",
        action="store_true",
        help="test whether the explicit cubic branch extends through hbar^7",
    )
    parser.add_argument(
        "--u36",
        default="0",
        help="rational coordinate on an explicit cubic hbar^5 solution line",
    )
    parser.add_argument(
        "--seventh-line",
        action="store_true",
        help="analyze hbar^7 over the full explicit u36 solution line",
    )
    parser.add_argument(
        "--radical-analysis",
        action="store_true",
        help="compute the radical and minimal-prime dimensions in Singular",
    )
    parser.add_argument(
        "--print-radical-basis",
        action="store_true",
        help="print the six-generator reduced nonlinear core",
    )
    parser.add_argument(
        "--seventh-component-profile",
        action="store_true",
        help="profile the hbar^7 matrix over the full reduced hbar^5 component",
    )
    parser.add_argument(
        "--seventh-component-samples",
        type=int,
        default=0,
        help="sample this many reduced-component points over the finite field",
    )
    parser.add_argument(
        "--seventh-component-elimination",
        action="store_true",
        help="compute the global reduced hbar^7 consistency ideal modulo prime",
    )
    args = parser.parse_args()
    prime = None if args.exact_cubic else args.prime
    (
        equations,
        variable_count,
        field,
        S,
        T,
        third,
        constant,
        nonconstant,
        corrections,
        s4_monomials,
        t4_monomials,
    ) = construct_projected_system(
        prime,
        args.a,
        args.tau,
    )
    if args.exact_cubic:
        numerator, denominator = args.u36.split("/", 1) if "/" in args.u36 else (
            args.u36,
            "1",
        )
        branch_u36 = field(int(numerator)) / field(int(denominator))
        verify_explicit_cubic_branch(
            equations,
            variable_count,
            field,
            S,
            T,
            third,
            constant,
            nonconstant,
            corrections,
            s4_monomials,
            t4_monomials,
            seventh_order=args.seventh_order,
            branch_u36=branch_u36,
        )
        if args.seventh_line:
            analyze_seventh_order_solution_line(
                field,
                S,
                T,
                third,
                constant,
                nonconstant,
                corrections,
                s4_monomials,
                t4_monomials,
            )
    if args.seventh_component_profile:
        component_a = field.unit if args.exact_cubic else field(args.a)
        profile_seventh_order_reduced_component(
            field,
            component_a,
            S,
            T,
            third,
            constant,
            nonconstant,
            corrections,
            s4_monomials,
            t4_monomials,
        )
    if args.seventh_component_samples:
        if args.exact_cubic:
            parser.error("--seventh-component-samples requires a finite field")
        sample_seventh_order_reduced_component(
            args.seventh_component_samples,
            field,
            field(args.a),
            S,
            T,
            third,
            constant,
            nonconstant,
            corrections,
            s4_monomials,
            t4_monomials,
        )
    if args.seventh_component_elimination:
        component_a = field.unit if args.exact_cubic else field(args.a)
        eliminate_seventh_order_reduced_component(
            field,
            component_a,
            S,
            T,
            third,
            constant,
            nonconstant,
            corrections,
            s4_monomials,
            t4_monomials,
        )
    if (
        args.groebner
        or args.print_basis
        or args.radical_analysis
        or args.print_radical_basis
    ):
        run_singular(
            equations,
            variable_count,
            prime,
            field,
            print_basis=args.print_basis,
            radical_analysis=(
                args.radical_analysis or args.print_radical_basis
            ),
            print_radical_basis=args.print_radical_basis,
        )


if __name__ == "__main__":
    main()

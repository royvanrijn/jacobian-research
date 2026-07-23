#!/usr/bin/env python3
"""Explore the five-dimensional residual odd first-correction space.

At u=0 every vector in this space satisfies the relaxed third-order
equation.  A nonzero branch scale additionally requires the cubic
second-order particular contribution to lie in the span of the 42 linear
kernel couplings.  After quotienting by the fixed third differential this is
a small rank test on projective four-space.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import (
    sdm_irref,
    sdm_nullspace_from_rref,
)

from explore_degree_five_a2_subprincipal import (
    add,
    filtered_monomials,
    poisson,
    scale,
)
from explore_degree_five_quantum_residue import column_rank
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    essential_problem,
    linear_combination,
    solve_many_particular,
    split_correction,
)
from explore_rank_two_odd_support_three_curves import (
    project_to_third_cokernel,
)


RESIDUAL_VECTORS = (
    {0: QQ(9) / QQ(5), 9: QQ.one},
    {0: QQ(26) / QQ(15), 10: QQ.one},
    {0: -QQ(37696) / QQ(1575), 23: QQ.one},
    {0: -QQ(142216) / QQ(4725), 24: QQ.one},
    {0: QQ(205276828) / QQ(165375), 35: QQ.one},
)


def add_sparse(left, right, coefficient):
    result = dict(left)
    for index, value in right.items():
        result[index] = result.get(index, QQ.zero) + coefficient * value
        if not result[index]:
            result.pop(index)
    return result


def residual_scale_problem():
    """Return projected linear coupling columns and cubic u-column terms."""

    S, T, _, pairs, _, _ = essential_problem()
    projected_couplings, _, project = project_to_third_cokernel()
    residual_pairs = [
        linear_combination(pairs, vector, QQ)
        for vector in RESIDUAL_VECTORS
    ]
    residual_couplings = []
    for vector in RESIDUAL_VECTORS:
        columns = []
        for kernel_index in range(42):
            column = {}
            for essential_index, coefficient in vector.items():
                column = add_sparse(
                    column,
                    projected_couplings[essential_index][kernel_index],
                    coefficient,
                )
            columns.append(column)
        residual_couplings.append(columns)

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    second_columns = [
        poisson({monomial: QQ.one}, T)
        for monomial in s2_monomials
    ]
    second_columns += [
        poisson(S, {monomial: QQ.one})
        for monomial in t2_monomials
    ]
    quadratic_indices = []
    quadratic_right_sides = []
    for left in range(5):
        for right in range(left, 5):
            left_s, left_t = residual_pairs[left]
            right_s, right_t = residual_pairs[right]
            if left == right:
                bracket = poisson(left_s, left_t)
            else:
                bracket = add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            quadratic_indices.append((left, right))
            quadratic_right_sides.append(scale(bracket, -QQ.one))
    quadratic_vectors = solve_many_particular(
        second_columns,
        quadratic_right_sides,
        QQ,
    )
    quadratic_pairs = [
        split_correction(
            vector,
            s2_monomials,
            t2_monomials,
        )
        for vector in quadratic_vectors
    ]

    cubic_terms = {}
    for first_index, first_pair in enumerate(residual_pairs):
        for quadratic_index, quadratic_pair in zip(
            quadratic_indices,
            quadratic_pairs,
        ):
            cubic_terms[(first_index, *quadratic_index)] = project(
                coupling(first_pair, quadratic_pair)
            )
    return residual_couplings, cubic_terms


def _column_space_coordinates(span_columns, extra_vectors):
    """Choose an exact basis of a sparse column space and project to it."""

    ambient = sorted(
        set().union(
            *(set(column) for column in span_columns),
            *(set(vector) for vector in extra_vectors),
        )
    )
    rows = {
        row: {
            column_index: value
            for column_index, column in enumerate(span_columns)
            if (value := column.get(coordinate, QQ.zero))
        }
        for row, coordinate in enumerate(ambient)
    }
    _, pivot_columns, _ = sdm_irref(rows)
    basis = [span_columns[index] for index in pivot_columns]

    ambient_index = {
        coordinate: index for index, coordinate in enumerate(ambient)
    }
    transpose = {
        column_index: {
            ambient_index[coordinate]: value
            for coordinate, value in column.items()
        }
        for column_index, column in enumerate(basis)
    }
    _, pivot_rows, _ = sdm_irref(transpose)
    row_coordinates = [ambient[index] for index in pivot_rows]
    matrix = sp.Matrix(
        [
            [
                sp.Rational(
                    column.get(coordinate, QQ.zero).numerator,
                    column.get(coordinate, QQ.zero).denominator,
                )
                for column in basis
            ]
            for coordinate in row_coordinates
        ]
    )
    inverse = matrix.inv()

    def coordinates(vector):
        values = sp.Matrix(
            [
                sp.Rational(
                    vector.get(coordinate, QQ.zero).numerator,
                    vector.get(coordinate, QQ.zero).denominator,
                )
                for coordinate in row_coordinates
            ]
        )
        return tuple(QQ.convert(value) for value in inverse * values)

    return basis, coordinates


def exact_scale_locus():
    """Prove the nonzero-scale locus is a hyperplane union a quadric."""

    couplings, cubic_terms = residual_scale_problem()
    all_couplings = [
        column for row in couplings for column in row
    ]
    all_vectors = all_couplings + list(cubic_terms.values())

    # First project the cubic u-column modulo the union of every linear
    # coupling image.  This union has dimension 14.
    ambient = sorted(set().union(*(set(vector) for vector in all_vectors)))
    ambient_index = {
        coordinate: index for index, coordinate in enumerate(ambient)
    }
    transpose = {
        column_index: {
            ambient_index[coordinate]: value
            for coordinate, value in column.items()
        }
        for column_index, column in enumerate(all_couplings)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose)
    union_cokernel, _ = sdm_nullspace_from_rref(
        reduced,
        QQ.one,
        len(ambient),
        pivots,
        nonzero,
    )
    assert len(pivots) == 14

    incidence = {}
    for functional_index, functional in enumerate(union_cokernel):
        for coordinate, value in functional.items():
            incidence.setdefault(
                ambient[coordinate],
                [],
            ).append((functional_index, value))

    def project_union_cokernel(vector):
        result = {}
        for coordinate, coefficient in vector.items():
            for functional_index, weight in incidence.get(coordinate, ()):
                result[functional_index] = (
                    result.get(functional_index, QQ.zero)
                    + coefficient * weight
                )
        return {
            coordinate: value
            for coordinate, value in result.items()
            if value
        }

    z = sp.symbols("z0:5")
    obstruction_polys = {}
    for (first, left, right), vector in cubic_terms.items():
        projected = project_union_cokernel(vector)
        for coordinate, coefficient in projected.items():
            obstruction_polys[coordinate] = (
                obstruction_polys.get(coordinate, sp.S.Zero)
                + sp.Rational(
                    coefficient.numerator,
                    coefficient.denominator,
                )
                * z[first]
                * z[left]
                * z[right]
            )
    obstruction_polys = [
        sp.Poly(sp.expand(expression), *z, domain=sp.QQ)
        for expression in obstruction_polys.values()
        if expression
    ]
    assert len(obstruction_polys) == 49
    common = obstruction_polys[0]
    for polynomial in obstruction_polys[1:]:
        common = sp.gcd(common, polynomial)
    assert common.total_degree() == 3
    assert all(
        sp.quo(polynomial, common).total_degree() == 0
        for polynomial in obstruction_polys
    )

    factorization = sp.factor_list(common)[1]
    assert sorted(
        factor.total_degree() for factor, _ in factorization
    ) == [1, 2]
    linear = next(
        factor.monic()
        for factor, _ in factorization
        if factor.total_degree() == 1
    )
    quadric = next(
        factor.primitive()[1]
        for factor, _ in factorization
        if factor.total_degree() == 2
    )
    expected_linear = sp.Poly(
        z[2] + 2 * z[3] - sp.Rational(9838, 105) * z[4],
        *z,
        domain=sp.QQ,
    ).monic()
    assert linear == expected_linear

    # The whole coupling map depends on only two linear forms a and b.
    basis, union_coordinates = _column_space_coordinates(
        all_couplings,
        list(cubic_terms.values()),
    )
    assert len(basis) == 14
    coupling_coordinates = [
        [
            union_coordinates(couplings[index][kernel_index])
            for kernel_index in range(42)
        ]
        for index in range(5)
    ]
    cubic_coordinates = {
        index: union_coordinates(vector)
        for index, vector in cubic_terms.items()
    }
    a_form = (
        z[0]
        + sp.Rational(4, 3) * z[1]
        + sp.Rational(2648, 189) * z[3]
        - sp.Rational(580504, 735) * z[4]
    )
    b_form = linear.as_expr()
    expected_relations = (
        (QQ.one, QQ.zero),
        (QQ(4) / QQ(3), QQ.zero),
        (QQ.zero, QQ.one),
        (QQ(2648) / QQ(189), QQ(2)),
        (-QQ(580504) / QQ(735), -QQ(9838) / QQ(105)),
    )
    for index, (a_coefficient, b_coefficient) in enumerate(
        expected_relations
    ):
        for kernel_index in range(42):
            expected = tuple(
                a_coefficient
                * coupling_coordinates[0][kernel_index][row]
                + b_coefficient
                * coupling_coordinates[2][kernel_index][row]
                for row in range(14)
            )
            assert coupling_coordinates[index][kernel_index] == expected

    def left_nullspace(columns, field):
        rows = {
            index: {
                coordinate: value
                for coordinate, value in enumerate(column)
                if value
            }
            for index, column in enumerate(columns)
            if any(column)
        }
        reduced_rows, column_pivots, nonzero_columns = sdm_irref(rows)
        nullspace, _ = sdm_nullspace_from_rref(
            reduced_rows,
            field.one,
            14,
            column_pivots,
            nonzero_columns,
        )
        return len(column_pivots), nullspace

    # On b != 0 use t=a/b.  The three remaining quotient functionals of
    # the 14-dimensional union vanish identically on the cubic column.
    parameter = sp.Symbol("t")
    function_field = QQ.frac_field(parameter)
    parameter_value = function_field.from_sympy(parameter)
    pencil_columns = []
    for kernel_index in range(42):
        pencil_columns.append(
            tuple(
                function_field.convert(
                    coupling_coordinates[0][kernel_index][row]
                )
                * parameter_value
                + function_field.convert(
                    coupling_coordinates[2][kernel_index][row]
                )
                for row in range(14)
            )
        )
    pencil_rank, pencil_nullspace = left_nullspace(
        pencil_columns,
        function_field,
    )
    assert (pencil_rank, len(pencil_nullspace)) == (11, 3)
    for functional in pencil_nullspace:
        expression = sp.S.Zero
        for (
            first,
            left_index,
            right_index,
        ), coordinates in cubic_coordinates.items():
            coefficient = sum(
                (
                    functional.get(row, function_field.zero)
                    * function_field.convert(coordinates[row])
                    for row in range(14)
                ),
                function_field.zero,
            )
            expression += (
                function_field.to_sympy(coefficient)
                * z[first]
                * z[left_index]
                * z[right_index]
            )
        assert sp.cancel(
            expression.subs(parameter, a_form / b_form)
        ) == 0

    # On b=0 the pencil is the fixed a-component.  Its three quotient
    # functionals also kill the cubic column.  At a=b=0 the cubic column
    # itself vanishes, covering the rank-zero kernel plane.
    a_columns = [
        coupling_coordinates[0][kernel_index]
        for kernel_index in range(42)
    ]
    a_rank, a_nullspace = left_nullspace(a_columns, QQ)
    assert (a_rank, len(a_nullspace)) == (11, 3)
    z4_on_hyperplane = sp.Rational(105, 9838) * (
        z[2] + 2 * z[3]
    )
    cubic_coordinate_polys = []
    for row in range(14):
        cubic_coordinate_polys.append(
            sp.expand(
                sum(
                    sp.Rational(
                        coordinates[row].numerator,
                        coordinates[row].denominator,
                    )
                    * z[first]
                    * z[left_index]
                    * z[right_index]
                    for (
                        first,
                        left_index,
                        right_index,
                    ), coordinates in cubic_coordinates.items()
                )
            )
        )
    for functional in a_nullspace:
        expression = sum(
            sp.Rational(value.numerator, value.denominator)
            * cubic_coordinate_polys[row]
            for row, value in functional.items()
        )
        assert sp.expand(
            expression.subs(z[4], z4_on_hyperplane)
        ) == 0

    a_on_hyperplane = sp.expand(
        a_form.subs(z[4], z4_on_hyperplane)
    )
    z1_on_kernel = sp.solve(a_on_hyperplane, z[1])[0]
    for expression in cubic_coordinate_polys:
        assert sp.expand(
            expression.subs(
                {
                    z[4]: z4_on_hyperplane,
                    z[1]: z1_on_kernel,
                }
            )
        ) == 0

    return linear.as_expr(), sp.factor(quadric.as_expr())


def reduce_value(value, field):
    return (
        field(int(value.numerator))
        / field(int(value.denominator))
    )


def reduce_vector(vector, field):
    result = {}
    for coordinate, value in vector.items():
        reduced = reduce_value(value, field)
        if reduced:
            result[coordinate] = reduced
    return result


def projective_points(field, dimension):
    prime = field.mod
    for pivot in range(dimension):
        for tail in product(range(prime), repeat=dimension - pivot - 1):
            yield (
                (field.zero,) * pivot
                + (field.one,)
                + tuple(field(value) for value in tail)
            )


def combine_vectors(vectors, coefficients, field):
    result = {}
    for vector, coefficient in zip(vectors, coefficients):
        if not coefficient:
            continue
        for coordinate, value in vector.items():
            result[coordinate] = (
                result.get(coordinate, field.zero)
                + coefficient * value
            )
            if not result[coordinate]:
                result.pop(coordinate)
    return result


def finite_field_scale_locus(prime):
    field = GF(prime)
    residual_couplings, cubic_terms = residual_scale_problem()
    modular_couplings = [
        [reduce_vector(column, field) for column in row]
        for row in residual_couplings
    ]
    modular_cubics = {
        index: reduce_vector(column, field)
        for index, column in cubic_terms.items()
    }

    survivors = []
    rank_profile = Counter()
    for coordinates in projective_points(field, 5):
        coupling_columns = [
            combine_vectors(
                [
                    modular_couplings[first_index][kernel_index]
                    for first_index in range(5)
                ],
                coordinates,
                field,
            )
            for kernel_index in range(42)
        ]
        cubic_column = {}
        for (first, left, right), column in modular_cubics.items():
            coefficient = (
                coordinates[first]
                * coordinates[left]
                * coordinates[right]
            )
            if not coefficient:
                continue
            cubic_column = combine_vectors(
                (cubic_column, column),
                (field.one, coefficient),
                field,
            )
        rank = column_rank(coupling_columns)
        augmented_rank = column_rank(
            coupling_columns + [cubic_column]
        )
        rank_profile[(rank, augmented_rank)] += 1
        if rank == augmented_rank:
            survivors.append(coordinates)
    return rank_profile, survivors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int)
    arguments = parser.parse_args()
    linear, quadric = exact_scale_locus()
    print(
        "PASS: the exact nonzero-scale locus in residual P4 is the union "
        "of one hyperplane and one quadric"
    )
    print(f"hyperplane={linear}")
    print(f"quadric={quadric}")
    if arguments.prime is not None:
        profile, survivors = finite_field_scale_locus(arguments.prime)
        print(f"prime={arguments.prime}")
        print(f"rank profile={sorted(profile.items())}")
        print(f"nonzero-scale projective points={len(survivors)}")
        print(
            "support profile="
            f"{sorted(Counter(sum(bool(value) for value in point) for point in survivors).items())}"
        )
        print(
            "survivors="
            f"{[tuple(int(value) for value in point) for point in survivors[:50]]}"
        )
        if len(survivors) > 50:
            print("survivor list truncated to first 50 points")


if __name__ == "__main__":
    main()

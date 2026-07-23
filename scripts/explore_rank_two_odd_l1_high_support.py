#!/usr/bin/env python3
"""Compress the high-support L1-normal quantization family.

The exact L1 normal cone is a rank-two quadric in fourteen variables.  A
linear third-order relaxation cuts it to a nine-dimensional section.  On
that section one global second correction solves the zero-scale affine
equation, and the coupling pencil depends on two linear forms.  This script
constructs the cubic nonzero-scale obstruction modulo the union of that
pencil.
"""

from __future__ import annotations

import sympy as sp
from sympy.polys.domains import QQ
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
from explore_degree_five_quantum_residue import solve_affine
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    essential_problem,
    linear_combination,
    solve_many_particular,
    split_correction,
)
from explore_rank_two_odd_normal_cones import normal_cone
from explore_rank_two_odd_support_three_curves import (
    project_to_third_cokernel,
)


SECTION_BASIS = (
    {3: QQ.one},
    {6: QQ.one, 0: -QQ(735704) / QQ(20853)},
    {7: QQ.one, 1: -QQ(3682276) / QQ(46431)},
    {8: QQ.one},
    {9: QQ.one, 4: -QQ(1429) / QQ(84)},
    {10: QQ.one},
    {11: QQ.one, 5: QQ(8278821) / QQ(146930)},
    {12: QQ.one, 4: -QQ(7397) / QQ(210)},
    {13: QQ.one, 0: -QQ(195850553) / QQ(81095)},
)

UNIFORM_KERNEL_CORRECTION = {
    18: -QQ(710122059) / QQ(392),
    33: -QQ(922185) / QQ(56),
}


def section_vectors():
    """Return the nine exact essential-coordinate section vectors."""

    data = normal_cone("L1", QQ)
    vectors = []
    for section_vector in SECTION_BASIS:
        essential = {}
        for tangent_index, coefficient in section_vector.items():
            for normal_index, value in data["tangent_basis"][
                tangent_index
            ].items():
                essential_index = data["normal_variable_indices"][
                    normal_index
                ]
                essential[essential_index] = (
                    essential.get(essential_index, QQ.zero)
                    + coefficient * value
                )
        vectors.append(
            {
                index: value
                for index, value in essential.items()
                if value
            }
        )
    return vectors


def add_sparse(left, right, coefficient, field):
    result = dict(left)
    for index, value in right.items():
        result[index] = (
            result.get(index, field.zero) + coefficient * value
        )
        if not result[index]:
            result.pop(index)
    return result


def projected_section_data(vectors, field):
    projected_couplings, projected_right_sides, project = (
        project_to_third_cokernel()
    )
    couplings = []
    right_sides = []
    for vector in vectors:
        row = []
        for kernel_index in range(42):
            column = {}
            for essential_index, coefficient in vector.items():
                column = add_sparse(
                    column,
                    projected_couplings[essential_index][kernel_index],
                    coefficient,
                    field,
                )
            row.append(column)
        couplings.append(row)
        right_side = {}
        for essential_index, coefficient in vector.items():
            right_side = add_sparse(
                right_side,
                projected_right_sides[essential_index],
                coefficient,
                field,
            )
        right_sides.append(right_side)
    return couplings, right_sides, project


def verify_uniform_affine_solution(couplings, right_sides, field):
    columns = [
        {
            (section_index, coordinate): value
            for section_index in range(len(couplings))
            for coordinate, value in couplings[section_index][
                kernel_index
            ].items()
        }
        for kernel_index in range(42)
    ]
    right_side = {
        (section_index, coordinate): value
        for section_index, vector in enumerate(right_sides)
        for coordinate, value in vector.items()
    }
    particular, kernel, rank = solve_affine(
        columns,
        right_side,
        field,
    )
    assert particular == {
        index: field.convert(value)
        for index, value in UNIFORM_KERNEL_CORRECTION.items()
    }
    assert (rank, len(kernel)) == (11, 31)


def quadratic_particular_pairs(S, T, pairs, field):
    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    second_columns = [
        poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    second_columns += [
        poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    indices = []
    right_sides = []
    for left in range(len(pairs)):
        for right in range(left, len(pairs)):
            if left == right:
                bracket = poisson(pairs[left][0], pairs[left][1])
            else:
                bracket = add(
                    poisson(pairs[left][0], pairs[right][1]),
                    poisson(pairs[right][0], pairs[left][1]),
                )
            indices.append((left, right))
            right_sides.append(scale(bracket, -field.one))
    vectors = solve_many_particular(
        second_columns,
        right_sides,
        field,
    )
    return indices, [
        split_correction(
            vector,
            s2_monomials,
            t2_monomials,
        )
        for vector in vectors
    ]


def cubic_terms(pairs, quadratic_indices, quadratic_pairs, project):
    return {
        (first, left, right): project(
            coupling(pairs[first], quadratic_pair)
        )
        for first in range(len(pairs))
        for (left, right), quadratic_pair in zip(
            quadratic_indices,
            quadratic_pairs,
        )
    }


def union_cokernel(couplings, extra_vectors, field):
    union_columns = [
        column for row in couplings for column in row if column
    ]
    all_vectors = union_columns + list(extra_vectors)
    ambient = sorted(set().union(*(set(vector) for vector in all_vectors)))
    ambient_index = {
        coordinate: index for index, coordinate in enumerate(ambient)
    }
    transpose = {
        column: {
            ambient_index[coordinate]: value
            for coordinate, value in vector.items()
        }
        for column, vector in enumerate(union_columns)
    }
    reduced, pivots, nonzero = sdm_irref(transpose)
    functionals, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(ambient),
        pivots,
        nonzero,
    )
    incidence = {}
    for functional_index, functional in enumerate(functionals):
        for coordinate, value in functional.items():
            incidence.setdefault(
                ambient[coordinate],
                [],
            ).append((functional_index, value))

    def project(vector):
        result = {}
        for coordinate, coefficient in vector.items():
            for functional_index, weight in incidence.get(coordinate, ()):
                result[functional_index] = (
                    result.get(functional_index, field.zero)
                    + coefficient * weight
                )
        return {
            coordinate: value
            for coordinate, value in result.items()
            if value
        }

    return len(pivots), project


def outside_union_polynomials(couplings, cubics, field):
    union_rank, project = union_cokernel(
        couplings,
        cubics.values(),
        field,
    )
    assert union_rank == 12
    y = sp.symbols(f"z0:{len(couplings)}")
    expressions = {}
    for (first, left, right), vector in cubics.items():
        for coordinate, coefficient in project(vector).items():
            expressions[coordinate] = (
                expressions.get(coordinate, sp.S.Zero)
                + field.to_sympy(coefficient)
                * y[first]
                * y[left]
                * y[right]
            )
    polynomials = [
        sp.Poly(sp.expand(expression), *y, domain=field)
        for expression in expressions.values()
        if expression
    ]
    if not polynomials:
        return (0, sp.S.Zero, 0, [])
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = sp.gcd(common, polynomial)
    reduced = [
        sp.quo(polynomial, common)
        for polynomial in polynomials
    ]
    reduced_basis = sp.groebner(
        [polynomial.as_expr() for polynomial in reduced],
        *y,
        order="grevlex",
        domain=field,
    )
    return (
        len(polynomials),
        sp.factor(common.as_expr()),
        len(reduced_basis.polys),
        [
            sp.factor(polynomial.as_expr())
            for polynomial in reduced_basis.polys
        ],
    )


def fixed_b_polynomials(couplings, cubics, field):
    """Restrict z4=0 and project the cubic column modulo the fixed B-map."""

    restricted = {
        index: vector
        for index, vector in cubics.items()
        if 4 not in index
    }
    rank, project = union_cokernel(
        [couplings[5]],
        restricted.values(),
        field,
    )
    assert rank == 11
    z = sp.symbols("z0:8")
    expressions = {}
    for (first, left, right), vector in restricted.items():
        for coordinate, coefficient in project(vector).items():
            expressions[coordinate] = (
                expressions.get(coordinate, sp.S.Zero)
                + field.to_sympy(coefficient)
                * z[first]
                * z[left]
                * z[right]
            )
    polynomials = [
        sp.Poly(sp.expand(expression), *z, domain=field)
        for expression in expressions.values()
        if expression
    ]
    if not polynomials:
        return (0, sp.S.Zero, 0, [])
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = sp.gcd(common, polynomial)
    reduced = [
        sp.quo(polynomial, common)
        for polynomial in polynomials
    ]
    reduced_basis = sp.groebner(
        [polynomial.as_expr() for polynomial in reduced],
        *z,
        order="grevlex",
        domain=field,
    )
    return (
        len(polynomials),
        sp.factor(common.as_expr()),
        len(reduced_basis.polys),
        [
            sp.factor(polynomial.as_expr())
            for polynomial in reduced_basis.polys
        ],
    )


def rank_zero_boundary_polynomials(cubics, field):
    """Restrict z4=z5=0 and test whether the cubic column itself vanishes."""

    z = sp.symbols("z0:8")
    expressions = {}
    for (first, left, right), vector in cubics.items():
        if 4 in (first, left, right) or 5 in (first, left, right):
            continue
        for coordinate, coefficient in vector.items():
            expressions[coordinate] = (
                expressions.get(coordinate, sp.S.Zero)
                + field.to_sympy(coefficient)
                * z[first]
                * z[left]
                * z[right]
            )
    polynomials = [
        sp.Poly(sp.expand(expression), *z, domain=field)
        for expression in expressions.values()
        if expression
    ]
    return [
        polynomial
        for polynomial in polynomials
        if not polynomial.is_zero
    ]


def main() -> None:
    alpha_expression = sp.sqrt(-5)
    field = QQ.algebraic_field(alpha_expression)
    alpha = field.from_sympy(alpha_expression)
    root = -field(2) + alpha / field(15)
    ambient_vectors = section_vectors()
    branch_coordinates = (0, 1, 2, 4, 5, 6, 7, 8)
    vectors = []
    for coordinate in branch_coordinates:
        vector = {
            index: field.convert(value)
            for index, value in ambient_vectors[coordinate].items()
        }
        if coordinate == 5:
            for index, value in ambient_vectors[3].items():
                vector[index] = (
                    vector.get(index, field.zero)
                    + root * field.convert(value)
                )
        vectors.append(
            {
                index: value
                for index, value in vector.items()
                if value
            }
        )
    S, T, _, essential_pairs, _, _ = essential_problem()
    S = {
        monomial: field.convert(value)
        for monomial, value in S.items()
    }
    T = {
        monomial: field.convert(value)
        for monomial, value in T.items()
    }
    essential_pairs = [
        (
            {
                monomial: field.convert(value)
                for monomial, value in pair[0].items()
            },
            {
                monomial: field.convert(value)
                for monomial, value in pair[1].items()
            },
        )
        for pair in essential_pairs
    ]
    pairs = [
        linear_combination(essential_pairs, vector, field)
        for vector in vectors
    ]
    couplings, right_sides, project = projected_section_data(
        vectors,
        field,
    )
    verify_uniform_affine_solution(couplings, right_sides, field)

    # On the branch, only z4=y5 and z5=y6 enter the coupling pencil.
    assert all(
        not any(couplings[index])
        for index in (0, 1, 2, 3, 6, 7)
    )
    flattened = [
        {
            (kernel, coordinate): value
            for kernel, column in enumerate(row)
            for coordinate, value in column.items()
        }
        for row in couplings
    ]
    indices, quadratic_pairs = quadratic_particular_pairs(
        S,
        T,
        pairs,
        field,
    )
    cubics = cubic_terms(pairs, indices, quadratic_pairs, project)
    profile = outside_union_polynomials(couplings, cubics, field)
    fixed_profile = fixed_b_polynomials(couplings, cubics, field)
    rank_zero_polynomials = rank_zero_boundary_polynomials(
        cubics,
        field,
    )
    z4 = sp.Symbol("z4")
    assert profile[0] == 329
    assert sp.Poly(profile[1], z4, domain=field).monic() == sp.Poly(
        z4**3,
        z4,
        domain=field,
    )
    assert profile[2:] == (1, [sp.S.One])
    assert fixed_profile == (0, sp.S.Zero, 0, [])
    assert not rank_zero_polynomials

    print(
        "PASS: the high-support section has one global zero-scale "
        "correction and a two-form rank-11 coupling pencil"
    )
    print(
        "branch field QQ(sqrt(-5)); coupling forms z4 and z5; "
        "union rank 12"
    )
    print(
        f"outside-union cubic coordinates={profile[0]}; "
        f"common factor={profile[1]}"
    )
    print(
        f"reduced Groebner basis size={profile[2]}; "
        f"basis={profile[3]}"
    )
    print(
        f"fixed-B cokernel coordinates={fixed_profile[0]}; "
        f"common factor={fixed_profile[1]}"
    )
    print(
        f"fixed-B reduced Groebner basis size={fixed_profile[2]}; "
        f"basis={fixed_profile[3]}"
    )
    print(
        "rank-zero boundary cubic coordinates="
        f"{len(rank_zero_polynomials)}"
    )


if __name__ == "__main__":
    main()

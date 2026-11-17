#!/usr/bin/env python3
"""Test the excess-conormal map controlling the degree-18 Tor defect.

After cutting by two common coordinates on E_(6,6,6), the ambient target has
fifteen coordinates and every normalization branch has four parameters.  The
two-jets determine the fibers of the conormal modules of the three branches,
their pairwise intersections, and their triple intersection.  This script
constructs those fibers directly and tests the map

    E_13 + E_23 -> E_(12),3

between the corresponding excess-conormal spaces.  Surjectivity is the
residue-field version of the Tor-surjectivity criterion for the conductor
equalizer.  It also checks that the two pairwise conductor ideals on one
normalization sheet have a two-generator common complete-intersection core,
which yields the quartic generator in their intersection.
"""

import itertools

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from explore_degree18_conductor_equalizer import (
    branch_jets,
    column_space,
    exponent_basis,
    multiply_polynomials,
    polynomial_dict,
    roots,
    transverse_base_fiber,
    verify_transverse_base_coordinates,
)


ORDER = 2
TARGET_COUNT = 15


def independent_columns(matrix):
    """Return indices of a left-to-right independent set of columns."""
    if matrix.cols == 0:
        return ()
    _, pivots = DomainMatrix.from_Matrix(matrix).rref()
    return tuple(pivots)


def target_monomials():
    return tuple(
        exponent
        for total in (1, 2)
        for exponent in itertools.product(range(total + 1), repeat=TARGET_COUNT)
        if sum(exponent) == total
    )


def evaluation_matrix(branch, target_basis):
    """Evaluate positive-degree target two-jets on one branch."""
    variables, jets = branch
    branch_basis = exponent_basis(len(variables), ORDER)[1:]
    branch_index = {exponent: index for index, exponent in enumerate(branch_basis)}
    zero_exponent = (0,) * len(variables)
    one = {zero_exponent: sp.Integer(1)}
    generators = tuple(
        polynomial_dict(function, variables, ORDER) for function in jets[2:]
    )
    entries = {}
    for column, target_exponent in enumerate(target_basis):
        value = one
        for generator_index, exponent in enumerate(target_exponent):
            for _ in range(exponent):
                value = multiply_polynomials(
                    value, generators[generator_index], ORDER
                )
        for exponent, coefficient in value.items():
            if exponent != zero_exponent:
                entries[(branch_index[exponent], column)] = coefficient
    return sp.MutableSparseMatrix(
        len(branch_basis), len(target_basis), entries
    ).as_immutable()


def ideal_kernel(evaluation):
    kernel = evaluation.nullspace()
    return sp.Matrix.hstack(*kernel) if kernel else sp.zeros(evaluation.cols, 0)


def maximal_ideal_multiple_space(ideal, target_basis, variable_count):
    """Return the degree-at-most-two part of m times a truncated ideal."""
    linear_count = variable_count
    linear_parts = column_space(ideal[:linear_count, :])
    basis_index = {exponent: index for index, exponent in enumerate(target_basis)}
    columns = []
    for generator in range(linear_parts.cols):
        for variable in range(variable_count):
            column = sp.zeros(len(target_basis), 1)
            for source in range(variable_count):
                coefficient = linear_parts[source, generator]
                if not coefficient:
                    continue
                exponent = list(target_basis[source])
                exponent[variable] += 1
                column[basis_index[tuple(exponent)], 0] += coefficient
            columns.append(column)
    if not columns:
        return sp.zeros(len(target_basis), 0)
    return column_space(sp.Matrix.hstack(*columns))


class ConormalFiber:
    """The quotient I/mI and reduction maps in a two-jet presentation."""

    def __init__(self, ideal, target_basis, variable_count=TARGET_COUNT):
        self.ideal = column_space(ideal)
        self.multiple = maximal_ideal_multiple_space(
            self.ideal, target_basis, variable_count
        )
        assert column_space(self.ideal.row_join(self.multiple)).cols == self.ideal.cols

        if self.multiple.cols:
            reduced_rows, pivots = DomainMatrix.from_Matrix(
                self.multiple.T
            ).rref()
            self.reduction_rows = reduced_rows.to_Matrix()
            self.reduction_pivots = tuple(pivots)
        else:
            self.reduction_rows = sp.zeros(0, self.ideal.rows)
            self.reduction_pivots = ()

        reduced_ideal = self.reduce(self.ideal)
        chosen = independent_columns(reduced_ideal)
        self.representatives = self.ideal[:, chosen]
        self.basis = reduced_ideal[:, chosen]
        row_pivots = independent_columns(self.basis.T)
        self.coordinate_rows = row_pivots
        square = self.basis[list(row_pivots), :]
        self.coordinate_inverse = square.inv()

    @property
    def dimension(self):
        return self.basis.cols

    def reduce(self, vectors):
        if not self.reduction_pivots:
            return vectors
        pivot_values = vectors[list(self.reduction_pivots), :]
        return vectors - self.reduction_rows.T * pivot_values

    def coordinates(self, representatives):
        reduced = self.reduce(representatives)
        coordinates = self.coordinate_inverse * reduced[list(self.coordinate_rows), :]
        assert self.basis * coordinates == reduced
        return coordinates


def kernel_basis(matrix):
    kernel = matrix.nullspace()
    return sp.Matrix.hstack(*kernel) if kernel else sp.zeros(matrix.cols, 0)


def main():
    verify_transverse_base_coordinates()
    branches = [
        transverse_base_fiber(branch_jets(root, ORDER), index, ORDER)
        for index, root in enumerate(roots)
    ]
    target_basis = target_monomials()
    evaluations = [evaluation_matrix(branch, target_basis) for branch in branches]
    branch_ideals = [ideal_kernel(evaluation) for evaluation in evaluations]
    for evaluation, ideal in zip(evaluations, branch_ideals):
        multiple = maximal_ideal_multiple_space(ideal, target_basis, TARGET_COUNT)
        assert evaluation * multiple == sp.zeros(evaluation.rows, multiple.cols)
    pair_ideals = {
        (left, right): column_space(
            branch_ideals[left].row_join(branch_ideals[right])
        )
        for left, right in itertools.combinations(range(3), 2)
    }
    triple_ideal = column_space(
        branch_ideals[0].row_join(branch_ideals[1]).row_join(branch_ideals[2])
    )

    branch_conormals = [ConormalFiber(ideal, target_basis) for ideal in branch_ideals]
    pair_conormals = {
        pair: ConormalFiber(ideal, target_basis) for pair, ideal in pair_ideals.items()
    }
    triple_conormal = ConormalFiber(triple_ideal, target_basis)

    def inclusion(source, target):
        return target.coordinates(source.representatives)

    pair_excess = {}
    for pair, conormal in pair_conormals.items():
        left, right = pair
        sum_map = inclusion(branch_conormals[left], conormal).row_join(
            inclusion(branch_conormals[right], conormal)
        )
        pair_excess[pair] = kernel_basis(sum_map)

    target_sum_map = inclusion(pair_conormals[(0, 1)], triple_conormal).row_join(
        inclusion(branch_conormals[2], triple_conormal)
    )
    target_excess = kernel_basis(target_sum_map)

    source_13 = pair_excess[(0, 2)]
    source_23 = pair_excess[(1, 2)]
    map_13 = sp.diag(
        inclusion(branch_conormals[0], pair_conormals[(0, 1)]),
        sp.eye(branch_conormals[2].dimension),
    ) * source_13
    map_23 = sp.diag(
        inclusion(branch_conormals[1], pair_conormals[(0, 1)]),
        sp.eye(branch_conormals[2].dimension),
    ) * source_23
    source_map = map_13.row_join(map_23)
    assert target_sum_map * source_map == sp.zeros(target_sum_map.rows, source_map.cols)
    target_rows = independent_columns(target_excess.T)
    target_inverse = target_excess[list(target_rows), :].inv()
    induced_coordinates = target_inverse * source_map[list(target_rows), :]
    assert target_excess * induced_coordinates == source_map
    induced_rank = induced_coordinates.rank()
    minor_columns = independent_columns(induced_coordinates)
    maximal_minor = induced_coordinates[:, list(minor_columns)].det()

    # On the third branch, pull back the first and second branch ideals.  The
    # resulting height-three ideals are the two pairwise conductor kernels.
    branch_basis = exponent_basis(4, ORDER)[1:]
    conductor_13 = evaluations[2] * branch_conormals[0].representatives
    conductor_23 = evaluations[2] * branch_conormals[1].representatives
    conductor_sum = column_space(conductor_13.row_join(conductor_23))
    conductor_conormal_13 = ConormalFiber(
        conductor_13, branch_basis, variable_count=4
    )
    conductor_conormal_23 = ConormalFiber(
        conductor_23, branch_basis, variable_count=4
    )
    conductor_sum_conormal = ConormalFiber(
        conductor_sum, branch_basis, variable_count=4
    )
    conductor_sum_map = inclusion(
        conductor_conormal_13, conductor_sum_conormal
    ).row_join(inclusion(conductor_conormal_23, conductor_sum_conormal))
    conductor_excess = kernel_basis(conductor_sum_map)
    left_excess_coordinates = conductor_excess[
        : conductor_conormal_13.dimension, :
    ]
    right_excess_coordinates = conductor_excess[
        conductor_conormal_13.dimension :, :
    ]
    left_projection_rank = left_excess_coordinates.rank()
    right_projection_rank = right_excess_coordinates.rank()
    left_common_generators = (
        conductor_conormal_13.representatives * left_excess_coordinates
    )
    right_common_generators = (
        conductor_conormal_23.representatives * right_excess_coordinates
    )
    common_linear_ranks = (
        left_common_generators[:4, :].rank(),
        right_common_generators[:4, :].rank(),
    )

    assert tuple(item.dimension for item in branch_conormals) == (11, 11, 11)
    assert tuple(
        pair_conormals[pair].dimension for pair in sorted(pair_conormals)
    ) == (14, 14, 14)
    assert triple_conormal.dimension == 15
    assert tuple(pair_excess[pair].cols for pair in sorted(pair_excess)) == (8, 8, 8)
    assert target_excess.cols == 10
    assert induced_rank == target_excess.cols
    assert maximal_minor != 0
    assert (
        conductor_conormal_13.dimension,
        conductor_conormal_23.dimension,
        conductor_sum_conormal.dimension,
    ) == (3, 3, 4)
    assert conductor_excess.cols == 2
    assert (left_projection_rank, right_projection_rank) == (2, 2)
    assert common_linear_ranks == (1, 1)

    print("branch conormal dimensions:", tuple(item.dimension for item in branch_conormals))
    print(
        "pair conormal dimensions:",
        tuple(pair_conormals[pair].dimension for pair in sorted(pair_conormals)),
    )
    print("triple conormal dimension:", triple_conormal.dimension)
    print(
        "pair excess dimensions:",
        tuple(pair_excess[pair].cols for pair in sorted(pair_excess)),
    )
    print("target excess dimension:", target_excess.cols)
    print("induced excess-map rank:", induced_rank)
    print("nonzero maximal minor:", maximal_minor)
    print(
        "one-sheet conductor conormal dimensions:",
        (
            conductor_conormal_13.dimension,
            conductor_conormal_23.dimension,
            conductor_sum_conormal.dimension,
        ),
    )
    print("one-sheet common-core dimension:", conductor_excess.cols)
    print(
        "common-core projection ranks:",
        (left_projection_rank, right_projection_rank),
    )
    print("common-core linear ranks:", common_linear_ranks)
    print("PASS: the residue excess-conormal map is surjective")


if __name__ == "__main__":
    main()

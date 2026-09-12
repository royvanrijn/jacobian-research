#!/usr/bin/env python3
"""Exact all-pole hbar^5 certificate on the kappa=-1 replacement chart."""

from __future__ import annotations

from collections import defaultdict

from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import (
    EXCEPTIONAL_PROFILE,
    add,
    degree_five_exceptional_family,
    laurent_monomials,
    scale,
    solve_affine,
)
from verify_degree_five_laurent_quantum_obstruction import (
    functional,
    split_coordinates,
)


SUPPORT = (
    (17, 1, 3),
    (14, 0, 2),
    (16, 0, 4),
    (16, 1, 2),
    (15, 0, 3),
    (15, 1, 1),
    (13, 0, 1),
)

COEFFICIENTS = (
    (1, 1),
    (35626, 1281),
    (10436077, 320250),
    (-388, 183),
    (-924622, 32025),
    (234, 61),
    (-2106, 61),
)


def monomial_bracket_functional(left, right, coefficients):
    """Evaluate Lambda({left,right}) for the exceptional Ore derivation."""

    i, j, k = left
    ii, jj, kk = right
    z_degree = k + kk - 1
    if z_degree < 0:
        return QQ.zero
    value = QQ.zero
    shifted = (i + ii + 2, j + jj, z_degree)
    shifted_coefficient = (
        k * (-2 * ii + 6 * jj) - (-2 * i + 6 * j) * kk
    )
    if shifted_coefficient and shifted in coefficients:
        value += shifted_coefficient * coefficients[shifted]
    unshifted = (i + ii, j + jj - 1, z_degree)
    unshifted_coefficient = 2 * (k * jj - j * kk)
    if unshifted_coefficient and unshifted in coefficients:
        value += unshifted_coefficient * coefficients[unshifted]
    return value


def main() -> None:
    field = QQ
    profile = EXCEPTIONAL_PROFILE
    coefficient_functional = {
        monomial: field(numerator, denominator)
        for monomial, (numerator, denominator) in zip(
            SUPPORT,
            COEFFICIENTS,
            strict=True,
        )
    }
    S, T = degree_five_exceptional_family(field, field.one)

    # The exceptional derivation shifts X-degree by at most two.  The d5
    # columns therefore shift by at most 15 and 16.  Only input degrees
    # -3 and above can meet Lambda's support; checking that finite band
    # proves annihilation on the whole Laurent union.
    assert min(x_degree for x_degree, _, _ in SUPPORT) == 13
    assert max(x_degree for x_degree, _, _ in T) == 13
    assert max(x_degree for x_degree, _, _ in S) == 14
    s4_monomials = laurent_monomials(24, 1, 3, 4)
    t4_monomials = laurent_monomials(20, 0, 3, 4)
    for monomial in s4_monomials:
        assert functional(
            profile.poisson({monomial: field.one}, T),
            coefficient_functional,
        ) == 0
    for monomial in t4_monomials:
        assert functional(
            profile.poisson(S, {monomial: field.one}),
            coefficient_functional,
        ) == 0

    # Quadratic terms meeting Lambda can use S2 degrees down to -13 and T2
    # degrees down to -17.  Inputs below -17 cannot affect the hbar^3
    # equations of output degree >= -1, because the two d3 shifts are at
    # most 15 and 16.
    pole_cutoff = 17
    s2_monomials = laurent_monomials(28, 3, pole_cutoff, 4)
    t2_monomials = laurent_monomials(24, 2, pole_cutoff, 4)
    split = len(s2_monomials)
    columns = [
        profile.poisson({monomial: field.one}, T)
        for monomial in s2_monomials
    ]
    columns += [
        profile.poisson(S, {monomial: field.one})
        for monomial in t2_monomials
    ]
    rhs = scale(
        profile.pi_power(S, T, 3),
        -field.one / field(24),
    )
    projected_columns = [
        {
            monomial: value
            for monomial, value in column.items()
            if monomial[0] >= -1
        }
        for column in columns
    ]
    projected_rhs = {
        monomial: value
        for monomial, value in rhs.items()
        if monomial[0] >= -1
    }
    base, kernel, rank = solve_affine(
        projected_columns,
        projected_rhs,
        field,
    )
    assert len(columns) == 5559
    assert rank == 3627
    assert len(kernel) == 1932

    bracket_matrix = defaultdict(list)
    for s_index, s_monomial in enumerate(s2_monomials):
        for t_index, t_monomial in enumerate(t2_monomials):
            value = monomial_bracket_functional(
                s_monomial,
                t_monomial,
                coefficient_functional,
            )
            if value:
                bracket_matrix[s_index].append((t_index, value))
    assert sum(map(len, bracket_matrix.values())) == 1147

    base_s, base_t = split_coordinates(base, split)
    kernel_pairs = [
        split_coordinates(vector, split)
        for vector in kernel
    ]

    def bilinear(s_part, t_part):
        return sum(
            s_value * matrix_value * t_part.get(t_index, field.zero)
            for s_index, s_value in s_part.items()
            for t_index, matrix_value in bracket_matrix.get(s_index, ())
        )

    left_images = []
    for s_part, _ in kernel_pairs:
        image = defaultdict(lambda: field.zero)
        for s_index, s_value in s_part.items():
            for t_index, matrix_value in bracket_matrix.get(s_index, ()):
                image[t_index] += s_value * matrix_value
        left_images.append(
            {
                t_index: value
                for t_index, value in image.items()
                if value
            }
        )
    t_incidence = defaultdict(list)
    for basis_index, (_, t_part) in enumerate(kernel_pairs):
        for t_index, value in t_part.items():
            t_incidence[t_index].append((basis_index, value))
    quadratic_matrix = [
        defaultdict(lambda: field.zero)
        for _ in kernel
    ]
    for left_index, image in enumerate(left_images):
        for t_index, left_value in image.items():
            for right_index, right_value in t_incidence.get(t_index, ()):
                quadratic_matrix[left_index][right_index] += (
                    left_value * right_value
                )
    assert sum(len(row) for row in quadratic_matrix) == 8
    for left_index, row in enumerate(quadratic_matrix):
        for right_index, value in row.items():
            assert (
                value
                + quadratic_matrix[right_index].get(left_index, field.zero)
                == 0
            )

    moyal_linear = []
    for monomial in s2_monomials:
        moyal_linear.append(
            functional(
                profile.pi_power({monomial: field.one}, T, 3),
                coefficient_functional,
            )
            / field(24)
        )
    for monomial in t2_monomials:
        moyal_linear.append(
            functional(
                profile.pi_power(S, {monomial: field.one}, 3),
                coefficient_functional,
            )
            / field(24)
        )
    for vector, (s_part, t_part) in zip(
        kernel,
        kernel_pairs,
        strict=True,
    ):
        derivative = sum(
            moyal_linear[column] * value
            for column, value in vector.items()
        )
        derivative += bilinear(base_s, t_part)
        derivative += bilinear(s_part, base_t)
        assert derivative == 0

    base_s_poly = {
        s2_monomials[column]: value
        for column, value in base_s.items()
    }
    base_t_poly = {
        t2_monomials[column]: value
        for column, value in base_t.items()
    }
    defect = profile.poisson(base_s_poly, base_t_poly)
    defect = add(
        defect,
        profile.pi_power(base_s_poly, T, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, base_t_poly, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        profile.pi_power(S, T, 5),
        field.one / field(1920),
    )
    obstruction_value = functional(defect, coefficient_functional)
    assert obstruction_value == field(40998496, 305)

    print("PASS: seven-term Lambda_-1 annihilates d5 on the full Laurent space")
    print("PASS: projected exceptional hbar^3 affine dimension is 1932 over Q")
    print("PASS: Lambda_-1(O5) is independent of every Laurent hbar^3 solution")
    print(f"PASS: Lambda_-1(O5)={obstruction_value} is nonzero")
    print("CONCLUSION: the kappa=-1 obstruction also survives localization at X")


if __name__ == "__main__":
    main()

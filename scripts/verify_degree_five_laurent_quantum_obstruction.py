#!/usr/bin/env python3
"""Exact-Q certificates that hbar^5 obstructions survive A[X^-1].

This is stronger than checking finitely many pole bands.  At the completed
degree-five seed ``(kappa,tau)=(0,1)``, an explicit 16-term functional
``Lambda`` has the following properties:

* it annihilates every bounded hbar^5 correction in
  ``Q[X,X^-1,Q,Z]``;
* its value on the hbar^5 defect is independent of every hbar^3 Laurent
  solution; and
* that value is the nonzero rational number
  ``-47547660815739/190658``.

Only a finite calculation is required.  Lambda is supported in X-degrees
14 through 19.  The hbar^5 linearized operator shifts the X-degree of an
input by at most 12, so negative-X correction monomials cannot meet its
support.  In the quadratic hbar^3 contribution, only S2 monomials of
X-degree at least -8 and T2 monomials of X-degree at least -12 can meet the
support.  All constraints on those coefficients which can be affected by a
Laurent solution occur in the projection of the hbar^3 equation to output
X-degree at least -1; inputs below degree -13 cannot enter those equations.

The projected affine space has 4065 variables, rank 2990, and dimension
1075.  The script proves exactly over Q that Lambda(O5) is constant on this
*larger* affine space.  It therefore applies a fortiori to every finite
Laurent hbar^3 solution, without a pole-order cutoff.

Besides the two frozen 16-term certificates, dynamic rational-seed mode
constructs the unique bounded 16-support period before applying the same
all-pole test.  A separate 15-term vertical certificate covers the rational
four-period common point (kappa,tau)=(0,-3).
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement

from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref, sdm_nullspace_from_rref

from explore_degree_five_quantum_residue import (
    GENERIC_PROFILE,
    add,
    degree_five_family,
    laurent_monomials,
    scale,
    solve_affine,
    third_order_family,
)


SUPPORT = (
    (19, 7, 0),
    (18, 2, 2),
    (14, 0, 1),
    (18, 4, 1),
    (17, 1, 2),
    (18, 6, 0),
    (16, 0, 2),
    (17, 3, 1),
    (19, 1, 3),
    (17, 5, 0),
    (16, 2, 1),
    (19, 3, 2),
    (15, 1, 1),
    (18, 0, 3),
    (19, 5, 1),
    (16, 4, 0),
)

COEFFICIENTS = (
    (1, 1),
    (-8034972745707, 26158277600),
    (-23945710924881, 4270739200),
    (12784957587, 373689680),
    (78951458085327, 104633110400),
    (-9991676, 3336515),
    (-532025132768811, 104633110400),
    (-214398955599, 3736896800),
    (-29318152220199, 13079138800),
    (77120461, 13346060),
    (-139051975887, 3736896800),
    (144065122713, 747379360),
    (-727690792977, 2135369600),
    (-364131907959447, 1464863545600),
    (-47876019, 3336515),
    (167077941, 16682575),
)

GENERIC_POINT_COEFFICIENTS = (
    (1, 1),
    (169871103312041400321, 279529273836772288),
    (2400468331172504613723, 1140935811578662400),
    (-289943964170372286, 15598731798927025),
    (-2517657172173452450637, 1996637670262659200),
    (-43648175856441, 356542441118332),
    (-63912093944830408405071, 27952927383677228800),
    (1646335131560893377, 10186926889095200),
    (-489296995398671811209, 174705796147982680),
    (-6450552927065511, 509346344454760),
    (617997613248231505677, 1996637670262659200),
    (60105473168763122293, 249579708782832400),
    (328464851230938900393, 570467905789331200),
    (-131396097305693207121123, 12229405730358787600),
    (-10332179057796907, 636682930568450),
    (-262082879373072357, 17827122055916600),
)

RATIONAL_COMMON_POINT_SUPPORT = (
    (14, 0, 1),
    (15, 1, 1),
    (16, 0, 2),
    (16, 2, 1),
    (16, 4, 0),
    (17, 1, 2),
    (17, 3, 1),
    (17, 5, 0),
    (18, 2, 2),
    (18, 4, 1),
    (18, 6, 0),
    (19, 1, 3),
    (19, 3, 2),
    (19, 5, 1),
    (20, 4, 2),
)

RATIONAL_COMMON_POINT_COEFFICIENTS = (
    (151258183, 2722734000),
    (296575901, 31505922000),
    (-12015101341, 1029193452000),
    (330761947, 330812181000),
    (-17176, 2531725875),
    (267797749, 441082908000),
    (1436809, 13502538000),
    (-8672, 1519035525),
    (3130373, 15752961000),
    (-781, 168781725),
    (-16, 217005075),
    (8496151, 42007896000),
    (-443, 168781725),
    (-16, 217005075),
    (-4, 72335025),
)


def functional(poly, coefficients):
    return sum(
        coefficients.get(monomial, QQ.zero) * value
        for monomial, value in poly.items()
    )


def split_coordinates(vector, split):
    return (
        {
            column: value
            for column, value in vector.items()
            if column < split
        },
        {
            column - split: value
            for column, value in vector.items()
            if column >= split
        },
    )


def monomial_bracket_functional(
    left,
    right,
    coefficients,
):
    """Evaluate Lambda({left,right}) without constructing the bracket."""

    i, j, k = left
    ii, jj, kk = right
    z_degree = k + kk - 1
    if z_degree < 0:
        return QQ.zero
    value = QQ.zero
    shifted = (i + ii + 1, j + jj, z_degree)
    shifted_coefficient = k * (3 * ii - 6 * jj) - (3 * i - 6 * j) * kk
    if shifted_coefficient and shifted in coefficients:
        value += shifted_coefficient * coefficients[shifted]
    unshifted = (i + ii, j + jj - 1, z_degree)
    unshifted_coefficient = 2 * (k * jj - j * kk)
    if unshifted_coefficient and unshifted in coefficients:
        value += unshifted_coefficient * coefficients[unshifted]
    return value


def bounded_period_functional(S, T, field):
    """Construct the unique normalized 16-support bounded period."""

    family = third_order_family(S, T, field)
    base_s, base_t = family.base
    constraints = []

    def restrict(poly):
        return {
            coordinate: poly.get(monomial, field.zero)
            for coordinate, monomial in enumerate(SUPPORT)
            if poly.get(monomial, field.zero)
        }

    for monomials, left in (
        (laurent_monomials(21, 1, 0, 3), True),
        (laurent_monomials(17, 0, 0, 3), False),
    ):
        for monomial in monomials:
            image = (
                GENERIC_PROFILE.poisson({monomial: field.one}, T)
                if left
                else GENERIC_PROFILE.poisson(S, {monomial: field.one})
            )
            row = restrict(image)
            if row:
                constraints.append(row)

    for direction_s, direction_t in family.kernel:
        variation = add(
            GENERIC_PROFILE.poisson(direction_s, base_t),
            GENERIC_PROFILE.poisson(base_s, direction_t),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(direction_s, T, 3),
            field.one / field(24),
        )
        variation = add(
            variation,
            GENERIC_PROFILE.pi_power(S, direction_t, 3),
            field.one / field(24),
        )
        row = restrict(variation)
        if row:
            constraints.append(row)

    for left, right in combinations_with_replacement(
        range(len(family.kernel)),
        2,
    ):
        left_s, left_t = family.kernel[left]
        right_s, right_t = family.kernel[right]
        variation = (
            GENERIC_PROFILE.poisson(left_s, left_t)
            if left == right
            else add(
                GENERIC_PROFILE.poisson(left_s, right_t),
                GENERIC_PROFILE.poisson(right_s, left_t),
            )
        )
        row = restrict(variation)
        if row:
            constraints.append(row)

    rows = {
        row_index: dict(row)
        for row_index, row in enumerate(constraints)
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(SUPPORT),
        pivots,
        nonzero,
    )
    if len(pivots) != 15 or len(kernel) != 1:
        raise ValueError(
            "the 16-support bounded period is not unique: "
            f"rank={len(pivots)}, kernel={len(kernel)}"
        )
    normalization_index = min(kernel[0])
    normalization = kernel[0][normalization_index]
    print(
        "BOUNDED_PERIOD_NORMALIZATION_INDEX="
        + str(normalization_index)
    )
    return {
        monomial: kernel[0].get(index, field.zero) / normalization
        for index, monomial in enumerate(SUPPORT)
        if kernel[0].get(index, field.zero)
    }


def parse_rational(value: str, field):
    parsed = Fraction(value)
    return field(parsed.numerator) / field(parsed.denominator)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--generic-point",
        action="store_true",
        help="verify (kappa,tau)=(1,1) instead of the known (0,1) sample",
    )
    parser.add_argument(
        "--rational-common-point",
        action="store_true",
        help="verify the vertical period at (kappa,tau)=(0,-3)",
    )
    parser.add_argument(
        "--kappa",
        help="dynamically construct the bounded period at this rational kappa",
    )
    parser.add_argument(
        "--tau",
        help="rational tau for --kappa (defaults to 1)",
    )
    parser.add_argument(
        "--profile-active-directions",
        action="store_true",
        help="print Laurent kernel directions which meet the period pairing",
    )
    parser.add_argument(
        "--compare-constraint-modules",
        action="store_true",
        help="build the full 16-coordinate Laurent period constraint module",
    )
    args = parser.parse_args()
    field = QQ
    profile = GENERIC_PROFILE
    dynamic = args.kappa is not None or args.tau is not None
    if sum((dynamic, args.generic_point, args.rational_common_point)) > 1:
        parser.error("choose only one alternate seed mode")
    if dynamic:
        kappa = parse_rational(args.kappa or "0", field)
        tau_value = parse_rational(args.tau or "1", field)
        if kappa in (-field(2), -field.one):
            parser.error("dynamic mode requires kappa != -2,-1")
        a_value = -(field.one + kappa) / (field(2) + kappa)
        expected_obstruction = None
        seed_label = f"(kappa,tau)=({kappa},{tau_value})"
        S, T = degree_five_family(field, a_value, tau_value)
        coefficient_functional = bounded_period_functional(
            S,
            T,
            field,
        )
        print(
            "DYNAMIC_BOUNDED_PERIOD_NONZERO_COEFFICIENTS="
            + str(len(coefficient_functional))
        )
        functional_support = tuple(coefficient_functional)
    elif args.rational_common_point:
        a_value = -field.one / field(2)
        tau_value = -field(3)
        expected_obstruction = field.one
        seed_label = "(kappa,tau)=(0,-3)"
        coefficient_functional = {
            monomial: field(numerator, denominator)
            for monomial, (numerator, denominator) in zip(
                RATIONAL_COMMON_POINT_SUPPORT,
                RATIONAL_COMMON_POINT_COEFFICIENTS,
                strict=True,
            )
        }
        functional_support = RATIONAL_COMMON_POINT_SUPPORT
        S, T = degree_five_family(
            field,
            a_value,
            tau_value,
        )
    elif args.generic_point:
        coefficient_data = GENERIC_POINT_COEFFICIENTS
        a_value = -field(2) / field(3)
        tau_value = field.one
        expected_obstruction = field(
            18408906052532094299482611,
            3338052203018715136,
        )
        seed_label = "(kappa,tau)=(1,1)"
    else:
        coefficient_data = COEFFICIENTS
        a_value = -field.one / field(2)
        tau_value = field.one
        expected_obstruction = field(-47547660815739, 190658)
        seed_label = "(kappa,tau)=(0,1)"
    if not dynamic:
        if not args.rational_common_point:
            coefficient_functional = {
                monomial: field(numerator, denominator)
                for monomial, (numerator, denominator) in zip(
                    SUPPORT,
                    coefficient_data,
                    strict=True,
                )
            }
            functional_support = SUPPORT
            S, T = degree_five_family(
                field,
                a_value,
                tau_value,
            )

    # The finite cutoff is derived from X-shift bounds, not guessed from a
    # pole-band search.
    assert min(x_degree for x_degree, _, _ in functional_support) == 14
    assert max(x_degree for x_degree, _, _ in T) == 10
    assert max(x_degree for x_degree, _, _ in S) == 11
    hbar5_s_shift = max(x_degree for x_degree, _, _ in T) + 1
    hbar5_t_shift = max(x_degree for x_degree, _, _ in S) + 1
    assert (hbar5_s_shift, hbar5_t_shift) == (11, 12)
    assert -1 + max(hbar5_s_shift, hbar5_t_shift) < 14

    # Every polynomial hbar^5 correction is killed directly.  The shift
    # inequality proves the same for every negative-X monomial, hence for
    # the full finite Laurent union.
    s4_monomials = laurent_monomials(21, 1, 0, 3)
    t4_monomials = laurent_monomials(17, 0, 0, 3)
    for monomial in s4_monomials:
        image = profile.poisson({monomial: field.one}, T)
        assert functional(image, coefficient_functional) == 0
    for monomial in t4_monomials:
        image = profile.poisson(S, {monomial: field.one})
        assert functional(image, coefficient_functional) == 0

    # Retain every hbar^2 variable which can influence either Lambda(O5) or
    # the high-output hbar^3 constraints on those variables.  Degree -13 is
    # the exact nuisance cutoff; lower inputs shift to output degree <= -2.
    pole_cutoff = 13
    s2_monomials = laurent_monomials(25, 3, pole_cutoff, 3)
    t2_monomials = laurent_monomials(21, 2, pole_cutoff, 3)
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
    assert len(columns) == 4065
    assert rank == 2990
    assert len(kernel) == 1075

    # Sparse matrix J represents Lambda({s,t}) on the ambient correction
    # monomial bases.  Only 3957 of the roughly four million possible pairs
    # are nonzero.
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
    bracket_nonzero = sum(map(len, bracket_matrix.values()))
    if not dynamic and not args.rational_common_point:
        assert bracket_nonzero == 3957

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

    # Compute A_ij=Lambda({s_i,t_j}) sparsely.  Constancy of the quadratic
    # term is exactly skew-symmetry of A.
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
    quadratic_nonzero = sum(len(row) for row in quadratic_matrix)
    if not dynamic and not args.rational_common_point:
        assert quadratic_nonzero == 45
    for left_index, row in enumerate(quadratic_matrix):
        for right_index, value in row.items():
            assert (
                value
                + quadratic_matrix[right_index].get(left_index, field.zero)
                == 0
            )

    # The linear derivative of Lambda(O5) along every kernel vector also
    # vanishes.  Precompute the Moyal-linear functional on ambient columns.
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
    active_linear_directions = []
    for basis_index, (vector, (s_part, t_part)) in enumerate(
        zip(
            kernel,
            kernel_pairs,
            strict=True,
        )
    ):
        moyal_part = sum(
            moyal_linear[column] * value
            for column, value in vector.items()
        )
        right_part = bilinear(base_s, t_part)
        left_part = bilinear(s_part, base_t)
        derivative = moyal_part + right_part + left_part
        if moyal_part or right_part or left_part:
            active_linear_directions.append(basis_index)
        assert derivative == 0

    base_s_poly = {
        s2_monomials[column]: value
        for column, value in base_s.items()
    }
    base_t_poly = {
        t2_monomials[column]: value
        for column, value in base_t.items()
    }
    if args.compare_constraint_modules:
        if tuple(functional_support) != SUPPORT:
            parser.error(
                "--compare-constraint-modules requires the standard "
                "16-term support"
            )
        support_index = {
            monomial: index for index, monomial in enumerate(SUPPORT)
        }

        def restrict_to_support(poly):
            return {
                support_index[monomial]: value
                for monomial, value in poly.items()
                if monomial in support_index and value
            }

        module_constraints = []
        module_types = []
        for monomials, left in (
            (s4_monomials, True),
            (t4_monomials, False),
        ):
            for monomial in monomials:
                image = (
                    profile.poisson({monomial: field.one}, T)
                    if left
                    else profile.poisson(S, {monomial: field.one})
                )
                row = restrict_to_support(image)
                if row:
                    module_constraints.append(row)
                    module_types.append("d5")

        # Store the bracket of every ambient S2/T2 monomial as a sparse
        # vector in the 16 functional coordinates.  This traverses the
        # ambient monomial pairs once; the 16 scalar induced matrices below
        # then reuse it.
        bracket_by_coordinate = [
            defaultdict(list) for _ in SUPPORT
        ]
        for s_index, left in enumerate(s2_monomials):
            i, j, k = left
            for t_index, right in enumerate(t2_monomials):
                ii, jj, kk = right
                z_degree = k + kk - 1
                if z_degree < 0:
                    continue
                shifted = (i + ii + 1, j + jj, z_degree)
                shifted_value = (
                    k * (3 * ii - 6 * jj)
                    - (3 * i - 6 * j) * kk
                )
                coordinate = support_index.get(shifted)
                if shifted_value and coordinate is not None:
                    bracket_by_coordinate[coordinate][s_index].append(
                        (t_index, field(shifted_value))
                    )
                unshifted = (i + ii, j + jj - 1, z_degree)
                unshifted_value = 2 * (k * jj - j * kk)
                coordinate = support_index.get(unshifted)
                if unshifted_value and coordinate is not None:
                    bracket_by_coordinate[coordinate][s_index].append(
                        (t_index, field(unshifted_value))
                    )

        moyal_vectors = []
        for monomial in s2_monomials:
            moyal_vectors.append(
                restrict_to_support(
                    scale(
                        profile.pi_power(
                            {monomial: field.one},
                            T,
                            3,
                        ),
                        field.one / field(24),
                    )
                )
            )
        for monomial in t2_monomials:
            moyal_vectors.append(
                restrict_to_support(
                    scale(
                        profile.pi_power(
                            S,
                            {monomial: field.one},
                            3,
                        ),
                        field.one / field(24),
                    )
                )
            )

        linear_rows = [defaultdict(lambda: field.zero) for _ in kernel]
        quadratic_rows = {}
        for coordinate, coordinate_bracket in enumerate(
            bracket_by_coordinate
        ):
            coordinate_left_images = []
            for s_part, _ in kernel_pairs:
                image = defaultdict(lambda: field.zero)
                for s_index, s_value in s_part.items():
                    for t_index, matrix_value in coordinate_bracket.get(
                        s_index,
                        (),
                    ):
                        image[t_index] += s_value * matrix_value
                coordinate_left_images.append(
                    {
                        t_index: value
                        for t_index, value in image.items()
                        if value
                    }
                )
            coordinate_matrix = [
                defaultdict(lambda: field.zero)
                for _ in kernel
            ]
            for left_index, image in enumerate(coordinate_left_images):
                for t_index, left_value in image.items():
                    for right_index, right_value in t_incidence.get(
                        t_index,
                        (),
                    ):
                        coordinate_matrix[left_index][right_index] += (
                            left_value * right_value
                        )
            candidate_pairs = {
                (
                    min(left_index, right_index),
                    max(left_index, right_index),
                )
                for left_index, row in enumerate(coordinate_matrix)
                for right_index, value in row.items()
                if value
            }
            for left_index, right_index in candidate_pairs:
                value = coordinate_matrix[left_index].get(
                    right_index,
                    field.zero,
                )
                if left_index != right_index:
                    value += coordinate_matrix[right_index].get(
                        left_index,
                        field.zero,
                    )
                if value:
                    quadratic_rows.setdefault(
                        (left_index, right_index),
                        {},
                    )[coordinate] = value

            for basis_index, (
                vector,
                (s_part, t_part),
            ) in enumerate(zip(kernel, kernel_pairs, strict=True)):
                value = sum(
                    moyal_vectors[column].get(
                        coordinate,
                        field.zero,
                    )
                    * coefficient
                    for column, coefficient in vector.items()
                )
                value += sum(
                    s_value
                    * matrix_value
                    * t_part.get(t_index, field.zero)
                    for s_index, s_value in base_s.items()
                    for t_index, matrix_value in coordinate_bracket.get(
                        s_index,
                        (),
                    )
                )
                value += sum(
                    s_value
                    * matrix_value
                    * base_t.get(t_index, field.zero)
                    for s_index, s_value in s_part.items()
                    for t_index, matrix_value in coordinate_bracket.get(
                        s_index,
                        (),
                    )
                )
                if value:
                    linear_rows[basis_index][coordinate] = value

        nonzero_linear_rows = [
            dict(row) for row in linear_rows if row
        ]
        nonzero_quadratic_rows = [
            row for _, row in sorted(quadratic_rows.items())
            if row
        ]
        module_constraints.extend(nonzero_linear_rows)
        module_types.extend(["linear"] * len(nonzero_linear_rows))
        module_constraints.extend(nonzero_quadratic_rows)
        module_types.extend(["quadratic"] * len(nonzero_quadratic_rows))
        module_rows = {
            index: dict(row)
            for index, row in enumerate(module_constraints)
        }
        reduced_module, module_pivots, module_nonzero = sdm_irref(
            module_rows
        )
        module_kernel, _ = sdm_nullspace_from_rref(
            reduced_module,
            field.one,
            len(SUPPORT),
            module_pivots,
            module_nonzero,
        )
        type_counts = {
            kind: module_types.count(kind)
            for kind in ("d5", "linear", "quadratic")
        }
        print(
            "LAURENT_CONSTRAINT_MODULE="
            f"conditions{len(module_constraints)},"
            f"rank{len(module_pivots)},"
            f"kernel{len(module_kernel)},"
            f"types{type_counts}",
        )
        if (
            len(module_constraints) != 41
            or len(module_pivots) != 15
            or len(module_kernel) != 1
            or type_counts
            != {"d5": 14, "linear": 3, "quadratic": 24}
        ):
            raise AssertionError(
                (
                    len(module_constraints),
                    len(module_pivots),
                    len(module_kernel),
                    type_counts,
                )
            )
        normalized_module = {
            coordinate: value / module_kernel[0][0]
            for coordinate, value in module_kernel[0].items()
        }
        for coordinate, monomial in enumerate(SUPPORT):
            expected = coefficient_functional.get(
                monomial,
                field.zero,
            )
            if normalized_module.get(coordinate, field.zero) != expected:
                raise AssertionError(
                    "Laurent and bounded period lines differ"
                )
        print("LAURENT_CONSTRAINT_MODULE_EQUALS_BOUNDED=1")

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
    if expected_obstruction is not None:
        assert obstruction_value == expected_obstruction

    print("PASS: explicit 16-term Lambda annihilates d5 on the full X-Laurent space")
    print("PASS: projected hbar^3 affine space has dimension 1075 over Q")
    print("PASS: Lambda(O5) is independent of every Laurent hbar^3 solution")
    print(f"LAURENT_BRACKET_MATRIX_NONZERO={bracket_nonzero}")
    print(f"LAURENT_QUADRATIC_MATRIX_NONZERO={quadratic_nonzero}")
    if args.profile_active_directions:
        active_quadratic_directions = sorted(
            {
                index
                for left_index, row in enumerate(quadratic_matrix)
                for index in (
                    [left_index, *row]
                    if row
                    else []
                )
            }
        )
        active_directions = sorted(
            set(active_linear_directions)
            | set(active_quadratic_directions)
        )
        print(
            "LAURENT_ACTIVE_LINEAR_DIRECTIONS="
            + ",".join(map(str, active_linear_directions))
        )
        print(
            "LAURENT_ACTIVE_QUADRATIC_DIRECTIONS="
            + ",".join(map(str, active_quadratic_directions))
        )
        print(
            "LAURENT_ACTIVE_DIRECTION_COUNT="
            + str(len(active_directions))
        )
        for basis_index in active_directions:
            s_part, t_part = kernel_pairs[basis_index]
            x_degrees = [
                monomials[column][0]
                for monomials, part in (
                    (s2_monomials, s_part),
                    (t2_monomials, t_part),
                )
                for column in part
            ]
            print(
                f"LAURENT_DIRECTION={basis_index},"
                f"S_TERMS={len(s_part)},T_TERMS={len(t_part)},"
                f"X_RANGE={min(x_degrees)}:{max(x_degrees)}"
            )
    print(f"LAURENT_PERIOD_VALUE={obstruction_value}")
    if obstruction_value:
        print(f"PASS: Lambda(O5)={obstruction_value} is nonzero")
        print(
            f"CONCLUSION: at {seed_label}, the hbar^5 obstruction "
            "does not vanish after localization at X"
        )
    else:
        print(
            f"CONCLUSION: at {seed_label}, this Laurent period vanishes; "
            "the obstruction test is inconclusive"
        )


if __name__ == "__main__":
    main()

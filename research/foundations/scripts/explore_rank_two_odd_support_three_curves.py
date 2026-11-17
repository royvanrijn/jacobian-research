#!/usr/bin/env python3
"""Collapse the positive-dimensional support-three charts at third order.

The third differential is fixed.  We quotient it out once, project the 42
second-order kernel couplings and the affine third-order right side to its
cokernel, and then perform exact row reduction on each coordinate plane.
This gives a linear necessary subspace for every line or conic chart of the
quadratic Maurer--Cartan locus.
"""

from __future__ import annotations

from collections import Counter

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import (
    sdm_irref,
    sdm_nullspace_from_rref,
)

from explore_degree_five_a2_subprincipal import add, pi_power, scale
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    essential_problem,
)
from explore_rank_two_odd_support_three import exact_support_three_audit
from verify_rank_two_odd_support_three_points import lower_lift_data


def project_to_third_cokernel():
    """Return the 38x42 coupling table and 38 affine classes modulo d_3."""

    S, T, _, pairs, _, _ = essential_problem()
    parity_pair, kernel_pairs, third_columns = lower_lift_data(S, T)

    basis_couplings = []
    basis_right_sides = []
    for pair in pairs:
        basis_couplings.append(
            [
                coupling(pair, kernel_pair)
                for kernel_pair in kernel_pairs
            ]
        )
        rhs = scale(
            add(
                pi_power(pair[0], T, 3),
                pi_power(S, pair[1], 3),
            ),
            -QQ.one / QQ(24),
        )
        basis_right_sides.append(
            add(rhs, coupling(pair, parity_pair), -QQ.one)
        )

    all_polys = list(third_columns)
    all_polys.extend(
        column
        for row in basis_couplings
        for column in row
    )
    all_polys.extend(basis_right_sides)
    monomials = sorted(set().union(*(set(poly) for poly in all_polys)))
    monomial_index = {
        monomial: index for index, monomial in enumerate(monomials)
    }

    # The rows below are the transpose of the third differential.  Its
    # nullspace is the exact dual cokernel basis.
    transpose_rows = {
        column_index: {
            monomial_index[monomial]: coefficient
            for monomial, coefficient in column.items()
        }
        for column_index, column in enumerate(third_columns)
        if column
    }
    reduced, pivots, nonzero = sdm_irref(transpose_rows)
    functionals, _ = sdm_nullspace_from_rref(
        reduced,
        QQ.one,
        len(monomials),
        pivots,
        nonzero,
    )
    assert len(pivots) == 1034
    assert len(functionals) == 2295

    incidence = {}
    for functional_index, functional in enumerate(functionals):
        for coordinate, coefficient in functional.items():
            incidence.setdefault(
                monomials[coordinate],
                [],
            ).append((functional_index, coefficient))

    def project(poly):
        result = {}
        for monomial, coefficient in poly.items():
            for functional_index, weight in incidence.get(monomial, ()):
                result[functional_index] = (
                    result.get(functional_index, QQ.zero)
                    + weight * coefficient
                )
        return {
            coordinate: coefficient
            for coordinate, coefficient in result.items()
            if coefficient
        }

    projected_couplings = [
        [project(column) for column in row]
        for row in basis_couplings
    ]
    projected_right_sides = [
        project(rhs) for rhs in basis_right_sides
    ]
    used_coordinates = set().union(
        *(
            set(column)
            for row in projected_couplings
            for column in row
        ),
        *(set(rhs) for rhs in projected_right_sides),
    )
    assert len(used_coordinates) == 920
    return projected_couplings, projected_right_sides, project


def necessary_direction_basis(
    triple,
    projected_couplings,
    projected_right_sides,
):
    """Return the exact local-coordinate kernel of the relaxed hbar^3 map."""

    columns = [
        projected_couplings[index][kernel_index]
        for index in triple
        for kernel_index in range(42)
        if projected_couplings[index][kernel_index]
    ]
    first_rhs = len(columns)
    all_columns = columns + [
        projected_right_sides[index] for index in triple
    ]
    coordinates = sorted(
        set().union(*(set(column) for column in all_columns))
    )
    rows = {
        row_index: {
            column_index: value
            for column_index, column in enumerate(all_columns)
            if (value := column.get(coordinate, QQ.zero))
        }
        for row_index, coordinate in enumerate(coordinates)
    }
    reduced, pivots, _ = sdm_irref(rows)
    rhs_pivots = {
        pivot for pivot in pivots if pivot >= first_rhs
    }
    free_rhs = [
        column
        for column in range(first_rhs, first_rhs + 3)
        if column not in rhs_pivots
    ]
    basis = []
    for free_column in free_rhs:
        vector = {free_column - first_rhs: QQ.one}
        for reduced_row, pivot in enumerate(pivots):
            if pivot < first_rhs:
                continue
            coefficient = reduced.get(reduced_row, {}).get(
                free_column,
                QQ.zero,
            )
            if coefficient:
                vector[pivot - first_rhs] = -coefficient
        basis.append(vector)
    return basis


def curve_polynomial(equations):
    """Return the primitive common line or conic in the affine chart."""

    r_symbol, s_symbol = sp.symbols("r s")
    common = sp.Poly(equations[0], s_symbol, r_symbol, domain=sp.QQ)
    for equation in equations[1:]:
        common = sp.gcd(
            common,
            sp.Poly(equation, s_symbol, r_symbol, domain=sp.QQ),
        )
    common = common.primitive()[1]
    assert common.total_degree() in (1, 2)
    return common


def homogenize_curve(polynomial):
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    r_symbol, s_symbol = sp.symbols("r s")
    degree = polynomial.total_degree()
    expression = sp.S.Zero
    for (s_degree, r_degree), coefficient in polynomial.terms():
        expression += (
            coefficient
            * x0 ** (degree - r_degree - s_degree)
            * x1 ** r_degree
            * x2 ** s_degree
        )
    assert sp.expand(
        expression.subs({x0: 1, x1: r_symbol, x2: s_symbol})
        - polynomial.as_expr()
    ) == 0
    return sp.Poly(expression, x0, x1, x2, domain=sp.QQ)


def vector_coordinates(vector):
    return tuple(vector.get(index, QQ.zero) for index in range(3))


def finite_candidates(triple, equations, basis):
    """Intersect a necessary projective subspace with the exact curve."""

    homogeneous = homogenize_curve(curve_polynomial(equations))
    x0, x1, x2 = homogeneous.gens
    if not basis:
        return [], False
    if len(basis) == 1:
        coordinates = vector_coordinates(basis[0])
        if any(not coordinate for coordinate in coordinates):
            return [], False
        value = homogeneous.as_expr().subs(
            dict(zip((x0, x1, x2), coordinates))
        )
        if value:
            return [], False
        return [
            {
                "triple": triple,
                "kind": "rational",
                "factor": sp.S.One,
                "coordinates": coordinates,
            }
        ], False

    assert len(basis) == 2
    parameter = sp.Symbol("t")
    first = vector_coordinates(basis[0])
    second = vector_coordinates(basis[1])
    coordinates = tuple(
        first[index] + parameter * second[index]
        for index in range(3)
    )
    intersection = sp.Poly(
        sp.expand(
            homogeneous.as_expr().subs(
                dict(zip((x0, x1, x2), coordinates))
            )
        ),
        parameter,
        domain=sp.QQ,
    )
    if intersection.is_zero:
        return [], True

    candidates = []
    for factor, _multiplicity in sp.factor_list(intersection)[1]:
        # A factor supported on a coordinate hyperplane belongs to a lower
        # support chart and is excluded from the exact-support-three torus.
        if any(
            sp.rem(
                sp.Poly(coordinate, parameter, domain=sp.QQ),
                factor,
            ).is_zero
            for coordinate in coordinates
        ):
            continue
        candidates.append(
            {
                "triple": triple,
                "kind": "finite",
                "factor": factor.monic().as_expr(),
                "coordinates": coordinates,
            }
        )

    # The parameter point at infinity is the second basis vector.
    if all(second):
        value = homogeneous.as_expr().subs(
            dict(zip((x0, x1, x2), second))
        )
        if value == 0:
            candidates.append(
                {
                    "triple": triple,
                    "kind": "infinity",
                    "factor": sp.S.One,
                    "coordinates": second,
                }
            )
    return candidates, False


def support_three_curve_candidates():
    projected_couplings, projected_right_sides, _ = (
        project_to_third_cokernel()
    )
    _, _, _, curves, _ = exact_support_three_audit()
    dimension_profile = Counter()
    candidates = []
    contained_necessary_lines = []
    for triple, equations in curves.items():
        basis = necessary_direction_basis(
            triple,
            projected_couplings,
            projected_right_sides,
        )
        dimension_profile[len(basis)] += 1
        local_candidates, contained = finite_candidates(
            triple,
            equations,
            basis,
        )
        candidates.extend(local_candidates)
        if contained:
            contained_necessary_lines.append(triple)
    return dimension_profile, candidates, contained_necessary_lines


def main() -> None:
    (
        dimension_profile,
        candidates,
        contained_necessary_lines,
    ) = support_three_curve_candidates()
    assert dimension_profile == Counter({0: 18, 1: 110, 2: 21})
    assert not contained_necessary_lines

    degree_profile = Counter()
    for candidate in candidates:
        if candidate["kind"] == "finite":
            degree = sp.Poly(
                candidate["factor"],
                sp.Symbol("t"),
                domain=sp.QQ,
            ).degree()
        else:
            degree = 1
        degree_profile[degree] += 1

    print(
        "PASS: exact relaxed hbar^3 necessary dimensions on the 149 "
        f"line/conic charts are {dict(sorted(dimension_profile.items()))}"
    )
    print(
        "PASS: no necessary projective line is contained in its quadratic "
        "curve"
    )
    print(
        "remaining exact-support-three closed-point candidates="
        f"{len(candidates)}"
    )
    print(f"candidate residue-degree profile={sorted(degree_profile.items())}")
    print(
        "candidate triples="
        f"{[candidate['triple'] for candidate in candidates]}"
    )


if __name__ == "__main__":
    main()

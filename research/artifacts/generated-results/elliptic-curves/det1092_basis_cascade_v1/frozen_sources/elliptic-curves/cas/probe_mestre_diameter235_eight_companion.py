#!/usr/bin/env python3
"""Local continuation and quotient audit at the new diameter-235 Mestre seed.

The all-prime post-census reconstruction finds eight affine companions at
the normalized roots ``(0,17,136,161,207,235)``.  This script follows the
transverse pair

    x1=-68/27-37*T/27,  x2=3067/459-5*T/27

without expanding a two-section residual.  It records exact tangent ranks,
full seven-row Hensel lifts, a characteristic-zero formal c1-slice, a finite
pair intersection, and mod-3 quotient ranks at small fibres.  It is local
evidence and a target-selection certificate, not a component reconstruction
or generic rank claim.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from mestre_root_tuples import SixRootMestreConstruction
from probe_mestre_two_section_local_continuation import (
    Field,
    FormalSeries,
    VARIABLES,
    rational_mod,
    residuals,
    residuals_from_jets,
    row_reduce,
    small_pade_models,
    solve_square,
    solve_square_over_q,
)
from screen_mestre_fermigier_two_section_escape import square_root
from screen_mestre_two_section_transverse_seeds import normalized_moduli
from search_mestre_root_tuple_scale import (
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
ROOTS = (0, 17, 136, 161, 207, 235)
MODULI = normalized_moduli(ROOTS)
SECTIONS = ((-Q(68, 27), -Q(37, 27)), (Q(3067, 459), -Q(5, 27)))
ALL_SECTIONS = (
    (-Q(68, 27), -Q(37, 27)),
    (-Q(68, 27), Q(37, 27)),
    (Q(3067, 459), -Q(5, 27)),
    (Q(3067, 459), Q(5, 27)),
    (Q(3947, 459), -Q(23, 27)),
    (Q(3947, 459), Q(23, 27)),
    (Q(5405, 459), -Q(17, 27)),
    (Q(5405, 459), Q(17, 27)),
)
HENSEL_PRIMES = (29, 41)
QUOTIENT_PARAMETERS = (Q(1), Q(2), Q(3), -Q(1))
# The order-32 root-chart recognition candidate in
# U=r6-235/17 and V=r3-8.  It is intentionally retained as a compact sparse
# polynomial rather than obtained by eliminating the incidence residuals.
ROOT_PLANE_QUARTIC = (
    (0, 1, 4_625_928), (0, 2, 17_506_872), (0, 3, -3_183_624),
    (0, 4, 353_736), (1, 0, -2_531_232), (1, 1, -17_221_340),
    (1, 2, 7_345_224), (1, 3, -1_002_252), (2, 0, 4_218_244),
    (2, 1, -4_959_818), (2, 2, 1_071_034), (3, 0, 1_041_556),
    (3, 1, -496_213), (4, 0, 83_521),
)


def rational_text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def truncated_product(left: list[Q], right: list[Q]) -> list[Q]:
    """Product of two univariate formal jets at the shared precision."""

    output = [Q(0)] * len(left)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right[: len(left) - left_index]):
            output[left_index + right_index] += left_value * right_value
    return output


def minimal_pade_models(
    coordinate_series: dict[str, list[str]], max_degree: int
) -> list[dict[str, object]]:
    """Return the smallest verified rational Padé model for each supplied jet."""

    models = []
    for name, rendered in coordinate_series.items():
        series = list(map(Q, rendered))
        found = None
        for total_degree in range(1, 2 * max_degree + 1):
            for numerator_degree in range(max_degree + 1):
                denominator_degree = total_degree - numerator_degree
                if not 1 <= denominator_degree <= max_degree:
                    continue
                equations = []
                right = []
                for degree in range(numerator_degree + 1, len(series)):
                    equations.append(
                        [
                            series[degree - offset] if degree >= offset else Q(0)
                            for offset in range(1, denominator_degree + 1)
                        ]
                    )
                    right.append(-series[degree])
                if len(equations) <= denominator_degree:
                    continue
                try:
                    denominator_tail = solve_square_over_q(
                        equations[:denominator_degree], right[:denominator_degree]
                    )
                except AssertionError:
                    continue
                if not all(
                    sum(left * right for left, right in zip(row, denominator_tail)) == target
                    for row, target in zip(equations, right)
                ):
                    continue
                numerator = [
                    series[degree]
                    + sum(
                        denominator_tail[offset - 1] * series[degree - offset]
                        for offset in range(1, min(degree, denominator_degree) + 1)
                    )
                    for degree in range(numerator_degree + 1)
                ]
                found = {
                    "coordinate": name,
                    "numerator_coefficients": [rational_text(value) for value in numerator],
                    "denominator_coefficients": ["1"]
                    + [rational_text(value) for value in denominator_tail],
                }
                break
            if found is not None:
                break
        if found is not None:
            models.append(found)
    return models


def low_bidegree_plane_search(
    coefficients: list[list[Q]], *, first: int, second: int, max_first: int = 8,
    max_second: int = 8, names: tuple[str, str] | None = None,
) -> dict[str, object]:
    """Look for formal plane relations without eliminating the incidence ideal.

    A displayed hit only means that the selected bidegree polynomial vanishes
    through the available jet order.  The absence of a hit is likewise a
    bounded formal exclusion, rather than a component theorem.
    """

    precision = len(coefficients[first])
    first_powers = [[Q(1)] + [Q(0)] * (precision - 1)]
    second_powers = [[Q(1)] + [Q(0)] * (precision - 1)]
    for _ in range(max_first):
        first_powers.append(truncated_product(first_powers[-1], coefficients[first]))
    for _ in range(max_second):
        second_powers.append(truncated_product(second_powers[-1], coefficients[second]))
    tested = []
    hits = []
    for first_degree in range(1, max_first + 1):
        for second_degree in range(1, max_second + 1):
            term_count = (first_degree + 1) * (second_degree + 1)
            if term_count > precision:
                continue
            columns = [
                truncated_product(first_powers[i], second_powers[j])
                for i in range(first_degree + 1)
                for j in range(second_degree + 1)
            ]
            rank, _ = row_reduce(
                [[column[degree] for column in columns] for degree in range(precision)],
                Field(),
            )
            record = {
                "first_degree": first_degree,
                "second_degree": second_degree,
                "term_count": term_count,
            }
            tested.append(record)
            if rank != term_count:
                hits.append(record)
    result = {
        "coordinates": list(names or (VARIABLES[first], VARIABLES[second])),
        "jet_order": precision - 1,
        "tested_bidegrees": tested,
        "formal_relation_hits": hits,
        "scope_limit": "a hit is only a jet relation; no hit is only a bounded low-bidegree exclusion",
    }


def split_root_jets(coefficients: list[list[Q]]) -> dict[str, list[Q]]:
    """Hensel-lift the four nonfixed normalized roots from the moduli jets."""

    order = len(coefficients[0]) - 1
    field = Field(series_order=order)
    c1, c2, c3, c4 = (FormalSeries.seed(field, row) for row in coefficients[:4])
    roots = {}
    for index, initial_root in enumerate(ROOTS[2:], start=3):
        root_coefficients = [Q(initial_root, ROOTS[1])] + [Q(0)] * order
        derivative = (
            4 * root_coefficients[0] ** 3
            + 3 * c1.value * root_coefficients[0] ** 2
            + 2 * c2.value * root_coefficients[0]
            + c3.value
        )
        if derivative == 0:
            raise AssertionError("a normalized seed root is not simple")
        for degree in range(1, order + 1):
            root = FormalSeries.seed(field, root_coefficients)
            value = root**4 + c1 * root**3 + c2 * root**2 + c3 * root + c4
            root_coefficients[degree] = -value.coefficients[degree] / derivative
            checked = FormalSeries.seed(field, root_coefficients)
            residual = checked**4 + c1 * checked**3 + c2 * checked**2 + c3 * checked + c4
            if residual.coefficients[degree]:
                raise AssertionError("split-root formal Hensel lift failed")
        roots[f"r{index}"] = root_coefficients
    return roots


def root_plane_quartic_residual(roots: dict[str, list[Q]]) -> list[Q]:
    """Evaluate the recognized quartic in the split-root formal chart."""

    order = len(roots["r3"]) - 1
    field = Field(series_order=order)
    u = FormalSeries.seed(field, roots["r6"]) - Q(235, 17)
    v = FormalSeries.seed(field, roots["r3"]) - Q(8)
    residual = FormalSeries.constant(field, 0)
    for u_degree, v_degree, coefficient in ROOT_PLANE_QUARTIC:
        residual += coefficient * u**u_degree * v**v_degree
    return list(residual.coefficients)


def root_plane_parameter_u(parameter: FormalSeries) -> FormalSeries:
    """The rational U-coordinate on the nodal-quartic normalization.

    The local seed is at ``parameter=-294``.  This is obtained by projecting
    the recognized plane quartic from its node
    ``(U,V)=(-235/17,-47/6)`` and rationalizing the resulting conic.
    """

    numerator = -45 * (parameter + 294) * (
        parameter**3 - 118 * parameter**2 - 2292 * parameter - 57416
    )
    denominator = 34 * (parameter - 66) * (parameter + 54) * (
        parameter**2 + 18 * parameter + 456
    )
    return numerator / denominator


def formal_root_plane_parameterization(order: int, pade_degree: int) -> dict[str, object]:
    """Re-solve the incidence branch using the quartic normalization parameter."""

    seed = (*MODULI, *SECTIONS[0], *SECTIONS[1])
    jacobian = jacobian_at(Field())
    columns = list(range(1, len(VARIABLES)))
    minor = [[row[index] for index in columns] for row in jacobian]
    if row_reduce(minor, Field())[0] != 7:
        raise AssertionError("c1 is not a characteristic-zero transverse coordinate")
    field = Field(series_order=order)
    parameter = FormalSeries.seed(field, [Q(-294), Q(1)] + [Q(0)] * (order - 1))
    c1 = -Q(739, 17) - 3 * root_plane_parameter_u(parameter)
    if c1.value != MODULI[0]:
        raise AssertionError("the root-plane parameter missed the base moduli point")
    coefficients = [[Q(0)] * (order + 1) for _ in VARIABLES]
    for index, value in enumerate(seed):
        coefficients[index][0] = value
    coefficients[0] = list(c1.coefficients)
    for degree in range(1, order + 1):
        values = residuals_from_jets([FormalSeries.seed(field, row) for row in coefficients])
        correction = solve_square_over_q(
            minor, [-value.coefficients[degree] for value in values]
        )
        for column, value in zip(columns, correction):
            coefficients[column][degree] = value
        checked = residuals_from_jets([FormalSeries.seed(field, row) for row in coefficients])
        if any(value.coefficients[degree] for value in checked):
            raise AssertionError("the root-plane parameter recursive solve failed")
    coordinate_series = {
        name: [str(value) for value in row]
        for name, row in zip(VARIABLES, coefficients)
    }
    return {
        "parameter": "p=-294+s",
        "c1_parameterization": "c1=-739/17-3*U(p), with U on the recognized nodal quartic",
        "order": order,
        "all_seven_residuals_vanish_through_order": True,
        "coordinate_series": coordinate_series,
        "recognized_low_degree_rational_series": small_pade_models(
            coordinate_series, max_degree=pade_degree
        ),
        "minimal_rational_series_models": minimal_pade_models(
            coordinate_series, max_degree=pade_degree
        ),
        "scope_limit": "formal recognition is not a global residual identity certificate",
    }


def jacobian_at(field: Field) -> list[list[int | Q]]:
    seed = (*MODULI, *SECTIONS[0], *SECTIONS[1])
    return [list(value.gradient) for value in residuals(seed, field)]


def full_hensel_lift(prime: int, precision: int) -> dict[str, object]:
    """Move c1 by p and lift all seven rows through p**precision."""

    seed = (*MODULI, *SECTIONS[0], *SECTIONS[1])
    coordinates = [rational_mod(value, prime) for value in seed]
    jacobian = jacobian_at(Field(prime))
    free = 0
    columns = [index for index in range(len(VARIABLES)) if index != free]
    minor = [[row[index] for index in columns] for row in jacobian]
    if row_reduce(minor, Field(prime))[0] != 7:
        raise AssertionError(f"c1 is not a transverse coordinate modulo {prime}")
    target = seed[free] + prime
    for exponent in range(1, precision):
        modulus = prime ** (exponent + 1)
        trial = list(coordinates)
        trial[free] = rational_mod(target, modulus)
        values = residuals(trial, Field(modulus, tangent=False))
        right = [(-int(value.value) // prime**exponent) % prime for value in values]
        correction = solve_square(minor, right, prime)
        for column, digit in zip(columns, correction):
            trial[column] = (trial[column] + prime**exponent * digit) % modulus
        coordinates = trial
    modulus = prime**precision
    if any(value.value for value in residuals(coordinates, Field(modulus, tangent=False))):
        raise AssertionError("the full rank-seven Hensel lift failed")
    return {
        "prime": prime,
        "precision": precision,
        "modulus": modulus,
        "free_coordinate": "c1",
        "all_seven_residuals_zero_mod_prime_power": True,
    }


def formal_c1_slice(
    order: int, *, include_coordinates: bool = False, pade_degree: int = 0,
    plane_relation_search: bool = False, root_chart_search: bool = False,
    root_plane_candidate: bool = False,
) -> dict[str, object]:
    """Solve all seven rows in Q[[t]] with c1=c1(0)+t."""

    seed = (*MODULI, *SECTIONS[0], *SECTIONS[1])
    jacobian = jacobian_at(Field())
    columns = list(range(1, len(VARIABLES)))
    minor = [[row[index] for index in columns] for row in jacobian]
    if row_reduce(minor, Field())[0] != 7:
        raise AssertionError("c1 is not a characteristic-zero transverse coordinate")
    coefficients = [[Q(0)] * (order + 1) for _ in VARIABLES]
    for index, value in enumerate(seed):
        coefficients[index][0] = value
    coefficients[0][1] = Q(1)
    field = Field(series_order=order)
    for degree in range(1, order + 1):
        values = residuals_from_jets([FormalSeries.seed(field, row) for row in coefficients])
        correction = solve_square_over_q(
            minor, [-value.coefficients[degree] for value in values]
        )
        for column, value in zip(columns, correction):
            coefficients[column][degree] = value
        checked = residuals_from_jets([FormalSeries.seed(field, row) for row in coefficients])
        if any(value.coefficients[degree] for value in checked):
            raise AssertionError("the formal rank-seven implicit solve failed")
    c2_expected = [MODULI[1], -Q(1621, 51), Q(13, 36)] + [Q(0)] * (order - 2)
    if coefficients[1] != c2_expected:
        raise AssertionError("the recognized c2(c1) jet changed")
    coordinate_series = {
        name: [str(value) for value in row]
        for name, row in zip(VARIABLES, coefficients)
    }
    result: dict[str, object] = {
        "parameter": "c1=-739/17+t",
        "order": order,
        "all_seven_residuals_vanish_through_order": True,
        "recognized_c2_jet": "c2=201815/289-(1621/51)*t+(13/36)*t^2",
    }
    if include_coordinates:
        result["coordinate_series"] = coordinate_series
    if pade_degree:
        result["recognized_low_degree_rational_series"] = small_pade_models(
            coordinate_series, max_degree=pade_degree
        )
    if plane_relation_search:
        result["low_bidegree_c3_c4_search"] = low_bidegree_plane_search(
            coefficients, first=2, second=3
        )
    if root_chart_search or root_plane_candidate:
        roots = split_root_jets(coefficients)
        root_plane_residual = root_plane_quartic_residual(roots)
        if any(root_plane_residual):
            raise AssertionError("the recognized split-root plane quartic failed")
    if root_plane_candidate:
        result["recognized_split_root_plane_quartic"] = {
            "variables": "U=r6-235/17, V=r3-8",
            "terms": [list(term) for term in ROOT_PLANE_QUARTIC],
            "vanishes_through_order": order,
            "scope_limit": "a finite formal-jet match is not a global component identity",
        }
    if root_chart_search:
        root_series = {name: [str(value) for value in row] for name, row in roots.items()}
        r6_expected = [Q(235, 17), -Q(1, 3)] + [Q(0)] * (order - 1)
        if roots["r6"] != r6_expected:
            raise AssertionError("the recognized r6(c1) jet changed")
        root_names = tuple(roots)
        root_searches = []
        for first in range(len(root_names)):
            for second in range(first + 1, len(root_names)):
                pair_coefficients = [roots[root_names[first]], roots[root_names[second]]]
                root_searches.append(
                    low_bidegree_plane_search(
                        pair_coefficients,
                        first=0,
                        second=1,
                        names=(root_names[first], root_names[second]),
                    )
                )
        result["split_root_chart"] = {
            "all_four_root_jets_vanish_in_the_quartic_through_order": order,
            "recognized_r6_jet": "r6=235/17-t/3, equivalently c1+3*r6=-34/17",
            "root_series": root_series,
            "other_roots_recognized_low_degree_rational_series": small_pade_models(
                {name: row for name, row in root_series.items() if name != "r6"},
                max_degree=8,
            ),
            "pair_low_bidegree_searches": root_searches,
        }
    return result


def finite_intersection() -> dict[str, object]:
    """Record the selected pair's finite affine meeting at the seed."""

    (a1, b1), (a2, b2) = SECTIONS
    parameter = (a2 - a1) / (b1 - b2)
    construction = SixRootMestreConstruction(tuple(Q(root, ROOTS[1]) for root in ROOTS))
    quartic = construction.primitive_quartic_coefficients(parameter)
    x_value = a1 + b1 * parameter
    y_value = square_root(quartic_value(quartic, x_value))
    if y_value is None or y_value == 0 or x_value != a2 + b2 * parameter:
        raise AssertionError("the selected finite intersection changed")
    return {
        "T": rational_text(parameter),
        "primitive_quartic_point": [rational_text(x_value), rational_text(y_value)],
        "scope_limit": "this is not a full section-intersection or Shioda calculation",
    }


def quotient_ranks() -> list[dict[str, object]]:
    """Compare visible and augmented exact mod-3 ranks at four fibres."""

    construction = SixRootMestreConstruction(tuple(Q(root, ROOTS[1]) for root in ROOTS))
    records = []
    for parameter in QUOTIENT_PARAMETERS:
        quartic = construction.primitive_quartic_coefficients(parameter)
        coefficients = construction.primitive_jacobian_coefficients(parameter)
        visible = tuple(
            quartic_point_to_jacobian(construction, parameter, point)
            for point in primitive_visible_points(construction, parameter)
        )
        affine = tuple(
            quartic_point_to_jacobian(
                construction,
                parameter,
                (intercept + slope * parameter, square_root(quartic_value(quartic, intercept + slope * parameter))),
            )
            for intercept, slope in ALL_SECTIONS
        )
        baseline = mod3_independence_certificate(coefficients, visible, prime_bound=251)
        augmented = mod3_independence_certificate(coefficients, (*visible, *affine), prime_bound=251)
        individual_ranks = [
            mod3_independence_certificate(
                coefficients, (*visible, companion), prime_bound=251
            )["combined_exact_rank_over_F3"]
            for companion in affine
        ]
        records.append(
            {
                "T": rational_text(parameter),
                "visible_mod3_rank": baseline["combined_exact_rank_over_F3"],
                "visible_plus_individual_companion_mod3_ranks": individual_ranks,
                "visible_plus_eight_companions_mod3_rank": augmented["combined_exact_rank_over_F3"],
            }
        )
    return records


def replay(
    *, order: int = 12, precision: int = 4, include_coordinates: bool = False,
    pade_degree: int = 0, plane_relation_search: bool = False,
    root_chart_search: bool = False, root_plane_candidate: bool = False,
    root_plane_parameterization: bool = False,
) -> dict[str, object]:
    seed = (*MODULI, *SECTIONS[0], *SECTIONS[1])
    values = residuals(seed, Field())
    if any(value.value for value in values):
        raise AssertionError("the selected pair is not an exact recursive residual zero")
    rank, pivots = row_reduce(jacobian_at(Field()), Field())
    if rank != 7:
        raise AssertionError(f"expected rank seven, obtained {rank}")
    modular = []
    for prime in HENSEL_PRIMES:
        rank_mod_prime, _ = row_reduce(jacobian_at(Field(prime)), Field(prime))
        if rank_mod_prime != 7:
            raise AssertionError(f"rank dropped modulo {prime}")
        modular.append({"prime": prime, "jacobian_rank": rank_mod_prime})
    quotient = quotient_ranks()
    if [(item["T"], item["visible_mod3_rank"], item["visible_plus_eight_companions_mod3_rank"]) for item in quotient] != [("1", 10, 10), ("2", 10, 11), ("3", 10, 11), ("-1", 10, 10)]:
        raise AssertionError("the finite-quotient target-selection profile changed")
    result = {
        "status": "diameter-235 eight-companion local continuation recorded",
        "roots": list(ROOTS),
        "selected_sections": [[rational_text(value) for value in section] for section in SECTIONS],
        "all_reconstructed_section_count": len(ALL_SECTIONS),
        "expanded_two_section_residual_materialized": False,
        "exact_jacobian_rank": rank,
        "exact_pivot_columns": [VARIABLES[index] for index in pivots],
        "finite_field_tangent_checks": modular,
        "hensel_continuations": [full_hensel_lift(prime, precision) for prime in HENSEL_PRIMES],
        "formal_c1_slice": formal_c1_slice(
            order,
            include_coordinates=include_coordinates,
            pade_degree=pade_degree,
            plane_relation_search=plane_relation_search,
            root_chart_search=root_chart_search,
            root_plane_candidate=root_plane_candidate,
        ),
        "selected_pair_finite_intersection": finite_intersection(),
        "finite_reduction_quotient_ranks": quotient,
        "conclusion": "the seed has a smooth rank-seven local two-section branch; at T=2 and T=3 every reconstructed companion separately raises the exact mod-3 quotient rank above the visible subgroup",
        "not_established": [
            "a rational parameterization or a global component identity",
            "generic independence, saturation, a Shioda Gram matrix, or a generic rank-14 claim",
        ],
    }
    if root_plane_parameterization:
        result["formal_root_plane_parameterization"] = formal_root_plane_parameterization(
            order, pade_degree or 8
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=12)
    parser.add_argument("--precision", type=int, default=4)
    parser.add_argument("--include-coordinate-series", action="store_true")
    parser.add_argument("--pade-degree", type=int, default=0)
    parser.add_argument("--plane-relation-search", action="store_true")
    parser.add_argument("--root-chart-search", action="store_true")
    parser.add_argument("--root-plane-candidate", action="store_true")
    parser.add_argument("--root-plane-parameterization", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        replay(
            order=args.order,
            precision=args.precision,
            include_coordinates=args.include_coordinate_series,
            pade_degree=args.pade_degree,
            plane_relation_search=args.plane_relation_search,
            root_chart_search=args.root_chart_search,
            root_plane_candidate=args.root_plane_candidate,
            root_plane_parameterization=args.root_plane_parameterization,
        ),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

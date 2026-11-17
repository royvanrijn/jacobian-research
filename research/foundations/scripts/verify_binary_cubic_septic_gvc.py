#!/usr/bin/env python3
"""Verify the complete binary ``(r, deg(P)) = (3, 7)`` GVC row.

The proof is in
``extended-geometry/BINARY_CUBIC_SEPTIC_GVC.md``.  This checker exhausts
the seven-factor Hall locus, constructs every rational Newton crossing from
the exact two-wing normal form, and checks the face radicals over ``QQ`` in
Singular.  The all-order promotion from the terminal common thresholds is
the written weight-defect argument, not a bounded-moment extrapolation.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

from verify_binary_quadratic_all_root_partitions_gvc import (
    all_moment_coefficients,
    compositions,
    has_matching,
    radical_equal,
)


Exponent = tuple[int, int]
ROOT_PARTITIONS = ((3,), (2, 1), (1, 1, 1))
MOMENT_BOUND = 10


def verify_hall_locus() -> None:
    """Exhaust Hall failure for three derivatives and seven factors."""

    for partition in ROOT_PARTITIONS:
        derivatives = tuple(
            direction
            for direction, multiplicity in enumerate(partition)
            for _ in range(multiplicity)
        )
        for counts in compositions(7, len(partition) + 1):
            annihilators = tuple(
                direction
                for direction, count in enumerate(counts[:-1])
                for _ in range(count)
            ) + (-1,) * counts[-1]
            expected_failure = any(
                counts[direction] >= 8 - multiplicity
                for direction, multiplicity in enumerate(partition)
            )
            assert has_matching(derivatives, annihilators) != expected_failure


def is_normal_operator_term(
    multiplicity: int, exponent: Exponent
) -> bool:
    """Return membership in the exact local two-wing normal support."""

    x_order, y_order = exponent
    total = x_order + y_order
    if total == 3:
        # The marked monomial and its strict homogeneous cofactor.
        return x_order >= multiplicity
    if total < 4:
        return False
    excess = total - 3
    return x_order < multiplicity or x_order > multiplicity + excess


def polynomial_terms(top_support: tuple[Exponent, ...]) -> tuple[Exponent, ...]:
    """The fixed degree-seven top support plus every lower monomial."""

    terms = list(top_support)
    terms.extend(
        (x_degree, y_degree)
        for x_degree in range(7)
        for y_degree in range(7 - x_degree)
    )
    return tuple(dict.fromkeys(terms))


def weight(exponent: Exponent, slope: Fraction) -> Fraction:
    return slope * exponent[0] + exponent[1]


def crossing_events(
    multiplicity: int,
    operator_base: Exponent,
    polynomial_base: Exponent,
    start: Fraction,
    stop: Fraction,
    top_support: tuple[Exponent, ...],
) -> tuple[Fraction, ...]:
    """Enumerate every crossing strictly between two Newton pivots."""

    events: set[Fraction] = set()

    # Only a term with smaller x-order can cross below the current operator
    # minimum as the slope increases.  Its y-order at a crossing is bounded
    # by the weight of the base at ``stop``.
    y_bound = int(weight(operator_base, stop)) + 1
    for x_order in range(operator_base[0]):
        for y_order in range(y_bound + 1):
            exponent = (x_order, y_order)
            if not is_normal_operator_term(multiplicity, exponent):
                continue
            slope = Fraction(
                operator_base[1] - y_order,
                x_order - operator_base[0],
            )
            if start < slope < stop:
                events.add(slope)

    # Only a term with larger x-degree can cross above the current polynomial
    # maximum.  The polynomial support is finite from the outset.
    for x_degree, y_degree in polynomial_terms(top_support):
        if x_degree <= polynomial_base[0]:
            continue
        slope = Fraction(
            polynomial_base[1] - y_degree,
            x_degree - polynomial_base[0],
        )
        if start < slope < stop:
            events.add(slope)

    return tuple(sorted(events))


def face_exponents(
    multiplicity: int,
    operator_base: Exponent,
    polynomial_base: Exponent,
    slope: Fraction,
    top_support: tuple[Exponent, ...],
) -> tuple[tuple[Exponent, ...], tuple[Exponent, ...]]:
    """Construct the complete equality supports at one rational slope."""

    operator_weight = weight(operator_base, slope)
    operator_terms = [operator_base]
    for x_order in range(int(operator_weight / slope) + 1):
        for y_order in range(int(operator_weight) + 1):
            exponent = (x_order, y_order)
            if (
                exponent != operator_base
                and is_normal_operator_term(multiplicity, exponent)
                and weight(exponent, slope) == operator_weight
            ):
                operator_terms.append(exponent)

    polynomial_weight = weight(polynomial_base, slope)
    polynomial_face = [polynomial_base]
    polynomial_face.extend(
        exponent
        for exponent in polynomial_terms(top_support)
        if exponent != polynomial_base
        and weight(exponent, slope) == polynomial_weight
    )
    return (
        tuple(dict.fromkeys(operator_terms)),
        tuple(dict.fromkeys(polynomial_face)),
    )


def face_data(
    label: str,
    multiplicity: int,
    operator_base: Exponent,
    polynomial_base: Exponent,
    slope: Fraction,
    top_support: tuple[Exponent, ...],
):
    operator_support, polynomial_support = face_exponents(
        multiplicity,
        operator_base,
        polynomial_base,
        slope,
        top_support,
    )
    suffix = str(slope).replace("/", "_")
    operator_variables = sp.symbols(
        f"{label}_{suffix}_o0:{len(operator_support) - 1}"
    )
    polynomial_variables = sp.symbols(
        f"{label}_{suffix}_p0:{len(polynomial_support) - 1}"
    )
    operator_map = dict(
        zip(operator_support[1:], operator_variables, strict=True)
    )
    polynomial_map = dict(
        zip(polynomial_support[1:], polynomial_variables, strict=True)
    )
    operator = {operator_base: sp.Integer(1), **operator_map}
    polynomial = {polynomial_base: sp.Integer(1), **polynomial_map}
    variables = (*operator_variables, *polynomial_variables)
    equations = all_moment_coefficients(
        polynomial, operator, MOMENT_BOUND
    )
    return variables, equations, operator_map, polynomial_map


def verify_scan(
    label: str,
    multiplicity: int,
    operator_base: Exponent,
    polynomial_base: Exponent,
    start: Fraction,
    stop: Fraction,
    top_support: tuple[Exponent, ...],
    expected_events: tuple[Fraction, ...],
    exceptional_radicals=None,
) -> None:
    """Audit one complete Newton interval and all of its face radicals."""

    actual_events = crossing_events(
        multiplicity,
        operator_base,
        polynomial_base,
        start,
        stop,
        top_support,
    )
    assert actual_events == expected_events, (label, actual_events)
    exceptional_radicals = exceptional_radicals or {}
    for slope in actual_events:
        variables, equations, operator_map, polynomial_map = face_data(
            label,
            multiplicity,
            operator_base,
            polynomial_base,
            slope,
            top_support,
        )
        builder = exceptional_radicals.get(slope)
        expected = (
            tuple(variables)
            if builder is None
            else tuple(builder(operator_map, polynomial_map))
        )
        radical_equal(tuple(variables), equations, expected)

    # The endpoint is a genuine common threshold.  Completeness of the
    # event list says every surviving term lies on the correct side of it.
    assert weight(operator_base, stop) == weight(polynomial_base, stop)


def triple_a_radical(operator, polynomial):
    return (
        polynomial[(4, 0)],
        polynomial[(3, 2)],
        operator[(0, 6)],
        operator[(1, 4)],
        operator[(2, 2)] * polynomial[(2, 4)],
    )


def triple_c_three_halves_radical(operator, polynomial):
    return (
        polynomial[(4, 1)],
        operator[(1, 3)] * polynomial[(2, 4)],
    )


def triple_c_two_radical(operator, polynomial):
    return (
        polynomial[(3, 1)],
        operator[(0, 6)],
        operator[(1, 4)] * polynomial[(2, 3)],
        operator[(2, 2)] * polynomial[(2, 3)],
        operator[(1, 4)] * polynomial[(1, 5)],
    )


def double_c_two_radical(operator, polynomial):
    return (
        polynomial[(3, 1)],
        polynomial[(2, 3)],
        operator[(0, 5)],
        operator[(1, 3)] * polynomial[(1, 5)],
    )


def verify_initial_charts() -> None:
    """Check all six Hall charts before a nontrivial Newton pivot."""

    F = Fraction
    verify_scan(
        "triple_B",
        3,
        (3, 0),
        (2, 5),
        F(1),
        F(5),
        ((2, 5), (0, 7)),
        tuple(
            map(
                F,
                (
                    "5/4",
                    "4/3",
                    "3/2",
                    "5/3",
                    "2",
                    "7/3",
                    "5/2",
                    "8/3",
                    "3",
                    "10/3",
                    "7/2",
                    "11/3",
                    "4",
                    "13/3",
                    "9/2",
                    "14/3",
                ),
            )
        ),
    )
    verify_scan(
        "triple_A",
        3,
        (3, 0),
        (1, 6),
        F(1),
        F(3),
        ((1, 6), (0, 7)),
        tuple(
            map(
                F,
                ("6/5", "5/4", "4/3", "3/2", "5/3", "2", "7/3", "5/2", "8/3"),
            )
        ),
        {F(2): triple_a_radical},
    )
    verify_scan(
        "triple_C",
        3,
        (3, 0),
        (0, 7),
        F(1),
        F(7, 3),
        ((0, 7),),
        tuple(
            map(
                F,
                ("7/6", "6/5", "5/4", "4/3", "7/5", "3/2", "5/3", "7/4", "2"),
            )
        ),
        {
            F(3, 2): triple_c_three_halves_radical,
            F(2): triple_c_two_radical,
        },
    )
    verify_scan(
        "double_A",
        2,
        (2, 1),
        (1, 6),
        F(1),
        F(5),
        ((1, 6), (0, 7)),
        tuple(
            map(
                F,
                ("6/5", "5/4", "4/3", "3/2", "5/3", "2", "5/2", "3", "7/2", "4", "9/2"),
            )
        ),
    )
    verify_scan(
        "double_C",
        2,
        (2, 1),
        (0, 7),
        F(1),
        F(3),
        ((0, 7),),
        tuple(
            map(
                F,
                ("7/6", "6/5", "5/4", "4/3", "7/5", "3/2", "5/3", "7/4", "2", "7/3", "5/2"),
            )
        ),
        {F(2): double_c_two_radical},
    )
    verify_scan(
        "simple",
        1,
        (1, 2),
        (0, 7),
        F(1),
        F(5),
        ((0, 7),),
        tuple(
            map(
                F,
                ("7/6", "6/5", "5/4", "4/3", "7/5", "3/2", "5/3", "7/4", "2", "7/3", "5/2", "3", "7/2", "4"),
            )
        ),
    )


def verify_pivot_tree() -> None:
    """Follow every axis of the three exceptional initial radicals."""

    F = Fraction
    scans = (
        (
            "triple_A_operator",
            3,
            (2, 2),
            (1, 6),
            F(2),
            F(4),
            ((1, 6), (0, 7)),
            (F(5, 2), F(3), F(7, 2)),
        ),
        (
            "triple_A_polynomial",
            3,
            (3, 0),
            (2, 4),
            F(2),
            F(4),
            ((1, 6), (0, 7)),
            (F(7, 3), F(5, 2), F(8, 3), F(3), F(10, 3), F(7, 2), F(11, 3)),
        ),
        (
            "triple_C_operator_13",
            3,
            (1, 3),
            (0, 7),
            F(3, 2),
            F(4),
            ((0, 7),),
            (F(5, 3), F(7, 4), F(2), F(7, 3), F(5, 2), F(3), F(7, 2)),
        ),
        (
            "triple_C_polynomial_24",
            3,
            (3, 0),
            (2, 4),
            F(3, 2),
            F(4),
            ((0, 7),),
            (F(5, 3), F(2), F(7, 3), F(5, 2), F(8, 3), F(3), F(10, 3), F(7, 2), F(11, 3)),
        ),
        (
            "triple_C_operator_22",
            3,
            (2, 2),
            (0, 7),
            F(2),
            F(5, 2),
            ((0, 7),),
            (F(7, 3),),
        ),
        (
            "triple_C_operator_14",
            3,
            (1, 4),
            (0, 7),
            F(2),
            F(3),
            ((0, 7),),
            (F(7, 3), F(5, 2)),
        ),
        (
            "triple_C_both",
            3,
            (2, 2),
            (1, 5),
            F(2),
            F(3),
            ((0, 7),),
            (F(5, 2),),
        ),
        (
            "triple_C_polynomial_15",
            3,
            (3, 0),
            (1, 5),
            F(2),
            F(5, 2),
            ((0, 7),),
            (F(7, 3),),
        ),
        (
            "triple_C_polynomial_23",
            3,
            (3, 0),
            (2, 3),
            F(2),
            F(3),
            ((0, 7),),
            (F(7, 3), F(5, 2), F(8, 3)),
        ),
        (
            "double_C_polynomial",
            2,
            (2, 1),
            (1, 5),
            F(2),
            F(4),
            ((0, 7),),
            (F(5, 2), F(3), F(7, 2)),
        ),
    )
    for scan in scans:
        verify_scan(*scan)

    # The operator axis of the double-root radical has the same exposed
    # supports as the already checked tail of the triple-root (XY^3,y^7)
    # scan; the extra X^2Y cofactor is strict for every slope > 2.
    double_events = crossing_events(
        2, (1, 3), (0, 7), F(2), F(4), ((0, 7),)
    )
    assert double_events == (F(7, 3), F(5, 2), F(3), F(7, 2))
    for slope in double_events:
        assert face_exponents(
            2, (1, 3), (0, 7), slope, ((0, 7),)
        ) == face_exponents(
            3, (1, 3), (0, 7), slope, ((0, 7),)
        )


def main() -> None:
    verify_hall_locus()
    verify_initial_charts()
    verify_pivot_tree()
    print("verified complete binary cubic-leading septic GVC row")


if __name__ == "__main__":
    main()

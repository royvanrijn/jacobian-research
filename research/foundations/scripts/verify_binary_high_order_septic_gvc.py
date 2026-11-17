#!/usr/bin/env python3
"""Verify the binary septic GVC rows with lowest order four through six.

The proof is in ``extended-geometry/BINARY_DEGREE_SEVEN_GVC.md``.  The
checker constructs every Hall chart and rational Newton crossing from the
two-wing normal form, verifies exact face radicals over ``QQ`` in Singular,
and follows every coordinate-axis primary to a common threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy as sp

from verify_binary_quadratic_all_root_partitions_gvc import (
    all_moment_coefficients,
    compositions,
    has_matching,
    radical_equal,
)


Exponent = tuple[int, int]
DEGREE = 7
MOMENT_BOUND = 10
ORDERS = (4, 5, 6)


@dataclass(frozen=True)
class AxisException:
    """The two free coordinates in a squarefree face radical ``(u*v)``."""

    operator_axis: Exponent
    polynomial_axis: Exponent


@dataclass(frozen=True)
class State:
    order: int
    multiplicity: int
    top_x_degree: int
    operator_base: Exponent
    polynomial_base: Exponent
    start: Fraction
    stop: Fraction


EXCEPTIONS: dict[tuple[int, int, int, Fraction], AxisException] = {
    # Order four: eight initial two-axis faces.
    (4, 2, 0, Fraction(2)): AxisException((1, 4), (1, 5)),
    (4, 3, 0, Fraction(3, 2)): AxisException((1, 4), (2, 4)),
    (4, 3, 1, Fraction(2)): AxisException((2, 3), (2, 4)),
    (4, 4, 0, Fraction(4, 3)): AxisException((1, 4), (3, 3)),
    (4, 4, 0, Fraction(3, 2)): AxisException((2, 3), (2, 4)),
    (4, 4, 0, Fraction(5, 3)): AxisException((1, 5), (3, 2)),
    (4, 4, 1, Fraction(3, 2)): AxisException((2, 3), (3, 3)),
    (4, 4, 2, Fraction(2)): AxisException((3, 2), (3, 3)),
    # Order five: seven initial two-axis faces.
    (5, 3, 0, Fraction(3, 2)): AxisException((1, 5), (2, 4)),
    (5, 4, 0, Fraction(4, 3)): AxisException((1, 5), (3, 3)),
    (5, 4, 1, Fraction(3, 2)): AxisException((2, 4), (3, 3)),
    (5, 5, 0, Fraction(5, 4)): AxisException((1, 5), (4, 2)),
    (5, 5, 0, Fraction(4, 3)): AxisException((2, 4), (3, 3)),
    (5, 5, 1, Fraction(4, 3)): AxisException((2, 4), (4, 2)),
    (5, 5, 2, Fraction(3, 2)): AxisException((3, 3), (4, 2)),
}


EXPECTED_CENSUS = {
    4: {
        "charts": 10,
        "initial_faces": 97,
        "exceptions": 8,
        "child_states": 16,
        "child_faces": 58,
    },
    5: {
        "charts": 15,
        "initial_faces": 112,
        "exceptions": 7,
        "child_states": 14,
        "child_faces": 40,
    },
    6: {
        "charts": 21,
        "initial_faces": 78,
        "exceptions": 0,
        "child_states": 0,
        "child_faces": 0,
    },
}


def integer_partitions(total: int, maximum: int | None = None):
    """Yield decreasing integer partitions of ``total``."""

    if total == 0:
        yield ()
        return
    maximum = total if maximum is None else min(maximum, total)
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


def verify_hall_loci() -> None:
    """Exhaust every root partition at orders four, five, and six."""

    for order in ORDERS:
        for partition in integer_partitions(order):
            derivatives = tuple(
                direction
                for direction, multiplicity in enumerate(partition)
                for _ in range(multiplicity)
            )
            for counts in compositions(DEGREE, len(partition) + 1):
                annihilators = tuple(
                    direction
                    for direction, count in enumerate(counts[:-1])
                    for _ in range(count)
                ) + (-1,) * counts[-1]
                expected_failure = any(
                    counts[direction] >= DEGREE + 1 - multiplicity
                    for direction, multiplicity in enumerate(partition)
                )
                assert (
                    has_matching(derivatives, annihilators)
                    != expected_failure
                )


def is_normal_operator_term(
    order: int, multiplicity: int, exponent: Exponent
) -> bool:
    """Membership in the exact marked-root two-wing normal support."""

    x_order, y_order = exponent
    total = x_order + y_order
    if total == order:
        return x_order >= multiplicity
    if total < order + 1:
        return False
    excess = total - order
    return (
        x_order < multiplicity
        or x_order > multiplicity + excess
    )


def polynomial_terms(top_x_degree: int) -> tuple[Exponent, ...]:
    """Return the degree-seven chart top and every lower monomial."""

    terms = [(x_degree, DEGREE - x_degree) for x_degree in range(top_x_degree + 1)]
    terms.extend(
        (x_degree, y_degree)
        for x_degree in range(DEGREE)
        for y_degree in range(DEGREE - x_degree)
    )
    return tuple(dict.fromkeys(terms))


def weight(exponent: Exponent, slope: Fraction) -> Fraction:
    return slope * exponent[0] + exponent[1]


def common_stop(
    operator_base: Exponent, polynomial_base: Exponent
) -> Fraction:
    x_gap = operator_base[0] - polynomial_base[0]
    assert x_gap > 0
    return Fraction(
        polynomial_base[1] - operator_base[1], x_gap
    )


def initial_state(
    order: int, multiplicity: int, top_x_degree: int
) -> State:
    operator_base = (multiplicity, order - multiplicity)
    polynomial_base = (top_x_degree, DEGREE - top_x_degree)
    return State(
        order,
        multiplicity,
        top_x_degree,
        operator_base,
        polynomial_base,
        Fraction(1),
        common_stop(operator_base, polynomial_base),
    )


def crossing_events(state: State) -> tuple[Fraction, ...]:
    """Enumerate every possible first crossing in one state interval."""

    events: set[Fraction] = set()
    operator_bound = int(weight(state.operator_base, state.stop)) + 1
    for x_order in range(state.operator_base[0]):
        for y_order in range(operator_bound + 1):
            exponent = (x_order, y_order)
            if not is_normal_operator_term(
                state.order, state.multiplicity, exponent
            ):
                continue
            slope = Fraction(
                state.operator_base[1] - y_order,
                x_order - state.operator_base[0],
            )
            if state.start < slope < state.stop:
                events.add(slope)

    for x_degree, y_degree in polynomial_terms(state.top_x_degree):
        if x_degree <= state.polynomial_base[0]:
            continue
        slope = Fraction(
            state.polynomial_base[1] - y_degree,
            x_degree - state.polynomial_base[0],
        )
        if state.start < slope < state.stop:
            events.add(slope)
    return tuple(sorted(events))


def face_data(label: str, state: State, slope: Fraction):
    """Construct one complete equality face and its first ten moments."""

    operator_weight = weight(state.operator_base, slope)
    operator_support = [state.operator_base]
    for x_order in range(int(operator_weight / slope) + 1):
        for y_order in range(int(operator_weight) + 1):
            exponent = (x_order, y_order)
            if (
                exponent != state.operator_base
                and is_normal_operator_term(
                    state.order, state.multiplicity, exponent
                )
                and weight(exponent, slope) == operator_weight
            ):
                operator_support.append(exponent)

    polynomial_weight = weight(state.polynomial_base, slope)
    polynomial_support = [state.polynomial_base]
    polynomial_support.extend(
        exponent
        for exponent in polynomial_terms(state.top_x_degree)
        if exponent != state.polynomial_base
        and weight(exponent, slope) == polynomial_weight
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
    operator = {state.operator_base: sp.Integer(1), **operator_map}
    polynomial = {
        state.polynomial_base: sp.Integer(1),
        **polynomial_map,
    }
    variables = (*operator_variables, *polynomial_variables)
    equations = all_moment_coefficients(
        polynomial, operator, MOMENT_BOUND
    )
    return variables, equations, operator_map, polynomial_map


def verify_face(
    label: str,
    state: State,
    slope: Fraction,
    exception: AxisException | None,
) -> None:
    variables, equations, operator_map, polynomial_map = face_data(
        label, state, slope
    )
    if exception is None:
        expected = tuple(variables)
    else:
        operator_axis = operator_map[exception.operator_axis]
        polynomial_axis = polynomial_map[exception.polynomial_axis]
        expected = tuple(
            variable
            for variable in variables
            if variable not in (operator_axis, polynomial_axis)
        ) + (operator_axis * polynomial_axis,)
    radical_equal(tuple(variables), equations, expected)


def child_states(
    state: State, slope: Fraction, exception: AxisException
) -> tuple[State, State]:
    """Generate the operator and polynomial axes of one face primary."""

    operator_stop = common_stop(
        exception.operator_axis, state.polynomial_base
    )
    polynomial_stop = common_stop(
        state.operator_base, exception.polynomial_axis
    )
    assert slope < operator_stop and slope < polynomial_stop

    # Each pivot strictly decreases the marked x-gap.  This is the
    # combinatorial no-cycle invariant behind the finite branch tree.
    old_gap = state.operator_base[0] - state.polynomial_base[0]
    assert (
        exception.operator_axis[0] - state.polynomial_base[0]
        < old_gap
    )
    assert (
        state.operator_base[0] - exception.polynomial_axis[0]
        < old_gap
    )

    return (
        State(
            state.order,
            state.multiplicity,
            state.top_x_degree,
            exception.operator_axis,
            state.polynomial_base,
            slope,
            operator_stop,
        ),
        State(
            state.order,
            state.multiplicity,
            state.top_x_degree,
            state.operator_base,
            exception.polynomial_axis,
            slope,
            polynomial_stop,
        ),
    )


def verify_order(order: int) -> None:
    census = {
        "charts": 0,
        "initial_faces": 0,
        "exceptions": 0,
        "child_states": 0,
        "child_faces": 0,
    }
    children: list[tuple[str, State]] = []

    for multiplicity in range(1, order + 1):
        for top_x_degree in range(multiplicity):
            census["charts"] += 1
            state = initial_state(order, multiplicity, top_x_degree)
            assert state.start < state.stop
            for slope in crossing_events(state):
                census["initial_faces"] += 1
                key = (order, multiplicity, top_x_degree, slope)
                exception = EXCEPTIONS.get(key)
                label = f"r{order}_e{multiplicity}_t{top_x_degree}"
                verify_face(label, state, slope, exception)
                if exception is not None:
                    census["exceptions"] += 1
                    slope_label = str(slope).replace("/", "_")
                    for branch, child in zip(
                        ("operator", "polynomial"),
                        child_states(state, slope, exception),
                        strict=True,
                    ):
                        children.append(
                            (f"{label}_{slope_label}_{branch}", child)
                        )

            assert weight(state.operator_base, state.stop) == weight(
                state.polynomial_base, state.stop
            )

    census["child_states"] = len(children)
    for label, state in children:
        for slope in crossing_events(state):
            census["child_faces"] += 1
            verify_face(label, state, slope, None)
        assert weight(state.operator_base, state.stop) == weight(
            state.polynomial_base, state.stop
        )

    assert census == EXPECTED_CENSUS[order], (order, census)


def main() -> None:
    verify_hall_loci()
    for order in ORDERS:
        verify_order(order)
    print(
        "verified septic rows r=4,5,6; prior rows complete binary GVC "
        "through polynomial degree seven"
    )


if __name__ == "__main__":
    main()

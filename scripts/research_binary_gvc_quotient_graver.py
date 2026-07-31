#!/usr/bin/env python3
"""Exact projected-scroll obstruction in the binary GVC ghost quotient.

The low-digit carry atoms have support at most three on one side or
two-plus-two across the two sides.  This does not imply that their
high-digit quotient fibers are connected by the corresponding circuit
moves after inactive channels have been projected away.

The first obstruction is the color-homogeneous partition identity

    R_3 B_1 B_2 = R_0 B_3^2.

Both monomials have one red selection, two blue selections, and total
level six.  On the projected support

    R = {0, 3},  B = {1, 2, 3},

they are the only two points of the fiber.  Their difference has support
five, so no rank-three circuit move connects them.  With the missing red
levels restored, the standard rational-normal-scroll path is

    R_3 B_1 B_2,
    R_2 B_2^2,
    R_0 B_3^2.

Alternatively, restoring ``B_0`` gives

    R_3 B_1 B_2,
    R_3 B_0 B_3,
    R_0 B_3^2.

This is an exact combinatorial obstruction to a circuit-only quotient
peeling lemma.  It is not a GVC counterexample.

The script also replays the two dominance-minimal exceptional scroll
fibers from Bogart--Hemmecke--Petrovic.  The ``S(6)`` witness is
one-color.  On the support of the ``S(5,4)`` witness, the fiber has
exactly four states, with circuit-component sizes three and one.  Its
three repeated-ray factorial signatures are ``(2,2)``, ``(2,1,1)``,
and ``(1,1,1,1)``.
"""

from __future__ import annotations

from collections import defaultdict, deque
from itertools import product


Vector = tuple[int, ...]


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def fiber_states(
    red_levels: tuple[int, ...],
    blue_levels: tuple[int, ...],
    red_count: int,
    blue_count: int,
    weight: int,
) -> tuple[Vector, ...]:
    states = []
    levels = red_levels + blue_levels
    for red in compositions(red_count, len(red_levels)):
        for blue in compositions(blue_count, len(blue_levels)):
            state = red + blue
            if sum(entry * level for entry, level in zip(state, levels)) == weight:
                states.append(state)
    return tuple(states)


def support_size(difference: Vector) -> int:
    return sum(entry != 0 for entry in difference)


def circuit_graph(states: tuple[Vector, ...]) -> dict[Vector, tuple[Vector, ...]]:
    graph: dict[Vector, list[Vector]] = {state: [] for state in states}
    for index, left in enumerate(states):
        for right in states[index + 1 :]:
            difference = tuple(a - b for a, b in zip(left, right))
            if support_size(difference) <= 4:
                graph[left].append(right)
                graph[right].append(left)
    return {state: tuple(neighbors) for state, neighbors in graph.items()}


def shortest_path(
    graph: dict[Vector, tuple[Vector, ...]],
    start: Vector,
    finish: Vector,
) -> tuple[Vector, ...] | None:
    queue = deque([start])
    previous: dict[Vector, Vector | None] = {start: None}
    while queue:
        state = queue.popleft()
        if state == finish:
            path = []
            current: Vector | None = state
            while current is not None:
                path.append(current)
                current = previous[current]
            return tuple(reversed(path))
        for neighbor in graph[state]:
            if neighbor not in previous:
                previous[neighbor] = state
                queue.append(neighbor)
    return None


def is_graver_primitive(
    difference: Vector,
    levels: tuple[int, ...],
    colors: tuple[int, ...],
) -> bool:
    """Test absence of a nonzero proper conformal colored subidentity."""

    positive = tuple(i for i, entry in enumerate(difference) if entry > 0)
    negative = tuple(i for i, entry in enumerate(difference) if entry < 0)
    positive_signatures: dict[
        tuple[tuple[int, int], int],
        list[tuple[int, ...]],
    ] = defaultdict(list)

    for values in product(
        *(range(difference[index] + 1) for index in positive)
    ):
        if not any(values):
            continue
        color_counts = [0, 0]
        weight = 0
        for index, value in zip(positive, values):
            color_counts[colors[index]] += value
            weight += value * levels[index]
        positive_signatures[(tuple(color_counts), weight)].append(values)

    for values in product(
        *(range(-difference[index] + 1) for index in negative)
    ):
        if not any(values):
            continue
        color_counts = [0, 0]
        weight = 0
        for index, value in zip(negative, values):
            color_counts[colors[index]] += value
            weight += value * levels[index]
        signature = (tuple(color_counts), weight)
        for positive_values in positive_signatures.get(signature, ()):
            full_positive = all(
                value == difference[index]
                for index, value in zip(positive, positive_values)
            )
            full_negative = all(
                value == -difference[index]
                for index, value in zip(negative, values)
            )
            if not (full_positive and full_negative):
                return False
    return True


def first_five_support_primitive(
    maximum_level: int,
) -> tuple[int, int, Vector] | None:
    """Find the first degree-three two-color primitive beyond circuits."""

    for radial_level in range(1, maximum_level + 1):
        levels = tuple(range(radial_level + 1)) * 2
        colors = (
            (0,) * (radial_level + 1)
            + (1,) * (radial_level + 1)
        )
        for red_count, blue_count in ((1, 2), (2, 1)):
            groups: dict[int, list[Vector]] = defaultdict(list)
            for red in compositions(red_count, radial_level + 1):
                for blue in compositions(blue_count, radial_level + 1):
                    state = red + blue
                    weight = sum(
                        entry * level
                        for entry, level in zip(state, levels)
                    )
                    groups[weight].append(state)
            for weight in sorted(groups):
                states = groups[weight]
                for left_index, left in enumerate(states):
                    for right in states[left_index + 1 :]:
                        difference = tuple(
                            a - b for a, b in zip(left, right)
                        )
                        if support_size(difference) <= 4:
                            continue
                        if any(a and b for a, b in zip(left, right)):
                            continue
                        if is_graver_primitive(difference, levels, colors):
                            return radial_level, weight, difference
    return None


def connected_components(
    graph: dict[Vector, tuple[Vector, ...]],
) -> tuple[tuple[Vector, ...], ...]:
    unseen = set(graph)
    components = []
    while unseen:
        start = min(unseen)
        queue = deque([start])
        unseen.remove(start)
        component = []
        while queue:
            state = queue.popleft()
            component.append(state)
            for neighbor in graph[state]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components, key=lambda component: (len(component), component)))


def positive_partition(state: Vector) -> tuple[int, ...]:
    return tuple(sorted((entry for entry in state if entry), reverse=True))


def entropy_base(partition: tuple[int, ...]) -> int:
    """Return product k^k controlling the factorial ray asymptotic."""

    answer = 1
    for part in partition:
        answer *= part**part
    return answer


def verify() -> None:
    projected_red = (0, 3)
    blue = (1, 2, 3)
    projected_states = fiber_states(projected_red, blue, 1, 2, 6)
    projected_start = (0, 1, 1, 1, 0)
    projected_finish = (1, 0, 0, 0, 2)
    assert projected_states == (projected_start, projected_finish)
    projected_graph = circuit_graph(projected_states)
    assert shortest_path(
        projected_graph,
        projected_start,
        projected_finish,
    ) is None

    projected_difference = tuple(
        a - b for a, b in zip(projected_start, projected_finish)
    )
    assert support_size(projected_difference) == 5
    assert is_graver_primitive(
        projected_difference,
        projected_red + blue,
        (0, 0, 1, 1, 1),
    )

    completed_red = (0, 1, 2, 3)
    completed_states = fiber_states(completed_red, blue, 1, 2, 6)
    completed_start = (0, 0, 0, 1, 1, 1, 0)
    completed_finish = (1, 0, 0, 0, 0, 0, 2)
    expected_path = (
        completed_start,
        (0, 0, 1, 0, 0, 2, 0),
        completed_finish,
    )
    completed_path = shortest_path(
        circuit_graph(completed_states),
        completed_start,
        completed_finish,
    )
    assert completed_path == expected_path
    assert all(
        support_size(
            tuple(a - b for a, b in zip(left, right))
        )
        == 4
        for left, right in zip(expected_path, expected_path[1:])
    )

    completed_blue = (0, 1, 2, 3)
    blue_completed_states = fiber_states(
        projected_red,
        completed_blue,
        1,
        2,
        6,
    )
    blue_completed_start = (0, 1, 0, 1, 1, 0)
    blue_completed_finish = (1, 0, 0, 0, 0, 2)
    expected_blue_path = (
        blue_completed_start,
        (0, 1, 1, 0, 0, 1),
        blue_completed_finish,
    )
    blue_completed_path = shortest_path(
        circuit_graph(blue_completed_states),
        blue_completed_start,
        blue_completed_finish,
    )
    assert blue_completed_path == expected_blue_path
    assert all(
        support_size(
            tuple(a - b for a, b in zip(left, right))
        )
        == 4
        for left, right in zip(
            expected_blue_path,
            expected_blue_path[1:],
        )
    )

    first = first_five_support_primitive(6)
    expected_difference = (
        -1,
        0,
        0,
        1,
        2,
        -1,
        -1,
        0,
    )
    assert first == (3, 3, expected_difference)

    # Bogart--Hemmecke--Petrovic's S(6) primitive non-UGB witness.
    s6_levels = tuple(range(7))
    s6_difference = (1, -1, 1, -1, -1, 0, 1)
    s6_positive = (1, 0, 1, 0, 0, 0, 1)
    s6_negative = (0, 1, 0, 1, 1, 0, 0)
    s6_auxiliary = (
        (1, 0, 0, 0, 2, 0, 0),
        (0, 2, 0, 0, 0, 0, 1),
        (0, 0, 1, 2, 0, 0, 0),
    )
    assert is_graver_primitive(
        s6_difference,
        s6_levels,
        (0,) * 7,
    )
    assert all(
        sum(state) == 3
        and sum(level * entry for level, entry in zip(s6_levels, state)) == 8
        for state in (s6_positive, s6_negative, *s6_auxiliary)
    )
    assert s6_difference == tuple(
        sum(state[index] for state in s6_auxiliary)
        - 3 * s6_negative[index]
        for index in range(7)
    )

    # The projected support of the genuinely two-color S(5,4) witness.
    s54_red = (0, 1, 4, 5)
    s54_blue = (0, 1, 4)
    s54_states = fiber_states(s54_red, s54_blue, 2, 2, 8)
    s54_auxiliary_left = (0, 0, 2, 0, 2, 0, 0)
    s54_negative = (0, 1, 0, 1, 0, 2, 0)
    s54_positive = (1, 0, 1, 0, 1, 0, 1)
    s54_auxiliary_right = (2, 0, 0, 0, 0, 0, 2)
    assert s54_states == (
        s54_auxiliary_left,
        s54_negative,
        s54_positive,
        s54_auxiliary_right,
    )
    s54_difference = tuple(
        left - right
        for left, right in zip(s54_positive, s54_negative)
    )
    assert is_graver_primitive(
        s54_difference,
        s54_red + s54_blue,
        (0,) * 4 + (1,) * 3,
    )
    assert all(
        2 * s54_positive[index]
        == s54_auxiliary_left[index] + s54_auxiliary_right[index]
        for index in range(7)
    )
    s54_components = connected_components(circuit_graph(s54_states))
    assert tuple(len(component) for component in s54_components) == (1, 3)
    assert s54_components[0] == (s54_negative,)

    s54_partitions = {
        positive_partition(state)
        for state in s54_states
    }
    assert s54_partitions == {
        (2, 2),
        (2, 1, 1),
        (1, 1, 1, 1),
    }
    assert {
        partition: entropy_base(partition)
        for partition in s54_partitions
    } == {
        (2, 2): 16,
        (2, 1, 1): 4,
        (1, 1, 1, 1): 1,
    }

    print(
        "PASS projected quotient obstruction: "
        "R3*B1*B2 = R0*B3^2 is, up to level reversal, "
        "the first five-support two-color primitive"
    )
    print(
        "projected fiber: two states and no support-at-most-four "
        "circuit path"
    )
    print(
        "completed-scroll path: "
        "R3*B1*B2 -> R2*B2^2 -> R0*B3^2"
    )
    print(
        "reversal-completed path: "
        "R3*B1*B2 -> R3*B0*B3 -> R0*B3^2"
    )
    print(
        "minimal exceptional scrolls: S(6) one-color witness replayed; "
        "projected S(5,4) fiber has four states and circuit components "
        "of sizes 1+3"
    )
    print(
        "S(5,4) repeated-ray factorial signatures: "
        "(2,2)->16, (2,1,1)->4, (1,1,1,1)->1"
    )
    print(
        "STATUS: the isolated ghost atoms remain terminal, but "
        "circuit-only high-quotient peeling is false on projected "
        "supports; projected Graver blocks are the refined gap"
    )


if __name__ == "__main__":
    verify()

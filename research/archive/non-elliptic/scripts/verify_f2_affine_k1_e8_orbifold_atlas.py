#!/usr/bin/env python3
"""Enumerate the simple-inertia orbifold quotients of the E8 cusp group."""

from __future__ import annotations

from collections import Counter

from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_r, free_group

from verify_f2_affine_k1_e8_monodromy import (
    Permutation,
    compose,
    generated_group,
    inverse,
    power,
)


def regular_orbifold_group() -> tuple[list[Permutation], int, int, int, int]:
    """Return the regular action of <a,b | a^3=b^5, (a^-1 b^2)^2=1>."""

    free, a_word, b_word = free_group("a,b")
    presentation = FpGroup(
        free,
        [a_word**3 * b_word**-5, (a_word**-1 * b_word**2) ** 2],
    )
    cosets = coset_enumeration_r(presentation, [], max_cosets=10_000)
    cosets.compress()
    cosets.standardize()
    assert cosets.n == 240
    a_image = tuple(row[0] for row in cosets.table)
    b_image = tuple(row[2] for row in cosets.table)
    group = list(generated_group((a_image, b_image)))
    identity = tuple(range(240))
    group.sort(key=lambda element: (element != identity, element))
    index = {element: position for position, element in enumerate(group)}
    meridian = compose(inverse(a_image), power(b_image, 2))
    central = power(a_image, 3)
    assert power(a_image, 3) == power(b_image, 5)
    assert power(meridian, 2) == identity
    assert power(central, 4) == identity
    assert power(central, 2) != identity
    assert compose(central, a_image) == compose(a_image, central)
    assert compose(central, b_image) == compose(b_image, central)
    return group, index[a_image], index[b_image], index[meridian], index[central]


def multiplication_data(
    group: list[Permutation],
) -> tuple[list[list[int]], list[int]]:
    """Build the multiplication and inverse tables from a regular action."""

    index = {element: position for position, element in enumerate(group)}
    point_to_index = {element[0]: position for position, element in enumerate(group)}
    assert len(point_to_index) == len(group)
    multiplication = [
        [point_to_index[left[right[0]]] for right in group] for left in group
    ]
    inverses = [index[inverse(element)] for element in group]
    return multiplication, inverses


def subgroup_conjugacy_classes(
    multiplication: list[list[int]], inverses: list[int]
) -> dict[tuple[int, ...], tuple[int, ...]]:
    """Enumerate subgroup classes, retaining one short generating tuple."""

    order = len(multiplication)

    def closure(generators: tuple[int, ...]) -> frozenset[int]:
        subgroup = {0}
        frontier = [0]
        moves = generators + tuple(inverses[element] for element in generators)
        while frontier:
            element = frontier.pop()
            for move in moves:
                product = multiplication[element][move]
                if product not in subgroup:
                    subgroup.add(product)
                    frontier.append(product)
        return frozenset(subgroup)

    def conjugate(subgroup: frozenset[int], element: int) -> tuple[int, ...]:
        return tuple(
            sorted(
                multiplication[multiplication[element][member]][inverses[element]]
                for member in subgroup
            )
        )

    def canonical(subgroup: frozenset[int]) -> tuple[int, ...]:
        return min(conjugate(subgroup, element) for element in range(order))

    trivial = (0,)
    classes: dict[tuple[int, ...], tuple[int, ...]] = {trivial: ()}
    frontier = [trivial]
    while frontier:
        representative = frontier.pop()
        subgroup = frozenset(representative)
        generators = classes[representative]
        for element in range(1, order):
            if element in subgroup:
                continue
            enlarged = closure(generators + (element,))
            key = canonical(enlarged)
            if key not in classes:
                classes[key] = generators + (element,)
                frontier.append(key)
    return classes


def fixed_cosets(
    subgroup: tuple[int, ...],
    meridian: int,
    multiplication: list[list[int]],
    inverses: list[int],
) -> int:
    """Count cosets fixed by left multiplication by the meridian."""

    subgroup_set = set(subgroup)
    conjugating_elements = sum(
        multiplication[multiplication[inverses[element]][meridian]][element]
        in subgroup_set
        for element in range(len(multiplication))
    )
    assert conjugating_elements % len(subgroup) == 0
    return conjugating_elements // len(subgroup)


def coset_action(
    subgroup: tuple[int, ...],
    element: int,
    multiplication: list[list[int]],
) -> Permutation:
    """Return left multiplication by an element on the left cosets of H."""

    subgroup_set = set(subgroup)
    unseen = set(range(len(multiplication)))
    cosets: list[frozenset[int]] = []
    while unseen:
        representative = min(unseen)
        coset = frozenset(
            multiplication[representative][member] for member in subgroup_set
        )
        cosets.append(coset)
        unseen.difference_update(coset)
    index = {coset: position for position, coset in enumerate(cosets)}
    action = []
    for coset in cosets:
        representative = next(iter(coset))
        image = frozenset(
            multiplication[multiplication[element][representative]][member]
            for member in subgroup_set
        )
        action.append(index[image])
    return tuple(action)


def peripheral_signature(
    subgroup: tuple[int, ...],
    meridian: int,
    central: int,
    multiplication: list[list[int]],
) -> tuple[int, int, int, int, int]:
    """Return fixed sheets, ramified f=1,2,4 rows, and unramified orbits."""

    meridian_action = coset_action(subgroup, meridian, multiplication)
    central_action = coset_action(subgroup, central, multiplication)
    degree = len(meridian_action)
    unseen = set(range(degree))
    ramified_f1 = 0
    ramified_f2 = 0
    ramified_f4 = 0
    unramified_orbits = 0
    while unseen:
        start = min(unseen)
        orbit = {start}
        frontier = [start]
        while frontier:
            sheet = frontier.pop()
            for move in (meridian_action, central_action):
                image = move[sheet]
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
        ramified = any(meridian_action[sheet] != sheet for sheet in orbit)
        if not ramified:
            unramified_orbits += 1
        elif len(orbit) == 2:
            ramified_f1 += 1
        elif len(orbit) == 4:
            ramified_f2 += 1
        else:
            assert len(orbit) == 8
            ramified_f4 += 1
    fixed = sum(
        meridian_action[sheet] == sheet for sheet in range(degree)
    )
    assert degree == (
        fixed + 2 * ramified_f1 + 4 * ramified_f2 + 8 * ramified_f4
    )
    return fixed, ramified_f1, ramified_f2, ramified_f4, unramified_orbits


def main() -> None:
    group, _, _, meridian, central = regular_orbifold_group()
    multiplication, inverses = multiplication_data(group)
    classes = subgroup_conjugacy_classes(multiplication, inverses)
    assert len(classes) == 30
    assert dict(sorted(Counter(map(len, classes)).items())) == {
        1: 1,
        2: 2,
        3: 1,
        4: 3,
        5: 1,
        6: 3,
        8: 3,
        10: 3,
        12: 3,
        16: 1,
        20: 3,
        24: 2,
        40: 1,
        48: 1,
        120: 1,
        240: 1,
    }
    rows = []
    for subgroup in classes:
        degree = len(group) // len(subgroup)
        fixed = fixed_cosets(subgroup, meridian, multiplication, inverses)
        assert (degree - fixed) % 2 == 0
        if fixed:
            signature = peripheral_signature(
                subgroup, meridian, central, multiplication
            )
            assert signature[0] == fixed
            ramified_weight = signature[1] + 2 * signature[2] + 4 * signature[3]
            row_count = signature[1] + signature[2] + signature[3]
            minimal_maximal = (
                7 * degree - 62 + 4 * row_count - 8 * ramified_weight
            )
            rows.append(
                (
                    degree,
                    len(subgroup),
                    central in subgroup,
                    *signature,
                    minimal_maximal,
                )
            )
    expected_rows = [
        (1, 240, True, 1, 0, 0, 0, 1, -55),
        (5, 48, True, 1, 2, 0, 0, 1, -35),
        (6, 40, True, 2, 2, 0, 0, 2, -28),
        (10, 24, True, 2, 4, 0, 0, 2, -8),
        (12, 20, False, 4, 0, 2, 0, 2, -2),
        (15, 16, True, 3, 6, 0, 0, 3, 19),
        (20, 12, False, 4, 0, 4, 0, 2, 30),
        (24, 10, False, 4, 0, 1, 2, 1, 38),
        (24, 10, False, 4, 0, 1, 2, 1, 38),
        (30, 8, False, 4, 1, 6, 0, 2, 72),
        (30, 8, True, 2, 14, 0, 0, 2, 92),
        (40, 6, False, 4, 0, 1, 4, 1, 94),
        (40, 6, False, 4, 0, 1, 4, 1, 94),
        (60, 4, False, 4, 0, 14, 0, 2, 190),
        (120, 2, False, 4, 0, 1, 14, 1, 374),
    ]
    assert sorted(rows) == expected_rows
    assert {
        row[0] for row in rows if row[0] >= 6
    } == {6, 10, 12, 15, 20, 24, 30, 40, 60, 120}
    print(
        "PASS: the order-240 E8 simple-inertia orbifold group has 30 "
        "subgroup classes and exactly 13 fixed-sheet F2 coset actions in "
        "degrees 6,10,12,15,20,24,30,40,60,120; peripheral rows have "
        "(e,f)=(2,f) with f in {1,2,4}, and the minimal maximal-contact "
        "squarefree ledger is negative only in degrees 6,10,12"
    )


if __name__ == "__main__":
    main()

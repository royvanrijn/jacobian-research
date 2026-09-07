#!/usr/bin/env python3
"""Verify the six-sheet filling obstruction for the k=1 A4+A2+A1 row."""

from __future__ import annotations

from collections.abc import Iterable
from itertools import combinations, permutations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


Permutation = tuple[int, ...]
Word = tuple[int, ...]


# Exact Tietze-reduced Zariski--van Kampen presentation at
# (a,b,c,d)=(-3,5/2,-5,-5).  The three generators remain geometric
# meridians.  Sage writes words using a right permutation action.
RELATORS: tuple[Word, ...] = (
    (3, -1, -3, -1, 3, 1),
    (3, 1, 2, 3, -1, -3, -2, -1),
    (1, 2, 1, 2, 1, -2, -1, -2, -1, -2),
)

# The unique simultaneous-conjugacy class of transitive degree-six actions
# in which all three geometric meridians are 3-cycles.
ACTION: tuple[Permutation, ...] = (
    (1, 2, 0, 3, 4, 5),
    (3, 1, 2, 4, 0, 5),
    (2, 1, 5, 3, 4, 0),
)


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(value: Permutation) -> Permutation:
    result = [0] * len(value)
    for index, image in enumerate(value):
        result[image] = index
    return tuple(result)


def cycle_type(value: Permutation) -> tuple[int, ...]:
    unseen = set(range(len(value)))
    lengths: list[int] = []
    while unseen:
        point = start = min(unseen)
        length = 0
        while point in unseen:
            unseen.remove(point)
            length += 1
            point = value[point]
        assert point == start
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def three_cycles() -> tuple[Permutation, ...]:
    result: list[Permutation] = []
    for support in combinations(range(6), 3):
        first, second, third = support
        for cycle in ((first, second, third), (first, third, second)):
            image = list(range(6))
            image[cycle[0]] = cycle[1]
            image[cycle[1]] = cycle[2]
            image[cycle[2]] = cycle[0]
            result.append(tuple(image))
    assert len(result) == 40
    return tuple(result)


def generated_group(images: tuple[Permutation, ...]) -> set[Permutation]:
    identity = tuple(range(len(images[0])))
    generators = images + tuple(inverse(image) for image in images)
    result = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(current, generator)
            if candidate not in result:
                result.add(candidate)
                frontier.append(candidate)
    return result


def evaluate_right_word(images: tuple[Permutation, ...], word: Word) -> Permutation:
    result = tuple(range(len(images[0])))
    inverses = tuple(inverse(image) for image in images)
    for letter in word:
        image = images[letter - 1] if letter > 0 else inverses[-letter - 1]
        result = compose(result, image)
    return result


def conjugate(image: Permutation, change: Permutation) -> Permutation:
    return compose(compose(change, image), inverse(change))


def cubic_actions() -> tuple[tuple[Permutation, ...], ...]:
    cycles = three_cycles()
    first = ACTION[0]
    identity = tuple(range(6))
    solutions: list[tuple[Permutation, ...]] = []
    for second, third in product(cycles, repeat=2):
        images = (first, second, third)
        if not all(
            evaluate_right_word(images, relation) == identity
            for relation in RELATORS
        ):
            continue
        group = generated_group(images)
        if {image[0] for image in group} == set(range(6)):
            solutions.append(images)
    return tuple(solutions)


def lifted_path(start: int, word: Word) -> tuple[int, tuple[int, ...]]:
    """Lift a Sage right-action word and return endpoint and cellular chain."""

    degree = len(ACTION[0])
    # A positive edge follows inverse(ACTION[j]) when sheet labels are
    # written on the left.  This is the left-action realization of Sage's
    # right-action tuple convention.
    positive = tuple(inverse(image) for image in ACTION)
    point = start
    chain = [0] * (len(ACTION) * degree)
    for letter in word:
        generator = abs(letter) - 1
        if letter > 0:
            chain[generator * degree + point] += 1
            point = positive[generator][point]
        else:
            previous = ACTION[generator][point]
            chain[generator * degree + previous] -= 1
            point = previous
    return point, tuple(chain)


def spanning_tree_chords() -> tuple[int, ...]:
    degree = len(ACTION[0])
    positive = tuple(inverse(image) for image in ACTION)
    parent = list(range(degree))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> bool:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return False
        parent[right_root] = left_root
        return True

    tree: list[int] = []
    for generator, image in enumerate(positive):
        for source, target in enumerate(image):
            edge = generator * degree + source
            if source != target and union(source, target):
                tree.append(edge)
    assert len(tree) == degree - 1
    return tuple(
        edge
        for edge in range(len(ACTION) * degree)
        if edge not in set(tree)
    )


def smith_data(
    relations: Iterable[tuple[int, ...]], chords: tuple[int, ...]
) -> tuple[tuple[int, ...], int, tuple[int, ...]]:
    columns = tuple(relations)
    matrix = sp.Matrix([[column[edge] for column in columns] for edge in chords])
    smith = smith_normal_form(matrix, domain=ZZ)
    diagonal = tuple(
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
    )
    rank = sum(value != 0 for value in diagonal)
    torsion = tuple(value for value in diagonal if value > 1)
    return diagonal, len(chords) - rank, torsion


def exponent_sum_matrix() -> sp.Matrix:
    return sp.Matrix(
        [
            [
                sum(
                    1 if letter == generator else -1 if letter == -generator else 0
                    for letter in relation
                )
                for relation in RELATORS
            ]
            for generator in (1, 2, 3)
        ]
    )


def main() -> None:
    identity = tuple(range(6))
    assert all(evaluate_right_word(ACTION, relation) == identity for relation in RELATORS)
    assert all(cycle_type(image) == (3, 1, 1, 1) for image in ACTION)
    group = generated_group(ACTION)
    assert len(group) == 360
    assert {image[0] for image in group} == set(range(6))

    solutions = cubic_actions()
    assert len(solutions) == 18
    assert ACTION in solutions
    assert {len(generated_group(images)) for images in solutions} == {360}
    centralizer = tuple(
        change
        for change in permutations(range(6))
        if compose(change, ACTION[0]) == compose(ACTION[0], change)
    )

    def canonical(images: tuple[Permutation, ...]) -> tuple[int, ...]:
        return min(
            tuple(
                value
                for image in tuple(conjugate(item, change) for item in images)
                for value in image
            )
            for change in centralizer
        )

    assert len({canonical(images) for images in solutions}) == 1

    base_smith = smith_normal_form(exponent_sum_matrix(), domain=ZZ)
    assert tuple(abs(int(base_smith[i, i])) for i in range(3)) == (1, 1, 0)

    relator_lifts: list[tuple[int, ...]] = []
    for relation in RELATORS:
        for sheet in range(6):
            endpoint, chain = lifted_path(sheet, relation)
            assert endpoint == sheet
            relator_lifts.append(chain)

    chords = spanning_tree_chords()
    assert len(chords) == 3 * 6 - 6 + 1 == 13
    unfilled = smith_data(relator_lifts, chords)
    assert unfilled == ((1,) * 11 + (0, 0), 2, ())

    all_fixed_lifts: list[tuple[int, ...]] = []
    for generator, image in enumerate(ACTION, start=1):
        fixed_lifts = []
        for sheet in range(6):
            if image[sheet] != sheet:
                continue
            endpoint, chain = lifted_path(sheet, (generator,))
            assert endpoint == sheet
            fixed_lifts.append(chain)
        assert len(fixed_lifts) == 3
        all_fixed_lifts.extend(fixed_lifts)
        filled = smith_data(relator_lifts + fixed_lifts, chords)
        assert filled == ((1,) * 12 + (0,), 1, ())

    # Adding all conjugate meridian representatives produces no extra kill:
    # the same primitive free class survives.
    fully_filled = smith_data(relator_lifts + all_fixed_lifts, chords)
    assert fully_filled == ((1,) * 12 + (0,), 1, ())

    print(
        "PASS: the unique cubic degree-six A4+A2+A1 action has image A6; "
        "its cover has H_1=Z^2, and filling the three fixed lifts of any "
        "geometric meridian (or all nine representatives) leaves H_1=Z"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env sage-python
"""Certify the generic F2 ``r=6`` cusp-pair braid obstruction.

This verifier requires SageMath with the optional SIROCCO package.  Unlike
``explore_f2_r8_cusp_braid.py``, every continued root path used here is
certified by SIROCCO.  The witness is a rational point of the main
determinant-open cusp surface.
"""

import argparse
from itertools import product

from sage.all import PolynomialRing, QQ
from sage.schemes.curves.zariski_vankampen import braid_monodromy


def implicit_pair(coefficients, field=QQ):
    """Return the two exact implicit quintics for one coefficient row."""

    ring = PolynomialRing(field, names=("P", "Q", "t"))
    P, Q, t = ring.gens()

    a, b, c, d, A, E, C, D, p0, q0 = map(field, coefficients)

    p1 = t**3 + a * t
    q1 = t**5 + b * t**4 + c * t**2 + d * t
    p2 = t**3 + A * t + p0
    q2 = t**5 + b * t**4 + E * t**3 + C * t**2 + D * t + q0
    first = (P - p1).resultant(Q - q1, t)
    second = (P - p2).resultant(Q - q2, t)

    target = PolynomialRing(field, names=("P", "Q"))
    first = target(first)
    second = target(second)
    assert first.total_degree() == second.total_degree() == 5
    assert first.degree(target.gen(1)) == second.degree(target.gen(1)) == 3
    for equation in (first, second):
        factors = list(equation.factor())
        assert len(factors) == 1 and factors[0][1] == 1
    return first, second


def transposition_audit(
    words,
    component_by_strand,
    label,
    expect_equal_disjoint=False,
):
    """Enumerate the transposition representations fixed by all braids."""

    identity = tuple(range(6))
    transpositions = []
    for left in range(6):
        for right in range(left + 1, 6):
            permutation = list(identity)
            permutation[left], permutation[right] = right, left
            transpositions.append(tuple(permutation))
    lookup = {permutation: index for index, permutation in enumerate(transpositions)}

    def compose(left, right):
        return tuple(left[right[index]] for index in range(6))

    conjugate = [[0] * 15 for _ in range(15)]
    for outer_index, outer in enumerate(transpositions):
        for inner_index, inner in enumerate(transpositions):
            conjugate[outer_index][inner_index] = lookup[
                compose(outer, compose(inner, outer))
            ]

    def act(values, word):
        result = list(values)
        for letter in word:
            position = abs(letter) - 1
            left, right = result[position], result[position + 1]
            if letter > 0:
                result[position] = right
                result[position + 1] = conjugate[right][left]
            else:
                result[position] = conjugate[left][right]
                result[position + 1] = left
        return tuple(result)

    def transitive(values):
        graph = [set() for _ in range(6)]
        for value in values:
            moved = [
                sheet
                for sheet in range(6)
                if transpositions[value][sheet] != sheet
            ]
            graph[moved[0]].add(moved[1])
            graph[moved[1]].add(moved[0])
        reached = {0}
        frontier = [0]
        while frontier:
            sheet = frontier.pop()
            for neighbor in graph[sheet]:
                if neighbor not in reached:
                    reached.add(neighbor)
                    frontier.append(neighbor)
        return len(reached) == 6

    survivors = []
    # Simultaneous conjugacy fixes the first meridian at (0,1).
    for tail in product(range(15), repeat=5):
        values = (0,) + tail
        if all(act(values, word) == values for word in words):
            survivors.append(values)

    first_component = component_by_strand[0]
    other_component = 1 - first_component
    expected = []
    for other in range(15):
        first_edge = set(index for index in range(6) if transpositions[0][index] != index)
        other_edge = set(
            index for index in range(6) if transpositions[other][index] != index
        )
        if other != 0 and first_edge & other_edge:
            continue
        expected.append(
            tuple(
                0 if component_by_strand[index] == first_component else other
                for index in range(6)
            )
        )

    if expect_equal_disjoint:
        assert sorted(survivors) == sorted(expected)
        assert len(survivors) == 7
    transitive_count = sum(transitive(values) for values in survivors)
    assert transitive_count == 0
    print(
        f"certified {label} braid: {len(words)} braids, "
        f"{len(survivors)} assignments, {transitive_count} transitive"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--witness",
        choices=("all", "r=6 generic A2", "r=6 E6-I", "r=6 E6-II"),
        default="all",
    )
    arguments = parser.parse_args()
    witnesses = (
        (
            "r=6 generic A2",
            (
                -3,
                QQ(35) / 12,
                QQ(-325) / 54,
                QQ(-125) / 27,
                -2,
                QQ(5) / 3,
                QQ(-205) / 54,
                QQ(-170) / 27,
                -1,
                QQ(-365) / 162,
            ),
            True,
        ),
        (
            "r=6 E6-I",
            (
                0,
                1,
                0,
                0,
                QQ(-36) / 25,
                QQ(-12) / 5,
                QQ(-24) / 25,
                QQ(48) / 25,
                QQ(72) / 125,
                QQ(-288) / 625,
            ),
            False,
        ),
        (
            "r=6 E6-II",
            (
                0,
                1,
                0,
                0,
                QQ(-32) / 25,
                QQ(-32) / 15,
                QQ(-128) / 225,
                QQ(2048) / 1125,
                QQ(256) / 375,
                QQ(-2048) / 3375,
            ),
            False,
        ),
    )
    for label, coefficients, expect_equal_disjoint in witnesses:
        if arguments.witness != "all" and arguments.witness != label:
            continue
        first, second = implicit_pair(coefficients)
        braids, component_by_strand, vertical, degree = braid_monodromy(
            first * second,
            arrangement=(first, second),
        )
        assert degree == 6
        if label == "r=6 generic A2":
            assert len(braids) == 17
        assert vertical == {}
        assert sorted(component_by_strand.values()) == [0, 0, 0, 1, 1, 1]
        words = [tuple(braid.Tietze()) for braid in braids]
        transposition_audit(
            words,
            component_by_strand,
            label,
            expect_equal_disjoint,
        )


if __name__ == "__main__":
    main()

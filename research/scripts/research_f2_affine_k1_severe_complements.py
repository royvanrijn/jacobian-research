#!/usr/bin/env python3
"""Explore complement groups on the remaining F2 ``k=1`` strata.

Run with SageMath and its optional ``sirocco`` library.  This is an
exploratory calculator: the theorem checker pins any results promoted to
the mathematical status ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations, product

from sage.all import PolynomialRing, QQ
from sage.schemes.curves.zariski_vankampen import braid_monodromy, fundamental_group


@dataclass(frozen=True)
class Witness:
    label: str
    parameters: tuple[object, object, object, object]


def witnesses():
    half = QQ(1) / 2
    quarter = QQ(1) / 4
    return (
        Witness("A2+A3+A1", (-3, -half, 3, -9)),
        Witness("A4+2A1", (-3, -9 * quarter, QQ(13) / 2, -9)),
        Witness("A4+A3", (-3, 10, 80, -205)),
        Witness("A6+A1", (-3, -QQ(7) / 2, -1, 11)),
        Witness("A4+A2+A1", (-3, QQ(5) / 2, -5, -5)),
        Witness("D4+A1", (1, 1, 0, -2)),
        Witness("D4+A2", (-3, 2, 0, 3)),
        Witness("D5+A1", (-3, 0, -2, -9)),
        Witness("D6", (1, 1, -1, -3)),
        Witness("D5+A2", (-3, -2, 4, -5)),
        Witness("E7-I3", (-3, -1, 1, -7)),
        Witness("E7-A4", (-3, 2, -8, -13)),
    )


def implicit_polynomial(parameters):
    ambient = PolynomialRing(QQ, names=("P", "Q", "t"))
    P3, Q3, t = ambient.gens()
    a, b, c, d = parameters
    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    resultant = (P3 - p).resultant(Q3 - q, t)

    plane = PolynomialRing(QQ, names=("P", "Q"))
    return plane(resultant)


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def evaluate_word(images, word):
    result = tuple(range(len(images[0])))
    inverses = tuple(inverse(image) for image in images)
    for letter in word:
        image = images[letter - 1] if letter > 0 else inverses[-letter - 1]
        result = compose(result, image)
    return result


def three_cycles(degree=6):
    result = []
    for support in combinations(range(degree), 3):
        first, second, third = support
        for oriented in ((first, second, third), (first, third, second)):
            image = list(range(degree))
            image[oriented[0]] = oriented[1]
            image[oriented[1]] = oriented[2]
            image[oriented[2]] = oriented[0]
            result.append(tuple(image))
    return tuple(result)


def generated_group(images):
    identity = tuple(range(len(images[0])))
    generators = tuple(images) + tuple(inverse(image) for image in images)
    elements = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(current, generator)
            if candidate not in elements:
                elements.add(candidate)
                frontier.append(candidate)
    return elements


def transitive(images):
    group = generated_group(images)
    return {element[0] for element in group} == set(range(len(images[0])))


def conjugate(image, change):
    return compose(compose(change, image), inverse(change))


def cubic_degree_six_actions(raw_group):
    """Enumerate degree-six actions with every meridian a 3-cycle."""

    cycles = three_cycles()
    first = cycles[0]
    identity = tuple(range(6))
    relations = tuple(
        tuple(relation.Tietze())
        for relation in raw_group.relations()
        if tuple(relation.Tietze())
    )
    centralizer = tuple(
        image
        for image in permutations(range(6))
        if compose(image, first) == compose(first, image)
    )
    solutions = []
    for second, third in product(cycles, repeat=2):
        images = (first, second, third)
        if not all(evaluate_word(images, relation) == identity for relation in relations):
            continue
        if transitive(images):
            solutions.append(images)

    def canonical(images):
        return min(
            tuple(value for image in changed for value in image)
            for change in centralizer
            for changed in [tuple(conjugate(image, change) for image in images)]
        )

    representatives = {}
    for images in solutions:
        representatives.setdefault(canonical(images), images)
    return solutions, tuple(representatives.values())


def main() -> None:
    noncyclic = {"A4+A3", "A6+A1", "A4+A2+A1", "D5+A2"}
    for witness in witnesses():
        polynomial = implicit_polynomial(witness.parameters)
        print(f"BEGIN {witness.label}", flush=True)
        print(f"parameters={witness.parameters}", flush=True)
        print(f"polynomial={polynomial}", flush=True)
        group = fundamental_group(polynomial, simplified=True, projective=False)
        print(f"group={group}", flush=True)
        print(f"ngens={group.ngens()}", flush=True)
        print(
            "relations="
            + repr(tuple(tuple(relation.Tietze()) for relation in group.relations())),
            flush=True,
        )
        print(f"abelian={tuple(group.abelian_invariants())}", flush=True)
        if witness.label in noncyclic:
            braid_data = braid_monodromy(polynomial)
            # Sage 9.5 returns the braid list directly; current Sage also
            # returns component, vertical-line, and projection metadata.
            braids = braid_data[0] if (
                isinstance(braid_data, tuple) and len(braid_data) == 4
            ) else braid_data
            print(
                "braids="
                + repr(tuple(tuple(braid.Tietze()) for braid in braids)),
                flush=True,
            )
            raw_group = fundamental_group(
                polynomial, simplified=False, projective=False
            )
            print(f"raw_ngens={raw_group.ngens()}", flush=True)
            print(
                "raw_relations="
                + repr(
                    tuple(
                        tuple(relation.Tietze())
                        for relation in raw_group.relations()
                    )
                ),
                flush=True,
            )
            solutions, representatives = cubic_degree_six_actions(raw_group)
            print(f"cubic_labeled_with_first_fixed={len(solutions)}", flush=True)
            print(f"cubic_conjugacy_classes={len(representatives)}", flush=True)
            for index, images in enumerate(representatives):
                print(
                    f"cubic_class_{index}="
                    f"order:{len(generated_group(images))},images:{images}",
                    flush=True,
                )
        print(f"END {witness.label}", flush=True)


if __name__ == "__main__":
    main()

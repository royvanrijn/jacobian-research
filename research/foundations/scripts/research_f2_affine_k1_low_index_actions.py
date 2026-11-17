#!/usr/bin/env sage
"""Enumerate low-index actions of the four noncyclic F2 k=1 complements.

This is an exploratory Sage/SIROCCO calculation.  It combines the strict
geometric-degree ceiling supplied by Makar-Limanov's Newton-polytope bound
with the exact Zariski--van Kampen presentations computed by
``research_f2_affine_k1_severe_complements.py``.

For every conjugacy class of subgroups through the requested index, the
script records the transitive coset action, the common geometric-meridian
cycle type, and the orders of the composition factors of its finite image.
The terminal F2 residue cover forces A6 to occur as a section (a subgroup
quotient) of the global monodromy image.  This implies nonsolvability and
divisibility of the image order by 360, but it does not force A6 to be a
composition factor.  The script therefore reports composition factors only
as diagnostics and uses the two valid weaker conditions as its preliminary
filter.  Passing that filter is not an exact A6-section test, a cover, or a
Keller-map construction.
"""

from __future__ import annotations

import argparse
from functools import reduce

import sage.all  # Initialize Sage before importing the libgap Cython module.
from sage.libs.gap.libgap import libgap

from research_f2_affine_k1_severe_complements import (
    fundamental_group,
    implicit_polynomial,
    witnesses,
)


NONCYCLIC_LABELS = ("A4+A3", "A6+A1", "A4+A2+A1", "D5+A2")


def cycle_type(permutation, degree: int) -> tuple[int, ...]:
    images = tuple(int(value) - 1 for value in libgap.ListPerm(permutation, degree))
    seen: set[int] = set()
    cycles: list[int] = []
    for start in range(degree):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = images[current]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def composition_factor_orders(group) -> tuple[int, ...]:
    series = tuple(libgap.CompositionSeries(group))
    orders = tuple(int(libgap.Size(term)) for term in series)
    return tuple(orders[index] // orders[index + 1] for index in range(len(orders) - 1))


def evaluate_word(word: tuple[int, ...], generators, group):
    return reduce(
        lambda current, letter: current
        * (generators[letter - 1] if letter > 0 else generators[-letter - 1] ** -1),
        word,
        libgap.One(group),
    )


def enumerate_label(label: str, maximum_index: int) -> None:
    witness = next(row for row in witnesses() if row.label == label)
    polynomial = implicit_polynomial(witness.parameters)
    raw_group = fundamental_group(polynomial, simplified=False, projective=False)

    simplification = raw_group.simplification_isomorphism()
    group = simplification.codomain()
    meridian_words = tuple(
        tuple(int(letter) for letter in simplification(generator).Tietze())
        for generator in raw_group.gens()
    )

    gap_group = group.gap()
    gap_generators = tuple(libgap.GeneratorsOfGroup(gap_group))
    meridians = tuple(
        evaluate_word(word, gap_generators, gap_group) for word in meridian_words
    )

    subgroups = tuple(gap_group.LowIndexSubgroupsFpGroup(maximum_index))
    allowed_degrees = {6, *range(8, maximum_index + 1)}
    rows = []
    for subgroup in subgroups:
        degree = int(libgap.Index(gap_group, subgroup))
        if degree not in allowed_degrees:
            continue

        cosets = libgap.RightCosets(gap_group, subgroup)
        action = libgap.ActionHomomorphism(gap_group, cosets, libgap.OnRight)
        permutations = tuple(libgap.Image(action, meridian) for meridian in meridians)
        types = tuple(cycle_type(permutation, degree) for permutation in permutations)
        if len(set(types)) != 1:
            raise AssertionError(f"geometric meridians are not conjugate in {label}: {types}")

        fixed_sheets = types[0].count(1)
        image = libgap.Image(action)
        factor_orders = composition_factor_orders(image)
        image_order = int(libgap.Size(image))
        is_nonsolvable = not bool(libgap.IsSolvableGroup(image))
        rows.append(
            {
                "degree": degree,
                "cycle_type": types[0],
                "fixed_sheets": fixed_sheets,
                "image_order": image_order,
                "composition_factor_orders": factor_orders,
                "has_affine_fixed_sheet": fixed_sheets > 0,
                "is_ramified": fixed_sheets < degree,
                "is_nonsolvable": is_nonsolvable,
                "order_divisible_by_360": image_order % 360 == 0,
                "has_A6_composition_factor_diagnostic": 360 in factor_orders,
            }
        )

    print(f"BEGIN {label}")
    print(f"simplified_group={group}")
    print(f"meridian_words={meridian_words}")
    print(f"subgroup_classes_through_{maximum_index}={len(subgroups)}")
    for row in sorted(
        rows,
        key=lambda item: (
            item["degree"],
            item["cycle_type"],
            item["image_order"],
            item["composition_factor_orders"],
        ),
    ):
        print(row)

    survivors = tuple(
        row
        for row in rows
        if row["has_affine_fixed_sheet"]
        and row["is_ramified"]
        and row["is_nonsolvable"]
        and row["order_divisible_by_360"]
    )
    print(f"weaker_necessary_filter_survivors={len(survivors)}")
    for row in survivors:
        print(f"SURVIVOR {row}")
    print(f"END {label}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-index", type=int, default=28)
    parser.add_argument(
        "--label",
        action="append",
        choices=NONCYCLIC_LABELS,
        help="restrict to one or more severe complement labels",
    )
    arguments = parser.parse_args()

    labels = tuple(arguments.label) if arguments.label else NONCYCLIC_LABELS
    for label in labels:
        enumerate_label(label, arguments.maximum_index)


if __name__ == "__main__":
    main()

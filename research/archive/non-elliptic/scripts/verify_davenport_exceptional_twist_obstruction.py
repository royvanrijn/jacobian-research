#!/usr/bin/env python3
"""Exact normalizer audit for the degree-seven Davenport action."""

from __future__ import annotations

import sys
from itertools import permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    gl32,
    inverses,
    line_permutation,
    point_permutation,
)


Permutation = tuple[int, ...]
identity = tuple(range(7))


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(7))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * 7
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def conjugate(
    sheet_permutation: Permutation,
    group_element: Permutation,
) -> Permutation:
    return compose(
        compose(sheet_permutation, group_element),
        inverse(sheet_permutation),
    )


matrices = gl32()
matrix_inverses = inverses(matrices)
point_group = {point_permutation(matrix) for matrix in matrices}
line_group = {
    line_permutation(matrix, matrix_inverses[matrix])
    for matrix in matrices
}
assert len(point_group) == 168
assert point_group == line_group
print("PASS: point and line actions have the same degree-seven image")


def generated_subgroup(generators: list[Permutation]) -> set[Permutation]:
    subgroup = {identity}
    changed = True
    while changed:
        changed = False
        for element in tuple(subgroup):
            for generator in generators:
                for product in (
                    compose(element, generator),
                    compose(generator, element),
                ):
                    if product not in subgroup:
                        subgroup.add(product)
                        changed = True
    return subgroup


generators: list[Permutation] = []
generated = {identity}
for element in point_group:
    if element not in generated:
        generators.append(element)
        generated = generated_subgroup(generators)
    if generated == point_group:
        break
assert generated == point_group

normalizer = {
    sheet_permutation
    for sheet_permutation in permutations(range(7))
    if all(
        conjugate(sheet_permutation, generator) in point_group
        for generator in generators
    )
}
assert len(normalizer) == 168
assert normalizer == point_group
print("PASS: N_{S7}(GL_3(F_2))=GL_3(F_2)")


fixed_point_counts = {
    sum(index == image for index, image in enumerate(element))
    for element in point_group
}
assert fixed_point_counts == {0, 1, 3, 7}
print("PASS: fixed-point counts {0,1,3,7} exclude exceptionality")
print("PASS Davenport exceptional-twist obstruction")

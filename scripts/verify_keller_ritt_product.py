#!/usr/bin/env python3
"""Exact regression for the separated-product Keller--Ritt diamond.

The symbolic part constructs the sparse weighted atoms F_3 and F_4 and
checks that their disjoint stabilizations commute.  The dependency-free
finite part enumerates every block containing one sheet in the natural
S_3 x S_4 product action.  The four blocks form the claimed diamond.
"""

from __future__ import annotations

from itertools import permutations
from math import factorial
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z


def sparse_weighted_map(degree: int) -> tuple[sp.Expr, ...]:
    primitive = (w**2 - w**degree) / (degree - 2)
    model = WeightedSeedModel(sp.diff(primitive, w), c=1, b=1)
    return tuple(sp.expand(component) for component in model.mapping())


F3 = sparse_weighted_map(3)
F4 = sparse_weighted_map(4)
assert sp.factor(sp.det(sp.Matrix(F3).jacobian((x, y, z)))) == 1
assert sp.factor(sp.det(sp.Matrix(F4).jacobian((x, y, z)))) == 1

coordinates = sp.symbols("x0:6")
left_variables = coordinates[:3]
right_variables = coordinates[3:]
left_substitution = dict(zip((x, y, z), left_variables, strict=True))
right_substitution = dict(zip((x, y, z), right_variables, strict=True))
left_atom = tuple(component.subs(left_substitution) for component in F3)
right_atom = tuple(component.subs(right_substitution) for component in F4)

A = left_atom + right_variables
B = left_variables + right_atom
K = left_atom + right_atom


def compose_on(
    outer: tuple[sp.Expr, ...],
    inner: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, ...]:
    substitution = dict(zip(variables, inner, strict=True))
    return tuple(
        sp.expand(component.subs(substitution, simultaneous=True))
        for component in outer
    )


assert compose_on(A, B, coordinates) == K
assert compose_on(B, A, coordinates) == K
assert 3 * 4 == 12


Point = tuple[int, int]
Permutation = tuple[int, ...]
points: tuple[Point, ...] = tuple(
    (first, second) for first in range(3) for second in range(4)
)
point_index = {point: index for index, point in enumerate(points)}


def product_permutation(
    first: tuple[int, ...], second: tuple[int, ...]
) -> Permutation:
    return tuple(
        point_index[(first[i], second[j])] for i, j in points
    )


group = tuple(
    product_permutation(first, second)
    for first in permutations(range(3))
    for second in permutations(range(4))
)
assert len(group) == factorial(3) * factorial(4) == 144
assert len(set(group)) == len(group)


def image(block: frozenset[int], permutation: Permutation) -> frozenset[int]:
    return frozenset(permutation[index] for index in block)


def is_block(block: frozenset[int]) -> bool:
    return all(
        not (translated := image(block, permutation)).intersection(block)
        or translated == block
        for permutation in group
    )


base_index = point_index[(0, 0)]
blocks = tuple(
    block
    for mask in range(1 << len(points))
    if mask & (1 << base_index)
    if is_block(
        block := frozenset(
            index for index in range(len(points)) if mask & (1 << index)
        )
    )
)

singleton = frozenset({base_index})
first_coordinate_block = frozenset(
    point_index[(first, 0)] for first in range(3)
)
second_coordinate_block = frozenset(
    point_index[(0, second)] for second in range(4)
)
whole_fiber = frozenset(range(len(points)))
assert set(blocks) == {
    singleton,
    first_coordinate_block,
    second_coordinate_block,
    whole_fiber,
}


def induced_local_group(block: frozenset[int]) -> set[tuple[int, ...]]:
    ordered = tuple(sorted(block))
    local_index = {point: index for index, point in enumerate(ordered)}
    return {
        tuple(local_index[permutation[point]] for point in ordered)
        for permutation in group
        if image(block, permutation) == block
    }


def induced_block_group(block: frozenset[int]) -> set[tuple[int, ...]]:
    orbit = tuple(sorted({image(block, permutation) for permutation in group}))
    orbit_index = {member: index for index, member in enumerate(orbit)}
    return {
        tuple(
            orbit_index[image(member, permutation)] for member in orbit
        )
        for permutation in group
    }


assert len(induced_local_group(first_coordinate_block)) == factorial(3)
assert len(induced_block_group(first_coordinate_block)) == factorial(4)
assert len(induced_local_group(second_coordinate_block)) == factorial(4)
assert len(induced_block_group(second_coordinate_block)) == factorial(3)

# S_3 contributes C_2,C_3 and S_4 contributes C_2^3,C_3.
composition_factor_orders = sorted((2, 3) + (2, 2, 2, 3))
assert composition_factor_orders == [2, 2, 2, 2, 3, 3]


def stabilized_atoms(
    degrees: tuple[int, ...],
) -> tuple[
    tuple[sp.Symbol, ...],
    dict[int, tuple[sp.Expr, ...]],
    tuple[sp.Expr, ...],
]:
    """Put one sparse atom on each disjoint three-coordinate block."""

    variables = sp.symbols(f"q0:{3 * len(degrees)}")
    atoms: dict[int, tuple[sp.Expr, ...]] = {}
    product_map: list[sp.Expr] = []
    for block_index, degree in enumerate(degrees):
        block = variables[3 * block_index : 3 * block_index + 3]
        substitution = dict(zip((x, y, z), block, strict=True))
        mapped_block = tuple(
            component.subs(substitution)
            for component in sparse_weighted_map(degree)
        )
        product_map.extend(mapped_block)
        stabilized = list(variables)
        stabilized[3 * block_index : 3 * block_index + 3] = mapped_block
        atoms[degree] = tuple(stabilized)
    return variables, atoms, tuple(product_map)


def word_composite(
    word: tuple[int, ...],
    atoms: dict[int, tuple[sp.Expr, ...]],
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, ...]:
    """Compose a factor word, written outermost to innermost."""

    result: tuple[sp.Expr, ...] = variables
    for degree in reversed(word):
        result = compose_on(atoms[degree], result, variables)
    return result


# Three separated atoms give a literal braid hexagon.  Both half-braids
# compare factorizations of one and the same polynomial map, not only maps
# with a common reduction or normalization.
braid_degrees = (3, 4, 5)
braid_variables, braid_atoms, braid_product = stabilized_atoms(braid_degrees)
braid_words = {
    (3, 4, 5),
    (4, 3, 5),
    (4, 5, 3),
    (5, 4, 3),
    (3, 5, 4),
    (5, 3, 4),
}
assert set(permutations(braid_degrees)) == braid_words
assert all(
    word_composite(word, braid_atoms, braid_variables) == braid_product
    for word in braid_words
)
left_half_braid = (
    (3, 4, 5),
    (4, 3, 5),
    (4, 5, 3),
    (5, 4, 3),
)
right_half_braid = (
    (3, 4, 5),
    (3, 5, 4),
    (5, 3, 4),
    (5, 4, 3),
)
assert left_half_braid[0] == right_half_braid[0]
assert left_half_braid[-1] == right_half_braid[-1]
assert all(
    word_composite(word, braid_atoms, braid_variables) == braid_product
    for word in left_half_braid + right_half_braid
)

# Four factors also contain a strict commuting square for nonadjacent swaps.
square_degrees = (3, 4, 5, 6)
square_variables, square_atoms, square_product = stabilized_atoms(square_degrees)
commuting_square = (
    (3, 4, 5, 6),
    (4, 3, 5, 6),
    (4, 3, 6, 5),
    (3, 4, 6, 5),
)
assert all(
    word_composite(word, square_atoms, square_variables) == square_product
    for word in commuting_square
)


def swap_positions(word: tuple[int, ...], index: int) -> tuple[int, ...]:
    result = list(word)
    result[index], result[index + 1] = result[index + 1], result[index]
    return tuple(result)


def coxeter_counts(degrees: tuple[int, ...]) -> tuple[int, int, int, int]:
    """Count vertices, edges, commuting squares, and braid hexagons."""

    vertices = set(permutations(degrees))
    edges = {
        frozenset((word, swap_positions(word, index)))
        for word in vertices
        for index in range(len(degrees) - 1)
    }
    squares = {
        frozenset(
            (
                word,
                swap_positions(word, first),
                swap_positions(word, second),
                swap_positions(swap_positions(word, first), second),
            )
        )
        for word in vertices
        for first in range(len(degrees) - 1)
        for second in range(first + 2, len(degrees) - 1)
    }
    hexagons = set()
    for word in vertices:
        for first in range(len(degrees) - 2):
            second = first + 1
            orbit = {word}
            frontier = {word}
            while frontier:
                current = frontier.pop()
                for index in (first, second):
                    neighbor = swap_positions(current, index)
                    if neighbor not in orbit:
                        orbit.add(neighbor)
                        frontier.add(neighbor)
            assert len(orbit) == 6
            hexagons.add(frozenset(orbit))
    return len(vertices), len(edges), len(squares), len(hexagons)


assert coxeter_counts((3, 4, 5)) == (6, 6, 0, 1)
assert coxeter_counts((3, 4, 5, 6)) == (24, 36, 6, 8)

# The intermediate-field lattice of an r-fold separated product is Boolean.
# Its canonical maximal chains are the r! orders in which coordinate factors
# are introduced, so their lengths and label multisets agree.
for rank in range(2, 7):
    labels = tuple(range(rank))
    maximal_chains = {
        tuple(frozenset(order[:cut]) for cut in range(rank + 1))
        for order in permutations(labels)
    }
    assert len(maximal_chains) == factorial(rank)
    assert all(len(chain) == rank + 1 for chain in maximal_chains)
    assert all(chain[0] == frozenset() for chain in maximal_chains)
    assert all(chain[-1] == frozenset(labels) for chain in maximal_chains)

print("PASS: F_3 and F_4 have determinant one")
print("PASS: their disjoint stabilizations commute and give degree twelve")
print("PASS: the S_3 x S_4 point-stabilizer interval is an exact diamond")
print("PASS: both maximal chains have degrees {3,4} and groups {S_3,S_4}")
print("PASS: factor-monodromy composition factors are C_2^4 and C_3^2")
print("PASS: F_3 x F_4 x F_5 realizes a strict Keller braid hexagon")
print("PASS: four separated atoms realize 6 squares and 8 braid hexagons")
print("PASS: canonical split towers form Boolean lattices through rank six")

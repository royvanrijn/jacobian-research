"""Combinatorial Ritt 2-complexes and their Dickson coefficient maps.

The one-skeleton is the usual graph of normalized complete decompositions:
an edge interchanges two adjacent coprime degrees.  Commuting squares and
Coxeter braid hexagons are retained as two-cells, so paths of elementary
Ritt moves are data until a two-cell identifies them.

This module deliberately separates the finite combinatorial complex from a
coefficient system on it.  ``dickson_vertex_factors`` supplies the
two-parameter Chebyshev/Dickson coefficient map used by the degree-thirty
braid certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from itertools import permutations
from math import comb, gcd, prod
from typing import Iterable

import sympy as sp


Word = tuple[int, ...]


class MoveType(str, Enum):
    """The two elementary cases in Ritt's second theorem."""

    POWER = "power"
    CHEBYSHEV = "chebyshev"


@dataclass(frozen=True)
class CoefficientSpace:
    """An affine coefficient chart attached to a cell of the complex."""

    name: str
    coordinates: tuple[str, ...]

    @property
    def dimension(self) -> int:
        return len(self.coordinates)


@dataclass(frozen=True)
class TransverseCompleteIntersection:
    """A monomial normal form for a transverse Artin slice.

    ``exponents=(a_1,...,a_r)`` denotes

    ``k[u_1,...,u_r]/(u_1**a_1,...,u_r**a_r)``.

    Besides its ordinary Hilbert vector, this records the conormal rank and
    residue-field Tor ranks forced by the regular sequence.  These are the
    first derived decorations needed on a nonreduced Ritt two-cell sector.
    """

    exponents: tuple[int, ...]

    def __post_init__(self) -> None:
        assert self.exponents
        assert all(exponent >= 2 for exponent in self.exponents)

    @property
    def length(self) -> int:
        return prod(self.exponents)

    @property
    def embedding_dimension(self) -> int:
        return len(self.exponents)

    @property
    def hilbert_vector(self) -> tuple[int, ...]:
        coefficients = [1]
        for exponent in self.exponents:
            updated = [0] * (len(coefficients) + exponent - 1)
            for left_degree, left_coefficient in enumerate(coefficients):
                for right_degree in range(exponent):
                    updated[left_degree + right_degree] += left_coefficient
            coefficients = updated
        return tuple(coefficients)

    @property
    def socle_degree(self) -> int:
        return sum(exponent - 1 for exponent in self.exponents)

    @property
    def conormal_rank(self) -> int:
        return len(self.exponents)

    @property
    def augmentation_ideal_length(self) -> int:
        """Length of the kernel of the reduction map to the closed point."""

        return self.length - 1

    @property
    def point_cotangent_homology_ranks(self) -> tuple[int, int]:
        """Ranks of ``H_0,H_1`` of the cotangent complex at the point.

        Every defining monomial has order at least two, so the Jacobian
        differential in the two-term complete-intersection cotangent complex
        vanishes after applying the augmentation.
        """

        rank = len(self.exponents)
        return rank, rank

    @property
    def residue_field_tor_ranks(self) -> tuple[int, ...]:
        """Return the Koszul Betti numbers over the smooth slice ambient."""

        rank = len(self.exponents)
        return tuple(comb(rank, degree) for degree in range(rank + 1))

    def intrinsic_tor_ranks(self, maximum_degree: int) -> tuple[int, ...]:
        """Return ``Tor^A_i(k,k)`` ranks through ``maximum_degree``.

        Here ``A`` is the transverse complete intersection and ``k`` its
        reduced closed point.  Since embedding dimension and codimension are
        both ``r``, Tate's complete-intersection series simplifies to
        ``1/(1-t)^r``.
        """

        assert maximum_degree >= 0
        rank = len(self.exponents)
        return tuple(
            comb(degree + rank - 1, rank - 1)
            for degree in range(maximum_degree + 1)
        )


@dataclass(frozen=True)
class BraidSectorDecoration:
    """Scheme data attached to one labelled half of a braid two-cell."""

    composite_omission: int
    prime_omission: int
    nilpotence_index: int
    conductor_power: int
    transverse_slice: TransverseCompleteIntersection


@dataclass(frozen=True)
class DecompositionVertex:
    """A normalized complete decomposition, ordered outer to inner."""

    word: Word
    coefficients: CoefficientSpace

    @property
    def cuts(self) -> tuple[int, ...]:
        """Return the proper outer-prefix degrees carried by this word."""

        return tuple(prod(self.word[:index]) for index in range(1, len(self.word)))


@dataclass(frozen=True)
class RittEdge:
    """An undirected adjacent interchange and its coefficient correspondence."""

    endpoints: tuple[Word, Word]
    swapped_index: int
    move_type: MoveType
    correspondence: CoefficientSpace

    def other(self, word: Word) -> Word:
        if word == self.endpoints[0]:
            return self.endpoints[1]
        if word == self.endpoints[1]:
            return self.endpoints[0]
        raise KeyError(f"{word} is not incident to {self.endpoints}")


@dataclass(frozen=True)
class RittTwoCell:
    """A relation between two edge paths with the same endpoints."""

    relation: str
    paths: tuple[tuple[Word, ...], tuple[Word, ...]]

    @property
    def boundary_vertices(self) -> frozenset[Word]:
        return frozenset(self.paths[0] + self.paths[1])


@dataclass(frozen=True)
class RittMoveComplex:
    """A finite two-complex of normalized decompositions and Ritt moves."""

    vertices: tuple[DecompositionVertex, ...]
    edges: tuple[RittEdge, ...]
    two_cells: tuple[RittTwoCell, ...]

    @property
    def words(self) -> frozenset[Word]:
        return frozenset(vertex.word for vertex in self.vertices)

    @property
    def euler_characteristic(self) -> int:
        return len(self.vertices) - len(self.edges) + len(self.two_cells)

    def edge_between(self, left: Word, right: Word) -> RittEdge:
        key = frozenset((left, right))
        for edge in self.edges:
            if frozenset(edge.endpoints) == key:
                return edge
        raise KeyError(f"no Ritt edge between {left} and {right}")

    def validate(self) -> None:
        """Check incidence, adjacency, and every declared two-cell boundary."""

        words = self.words
        assert len(words) == len(self.vertices)
        seen_edges: set[frozenset[Word]] = set()
        for edge in self.edges:
            left, right = edge.endpoints
            assert left in words and right in words and left != right
            assert len(left) == len(right)
            differences = [
                index
                for index, pair in enumerate(zip(left, right))
                if pair[0] != pair[1]
            ]
            assert differences == [edge.swapped_index, edge.swapped_index + 1]
            index = edge.swapped_index
            assert left[index] == right[index + 1]
            assert left[index + 1] == right[index]
            assert gcd(left[index], left[index + 1]) == 1
            key = frozenset(edge.endpoints)
            assert key not in seen_edges
            seen_edges.add(key)

        for cell in self.two_cells:
            first, second = cell.paths
            assert first[0] == second[0] and first[-1] == second[-1]
            for path in cell.paths:
                assert len(path) >= 2
                for left, right in zip(path, path[1:]):
                    self.edge_between(left, right)


def normalized_coefficient_space(word: Word) -> CoefficientSpace:
    """Return the affine chart for monic-original factors of ``word``."""

    coordinates = tuple(
        f"c{position}_{power}"
        for position, degree in enumerate(word)
        for power in range(1, degree)
    )
    return CoefficientSpace(f"D_{'_'.join(map(str, word))}", coordinates)


def _swap(word: Word, index: int) -> Word:
    target = list(word)
    target[index], target[index + 1] = target[index + 1], target[index]
    return tuple(target)


def permutation_ritt_complex(
    degrees: Iterable[int],
    move_type: MoveType = MoveType.CHEBYSHEV,
) -> RittMoveComplex:
    """Build the Coxeter 2-complex for distinct pairwise-coprime degrees.

    Its two-cells are the commuting squares ``s_i s_j=s_j s_i`` for
    ``|i-j|>1`` and the braid hexagons
    ``s_i s_(i+1) s_i=s_(i+1) s_i s_(i+1)``.
    """

    base = tuple(degrees)
    assert len(base) >= 2 and len(set(base)) == len(base)
    assert all(
        gcd(base[i], base[j]) == 1
        for i in range(len(base))
        for j in range(i)
    )
    words = tuple(sorted(set(permutations(base))))
    vertices = tuple(
        DecompositionVertex(word, normalized_coefficient_space(word)) for word in words
    )

    edge_records: dict[frozenset[Word], RittEdge] = {}
    for word in words:
        for index in range(len(base) - 1):
            target = _swap(word, index)
            endpoints = tuple(sorted((word, target)))
            key = frozenset(endpoints)
            edge_records[key] = RittEdge(
                endpoints=endpoints,
                swapped_index=index,
                move_type=move_type,
                correspondence=CoefficientSpace(
                    f"{move_type.value}_{endpoints[0]}_{endpoints[1]}",
                    ("t", "z"),
                ),
            )

    cell_records: dict[frozenset[frozenset[Word]], RittTwoCell] = {}

    def add_cell(relation: str, paths: tuple[tuple[Word, ...], tuple[Word, ...]]) -> None:
        boundary_edges = frozenset(
            frozenset((left, right))
            for path in paths
            for left, right in zip(path, path[1:])
        )
        cell_records.setdefault(boundary_edges, RittTwoCell(relation, paths))

    for word in words:
        for left_index in range(len(base) - 1):
            for right_index in range(left_index + 2, len(base) - 1):
                first = (
                    word,
                    _swap(word, left_index),
                    _swap(_swap(word, left_index), right_index),
                )
                second = (
                    word,
                    _swap(word, right_index),
                    _swap(_swap(word, right_index), left_index),
                )
                add_cell("commuting", (first, second))
        for index in range(len(base) - 2):
            first_one = _swap(word, index)
            first_two = _swap(first_one, index + 1)
            second_one = _swap(word, index + 1)
            second_two = _swap(second_one, index)
            paths = (
                (word, first_one, first_two, _swap(first_two, index)),
                (word, second_one, second_two, _swap(second_two, index + 1)),
            )
            add_cell("braid", paths)

    complex_ = RittMoveComplex(
        vertices, tuple(edge_records.values()), tuple(cell_records.values())
    )
    complex_.validate()
    return complex_


def symmetric_braid_complex(
    degrees: Iterable[int],
    move_type: MoveType = MoveType.CHEBYSHEV,
) -> RittMoveComplex:
    """Build the filled ``S_3`` braid hexagon for three coprime degrees."""

    base = tuple(degrees)
    assert len(base) == 3
    complex_ = permutation_ritt_complex(base, move_type)
    assert len(complex_.two_cells) == 1
    assert complex_.two_cells[0].relation == "braid"
    assert complex_.two_cells[0].boundary_vertices == complex_.words
    return complex_


def degree_thirty_braid_decorations() -> tuple[BraidSectorDecoration, ...]:
    """Return the exact three-sector decoration of the ``(2,3,5)`` cell."""

    return (
        BraidSectorDecoration(
            composite_omission=10,
            prime_omission=3,
            nilpotence_index=4,
            conductor_power=2,
            transverse_slice=TransverseCompleteIntersection((5,)),
        ),
        BraidSectorDecoration(
            composite_omission=15,
            prime_omission=2,
            nilpotence_index=3,
            conductor_power=2,
            transverse_slice=TransverseCompleteIntersection((2, 2)),
        ),
        BraidSectorDecoration(
            composite_omission=6,
            prime_omission=5,
            nilpotence_index=4,
            conductor_power=4,
            transverse_slice=TransverseCompleteIntersection((4, 2)),
        ),
    )


def dickson(degree: int, variable: sp.Expr, parameter: sp.Expr) -> sp.Expr:
    """Return the monic Dickson polynomial ``D_degree(variable, parameter)``."""

    assert degree >= 1
    if degree == 1:
        return variable
    previous, current = sp.Integer(2), variable
    for _ in range(2, degree + 1):
        previous, current = current, sp.expand(
            variable * current - parameter * previous
        )
    return current


def dickson_vertex_factors(
    word: Word,
    variable: sp.Symbol,
    translation: sp.Expr,
    parameter: sp.Expr,
) -> tuple[sp.Expr, ...]:
    """Map ``A^2_(translation,parameter)`` to one vertex coefficient chart.

    The returned factors are monic and original, ordered outer to inner.  The
    parameter powers encode

    ``D_m(D_n(x,z), z**n) = D_(mn)(x,z)``.
    """

    factors = []
    for index, degree in enumerate(word):
        suffix_degree = prod(word[index + 1 :])
        center = dickson(suffix_degree, translation, parameter)
        factor_parameter = parameter**suffix_degree
        factor = sp.expand(
            dickson(degree, variable + center, factor_parameter)
            - dickson(degree, center, factor_parameter)
        )
        assert sp.Poly(factor, variable).degree() == degree
        assert sp.Poly(factor, variable).LC() == 1
        assert factor.subs(variable, 0) == 0
        factors.append(factor)
    return tuple(factors)


def compose_factors(
    factors: tuple[sp.Expr, ...], variable: sp.Symbol
) -> sp.Expr:
    """Compose an outer-to-inner tuple of univariate factors."""

    assert factors
    result = factors[-1]
    for factor in reversed(factors[:-1]):
        result = sp.expand(factor.subs(variable, result))
    return result


def composition_differential(
    factors: tuple[sp.Expr, ...],
    variations: tuple[sp.Expr, ...],
    variable: sp.Symbol,
) -> sp.Expr:
    """Differentiate an outer-to-inner polynomial composition.

    This implements the tree-local formula

    ``sum_i (A_i' o f_i o B_i) * (dot(f_i) o B_i)``,

    where ``A_i`` and ``B_i`` are the compositions strictly outside and
    inside the ``i``-th factor.  Leading coefficients and constant terms are
    not treated specially here; callers choose the desired tangent slice by
    choosing ``variations``.
    """

    assert factors and len(factors) == len(variations)
    differential = sp.Integer(0)
    for index, (factor, variation) in enumerate(zip(factors, variations)):
        outer = (
            compose_factors(factors[:index], variable)
            if index
            else variable
        )
        inner = (
            compose_factors(factors[index + 1 :], variable)
            if index + 1 < len(factors)
            else variable
        )
        current_suffix = sp.expand(factor.subs(variable, inner))
        outer_derivative = sp.diff(outer, variable).subs(
            variable, current_suffix
        )
        differential += outer_derivative * variation.subs(variable, inner)
    return sp.expand(differential)

"""Lyashko--Looijenga signatures and prime-word Ritt relation graphs.

For a generically simply ramified normalized composition

``f = f_1 o ... o f_r``

the critical points contributed by ``f_i`` occur in packets of size equal to
the degree of the inner suffix.  The corresponding branch cycle is a product
of that many disjoint transpositions.  This module records those packets and
the finite prime-word/block-system census used by the degree-12 and degree-30
Hessian--Ritt component-completeness argument.

The combinatorics are a verifier for the written Hurwitz/Ritt proof.  They do
not replace the synchronization theorem or the tame multi-collision theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from math import prod


Word = tuple[int, ...]


@dataclass(frozen=True, order=True)
class CriticalValuePacket:
    """Generic branch values with one common transposition multiplicity."""

    transpositions: int
    branch_values: int

    @property
    def ramification_points(self) -> int:
        return self.transpositions * self.branch_values


@dataclass(frozen=True)
class LLSignature:
    """Generic LL branch-cycle signature of a normalized factor word."""

    degree: int
    packets: tuple[CriticalValuePacket, ...]

    @property
    def branch_value_count(self) -> int:
        return sum(packet.branch_values for packet in self.packets)

    @property
    def ramification_point_count(self) -> int:
        return sum(packet.ramification_points for packet in self.packets)

    @property
    def marked_dimension(self) -> int:
        """Dimension after quotienting the one-dimensional target scaling."""

        return self.branch_value_count - 1

    @property
    def expanded_transposition_multiplicities(self) -> tuple[int, ...]:
        return tuple(
            transpositions
            for packet in self.packets
            for transpositions in (packet.transpositions,) * packet.branch_values
        )


def ll_signature(word: Word) -> LLSignature:
    """Return the generic LL signature of an outer-to-inner factor word."""

    if not word or any(degree < 2 for degree in word):
        raise ValueError("factor degrees must all be at least two")
    packets = []
    for index, degree in enumerate(word):
        suffix_degree = prod(word[index + 1 :])
        packets.append(
            CriticalValuePacket(
                transpositions=suffix_degree,
                branch_values=degree - 1,
            )
        )
    signature = LLSignature(prod(word), tuple(packets))
    assert signature.ramification_point_count == signature.degree - 1
    return signature


@dataclass(frozen=True, order=True)
class PrimeOccurrence:
    """One labelled occurrence in the prime multiset of a total degree."""

    degree: int
    ordinal: int

    @property
    def label(self) -> str:
        return (
            str(self.degree)
            if self.ordinal == 1
            else f"{self.degree}_{self.ordinal}"
        )


OccurrenceWord = tuple[PrimeOccurrence, ...]
DirectedEdge = tuple[PrimeOccurrence, PrimeOccurrence]
UndirectedEdge = frozenset[PrimeOccurrence]


@dataclass(frozen=True)
class RelationGraphType:
    """The orientation data forced by a selection of prime refinements."""

    vertices: tuple[PrimeOccurrence, ...]
    two_way_edges: frozenset[UndirectedEdge]
    one_way_edges: frozenset[DirectedEdge]

    def strongly_connected_components(
        self,
    ) -> tuple[frozenset[PrimeOccurrence], ...]:
        adjacency = {vertex: set() for vertex in self.vertices}
        reverse = {vertex: set() for vertex in self.vertices}
        for edge in self.two_way_edges:
            left, right = tuple(edge)
            adjacency[left].add(right)
            adjacency[right].add(left)
            reverse[left].add(right)
            reverse[right].add(left)
        for left, right in self.one_way_edges:
            adjacency[left].add(right)
            reverse[right].add(left)

        visited: set[PrimeOccurrence] = set()
        finish_order: list[PrimeOccurrence] = []

        def finish(vertex: PrimeOccurrence) -> None:
            visited.add(vertex)
            for neighbor in sorted(adjacency[vertex]):
                if neighbor not in visited:
                    finish(neighbor)
            finish_order.append(vertex)

        for vertex in self.vertices:
            if vertex not in visited:
                finish(vertex)

        components = []
        visited.clear()
        for start in reversed(finish_order):
            if start in visited:
                continue
            stack = [start]
            visited.add(start)
            component = set()
            while stack:
                vertex = stack.pop()
                component.add(vertex)
                for neighbor in reverse[vertex]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            components.append(frozenset(component))
        return tuple(sorted(components, key=lambda component: min(component)))


def prime_factors(degree: int) -> tuple[int, ...]:
    """Return the prime factors with multiplicity."""

    if degree < 2:
        raise ValueError("degree must be at least two")
    factors = []
    divisor = 2
    remainder = degree
    while divisor * divisor <= remainder:
        while remainder % divisor == 0:
            factors.append(divisor)
            remainder //= divisor
        divisor += 1
    if remainder > 1:
        factors.append(remainder)
    return tuple(factors)


def prime_occurrences(degree: int) -> tuple[PrimeOccurrence, ...]:
    """Label repeated prime factors in their fixed relative order."""

    counts: dict[int, int] = {}
    result = []
    for factor in prime_factors(degree):
        counts[factor] = counts.get(factor, 0) + 1
        result.append(PrimeOccurrence(factor, counts[factor]))
    return tuple(result)


def prime_words(degree: int) -> tuple[OccurrenceWord, ...]:
    """Return prime words, identifying permutations of equal occurrences.

    Equal-degree factors cannot be exchanged by a coprime Ritt move.  We
    therefore keep their labelled occurrences in ordinal order.
    """

    occurrences = prime_occurrences(degree)
    words = {
        word
        for word in permutations(occurrences)
        if all(
            [
                occurrence.ordinal
                for occurrence in word
                if occurrence.degree == prime
            ]
            == sorted(
                occurrence.ordinal
                for occurrence in word
                if occurrence.degree == prime
            )
            for prime in set(prime_factors(degree))
        )
    }
    return tuple(sorted(words))


def outer_cuts(word: OccurrenceWord) -> frozenset[int]:
    return frozenset(
        prod(occurrence.degree for occurrence in word[:index])
        for index in range(1, len(word))
    )


def proper_outer_cuts(degree: int) -> tuple[int, ...]:
    return tuple(
        divisor
        for divisor in range(2, degree)
        if degree % divisor == 0
    )


def cut_carriers(
    degree: int, requested_cuts: tuple[int, ...]
) -> dict[int, tuple[OccurrenceWord, ...]]:
    words = prime_words(degree)
    carriers = {
        cut: tuple(word for word in words if cut in outer_cuts(word))
        for cut in requested_cuts
    }
    if any(not words_for_cut for words_for_cut in carriers.values()):
        raise ValueError("a requested cut has no prime refinement")
    return carriers


def relation_graph_type(words: tuple[OccurrenceWord, ...]) -> RelationGraphType:
    if not words:
        raise ValueError("at least one word is required")
    vertices = tuple(sorted(set(words[0])))
    if any(set(word) != set(vertices) for word in words):
        raise ValueError("all words must permute the same occurrences")

    orientations: set[DirectedEdge] = set()
    for word in words:
        for index, left in enumerate(word):
            for right in word[index + 1 :]:
                orientations.add((left, right))

    two_way = frozenset(
        frozenset((left, right))
        for left, right in orientations
        if (right, left) in orientations
    )
    one_way = frozenset(
        (left, right)
        for left, right in orientations
        if (right, left) not in orientations
    )
    return RelationGraphType(vertices, two_way, one_way)


def selected_relation_graph_types(
    degree: int, requested_cuts: tuple[int, ...]
) -> dict[RelationGraphType, int]:
    """Enumerate graph types over all choices of one carrier per cut."""

    carriers = cut_carriers(degree, requested_cuts)
    result: dict[RelationGraphType, int] = {}
    for selection in product(*(carriers[cut] for cut in requested_cuts)):
        graph = relation_graph_type(selection)
        result[graph] = result.get(graph, 0) + 1
    return result

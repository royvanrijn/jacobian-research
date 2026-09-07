"""Exact additive relation complexes on primitive unoriented vectors."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from typing import Mapping, Sequence

from .integer import (
    add,
    canonical_unoriented,
    content,
    rational_rank,
    row_basis_coordinates,
    subtract,
)


@dataclass(frozen=True)
class RelationComplex:
    """A finite, ordering-independent additive complex.

    Vertices are primitive unoriented rays.  A ternary relation is an
    unordered triple of ray indices for which some choice of signs sums to
    zero with unit content.  Non-unit identities ``a +/- b = m c`` are kept
    separately, including the multiplier and source/target roles.  Forgetting
    signs is essential: a ray has no preferred orientation, so the distinction
    between ``a+b=c`` and ``a-b=c`` is not invariant.  The representation is
    exact for the supplied finite vertex set; it makes no completeness claim
    outside that set.
    """

    vertices: tuple[tuple[int, ...], ...]
    ternary_relations: tuple[tuple[int, int, int], ...]
    scaled_relations: tuple[tuple[int, int, int, int], ...]
    additive_degrees: tuple[int, ...]
    divisibility_degrees: tuple[int, ...]
    metadata: tuple[Mapping[str, object], ...]
    canonical_digest: str

    def wl_colors(self, rounds: int = 4) -> tuple[str, ...]:
        """Return coordinate-free hypergraph color-refinement labels.

        Equal labels are necessary, not sufficient, for vertices to
        correspond under an isomorphism.  The resulting multiset is a safe
        isomorphism invariant; no complete canonical labelling is claimed.
        """

        if rounds < 0:
            raise ValueError("rounds must be nonnegative")
        colors = [
            sha256(_stable_metadata(item).encode()).hexdigest()
            for item in self.metadata
        ]
        incidence: list[list[tuple[str, int, int, int]]] = [
            [] for _ in self.vertices
        ]
        for left, middle, right in self.ternary_relations:
            incidence[left].append(("unit", 1, middle, right))
            incidence[middle].append(("unit", 1, left, right))
            incidence[right].append(("unit", 1, left, middle))
        for left, right, target, multiplier in self.scaled_relations:
            incidence[left].append(("source", multiplier, right, target))
            incidence[right].append(("source", multiplier, left, target))
            incidence[target].append(("target", multiplier, left, right))
        for _round in range(rounds):
            descriptions = []
            for vertex, neighbours in enumerate(incidence):
                edge_colors = sorted(
                    (
                        role,
                        multiplier,
                        *sorted((colors[left], colors[right])),
                    )
                    for role, multiplier, left, right in neighbours
                )
                descriptions.append(
                    json.dumps(
                        [colors[vertex], edge_colors],
                        separators=(",", ":"),
                    )
                )
            palette = {
                description: sha256(description.encode()).hexdigest()
                for description in sorted(set(descriptions))
            }
            colors = [palette[description] for description in descriptions]
        return tuple(colors)

    def wl_profile(self, rounds: int = 4) -> tuple[tuple[str, int], ...]:
        """Return sorted color multiplicities after refinement."""

        return tuple(sorted(Counter(self.wl_colors(rounds)).items()))

    def to_record(self, *, include_relations: bool = True) -> dict[str, object]:
        record: dict[str, object] = {
            "vertex_count": len(self.vertices),
            "vertices": [list(vertex) for vertex in self.vertices],
            "additive_degrees": list(self.additive_degrees),
            "divisibility_degrees": list(self.divisibility_degrees),
            "metadata": [dict(item) for item in self.metadata],
            "ternary_relation_count": len(self.ternary_relations),
            "scaled_relation_count": len(self.scaled_relations),
            "canonical_digest": self.canonical_digest,
        }
        if include_relations:
            record["ternary_relations"] = [list(edge) for edge in self.ternary_relations]
            record["scaled_relations"] = [list(edge) for edge in self.scaled_relations]
        return record


def _stable_metadata(metadata: Mapping[str, object]) -> str:
    return json.dumps(metadata, sort_keys=True, separators=(",", ":"), default=str)


def build_relation_complex(
    vectors: Sequence[Sequence[int]],
    metadata: Sequence[Mapping[str, object]] | None = None,
    *,
    pair_limit: int | None = None,
) -> RelationComplex:
    """Build all visible ``a+b=c`` and ``a-b=c`` relations in a declared box.

    ``pair_limit`` is an explicit computational boundary.  When supplied,
    only pairs among the first ``pair_limit`` canonical vertices are tested,
    while the result still records all supplied vertices.
    """

    canonical = sorted({canonical_unoriented(vector) for vector in vectors})
    if not canonical:
        raise ValueError("a relation complex needs at least one nonzero vector")
    width = len(canonical[0])
    if any(len(vector) != width for vector in canonical):
        raise ValueError("vector widths differ")
    metadata_by_vertex: dict[tuple[int, ...], Mapping[str, object]] = {}
    if metadata is not None:
        if len(metadata) != len(vectors):
            raise ValueError("metadata and vector counts differ")
        for vector, item in zip(vectors, metadata):
            key = canonical_unoriented(vector)
            previous = metadata_by_vertex.get(key)
            if previous is not None and _stable_metadata(previous) != _stable_metadata(item):
                raise ValueError("opposite/duplicate vectors have inconsistent metadata")
            metadata_by_vertex[key] = dict(item)
    ordered_metadata = tuple(metadata_by_vertex.get(vertex, {}) for vertex in canonical)
    index = {vertex: position for position, vertex in enumerate(canonical)}
    stop = len(canonical) if pair_limit is None else min(int(pair_limit), len(canonical))
    if stop < 0:
        raise ValueError("pair_limit must be nonnegative")
    relations: set[tuple[int, int, int]] = set()
    scaled_relations: set[tuple[int, int, int, int]] = set()
    degrees = [0] * len(canonical)
    divisibility_degrees = [0] * len(canonical)
    for left in range(stop):
        a_vector = canonical[left]
        for right in range(left):
            b_vector = canonical[right]
            for result in (add(a_vector, b_vector), subtract(a_vector, b_vector)):
                if not any(result):
                    continue
                multiplier = content(result)
                target = index.get(canonical_unoriented(result))
                if target is None:
                    continue
                if multiplier == 1:
                    edge = tuple(sorted((right, left, target)))
                    if edge in relations:
                        continue
                    relations.add(edge)
                    for vertex in edge:
                        degrees[vertex] += 1
                else:
                    source_left, source_right = sorted((right, left))
                    edge = (source_left, source_right, target, multiplier)
                    if edge in scaled_relations:
                        continue
                    scaled_relations.add(edge)
                    divisibility_degrees[source_left] += 1
                    divisibility_degrees[source_right] += 1
                    divisibility_degrees[target] += 1
    ordered_relations = tuple(sorted(relations))
    ordered_scaled = tuple(sorted(scaled_relations))
    provisional = RelationComplex(
        vertices=tuple(canonical),
        ternary_relations=ordered_relations,
        scaled_relations=ordered_scaled,
        additive_degrees=tuple(degrees),
        divisibility_degrees=tuple(divisibility_degrees),
        metadata=ordered_metadata,
        canonical_digest="",
    )
    invariant_payload = {
        "vertex_count": len(canonical),
        "edge_count": len(ordered_relations),
        "scaled_edge_count": len(ordered_scaled),
        "degree_multiset": sorted(degrees),
        "divisibility_degree_multiset": sorted(divisibility_degrees),
        "scaled_multiplier_multiset": sorted(edge[3] for edge in ordered_scaled),
        "wl_profile_round_4": provisional.wl_profile(4),
        "pair_limit": stop,
    }
    digest = sha256(
        json.dumps(invariant_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return RelationComplex(
        vertices=provisional.vertices,
        ternary_relations=provisional.ternary_relations,
        scaled_relations=provisional.scaled_relations,
        additive_degrees=provisional.additive_degrees,
        divisibility_degrees=provisional.divisibility_degrees,
        metadata=provisional.metadata,
        canonical_digest=digest,
    )


def _ternary_signs(
    vertices: Sequence[Sequence[int]], edge: Sequence[int]
) -> tuple[int, int, int]:
    vectors = [vertices[index] for index in edge]
    for signs in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)):
        if all(
            sum(sign * int(vector[column]) for sign, vector in zip(signs, vectors))
            == 0
            for column in range(len(vectors[0]))
        ):
            return signs
    raise ArithmeticError("declared ternary edge has no signed unit relation")


def _fraction_determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    work = [list(row) for row in matrix]
    answer = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (index for index in range(column, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        value = work[column][column]
        answer *= value
        for index in range(column + 1, len(work)):
            if not work[index][column]:
                continue
            multiplier = work[index][column] / value
            for offset in range(column, len(work)):
                work[index][offset] -= multiplier * work[column][offset]
    return answer


def lift_relation_vertex_bijection(
    source: RelationComplex,
    target: RelationComplex,
    vertex_map: Sequence[int],
    *,
    maximum_sign_components: int = 16,
    _allow_injection: bool = False,
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Lift an unoriented hypergraph bijection to certified integral maps.

    The returned matrices act on row vectors on the right.  Ternary relations
    propagate the unknown ray orientations; remaining connected-component
    signs are exhaustively enumerated within the declared bound.  Every
    survivor is required to lie in ``GL(k,Z)`` and is replayed on all rays.
    """

    if not _allow_injection and len(source.vertices) != len(target.vertices):
        raise ValueError("relation complexes have different vertex counts")
    mapping = tuple(map(int, vertex_map))
    if len(mapping) != len(source.vertices) or len(set(mapping)) != len(mapping):
        raise ValueError("vertex_map is not an injection")
    if any(index < 0 or index >= len(target.vertices) for index in mapping):
        raise ValueError("vertex_map target is outside the relation complex")
    if not _allow_injection and set(mapping) != set(range(len(target.vertices))):
        raise ValueError("vertex_map is not a bijection")
    source_edges = set(source.ternary_relations)
    target_edges = set(target.ternary_relations)
    mapped_edges = {
        tuple(sorted(mapping[index] for index in edge)) for edge in source_edges
    }
    if (not _allow_injection and mapped_edges != target_edges) or (
        _allow_injection and not mapped_edges <= target_edges
    ):
        raise ValueError("vertex_map does not preserve ternary hyperedges")
    rank = rational_rank(source.vertices)
    mapped_target_vertices = tuple(target.vertices[index] for index in mapping)
    if rank != rational_rank(mapped_target_vertices):
        raise ValueError("mapped relation vertices have different vector ranks")
    if rank != len(source.vertices[0]):
        raise ValueError("source relation vertices must span their ambient lattice")
    if not _allow_injection:
        if rank != rational_rank(target.vertices):
            raise ValueError("relation complexes have different vector ranks")
        if rank != len(target.vertices[0]):
            raise ValueError("target relation vertices must span its ambient lattice")

    constraints: list[list[tuple[int, int]]] = [
        [] for _ in source.vertices
    ]
    target_edge_by_set = {frozenset(edge): edge for edge in target.ternary_relations}
    for source_edge in source.ternary_relations:
        source_signs = _ternary_signs(source.vertices, source_edge)
        mapped_set = frozenset(mapping[index] for index in source_edge)
        target_edge = target_edge_by_set[mapped_set]
        target_sign_tuple = _ternary_signs(target.vertices, target_edge)
        target_signs = dict(zip(target_edge, target_sign_tuple))
        ratios = {
            vertex: source_sign * target_signs[mapping[vertex]]
            for vertex, source_sign in zip(source_edge, source_signs)
        }
        anchor = source_edge[0]
        for vertex in source_edge[1:]:
            ratio = ratios[vertex] * ratios[anchor]
            constraints[anchor].append((vertex, ratio))
            constraints[vertex].append((anchor, ratio))

    relative: list[int | None] = [None] * len(source.vertices)
    components: list[tuple[int, ...]] = []
    for start in range(len(source.vertices)):
        if relative[start] is not None:
            continue
        relative[start] = 1
        stack = [start]
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbour, ratio in constraints[vertex]:
                proposed = int(relative[vertex]) * ratio
                if relative[neighbour] is None:
                    relative[neighbour] = proposed
                    stack.append(neighbour)
                elif relative[neighbour] != proposed:
                    raise ValueError("edge signs are inconsistent under vertex_map")
        components.append(tuple(component))
    if len(components) > maximum_sign_components:
        raise ValueError("relation sign-component bound exceeded")

    basis_indices = []
    basis_vectors = []
    for index, vector in enumerate(source.vertices):
        if rational_rank(tuple(basis_vectors) + (vector,)) == len(basis_vectors) + 1:
            basis_indices.append(index)
            basis_vectors.append(vector)
            if len(basis_vectors) == rank:
                break
    source_basis = tuple(basis_vectors)
    source_transpose = tuple(tuple(column) for column in zip(*source_basis))
    answers = set()
    for component_signs in product((-1, 1), repeat=len(components)):
        orientations = [int(value) for value in relative]
        for component, sign in zip(components, component_signs):
            for vertex in component:
                orientations[vertex] *= sign
        target_basis = tuple(
            tuple(
                orientations[index] * value
                for value in target.vertices[mapping[index]]
            )
            for index in basis_indices
        )
        target_columns = tuple(tuple(column) for column in zip(*target_basis))
        map_columns = row_basis_coordinates(
            target_columns, source_transpose, require_integral=False
        )
        target_width = len(target.vertices[0])
        matrix_fraction = tuple(
            tuple(map_columns[column][row] for column in range(target_width))
            for row in range(rank)
        )
        if any(value.denominator != 1 for row in matrix_fraction for value in row):
            continue
        if not _allow_injection and abs(_fraction_determinant(matrix_fraction)) != 1:
            continue
        matrix = tuple(tuple(int(value) for value in row) for row in matrix_fraction)
        if rational_rank(matrix) != rank:
            continue
        if any(
            canonical_unoriented(
                tuple(
                    sum(int(vector[row]) * matrix[row][column] for row in range(rank))
                    for column in range(target_width)
                )
            )
            != target.vertices[mapping[index]]
            for index, vector in enumerate(source.vertices)
        ):
            continue
        answers.add(matrix)
    return tuple(sorted(answers))


def lift_relation_vertex_injection(
    source: RelationComplex,
    target: RelationComplex,
    vertex_map: Sequence[int],
    *,
    maximum_sign_components: int = 16,
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Lift a full-rank relation subcomplex injection to integral maps.

    The source must span its coordinate lattice, while the target may be a
    larger cloud in a higher-rank displayed subgroup.  Returned matrices have
    shape ``source_rank x target_ambient_rank`` and are replayed exactly on
    every mapped ray.  This routine certifies integrality and full row rank;
    primitivity of the rectangular image is a separate Smith-form check.
    """

    return lift_relation_vertex_bijection(
        source,
        target,
        vertex_map,
        maximum_sign_components=maximum_sign_components,
        _allow_injection=True,
    )

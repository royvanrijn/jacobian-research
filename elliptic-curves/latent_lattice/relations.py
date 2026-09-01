"""Exact additive relation complexes on primitive unoriented vectors."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
import json
from typing import Mapping, Sequence

from .integer import add, canonical_unoriented, content, subtract


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

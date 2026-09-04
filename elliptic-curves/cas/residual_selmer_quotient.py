#!/usr/bin/env python3
"""Small, BNF-free bookkeeping for a residual 2-Selmer computation.

This module intentionally does *not* compute a class group or assert a Selmer
bound.  It is the exact GF(2) layer between a relation collector and local
descent code.  In particular it keeps a principal generator with each
relation, eliminates relation rows sparsely, and projects candidate global
squareclasses through the local/fingerprint/Mordell--Weil quotient as soon as
they are available.

The code is dependency-free so it can audit collector output independently of
Sage, PARI, Magma, or a BNF computation.  Field-specific code must establish
the ideal valuations and the local images before constructing these records.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping


class F2Error(ValueError):
    """An invalid packed GF(2) vector or unsupported certification claim."""


def _check_mask(mask: int, dimension: int, name: str) -> int:
    if dimension < 0 or mask < 0 or mask >> dimension:
        raise F2Error(f"{name} is not a {dimension}-dimensional GF(2) mask")
    return mask


def mask_support(mask: int) -> tuple[int, ...]:
    """Return the set-bit coordinates of a packed GF(2) vector."""

    return tuple(index for index in range(mask.bit_length()) if mask >> index & 1)


def independent_masks(masks: Iterable[tuple[str, int]]) -> tuple[int, tuple[str, ...]]:
    """Return the GF(2) rank and a deterministic independent label subset."""

    pivots: dict[int, int] = {}
    labels: list[str] = []
    for label, mask in masks:
        reduced = mask
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot in pivots:
                reduced ^= pivots[pivot]
            else:
                pivots[pivot] = reduced
                labels.append(label)
                break
    return len(pivots), tuple(labels)


def f2_row_basis(rows: Iterable[Iterable[int]], width: int) -> list[list[int]]:
    """Return a deterministic independent basis of row vectors over GF(2)."""

    if width < 0:
        raise F2Error("binary-vector width must be nonnegative")
    pivots: dict[int, int] = {}
    for row_index, row in enumerate(rows):
        values = [int(value) for value in row]
        if len(values) != width:
            raise F2Error(f"binary row {row_index} has width {len(values)}, not {width}")
        if any(value not in (0, 1) for value in values):
            raise F2Error(f"binary row {row_index} contains a non-binary entry")
        mask = sum(value << index for index, value in enumerate(values))
        reduced = mask
        while reduced:
            pivot = reduced.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = reduced
                break
            reduced ^= previous
    return [
        [int(mask >> index & 1) for index in range(width)]
        for _pivot, mask in sorted(pivots.items(), reverse=True)
    ]


def f2_rank_rows(rows: Iterable[Iterable[int]], width: int) -> int:
    return len(f2_row_basis(rows, width))


def f2_dot(left: Iterable[int], right: Iterable[int]) -> int:
    left_values = [int(value) for value in left]
    right_values = [int(value) for value in right]
    if len(left_values) != len(right_values):
        raise F2Error("binary dot-product vectors have different widths")
    if any(value not in (0, 1) for value in [*left_values, *right_values]):
        raise F2Error("binary dot-product vector contains a non-binary entry")
    return sum(a * b for a, b in zip(left_values, right_values)) & 1


def f2_nullspace_basis(rows: Iterable[Iterable[int]], width: int) -> list[list[int]]:
    """Return a basis for vectors annihilated by the supplied row matrix."""

    masks = [sum(bit << index for index, bit in enumerate(row))
             for row in f2_row_basis(rows, width)]
    # The basis is in echelon form with distinct highest pivots.  Each pivot
    # equation only contains lower-numbered coordinates, so solve in ascending
    # pivot order after choosing one free coordinate.
    pivot_rows = {mask.bit_length() - 1: mask for mask in masks}
    pivot_columns = set(pivot_rows)
    free_columns = [index for index in range(width) if index not in pivot_columns]
    result = []
    for free in free_columns:
        vector_mask = 1 << free
        for pivot in sorted(pivot_rows):
            row_without_pivot = pivot_rows[pivot] ^ (1 << pivot)
            if (row_without_pivot & vector_mask).bit_count() & 1:
                vector_mask |= 1 << pivot
        result.append(
            [int(vector_mask >> index & 1) for index in range(width)]
        )
    return result


def f2_standard_complement(
    subspace_rows: Iterable[Iterable[int]], width: int
) -> tuple[list[list[int]], list[list[int]]]:
    """Return a basis of a subspace and a standard-coordinate complement."""

    subspace = f2_row_basis(subspace_rows, width)
    combined = list(subspace)
    complement = []
    rank = len(subspace)
    for index in range(width):
        unit = [int(column == index) for column in range(width)]
        new_rank = f2_rank_rows([*combined, unit], width)
        if new_rank > rank:
            complement.append(unit)
            combined.append(unit)
            rank = new_rank
        if rank == width:
            break
    if rank != width:
        raise F2Error("failed to extend the binary subspace to the ambient space")
    return subspace, complement


def build_relative_local_condition_matrix(
    *,
    ambient_dimension: int,
    known_mw_rows: Iterable[Iterable[int]],
    places: Iterable[Mapping[str, object]],
    maximum_cut_size: int = 8,
    maximum_cut_subsets: int = 1_000_000,
) -> Mapping[str, object]:
    """Quotient known Mordell--Weil directions before applying local conditions.

    The ambient coordinates must already incorporate the global norm
    condition. Each place supplies a basis of its allowed subspace in those
    coordinates. The resulting matrix acts only on a complement to the known
    Mordell--Weil space, so its kernel is precisely the unexplained excess
    inside the supplied ambient envelope.
    """

    if ambient_dimension < 0 or maximum_cut_size < 0 or maximum_cut_subsets < 0:
        raise F2Error("relative dimensions and cut bounds must be nonnegative")
    known_input = [[int(value) for value in row] for row in known_mw_rows]
    if any(len(row) != ambient_dimension for row in known_input):
        raise F2Error("a known Mordell--Weil row has the wrong ambient width")
    known_basis = f2_row_basis(known_input, ambient_dimension)
    if len(known_basis) != len(known_input):
        raise F2Error("known Mordell--Weil Kummer rows are not independent")
    _known, complement = f2_standard_complement(known_basis, ambient_dimension)
    residual_dimension = len(complement)

    place_records = []
    for place_index, place in enumerate(places):
        if not isinstance(place, Mapping):
            raise F2Error(f"place {place_index} is not an object")
        label = str(place.get("place", f"place-{place_index}"))
        allowed_input = place.get("allowed_subspace_basis", ())
        if isinstance(allowed_input, (str, bytes, Mapping)) or not isinstance(
            allowed_input, Iterable
        ):
            raise F2Error(f"place {label} has no allowed-subspace basis")
        allowed = f2_row_basis(allowed_input, ambient_dimension)
        annihilator = f2_nullspace_basis(allowed, ambient_dimension)
        for known_index, known_row in enumerate(known_basis):
            if any(f2_dot(functional, known_row) for functional in annihilator):
                raise F2Error(
                    f"known Mordell--Weil row {known_index} violates place {label}"
                )
        residual_rows = f2_row_basis(
            (
                [f2_dot(functional, vector) for vector in complement]
                for functional in annihilator
            ),
            residual_dimension,
        )
        place_records.append(
            {
                "place": label,
                "allowed_subspace_dimension_in_global_ambient": len(allowed),
                "condition_rank_on_global_ambient": len(annihilator),
                "condition_rows_on_mw_quotient": residual_rows,
                "condition_rank_on_mw_quotient": len(residual_rows),
                "kernel_dimension_on_mw_quotient_for_this_place_alone": (
                    residual_dimension - len(residual_rows)
                ),
            }
        )

    def stacked_rows(indices: Iterable[int]) -> list[list[int]]:
        return [
            row
            for index in indices
            for row in place_records[index]["condition_rows_on_mw_quotient"]
        ]

    all_indices = tuple(range(len(place_records)))
    full_rank = f2_rank_rows(stacked_rows(all_indices), residual_dimension)
    full_kernel_dimension = residual_dimension - full_rank
    full_rows = f2_row_basis(stacked_rows(all_indices), residual_dimension)
    kernel_basis = f2_nullspace_basis(full_rows, residual_dimension)

    def quotient_vector_in_ambient(vector: Iterable[int]) -> list[int]:
        coefficients = [int(value) & 1 for value in vector]
        if len(coefficients) != residual_dimension:
            raise F2Error("relative vector has the wrong width")
        return [
            sum(
                coefficient * complement[index][coordinate]
                for index, coefficient in enumerate(coefficients)
            )
            & 1
            for coordinate in range(ambient_dimension)
        ]

    greedy_order = []
    chosen: list[int] = []
    remaining = list(all_indices)
    current_rank = 0
    while remaining:
        scored = []
        for index in remaining:
            candidate_rank = f2_rank_rows(
                stacked_rows([*chosen, index]), residual_dimension
            )
            scored.append((candidate_rank - current_rank, -index, candidate_rank, index))
        gain, _negative_index, candidate_rank, selected = max(scored)
        chosen.append(selected)
        remaining.remove(selected)
        current_rank = candidate_rank
        greedy_order.append(
            {
                "step": len(chosen),
                "place": place_records[selected]["place"],
                "rank_gain": gain,
                "cumulative_rank": current_rank,
                "remaining_kernel_dimension": residual_dimension - current_rank,
            }
        )

    annihilating_prefix = None
    if full_kernel_dimension == 0:
        for row in greedy_order:
            if row["remaining_kernel_dimension"] == 0:
                annihilating_prefix = {
                    "size": row["step"],
                    "places": [
                        entry["place"] for entry in greedy_order[: int(row["step"])]
                    ],
                }
                break

    delete_one = []
    for omitted in all_indices:
        rank = f2_rank_rows(
            stacked_rows(index for index in all_indices if index != omitted),
            residual_dimension,
        )
        delete_one.append(
            {
                "deleted_place": place_records[omitted]["place"],
                "matrix_rank": rank,
                "rank_drop": full_rank - rank,
                "kernel_dimension": residual_dimension - rank,
            }
        )

    pairwise = []
    for left, right in combinations(all_indices, 2):
        left_rank = int(place_records[left]["condition_rank_on_mw_quotient"])
        right_rank = int(place_records[right]["condition_rank_on_mw_quotient"])
        joint_rank = f2_rank_rows(stacked_rows((left, right)), residual_dimension)
        pairwise.append(
            {
                "left_place": place_records[left]["place"],
                "right_place": place_records[right]["place"],
                "joint_condition_rank": joint_rank,
                "condition_rowspace_intersection_dimension": (
                    left_rank + right_rank - joint_rank
                ),
                "joint_kernel_dimension": residual_dimension - joint_rank,
            }
        )

    minimum_cut = None
    cut_search_limit = min(maximum_cut_size, len(place_records))
    cut_search_complete_through = -1
    cut_subsets_examined = 0
    cut_search_truncated = False
    if full_kernel_dimension == 0:
        for size in range(cut_search_limit + 1):
            completed_size = True
            for subset in combinations(all_indices, size):
                if cut_subsets_examined >= maximum_cut_subsets:
                    completed_size = False
                    cut_search_truncated = True
                    break
                cut_subsets_examined += 1
                if f2_rank_rows(stacked_rows(subset), residual_dimension) == residual_dimension:
                    minimum_cut = {
                        "size": size,
                        "places": [place_records[index]["place"] for index in subset],
                        "minimality_proved": True,
                    }
                    break
            if minimum_cut is not None:
                cut_search_complete_through = size
                break
            if not completed_size:
                break
            cut_search_complete_through = size
    else:
        # Rank is monotone under adding places, so the full-set failure proves
        # that no annihilating subset exists among the supplied places.
        cut_search_complete_through = cut_search_limit

    return {
        "ambient_norm_square_dimension": ambient_dimension,
        "known_mw_kummer_dimension": len(known_basis),
        "mw_quotient_ambient_dimension": residual_dimension,
        "known_mw_basis_in_ambient_coordinates": known_basis,
        "mw_quotient_complement_basis_in_ambient_coordinates": complement,
        "places": place_records,
        "full_relative_local_condition_matrix_rows": full_rows,
        "full_relative_local_condition_matrix_rank": full_rank,
        "unexplained_selmer_excess_kernel_dimension": full_kernel_dimension,
        "unexplained_kernel_basis_in_mw_quotient_coordinates": kernel_basis,
        "unexplained_kernel_basis_in_global_ambient_coordinates": [
            quotient_vector_in_ambient(vector) for vector in kernel_basis
        ],
        "relative_rank_nullity_verified": full_rank + full_kernel_dimension
        == residual_dimension,
        "greedy_place_order": greedy_order,
        "greedy_annihilating_prefix": annihilating_prefix,
        "delete_one_place_ranks": delete_one,
        "pairwise_place_condition_intersections": pairwise,
        "minimum_annihilating_place_cut": minimum_cut,
        "minimum_cut_search_complete_through_size": cut_search_complete_through,
        "minimum_cut_search_requested_through_size": cut_search_limit,
        "minimum_cut_subsets_examined": cut_subsets_examined,
        "minimum_cut_search_truncated_by_subset_budget": cut_search_truncated,
        "annihilating_place_cut_impossible_for_supplied_places": (
            full_kernel_dimension != 0
        ),
        "all_residual_candidates_annihilated": full_kernel_dimension == 0,
    }


@dataclass(frozen=True)
class PrincipalRelation:
    """An exact principal ideal relation, stored with its actual generator.

    ``generator`` should be a reproducible expression (for example power-basis
    coordinates in the maximal-order basis), not merely a valuation vector.
    ``ideal_valuations`` is its valuation parity vector in the collector's
    explicitly ordered prime-ideal factor base.
    """

    label: str
    generator: str
    ideal_valuations: int


@dataclass(frozen=True)
class PrincipalDependency:
    """A mod-two dependency and its compact product of principal generators."""

    relation_labels: tuple[str, ...]
    generator_product: tuple[str, ...]


class SparseF2Relations:
    """Incremental sparse elimination which preserves dependency witnesses."""

    def __init__(self, ideal_dimension: int) -> None:
        self.ideal_dimension = ideal_dimension
        self._pivots: dict[int, tuple[int, int]] = {}
        self._relations: list[PrincipalRelation] = []

    @property
    def rank(self) -> int:
        return len(self._pivots)

    @property
    def relations(self) -> tuple[PrincipalRelation, ...]:
        return tuple(self._relations)

    def add(self, relation: PrincipalRelation) -> PrincipalDependency | None:
        """Add one relation; return a generator witness if it is dependent."""

        row = _check_mask(
            relation.ideal_valuations, self.ideal_dimension, "ideal_valuations"
        )
        relation_index = len(self._relations)
        self._relations.append(relation)
        combination = 1 << relation_index

        while row:
            pivot = row.bit_length() - 1
            previous = self._pivots.get(pivot)
            if previous is None:
                self._pivots[pivot] = (row, combination)
                return None
            previous_row, previous_combination = previous
            row ^= previous_row
            combination ^= previous_combination

        indices = mask_support(combination)
        chosen = tuple(self._relations[index] for index in indices)
        return PrincipalDependency(
            relation_labels=tuple(item.label for item in chosen),
            generator_product=tuple(item.generator for item in chosen),
        )


@dataclass(frozen=True)
class SquareclassImage:
    """Images of one global squareclass in the early-quotient targets."""

    label: str
    generator: str
    local: int
    fingerprint: int


@dataclass(frozen=True)
class QuotientImage:
    """The reduced signature of a global squareclass modulo known MW images."""

    label: str
    generator: str
    raw_signature: int
    residual_signature: int

    @property
    def killed_by_known_mw(self) -> bool:
        return self.residual_signature == 0


class EarlyQuotient:
    """Project images to local/fingerprint space modulo known MW directions.

    The two target spaces are deliberately concatenated, rather than treated
    as interchangeable.  A zero residual signature only says that this
    *chosen* target fails to distinguish a candidate from known MW;
    it is not a global-square or Selmer conclusion.
    """

    def __init__(
        self,
        *,
        local_dimension: int,
        fingerprint_dimension: int,
        known_mw_images: Iterable[SquareclassImage] = (),
    ) -> None:
        if local_dimension < 0 or fingerprint_dimension < 0:
            raise F2Error("target dimensions must be nonnegative")
        self.local_dimension = local_dimension
        self.fingerprint_dimension = fingerprint_dimension
        self.dimension = local_dimension + fingerprint_dimension
        self._pivots: dict[int, int] = {}
        self.known_mw_images = tuple(known_mw_images)
        for image in self.known_mw_images:
            self._insert(self.pack(image))

    def pack(self, image: SquareclassImage) -> int:
        local = _check_mask(image.local, self.local_dimension, "local image")
        fingerprint = _check_mask(
            image.fingerprint, self.fingerprint_dimension, "fingerprint image")
        return local | (fingerprint << self.local_dimension)

    def _insert(self, row: int) -> bool:
        row = _check_mask(row, self.dimension, "quotient row")
        while row:
            pivot = row.bit_length() - 1
            if pivot in self._pivots:
                row ^= self._pivots[pivot]
            else:
                self._pivots[pivot] = row
                return True
        return False

    @property
    def known_mw_rank(self) -> int:
        return len(self._pivots)

    def reduce(self, signature: int) -> int:
        signature = _check_mask(signature, self.dimension, "quotient signature")
        # Membership testing may stop at the first free coordinate, but a
        # quotient projection must eliminate every pivot. Otherwise two
        # representatives of the same coset can acquire independent residues.
        for pivot in sorted(self._pivots, reverse=True):
            if signature >> pivot & 1:
                signature ^= self._pivots[pivot]
        return signature

    def image(self, candidate: SquareclassImage) -> QuotientImage:
        raw = self.pack(candidate)
        return QuotientImage(
            label=candidate.label,
            generator=candidate.generator,
            raw_signature=raw,
            residual_signature=self.reduce(raw),
        )


@dataclass(frozen=True)
class ClassQuotientCertification:
    """An explicitly scoped upper bound for the remaining mod-2 S-class part."""

    method: str
    remaining_dimension_upper_bound: int | None
    hypothesis: str | None = None

    def status(self) -> str:
        if self.remaining_dimension_upper_bound is None:
            return "UNCERTIFIED_RELATION_STABILIZATION"
        if self.remaining_dimension_upper_bound < 0:
            raise F2Error("remaining mod-2 class quotient bound must be nonnegative")
        return "CERTIFIED_UNDER_HYPOTHESIS" if self.hypothesis else "CERTIFIED"

    def require_valid_method(self) -> None:
        method = self.method.strip().lower()
        if self.remaining_dimension_upper_bound is not None and (
            not method or method == "none" or "stabili" in method
        ):
            raise F2Error("relation-rank stabilization is not a completeness certificate")


def certification_record(certification: ClassQuotientCertification) -> Mapping[str, object]:
    """Return a JSON-ready record while enforcing the no-stabilization rule."""

    certification.require_valid_method()
    return {
        "status": certification.status(),
        "method": certification.method,
        "remaining_mod2_s_class_dimension_upper_bound": (
            certification.remaining_dimension_upper_bound
        ),
        "hypothesis": certification.hypothesis,
        "interpretation": (
            "This only certifies the stated S-class quotient bound; a separate "
            "local-descent argument is required for a 2-Selmer or rank bound."
        ),
    }


def _integer(value: object, name: str) -> int:
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        return int(value, 0)
    raise F2Error(f"{name} must be an integer or a base-prefixed integer string")


def _image_from_json(record: Mapping[str, object]) -> SquareclassImage:
    try:
        return SquareclassImage(
            label=str(record["label"]),
            generator=str(record["generator"]),
            local=_integer(record["local"], "local"),
            fingerprint=_integer(record["fingerprint"], "fingerprint"),
        )
    except KeyError as error:
        raise F2Error(f"squareclass image is missing {error.args[0]}") from error


def audit_manifest(manifest: Mapping[str, object]) -> Mapping[str, object]:
    """Audit a JSON early-quotient manifest and return JSON-ready output.

    This is deliberately a transport layer: the field-specific local-square
    calculations and the Kummer images are inputs, whose exact generators stay
    visible in the output.
    """

    try:
        local_dimension = _integer(manifest["local_dimension"], "local_dimension")
        fingerprint_dimension = _integer(
            manifest["fingerprint_dimension"], "fingerprint_dimension"
        )
    except KeyError as error:
        raise F2Error(f"manifest is missing {error.args[0]}") from error

    known = tuple(
        _image_from_json(record)
        for record in manifest.get("known_mw_images", ())
        if isinstance(record, Mapping)
    )
    if len(known) != len(manifest.get("known_mw_images", ())):
        raise F2Error("known_mw_images must contain objects")
    quotient = EarlyQuotient(
        local_dimension=local_dimension,
        fingerprint_dimension=fingerprint_dimension,
        known_mw_images=known,
    )
    candidates = tuple(
        _image_from_json(record)
        for record in manifest.get("candidate_images", ())
        if isinstance(record, Mapping)
    )
    if len(candidates) != len(manifest.get("candidate_images", ())):
        raise F2Error("candidate_images must contain objects")
    images = tuple(quotient.image(candidate) for candidate in candidates)
    residual_rank, independent_labels = independent_masks(
        (image.label, image.residual_signature) for image in images
    )

    certification_input = manifest.get("class_quotient_certification", {})
    if not isinstance(certification_input, Mapping):
        raise F2Error("class_quotient_certification must be an object")
    certification = ClassQuotientCertification(
        method=str(certification_input.get("method", "none")),
        remaining_dimension_upper_bound=(
            None
            if certification_input.get("remaining_dimension_upper_bound") is None
            else _integer(
                certification_input["remaining_dimension_upper_bound"],
                "remaining_dimension_upper_bound",
            )
        ),
        hypothesis=(
            None
            if certification_input.get("hypothesis") is None
            else str(certification_input["hypothesis"])
        ),
    )
    return {
        "protocol": "BNFFREE2SEL-v1",
        "known_mw_target_rank": quotient.known_mw_rank,
        "candidate_residual_rank": residual_rank,
        "independent_candidate_labels": independent_labels,
        "candidate_images": [
            {
                "label": image.label,
                "generator": image.generator,
                "raw_signature": image.raw_signature,
                "residual_signature": image.residual_signature,
                "residual_support": mask_support(image.residual_signature),
                "killed_by_known_mw_in_this_target": image.killed_by_known_mw,
            }
            for image in images
        ],
        "class_quotient_certification": certification_record(certification),
        "status": (
            "BOOKKEEPING_ONLY: no 2-Selmer dimension, Cassels--Tate result, "
            "or Mordell--Weil rank follows from this audit alone."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="JSON image manifest")
    parser.add_argument("--output", type=Path, required=True, help="JSON audit output")
    args = parser.parse_args()
    manifest = json.loads(args.input.read_text())
    if not isinstance(manifest, Mapping):
        raise F2Error("top-level JSON value must be an object")
    output = audit_manifest(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"BNFFREE2SEL|status=WROTE|output={args.output}")


if __name__ == "__main__":
    main()

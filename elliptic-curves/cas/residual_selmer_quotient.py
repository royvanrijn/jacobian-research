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
    *chosen* faithful target fails to distinguish a candidate from known MW;
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
        while signature:
            pivot = signature.bit_length() - 1
            if pivot not in self._pivots:
                break
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

#!/usr/bin/env python3
"""Finite combinatorial atlas for plane finite-normalization signatures.

For one irreducible target nonproperness curve C, a canonical finite cover
of geometric degree d has

    d = sum_boundary(e_i f_i) + sum_affine(a_j),

where every affine row is unramified and the affine contribution is
positive.  A boundary residue map of degree f_i has s_i punctures with
1 <= s_i <= f_i and Riemann--Hurwitz forces f_i+s_i-2 ramification on its
affine normalization.

This module enumerates the finite set of coarse signatures at fixed d.  It
does not assert existence of a cover or a Keller map.  Its purpose is the
same as a support atlas before coefficient elimination: reject impossible
combinatorics first and expose the exact geometric certificates still needed.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from functools import lru_cache


@lru_cache(maxsize=None)
def integer_partitions(total: int, minimum: int = 1) -> tuple[tuple[int, ...], ...]:
    """Return nondecreasing positive partitions of ``total``."""

    if total < 0 or minimum <= 0:
        return ()
    if total == 0:
        return ((),)
    result: list[tuple[int, ...]] = []
    for first in range(minimum, total + 1):
        for tail in integer_partitions(total - first, first):
            result.append((first, *tail))
    return tuple(result)


@dataclass(frozen=True, order=True)
class BoundaryRow:
    """One boundary valuation row above the chosen target curve."""

    ramification_index: int
    residue_degree: int
    punctures: int

    def __post_init__(self) -> None:
        if self.ramification_index <= 0:
            raise ValueError("ramification index must be positive")
        if self.residue_degree <= 0:
            raise ValueError("residue degree must be positive")
        if not 1 <= self.punctures <= self.residue_degree:
            raise ValueError("punctures must lie between one and residue degree")

    @property
    def degree_contribution(self) -> int:
        return self.ramification_index * self.residue_degree

    @property
    def forced_affine_residue_ramification(self) -> int:
        return self.residue_degree + self.punctures - 2

    @property
    def puncture_ramification(self) -> int:
        return self.residue_degree - self.punctures

    @property
    def total_residue_ramification(self) -> int:
        return 2 * self.residue_degree - 2

    @property
    def is_transversely_ramified(self) -> bool:
        return self.ramification_index > 1

    @property
    def is_residue_immersive_compatible(self) -> bool:
        return self.forced_affine_residue_ramification == 0


@dataclass(frozen=True)
class TargetNormalizationSignature:
    """Coarse generic and puncture data over one nonproperness component."""

    geometric_degree: int
    boundary_rows: tuple[BoundaryRow, ...]
    affine_residue_degrees: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.geometric_degree <= 0:
            raise ValueError("geometric degree must be positive")
        if not self.boundary_rows:
            raise ValueError("a nonproperness component needs a boundary row")
        if not self.affine_residue_degrees:
            raise ValueError("every target curve needs an affine sheet")
        if any(degree <= 0 for degree in self.affine_residue_degrees):
            raise ValueError("affine residue degrees must be positive")
        if self.total_degree != self.geometric_degree:
            raise ValueError("boundary and affine rows must exhaust the degree")

    @property
    def boundary_degree(self) -> int:
        return sum(row.degree_contribution for row in self.boundary_rows)

    @property
    def affine_degree(self) -> int:
        return sum(self.affine_residue_degrees)

    @property
    def total_degree(self) -> int:
        return self.boundary_degree + self.affine_degree

    @property
    def has_transverse_ramification(self) -> bool:
        return any(row.is_transversely_ramified for row in self.boundary_rows)

    @property
    def residue_immersive_compatible(self) -> bool:
        return all(
            row.is_residue_immersive_compatible for row in self.boundary_rows
        )

    @property
    def total_punctures(self) -> int:
        return sum(row.punctures for row in self.boundary_rows)

    @property
    def residual_ramification_cost(self) -> int:
        return sum(
            row.forced_affine_residue_ramification
            for row in self.boundary_rows
        )

    @property
    def pareto_complexity(self) -> tuple[int, int, int, int, int]:
        """Canonical local complexity coordinates for bounded atlases."""

        return (
            self.geometric_degree,
            len(self.boundary_rows),
            len(self.affine_residue_degrees),
            self.total_punctures,
            self.residual_ramification_cost,
        )

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result.update(
            boundary_degree=self.boundary_degree,
            affine_degree=self.affine_degree,
            has_transverse_ramification=self.has_transverse_ramification,
            residue_immersive_compatible=self.residue_immersive_compatible,
            residual_ramification_cost=self.residual_ramification_cost,
            pareto_complexity=self.pareto_complexity,
        )
        for row in result["boundary_rows"]:
            residue_degree = row["residue_degree"]
            punctures = row["punctures"]
            row.update(
                puncture_ramification=residue_degree - punctures,
                forced_affine_residue_ramification=(
                    residue_degree + punctures - 2
                ),
                total_residue_ramification=2 * residue_degree - 2,
            )
        return result


def _boundary_row_types(maximum_contribution: int) -> tuple[BoundaryRow, ...]:
    rows: list[BoundaryRow] = []
    for ramification_index in range(1, maximum_contribution + 1):
        for residue_degree in range(
            1, maximum_contribution // ramification_index + 1
        ):
            for punctures in range(1, residue_degree + 1):
                rows.append(
                    BoundaryRow(
                        ramification_index,
                        residue_degree,
                        punctures,
                    )
                )
    return tuple(sorted(rows))


@lru_cache(maxsize=None)
def _boundary_multisets(
    row_types: tuple[BoundaryRow, ...],
    maximum_total: int,
    start: int = 0,
) -> tuple[tuple[BoundaryRow, ...], ...]:
    """Enumerate nonempty row multisets of bounded degree contribution."""

    result: list[tuple[BoundaryRow, ...]] = []
    for index in range(start, len(row_types)):
        row = row_types[index]
        contribution = row.degree_contribution
        if contribution > maximum_total:
            continue
        result.append((row,))
        for tail in _boundary_multisets(
            row_types,
            maximum_total - contribution,
            index,
        ):
            result.append((row, *tail))
    return tuple(result)


def enumerate_target_signatures(
    geometric_degree: int,
    *,
    require_transverse_ramification: bool = False,
    require_residue_immersion: bool = False,
) -> tuple[TargetNormalizationSignature, ...]:
    """Enumerate every coarse signature over one target component."""

    if geometric_degree <= 1:
        return ()
    row_types = _boundary_row_types(geometric_degree - 1)
    signatures: list[TargetNormalizationSignature] = []
    for boundary_rows in _boundary_multisets(
        row_types,
        geometric_degree - 1,
    ):
        boundary_degree = sum(
            row.degree_contribution for row in boundary_rows
        )
        affine_degree = geometric_degree - boundary_degree
        for affine_partition in integer_partitions(affine_degree):
            signature = TargetNormalizationSignature(
                geometric_degree,
                boundary_rows,
                affine_partition,
            )
            if (
                require_transverse_ramification
                and not signature.has_transverse_ramification
            ):
                continue
            if (
                require_residue_immersion
                and not signature.residue_immersive_compatible
            ):
                continue
            signatures.append(signature)
    return tuple(
        sorted(
            set(signatures),
            key=lambda signature: (
                signature.pareto_complexity,
                signature.boundary_rows,
                signature.affine_residue_degrees,
            ),
        )
    )


def pareto_minimal_signatures(
    signatures: tuple[TargetNormalizationSignature, ...],
) -> tuple[TargetNormalizationSignature, ...]:
    """Keep signatures undominated in the four non-degree coordinates."""

    def dominates(
        left: TargetNormalizationSignature,
        right: TargetNormalizationSignature,
    ) -> bool:
        left_complexity = left.pareto_complexity[1:]
        right_complexity = right.pareto_complexity[1:]
        return all(
            x <= y for x, y in zip(left_complexity, right_complexity)
        ) and any(
            x < y for x, y in zip(left_complexity, right_complexity)
        )

    return tuple(
        signature
        for signature in signatures
        if not any(
            dominates(other, signature)
            for other in signatures
            if other != signature
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("degree", type=int)
    parser.add_argument("--ramified", action="store_true")
    parser.add_argument("--immersive", action="store_true")
    parser.add_argument("--pareto", action="store_true")
    args = parser.parse_args()
    signatures = enumerate_target_signatures(
        args.degree,
        require_transverse_ramification=args.ramified,
        require_residue_immersion=args.immersive,
    )
    if args.pareto:
        signatures = pareto_minimal_signatures(signatures)
    print(
        json.dumps(
            [signature.as_dict() for signature in signatures],
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Parameter-independent API for skew-height searches on Nagao rank-13 curves.

The first pinned application lives in ``search_nagao_u42_skew_height.py``.
This module exposes its low-level exact search machinery through target and
checkpoint objects, so the same boxes and Mobius-chart construction can be
applied to another rational base-change parameter without copying the engine.

In particular, ``load_rank17_target`` accepts the exact rank-17 frontier
certificate, while ``classify_uniform_checkpoint`` turns a fresh uniform
``hyperellratpoints`` result into the chart centres needed for the next stage.
No search is launched merely by importing or calling the planning helpers.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Sequence

from nagao_1994 import (
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_value,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from nagao_linear_sections import omitted_companion_sections, point_on_short_curve
from search_extra_points import signless_quartic_points
from search_nagao_u42_skew_height import (
    MOBIUS_HEIGHT,
    SEARCH_BOXES,
    SearchBox,
    centered_unimodular_matrix,
    discover_relations,
    exact_linear_combination,
    map_chart_point,
    run_mobius_charts,
    run_skew_box,
    short_add,
    short_multiply,
    transform_binary_quartic,
)


Q = Fraction


@dataclass(frozen=True)
class Rank17Target:
    parameter_u: Fraction
    parameter_t: Fraction
    quartic_coefficients: tuple[Fraction, ...]
    jacobian_coefficients: tuple[Fraction, ...]
    saturated_basis: tuple[tuple[Fraction, Fraction], ...]
    certified_rank_lower_bound: int
    conductor: int
    log_conductor: str
    root_number: int
    certificate_source: str

    @property
    def identifier(self) -> str:
        numerator = self.parameter_u.numerator
        denominator = self.parameter_u.denominator
        return f"nagao-u-{numerator}-{denominator}"


@dataclass(frozen=True)
class UniformCheckpoint:
    target: Rank17Target
    height_bound: int
    raw_signed_points: tuple[tuple[Fraction, Fraction], ...]
    signless_points: tuple[tuple[Fraction, Fraction], ...]
    unexpected_points: tuple[tuple[Fraction, Fraction], ...]
    displayed_x_count: int
    companion_x_count: int
    zero_ordinate_count: int

    @property
    def chart_centers(self) -> tuple[Fraction, ...]:
        return tuple(point[0] for point in self.unexpected_points)


@dataclass(frozen=True)
class CheckpointReference:
    parameter_u: Fraction
    height_bound: int
    signed_point_count: int
    unexpected_x_count: int
    distinct_x_count: int
    unexpected_point_sha256: str
    stable_pool_numerical_rank: int | None
    source: str


def load_rank17_target(
    certificate_path: Path, parameter_u: Fraction
) -> Rank17Target:
    """Load and exactly validate one target from the frontier certificate."""

    parameter_u = Q(parameter_u)
    data = json.loads(certificate_path.read_text())
    matches = [
        record
        for record in data["certificates"]
        if Q(record["parameter_u"]) == parameter_u
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one certificate for u={parameter_u}, found {len(matches)}"
        )
    record = matches[0]
    parameter_t = rank13_base_parameter(parameter_u)
    if Q(record["parameter_t"]) != parameter_t:
        raise AssertionError("the certificate has an inconsistent base parameter")
    quartic_coefficients = primitive_quartic_coefficients(
        RANK13_CONSTRUCTION, parameter_t
    )
    jacobian_coefficients = rank13_base_changed_short_jacobian_coefficients(
        parameter_u
    )
    stored_coefficients = tuple(
        Q(value) for value in record["short_weierstrass_coefficients"]
    )
    if stored_coefficients != jacobian_coefficients:
        raise AssertionError("the certificate has inconsistent curve coefficients")
    basis = tuple(
        (Q(point["jacobian_x"]), Q(point["jacobian_y"]))
        for point in record["saturated_basis"]
    )
    if any(not point_on_short_curve(jacobian_coefficients, point) for point in basis):
        raise AssertionError("a certificate basis point is off the exact curve")
    certified_lower_bound = int(
        record["finite_reduction_certificate"]["certified_algebraic_rank_lower_bound"]
    )
    if len(basis) != certified_lower_bound:
        raise AssertionError("basis length and certified lower bound disagree")
    return Rank17Target(
        parameter_u=parameter_u,
        parameter_t=parameter_t,
        quartic_coefficients=quartic_coefficients,
        jacobian_coefficients=jacobian_coefficients,
        saturated_basis=basis,
        certified_rank_lower_bound=certified_lower_bound,
        conductor=int(record["conductor"]),
        log_conductor=str(record["log_conductor"]),
        root_number=int(record["root_number"]),
        certificate_source=certificate_path.name,
    )


def checkpoint_reference(
    rank_gain_path: Path,
    parameter_u: Fraction,
    *,
    height_bound: int = 1_000_000,
) -> CheckpointReference:
    """Read the compact yield/rank summary for an already-run uniform box."""

    parameter_u = Q(parameter_u)
    data = json.loads(rank_gain_path.read_text())
    matches: list[dict[str, Any]] = []
    for section_name in ("final_box", "escalation_box"):
        section = data.get(section_name, {})
        for record in section.get("records", []):
            if (
                Q(record["parameter_u"]) == parameter_u
                and int(record["quartic_naive_height_bound"]) == height_bound
            ):
                matches.append(record)
    if len(matches) != 1:
        raise ValueError(
            f"expected one H={height_bound} reference for u={parameter_u}, "
            f"found {len(matches)}"
        )
    record = matches[0]
    if not record["all_returned_points_checked_exactly"]:
        raise AssertionError("the reference point search was not exactly checked")
    stable_rank = record.get("stable_pool_numerical_rank")
    return CheckpointReference(
        parameter_u=parameter_u,
        height_bound=height_bound,
        signed_point_count=int(record["signed_points_found"]),
        unexpected_x_count=int(record["unexpected_nonzero_quartic_x_values"]),
        distinct_x_count=int(record["distinct_quartic_x_values"]),
        unexpected_point_sha256=str(record["unexpected_point_sha256"]),
        stable_pool_numerical_rank=(
            int(stable_rank) if stable_rank is not None else None
        ),
        source=rank_gain_path.name,
    )


def classify_uniform_checkpoint(
    target: Rank17Target,
    raw_points: Sequence[tuple[Fraction, Fraction]],
    *,
    height_bound: int,
) -> UniformCheckpoint:
    """Check and classify a raw uniform-search result exactly.

    The 13 displayed sections and the five omitted linear companion sections
    are excluded from ``unexpected_points``.  One sign from each ``z/-z`` pair
    is retained, matching the later Jacobian-image search convention.
    """

    if height_bound <= 0:
        raise ValueError("the checkpoint height must be positive")
    raw = tuple((Q(x_value), Q(z_value)) for x_value, z_value in raw_points)
    if any(
        z_value**2 != quartic_value(target.quartic_coefficients, x_value)
        for x_value, z_value in raw
    ):
        raise AssertionError("a raw checkpoint point is off the exact quartic")
    signless = signless_quartic_points(raw)
    displayed_x = {
        point[0] for point in rank13_known_quartic_points(target.parameter_t)
    }
    companion_x = {
        section.point(target.parameter_t)[0]
        for section in omitted_companion_sections()
    }
    unexpected = tuple(
        point
        for point in signless
        if point[1] != 0
        and point[0] not in displayed_x
        and point[0] not in companion_x
    )
    return UniformCheckpoint(
        target=target,
        height_bound=height_bound,
        raw_signed_points=raw,
        signless_points=signless,
        unexpected_points=unexpected,
        displayed_x_count=sum(point[0] in displayed_x for point in signless),
        companion_x_count=sum(point[0] in companion_x for point in signless),
        zero_ordinate_count=sum(point[1] == 0 for point in signless),
    )


def build_mobius_chart_plan(
    checkpoint: UniformCheckpoint,
    *,
    shifts: Sequence[int] = (0, -1),
) -> tuple[tuple[str, tuple[int, int, int, int]], ...]:
    """Build deterministic determinant-one charts around checkpoint extras."""

    if not shifts or len(set(shifts)) != len(shifts):
        raise ValueError("chart shifts must be nonempty and distinct")
    charts = []
    for index, center in enumerate(checkpoint.chart_centers):
        for shift in shifts:
            charts.append(
                (
                    f"unexpected_{index:03d}_shift_{shift}",
                    centered_unimodular_matrix(center, int(shift)),
                )
            )
    return tuple(charts)


def target_plan(
    target: Rank17Target,
    reference: CheckpointReference,
    *,
    search_boxes: Sequence[SearchBox] = SEARCH_BOXES,
    chart_height: int = MOBIUS_HEIGHT,
    chart_shifts: Sequence[int] = (0, -1),
) -> dict[str, Any]:
    """Return a no-compute plan for the next parameterized search."""

    if target.parameter_u != reference.parameter_u:
        raise ValueError("the target and checkpoint reference disagree")
    return {
        "target": target.identifier,
        "parameter_u": str(target.parameter_u),
        "parameter_t": str(target.parameter_t),
        "certified_rank_lower_bound": target.certified_rank_lower_bound,
        "log_conductor": target.log_conductor,
        "uniform_checkpoint": {
            "height_bound": reference.height_bound,
            "expected_unexpected_x_count": reference.unexpected_x_count,
            "expected_distinct_x_count": reference.distinct_x_count,
            "expected_stable_numerical_rank": reference.stable_pool_numerical_rank,
        },
        "skew_boxes": [
            {
                "id": item.identifier,
                "numerator_bound": item.numerator_bound,
                "denominator_lower": item.denominator_lower,
                "denominator_upper": item.denominator_upper,
            }
            for item in search_boxes
        ],
        "mobius": {
            "chart_height": chart_height,
            "shifts": list(chart_shifts),
            "expected_chart_count_after_checkpoint": (
                reference.unexpected_x_count * len(chart_shifts)
            ),
        },
        "launches_search": False,
    }


__all__ = [
    "CheckpointReference",
    "MOBIUS_HEIGHT",
    "Rank17Target",
    "SEARCH_BOXES",
    "SearchBox",
    "UniformCheckpoint",
    "build_mobius_chart_plan",
    "centered_unimodular_matrix",
    "checkpoint_reference",
    "classify_uniform_checkpoint",
    "discover_relations",
    "exact_linear_combination",
    "load_rank17_target",
    "map_chart_point",
    "run_mobius_charts",
    "run_skew_box",
    "short_add",
    "short_multiply",
    "target_plan",
    "transform_binary_quartic",
]

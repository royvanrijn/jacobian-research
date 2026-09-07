#!/usr/bin/env python3
"""Independent proof and resource gates for production rank searches.

The statistical laboratories may impose stronger blinding/phase barriers on
their own outcomes.  Those experiment-specific barriers are not production
exclusion theorems.  In production, only a certified rank/Selmer upper bound
below the target blocks a point search.  Missing or conditional descent data
remain scheduling information.  Every search must still have an explicit
finite resource budget, and every point-based rank claim must be certified
independently of the descent state.
"""

from __future__ import annotations

from math import isfinite
from typing import Mapping, Sequence


SCHEMA = "elliptic-curves.production-rank-search-gates.v1"
EXCLUDED_STATUS = "EXCLUDED_BY_CERTIFIED_UPPER_BOUND"
OPEN_STATUS = "NOT_EXCLUDED_BY_CERTIFIED_UPPER_BOUND"
SEARCH_OPEN_STATUS = "OPEN_BOUNDED_SEARCH_BUDGET"
SEARCH_BLOCKED_STATUS = "BLOCKED_BY_PROVED_TARGET_EXCLUSION"
LOWER_BOUND_SUCCESS_STATUS = "PASS_CERTIFIED_POINT_LOWER_BOUND_TARGET_MET"
LOWER_BOUND_PARTIAL_STATUS = "PASS_CERTIFIED_POINT_LOWER_BOUND_BELOW_TARGET"


class ProductionSearchGateError(ValueError):
    """A proposed proof, budget, or lower-bound record is malformed."""


def _positive_search_limits(limits: Mapping[str, object]) -> dict[str, int | float]:
    if not isinstance(limits, Mapping) or not limits:
        raise ProductionSearchGateError("production search requires explicit limits")
    checked: dict[str, int | float] = {}
    for name, value in limits.items():
        if (
            not isinstance(name, str)
            or not name
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or value <= 0
            or (isinstance(value, float) and not isfinite(value))
        ):
            raise ProductionSearchGateError(
                "every production search limit must be named, numeric, finite, and positive"
            )
        checked[name] = value
    return checked


def production_gate_record(
    *,
    target_rank: int,
    search_limits: Mapping[str, object],
    certified_rank_upper_bound: int | None = None,
    upper_bound_kind: str | None = None,
    upper_bound_evidence: str | None = None,
    scheduling_information: Sequence[Mapping[str, object]] = (),
) -> dict[str, object]:
    """Build separate mathematical-exclusion and search-budget decisions.

    ``certified_rank_upper_bound`` must be an unconditional theorem-level
    upper bound.  Conditional computations, incomplete descents, timeouts,
    scores, and local signatures belong only in ``scheduling_information``.
    """

    if not isinstance(target_rank, int) or target_rank < 0:
        raise ProductionSearchGateError("target rank must be a nonnegative integer")
    limits = _positive_search_limits(search_limits)
    if certified_rank_upper_bound is None:
        if upper_bound_kind is not None or upper_bound_evidence is not None:
            raise ProductionSearchGateError(
                "upper-bound provenance was supplied without a certified upper bound"
            )
    else:
        if not isinstance(certified_rank_upper_bound, int) or certified_rank_upper_bound < 0:
            raise ProductionSearchGateError(
                "certified rank upper bound must be a nonnegative integer"
            )
        if not isinstance(upper_bound_kind, str) or not upper_bound_kind:
            raise ProductionSearchGateError("certified upper bound lacks its proof kind")
        if not isinstance(upper_bound_evidence, str) or not upper_bound_evidence:
            raise ProductionSearchGateError("certified upper bound lacks evidence")

    scheduling = []
    for item in scheduling_information:
        if not isinstance(item, Mapping) or not item:
            raise ProductionSearchGateError("scheduling information is malformed")
        record = dict(item)
        record["used_as_mathematical_exclusion"] = False
        scheduling.append(record)

    excluded = (
        certified_rank_upper_bound is not None
        and certified_rank_upper_bound < target_rank
    )
    return {
        "schema": SCHEMA,
        "target_rank": target_rank,
        "proof_gate": {
            "status": EXCLUDED_STATUS if excluded else OPEN_STATUS,
            "certified_rank_upper_bound": certified_rank_upper_bound,
            "upper_bound_kind": upper_bound_kind,
            "evidence": upper_bound_evidence,
            "target_mathematically_excluded": excluded,
            "decision": (
                "target excluded by an unconditional certified upper bound"
                if excluded
                else "no unconditional certified upper bound excludes the target"
            ),
        },
        "search_budget_gate": {
            "status": SEARCH_BLOCKED_STATUS if excluded else SEARCH_OPEN_STATUS,
            "bounded_search_authorized": not excluded,
            "limits": limits,
            "descent_completion_required": False,
            "rank_or_selmer_claims_from_search_misses": False,
        },
        "scheduling_information": scheduling,
        "claim_boundary": (
            "Incomplete or conditional arithmetic may prioritize work but cannot "
            "exclude the target. A bounded miss is not a rank upper bound."
        ),
    }


def certified_point_lower_bound_record(
    *,
    certified_independent_rank: int,
    target_rank: int,
    curve_equations_verified: bool,
    independence_evidence: str,
) -> dict[str, object]:
    """Record an unconditional point-based rank lower bound.

    This certificate is deliberately independent of whether a descent is
    complete.  It does not claim exact rank or saturation.
    """

    if (
        not isinstance(certified_independent_rank, int)
        or certified_independent_rank < 0
        or not isinstance(target_rank, int)
        or target_rank < 0
    ):
        raise ProductionSearchGateError("point ranks must be nonnegative integers")
    if curve_equations_verified is not True:
        raise ProductionSearchGateError("point coordinates lack exact equation checks")
    if not isinstance(independence_evidence, str) or not independence_evidence:
        raise ProductionSearchGateError("point lower bound lacks independence evidence")
    success = certified_independent_rank >= target_rank
    return {
        "schema": SCHEMA,
        "status": (
            LOWER_BOUND_SUCCESS_STATUS if success else LOWER_BOUND_PARTIAL_STATUS
        ),
        "certified_rank_lower_bound": certified_independent_rank,
        "target_rank": target_rank,
        "target_lower_bound_met": success,
        "curve_equations_verified": True,
        "independence_evidence": independence_evidence,
        "descent_completion_required": False,
        "exact_rank_claimed": False,
        "saturation_claimed": False,
    }


__all__ = [
    "EXCLUDED_STATUS",
    "LOWER_BOUND_PARTIAL_STATUS",
    "LOWER_BOUND_SUCCESS_STATUS",
    "OPEN_STATUS",
    "ProductionSearchGateError",
    "SCHEMA",
    "SEARCH_BLOCKED_STATUS",
    "SEARCH_OPEN_STATUS",
    "certified_point_lower_bound_record",
    "production_gate_record",
]

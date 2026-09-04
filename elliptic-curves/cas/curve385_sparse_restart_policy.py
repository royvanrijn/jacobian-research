#!/usr/bin/env python3
"""Independent restart budgets for the curve-385 sparse rank-32 search."""

from __future__ import annotations

from typing import Any, Iterable


POLICY_ID = "curve385-sparse-independent-restart-budgets-v2"
RANK_CHANGING = "RANK_CHANGING"
SATURATION_ONLY = "SATURATION_ONLY"
DEFAULT_MAXIMUM_RANK_CHANGING_RESTARTS = 3
DEFAULT_MAXIMUM_SATURATION_ONLY_RESTARTS = 4

ACCOUNTING_FIELDS = (
    "rank_changing_group_change_count",
    "saturation_only_group_change_count",
)


def empty_accounting() -> dict[str, int]:
    return {field: 0 for field in ACCOUNTING_FIELDS}


def validate_accounting(accounting: dict[str, int]) -> None:
    """Reject legacy, incomplete, or malformed checkpoint counters."""

    if set(accounting) != set(ACCOUNTING_FIELDS):
        raise ValueError("restart accounting must contain exactly the v2 counters")
    if any(
        type(accounting[field]) is not int or accounting[field] < 0
        for field in ACCOUNTING_FIELDS
    ):
        raise ValueError("restart accounting counters must be nonnegative integers")


def classify_group_change(
    *,
    rank_before: int,
    rank_after: int,
    basis_before_sha256: str,
    basis_after_sha256: str,
    finite_index_saturation_event_count: int,
) -> str | None:
    """Classify a discovered-group transition into one budget category."""

    if rank_after < rank_before:
        raise ValueError("a discovered-group update cannot lower rank")
    basis_changed = basis_before_sha256 != basis_after_sha256
    if rank_after > rank_before:
        if not basis_changed:
            raise ValueError("a rank-changing update must change the basis")
        return RANK_CHANGING
    if basis_changed:
        if finite_index_saturation_event_count < 1:
            raise ValueError(
                "a same-rank basis change needs an exact finite-index saturation event"
            )
        return SATURATION_ONLY
    if finite_index_saturation_event_count:
        raise ValueError("a recorded saturation event did not change the basis")
    return None


def charge(
    accounting: dict[str, int],
    restart_kind: str,
) -> dict[str, int]:
    """Return updated accounting without mutating the checkpoint object."""

    validate_accounting(accounting)
    updated = dict(accounting)
    if restart_kind == RANK_CHANGING:
        field = "rank_changing_group_change_count"
    elif restart_kind == SATURATION_ONLY:
        field = "saturation_only_group_change_count"
    else:
        raise ValueError(f"unknown restart kind: {restart_kind!r}")
    updated[field] += 1
    return updated


def exceeded_budget(
    accounting: dict[str, int],
    *,
    maximum_rank_changing_restarts: int,
    maximum_saturation_only_restarts: int,
) -> str | None:
    """Return the exceeded category, if any."""

    validate_accounting(accounting)
    if (
        maximum_rank_changing_restarts < 0
        or maximum_saturation_only_restarts < 0
    ):
        raise ValueError("restart budgets must be nonnegative")
    if (
        accounting["rank_changing_group_change_count"]
        > maximum_rank_changing_restarts
    ):
        return RANK_CHANGING
    if (
        accounting["saturation_only_group_change_count"]
        > maximum_saturation_only_restarts
    ):
        return SATURATION_ONLY
    return None


def simulate_unit_rank_path(
    restart_kinds: Iterable[str],
    *,
    starting_rank: int = 29,
    target_rank: int = 32,
    maximum_rank_changing_restarts: int = DEFAULT_MAXIMUM_RANK_CHANGING_RESTARTS,
    maximum_saturation_only_restarts: int = DEFAULT_MAXIMUM_SATURATION_ONLY_RESTARTS,
) -> dict[str, Any]:
    """Replay budget accounting for a deterministic regression path.

    Each ``RANK_CHANGING`` event raises the rank by one.  Reaching the target
    succeeds before any further state is requested, exactly as in the search
    runner.
    """

    rank = int(starting_rank)
    accounting = empty_accounting()
    searched_state_count = 0
    for restart_kind in restart_kinds:
        searched_state_count += 1
        accounting = charge(accounting, restart_kind)
        if restart_kind == RANK_CHANGING:
            rank += 1
        if rank >= target_rank:
            return {
                "status": "TARGET_REACHED",
                "rank": rank,
                "searched_state_count": searched_state_count,
                "restart_accounting": accounting,
            }
        exceeded = exceeded_budget(
            accounting,
            maximum_rank_changing_restarts=maximum_rank_changing_restarts,
            maximum_saturation_only_restarts=maximum_saturation_only_restarts,
        )
        if exceeded is not None:
            return {
                "status": "RESTART_BUDGET_EXCEEDED",
                "exceeded": exceeded,
                "rank": rank,
                "searched_state_count": searched_state_count,
                "restart_accounting": accounting,
            }
    return {
        "status": "EVENT_SEQUENCE_ENDED",
        "rank": rank,
        "searched_state_count": searched_state_count,
        "restart_accounting": accounting,
    }


def policy_document() -> dict[str, Any]:
    return {
        "policy_id": POLICY_ID,
        "budgets": {
            "maximum_rank_changing_restarts": (
                DEFAULT_MAXIMUM_RANK_CHANGING_RESTARTS
            ),
            "maximum_saturation_only_restarts": (
                DEFAULT_MAXIMUM_SATURATION_ONLY_RESTARTS
            ),
            "combined_lattice_state_cap": None,
        },
        "classification": {
            RANK_CHANGING: (
                "rank_after exceeds rank_before; any simultaneous saturation "
                "events are recorded but do not also consume saturation-only budget"
            ),
            SATURATION_ONLY: (
                "rank is unchanged, the basis hash changes, and at least one exact "
                "finite-index saturation event is recorded"
            ),
            "NO_CHANGE": "rank and basis hash are unchanged; no restart is charged",
        },
        "target_precedence": (
            "a certified rank at least 32 succeeds immediately after classification, "
            "before another restart or any budget-stop decision"
        ),
        "safety_boundary": (
            "exhausting either category stops only the search campaign; it supplies "
            "no rank upper bound or saturation theorem"
        ),
    }

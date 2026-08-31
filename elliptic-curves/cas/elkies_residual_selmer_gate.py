#!/usr/bin/env python3
"""Exact policy gate between Elkies scoring and expensive point searches.

A Kummer signature, norm-one element, class-group envelope, or incomplete
descent is deliberately insufficient.  Authorization requires an
unconditional, completed 2-descent whose class-group completeness and every
local-solubility condition are part of the backend result.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence


SCHEMA = "elliptic-curves.elkies-2026-residual-2-selmer-gate.v1"
PASS_STATUS = "PASS_RANK32_RESIDUAL_2_SELMER_GATE"
REJECT_STATUS = "REJECT_RANK32_BY_RESIDUAL_2_SELMER"
INCOMPLETE_STATUS = "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN"


class ResidualSelmerGateError(ValueError):
    """A purported descent artifact cannot authorize expensive search."""


def gate_record(
    *,
    total_two_selmer_dimension: int,
    known_generic_rank: int = 17,
    target_rank: int = 32,
    two_torsion_dimension: int = 0,
) -> dict[str, object]:
    """Return the exact residual-dimension rejection decision.

    ``total_two_selmer_dimension`` includes rational 2-torsion, as usual for
    ``Sel_2(E/Q)``.  Quotienting by the known Mordell--Weil image removes only
    the ``known_generic_rank`` directions; rational 2-torsion is reported and
    retained unless it is zero.
    """

    values = (
        total_two_selmer_dimension,
        known_generic_rank,
        target_rank,
        two_torsion_dimension,
    )
    if any(not isinstance(value, int) or value < 0 for value in values):
        raise ResidualSelmerGateError("Selmer and rank dimensions must be nonnegative integers")
    if target_rank < known_generic_rank:
        raise ResidualSelmerGateError("target rank cannot be below the known generic rank")
    residual_dimension = total_two_selmer_dimension - known_generic_rank
    if residual_dimension < 0:
        raise ResidualSelmerGateError("the Selmer dimension contradicts the certified subgroup")
    required = target_rank - known_generic_rank
    survives = residual_dimension >= required
    return {
        "known_generic_rank": known_generic_rank,
        "target_rank": target_rank,
        "two_torsion_dimension": two_torsion_dimension,
        "total_two_selmer_dimension": total_two_selmer_dimension,
        "residual_two_selmer_quotient_dimension": residual_dimension,
        "required_residual_dimension": required,
        "status": PASS_STATUS if survives else REJECT_STATUS,
        "expensive_search_authorized": survives,
        "decision": (
            "survives exact residual 2-Selmer gate"
            if survives
            else "rank 32 is excluded by the exact residual 2-Selmer upper bound"
        ),
    }


def require_expensive_search_gate(
    path: Path,
    *,
    expected_parameter: str | None = None,
    expected_model: Sequence[object] | None = None,
) -> Mapping[str, object]:
    """Validate a completed unconditional gate artifact or raise."""

    document = json.loads(path.read_text())
    if document.get("schema") != SCHEMA:
        raise ResidualSelmerGateError("unknown residual 2-Selmer gate schema")
    if document.get("status") != PASS_STATUS:
        raise ResidualSelmerGateError(
            "expensive search is forbidden without a passing residual 2-Selmer gate"
        )
    backend = document.get("descent_backend")
    if not isinstance(backend, dict):
        raise ResidualSelmerGateError("missing descent-backend attestation")
    required_flags = (
        "unconditional",
        "class_group_completeness_completed",
        "all_local_solubility_conditions_completed",
    )
    if any(backend.get(flag) is not True for flag in required_flags):
        raise ResidualSelmerGateError(
            "the descent backend does not attest unconditional global/local completeness"
        )
    gate = document.get("gate")
    if not isinstance(gate, dict) or gate.get("expensive_search_authorized") is not True:
        raise ResidualSelmerGateError("the residual-dimension gate did not authorize search")
    if int(gate.get("residual_two_selmer_quotient_dimension", -1)) < int(
        gate.get("required_residual_dimension", 10**9)
    ):
        raise ResidualSelmerGateError("the stored residual dimension is below its threshold")
    if expected_parameter is not None and document.get("parameter") != expected_parameter:
        raise ResidualSelmerGateError("the descent gate belongs to a different fibre")
    if expected_model is not None and document.get("global_minimal_model") != [
        str(value) for value in expected_model
    ]:
        raise ResidualSelmerGateError("the descent gate belongs to a different minimal curve")
    return document


__all__ = [
    "INCOMPLETE_STATUS",
    "PASS_STATUS",
    "REJECT_STATUS",
    "ResidualSelmerGateError",
    "SCHEMA",
    "gate_record",
    "require_expensive_search_gate",
]

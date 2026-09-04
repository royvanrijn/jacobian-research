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
from typing import Any, Mapping, Sequence


SCHEMA = "elliptic-curves.elkies-2026-residual-2-selmer-gate.v1"
PASS_STATUS = "PASS_RANK32_RESIDUAL_2_SELMER_GATE"
REJECT_STATUS = "REJECT_RANK32_BY_RESIDUAL_2_SELMER"
INCOMPLETE_STATUS = "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN"


class ResidualSelmerGateError(ValueError):
    """A purported descent artifact cannot authorize expensive search."""


def pari_ellrank_total_two_selmer_dimension(
    *,
    rank_lower: int,
    rank_upper: int,
    cassels_pairing_rank: int,
    two_torsion_dimension: int,
) -> int:
    """Recover ``dim Sel_2(E/Q)`` from PARI ``ellrank`` output.

    For PARI's four-component result ``[r1, r2, s, L]``, ``s`` is the
    (even) rank of ``Sha(E)[2] / 2 Sha(E)[4]`` detected by the Cassels
    pairing; it is not the full dimension of ``Sha(E)[2]``.  PARI defines
    ``r2 = C - T - s``, where ``C`` is the 2-Selmer dimension and ``T`` is
    the rational 2-torsion dimension.  Hence ``C = r2 + T + s``.
    """

    values = (
        rank_lower,
        rank_upper,
        cassels_pairing_rank,
        two_torsion_dimension,
    )
    if any(not isinstance(value, int) or value < 0 for value in values):
        raise ResidualSelmerGateError(
            "PARI rank, Cassels-pairing, and torsion dimensions must be "
            "nonnegative integers"
        )
    if rank_lower > rank_upper:
        raise ResidualSelmerGateError("PARI ellrank returned a reversed rank interval")
    if cassels_pairing_rank % 2:
        raise ResidualSelmerGateError(
            "PARI ellrank returned an odd Cassels-pairing quotient rank"
        )
    return rank_upper + two_torsion_dimension + cassels_pairing_rank


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


def require_gate_for_specialization(
    path: Path, specialization: Mapping[str, Any]
) -> Mapping[str, object]:
    """Bind a passing gate to an exact q12o5867 specialization artifact."""

    parameter = specialization.get("parameter")
    minimal = specialization.get("global_minimal_specialization")
    if not isinstance(parameter, dict) or not isinstance(minimal, dict):
        raise ResidualSelmerGateError("specialization lacks exact parameter/model data")
    affine_value = parameter.get("affine_value")
    model = minimal.get("model")
    if not isinstance(affine_value, str) or not isinstance(model, list):
        raise ResidualSelmerGateError("specialization has malformed parameter/model data")
    return require_expensive_search_gate(
        path,
        expected_parameter=affine_value,
        expected_model=model,
    )


__all__ = [
    "INCOMPLETE_STATUS",
    "PASS_STATUS",
    "REJECT_STATUS",
    "ResidualSelmerGateError",
    "SCHEMA",
    "gate_record",
    "require_gate_for_specialization",
    "require_expensive_search_gate",
]

#!/usr/bin/env python3
"""Residual-Selmer policy gates for theorem claims and bounded point searches.

An exact rank/Selmer claim still requires a complete unconditional descent.
Bounded point search uses a weaker monotone policy: every proved residual
upper bound is accumulated, a candidate is rejected once that bound is below
the target, and missing class-group/BNF data remain ``no finite bound yet``.
The latter authorizes only an explicitly bounded search, never a Selmer or
rank assertion.  When an exact root-number/2-parity certificate is supplied,
upper bounds are rounded down and target thresholds up to the permitted
Selmer parity.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "elliptic-curves.elkies-2026-residual-2-selmer-gate.v2"
LEGACY_SCHEMA = "elliptic-curves.elkies-2026-residual-2-selmer-gate.v1"
PASS_STATUS = "PASS_RANK32_RESIDUAL_2_SELMER_GATE"
REJECT_STATUS = "REJECT_RANK32_BY_RESIDUAL_2_SELMER"
INCOMPLETE_STATUS = "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN"
OPEN_STATUS = "OPEN_MONOTONE_RESIDUAL_SELMER_SIEVE"


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


def _validate_monotone_stages(stages: object) -> int | None:
    if not isinstance(stages, list) or not stages:
        raise ResidualSelmerGateError("the monotone sieve has no evidence stages")
    previous: int | None = None
    final: int | None = None
    for index, stage in enumerate(stages):
        if not isinstance(stage, dict):
            raise ResidualSelmerGateError("a monotone sieve stage is malformed")
        bound = stage.get("residual_upper_bound")
        if bound is None:
            if stage.get("proof_status") != "NO_FINITE_UPPER_BOUND_YET":
                raise ResidualSelmerGateError("an infinite sieve stage has the wrong status")
            if previous is not None:
                raise ResidualSelmerGateError("a finite upper bound cannot revert to infinity")
            continue
        if not isinstance(bound, int) or bound < 0:
            raise ResidualSelmerGateError("a residual upper bound must be nonnegative")
        if stage.get("proof_status") != "PROVED_UPPER_BOUND":
            raise ResidualSelmerGateError("a finite sieve bound lacks proof status")
        if not stage.get("evidence"):
            raise ResidualSelmerGateError("a finite sieve bound lacks evidence provenance")
        if previous is not None and bound > previous:
            raise ResidualSelmerGateError("residual upper bounds are not monotone")
        previous = bound
        final = bound
    return final


def parity_adjusted_residual_thresholds(
    *,
    known_generic_rank: int,
    target_rank: int,
    residual_upper_bound: int | None,
    total_two_selmer_parity: int,
) -> tuple[int, int | None, int]:
    """Round residual thresholds using ``dim Sel_2(E/Q) mod 2``."""

    if total_two_selmer_parity not in (0, 1):
        raise ResidualSelmerGateError("the total 2-Selmer parity must be zero or one")
    residual_parity = (total_two_selmer_parity - known_generic_rank) & 1
    raw_required = target_rank - known_generic_rank
    required = raw_required + ((residual_parity - raw_required) & 1)
    effective_upper = residual_upper_bound
    if effective_upper is not None:
        effective_upper -= (effective_upper - residual_parity) & 1
    return required, effective_upper, residual_parity


def monotone_sieve_gate_record(
    *,
    stages: list[Mapping[str, object]],
    search_limits: Mapping[str, object],
    known_generic_rank: int = 17,
    target_rank: int = 32,
    total_two_selmer_parity: int | None = None,
    parity_evidence: str | None = None,
) -> dict[str, object]:
    """Authorize a bounded search unless proved upper bounds reject the fibre."""

    if target_rank < known_generic_rank:
        raise ResidualSelmerGateError("target rank cannot be below the known generic rank")
    if not isinstance(search_limits, Mapping) or not search_limits:
        raise ResidualSelmerGateError("bounded search authorization needs explicit limits")
    final_bound = _validate_monotone_stages(stages)
    raw_required = target_rank - known_generic_rank
    required = raw_required
    effective_bound = final_bound
    parity_record = None
    if total_two_selmer_parity is not None:
        if not isinstance(parity_evidence, str) or not parity_evidence:
            raise ResidualSelmerGateError(
                "a Selmer-parity refinement requires evidence provenance"
            )
        required, effective_bound, residual_parity = (
            parity_adjusted_residual_thresholds(
                known_generic_rank=known_generic_rank,
                target_rank=target_rank,
                residual_upper_bound=final_bound,
                total_two_selmer_parity=total_two_selmer_parity,
            )
        )
        parity_record = {
            "total_two_selmer_dimension_mod_2": total_two_selmer_parity,
            "residual_two_selmer_dimension_mod_2": residual_parity,
            "evidence": parity_evidence,
            "method": (
                "proved 2-Selmer parity over Q plus exact root number and "
                "rational 2-torsion audit"
            ),
        }
    elif parity_evidence is not None:
        raise ResidualSelmerGateError("parity evidence was supplied without a parity")
    rejected = effective_bound is not None and effective_bound < required
    return {
        "known_generic_rank": known_generic_rank,
        "target_rank": target_rank,
        "raw_required_residual_dimension": raw_required,
        "required_residual_dimension": required,
        "proved_residual_upper_bound": final_bound,
        "parity_adjusted_proved_residual_upper_bound": effective_bound,
        "selmer_parity": parity_record,
        "status": REJECT_STATUS if rejected else OPEN_STATUS,
        "expensive_search_authorized": False,
        "bounded_point_search_authorized": not rejected,
        "theorem_claim_authorized": False,
        "decision": (
            "rank target excluded by a proved monotone residual upper bound"
            if rejected
            else "not rejected; bounded point search may run within the recorded limits"
        ),
        "sieve": {
            "order": "nonincreasing proved residual upper bounds; null means infinity",
            "stages": [dict(stage) for stage in stages],
        },
        "search_authorization": {
            "kind": "bounded_point_search",
            "bounded": True,
            "limits": dict(search_limits),
            "rank_or_selmer_claims_from_search": False,
        },
    }


def require_expensive_search_gate(
    path: Path,
    *,
    expected_parameter: str | None = None,
    expected_model: Sequence[object] | None = None,
    requested_search_limits: Mapping[str, object] | None = None,
) -> Mapping[str, object]:
    """Validate an exact gate or an explicitly bounded monotone-sieve gate."""

    document = json.loads(path.read_text())
    if document.get("schema") not in (SCHEMA, LEGACY_SCHEMA):
        raise ResidualSelmerGateError("unknown residual 2-Selmer gate schema")
    status = document.get("status")
    if status not in (PASS_STATUS, OPEN_STATUS):
        raise ResidualSelmerGateError(
            "expensive search is forbidden by the residual 2-Selmer gate"
        )
    gate = document.get("gate")
    if not isinstance(gate, dict):
        raise ResidualSelmerGateError("missing residual-dimension gate")
    if status == PASS_STATUS:
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
        if gate.get("expensive_search_authorized") is not True:
            raise ResidualSelmerGateError("the exact residual gate did not authorize search")
        if int(gate.get("residual_two_selmer_quotient_dimension", -1)) < int(
            gate.get("required_residual_dimension", 10**9)
        ):
            raise ResidualSelmerGateError("the stored residual dimension is below its threshold")
    else:
        if gate.get("bounded_point_search_authorized") is not True:
            raise ResidualSelmerGateError("the monotone sieve did not authorize bounded search")
        if gate.get("theorem_claim_authorized") is not False:
            raise ResidualSelmerGateError("an incomplete sieve cannot authorize theorem claims")
        sieve = gate.get("sieve")
        if not isinstance(sieve, dict):
            raise ResidualSelmerGateError("missing monotone sieve record")
        final_bound = _validate_monotone_stages(sieve.get("stages"))
        if final_bound != gate.get("proved_residual_upper_bound"):
            raise ResidualSelmerGateError("the final monotone bound is inconsistent")
        parity_record = gate.get("selmer_parity")
        effective_bound = final_bound
        if parity_record is not None:
            if (
                not isinstance(parity_record, dict)
                or parity_record.get("total_two_selmer_dimension_mod_2") not in (0, 1)
                or not parity_record.get("evidence")
            ):
                raise ResidualSelmerGateError("the Selmer-parity record is malformed")
            known_rank = gate.get("known_generic_rank")
            target_rank = gate.get("target_rank")
            if not isinstance(known_rank, int) or not isinstance(target_rank, int):
                raise ResidualSelmerGateError("the parity-bound ranks are malformed")
            recomputed_required, effective_bound, residual_parity = (
                parity_adjusted_residual_thresholds(
                    known_generic_rank=known_rank,
                    target_rank=target_rank,
                    residual_upper_bound=final_bound,
                    total_two_selmer_parity=int(
                        parity_record["total_two_selmer_dimension_mod_2"]
                    ),
                )
            )
            if (
                recomputed_required != gate.get("required_residual_dimension")
                or residual_parity
                != parity_record.get("residual_two_selmer_dimension_mod_2")
                or effective_bound
                != gate.get("parity_adjusted_proved_residual_upper_bound")
            ):
                raise ResidualSelmerGateError(
                    "the Selmer-parity refinement is inconsistent"
                )
        required = gate.get("required_residual_dimension")
        if not isinstance(required, int) or required < 0:
            raise ResidualSelmerGateError("the residual target is malformed")
        if effective_bound is not None and effective_bound < required:
            raise ResidualSelmerGateError("a rejected residual upper bound authorized search")
        authorization = gate.get("search_authorization")
        if (
            not isinstance(authorization, dict)
            or authorization.get("kind") != "bounded_point_search"
            or authorization.get("bounded") is not True
            or not isinstance(authorization.get("limits"), dict)
            or not authorization["limits"]
            or authorization.get("rank_or_selmer_claims_from_search") is not False
        ):
            raise ResidualSelmerGateError("bounded search limits are missing or malformed")
        declared_limits = authorization["limits"]
        if not isinstance(requested_search_limits, Mapping) or not requested_search_limits:
            raise ResidualSelmerGateError(
                "an open sieve requires the entrypoint's requested search limits"
            )
        for name, requested in requested_search_limits.items():
            declared = declared_limits.get(name)
            if (
                isinstance(requested, bool)
                or isinstance(declared, bool)
                or not isinstance(requested, (int, float))
                or not isinstance(declared, (int, float))
                or requested <= 0
                or declared <= 0
            ):
                raise ResidualSelmerGateError(
                    f"bounded search limit {name!r} is missing or nonnumeric"
                )
            if requested > declared:
                raise ResidualSelmerGateError(
                    f"requested search limit {name!r} exceeds its authorization"
                )
    if expected_parameter is not None and document.get("parameter") != expected_parameter:
        raise ResidualSelmerGateError("the descent gate belongs to a different fibre")
    if expected_model is not None and document.get("global_minimal_model") != [
        str(value) for value in expected_model
    ]:
        raise ResidualSelmerGateError("the descent gate belongs to a different minimal curve")
    return document


def require_gate_for_specialization(
    path: Path,
    specialization: Mapping[str, Any],
    *,
    requested_search_limits: Mapping[str, object] | None = None,
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
        requested_search_limits=requested_search_limits,
    )


__all__ = [
    "INCOMPLETE_STATUS",
    "LEGACY_SCHEMA",
    "OPEN_STATUS",
    "PASS_STATUS",
    "REJECT_STATUS",
    "ResidualSelmerGateError",
    "SCHEMA",
    "gate_record",
    "monotone_sieve_gate_record",
    "parity_adjusted_residual_thresholds",
    "require_gate_for_specialization",
    "require_expensive_search_gate",
]

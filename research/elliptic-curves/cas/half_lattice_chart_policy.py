#!/usr/bin/env python3
"""State-bound policy for half-lattice chart-order heuristics.

The quartics used by the half-lattice searches are pointed birational models
of the same elliptic curve.  This module deliberately gives their ordering
scores no covering, Selmer, rank, or negative-search semantics.  It also makes
an ordering certificate valid only for the exact lattice presentation from
which the chart representatives and scores were recomputed.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Iterable, Sequence


POLICY_ID = "half-lattice-birational-chart-order-v1"

HEURISTIC_MEANINGS = {
    "legacy_half_lattice_depth": {
        "permitted": (
            "a fixed-presentation CVP-derived score for ordering birational "
            "point-search charts"
        ),
        "forbidden": [
            "arithmetic depth or filtration",
            "2-covering or torsor class",
            "Selmer-class dimension or structure",
            "rank or rank-jump evidence without exact returned points",
        ],
    },
    "old_deep_43": {
        "permitted": (
            "a historically calibrated list of 43 chart labels for one exact "
            "generic basis and height form"
        ),
        "forbidden": [
            "a basis-invariant distinguished arithmetic subset",
            "a complete set of locally soluble coverings",
            "a reason to exclude unsearched chart labels",
        ],
    },
    "quotient_hamming_weight": {
        "permitted": (
            "an enumeration order in one explicitly recorded quotient-coordinate basis"
        ),
        "forbidden": [
            "an invariant of the quotient group",
            "Selmer filtration or residual Selmer dimension",
            "evidence that higher-weight charts contain no points",
        ],
    },
}

MISS_INFERENCES = {
    "rational_point_absence_on_the_elliptic_curve": False,
    "absence_in_unsearched_or_differently_presented_charts": False,
    "rank_upper_bound": False,
    "mordell_weil_saturation": False,
    "selmer_group_or_selmer_class_structure": False,
    "local_insolubility_of_a_nontrivial_2_covering": False,
}

INVALIDATING_CHANGES = [
    "adjoining an independent point",
    "finite-index enlargement or saturation of the discovered subgroup",
    "adding, deleting, reordering, negating, or otherwise changing a basis generator",
    "changing the height Gram used for representative selection",
    "changing the embedded generic-coordinate rows",
    "changing the quotient-coordinate basis or complement",
    "changing the chart universe or its coordinate labels",
]


def canonical_hash(value: Any) -> str:
    """Return the SHA-256 of a canonical JSON representation."""

    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()
    return sha256(encoded).hexdigest()


def _rows(rows: Iterable[Iterable[Any]]) -> list[list[str]]:
    """Normalize exact/numerical matrix entries without rounding them again."""

    return [[str(value) for value in row] for row in rows]


def lattice_state(
    *,
    basis_records: Sequence[Any],
    height_gram_rows: Iterable[Iterable[Any]],
    generic_coordinate_rows: Iterable[Iterable[Any]],
    quotient_coordinate_rows: Iterable[Iterable[Any]],
    chart_universe_id: str,
) -> dict[str, Any]:
    """Build the complete presentation state to which an order may be bound."""

    if not basis_records:
        raise ValueError("a chart order requires a nonempty lattice basis")
    if not chart_universe_id:
        raise ValueError("a chart order requires an explicit chart-universe id")
    return {
        "basis_records": list(basis_records),
        "height_gram_rows": _rows(height_gram_rows),
        "generic_coordinate_rows": _rows(generic_coordinate_rows),
        "quotient_coordinate_rows": _rows(quotient_coordinate_rows),
        "chart_universe_id": chart_universe_id,
    }


def bind_ordering(
    *,
    basis_records: Sequence[Any],
    height_gram_rows: Iterable[Iterable[Any]],
    generic_coordinate_rows: Iterable[Iterable[Any]],
    quotient_coordinate_rows: Iterable[Iterable[Any]],
    chart_universe_id: str,
    ordered_chart_ids: Sequence[Any],
    heuristics: Sequence[str],
) -> dict[str, Any]:
    """Certify that an ordering was recomputed for one exact presentation.

    This certificate says nothing about the arithmetic content or expected
    success of any chart.  It only prevents a cached order from being reused
    after the basis, height form, quotient coordinates, or chart universe has
    changed.
    """

    unknown = sorted(set(heuristics) - set(HEURISTIC_MEANINGS))
    if unknown:
        raise ValueError(f"unknown half-lattice chart heuristic(s): {unknown}")
    if not heuristics:
        raise ValueError("at least one chart-order heuristic must be named")
    if len(set(map(str, ordered_chart_ids))) != len(ordered_chart_ids):
        raise ValueError("ordered chart ids must be unique")
    state = lattice_state(
        basis_records=basis_records,
        height_gram_rows=height_gram_rows,
        generic_coordinate_rows=generic_coordinate_rows,
        quotient_coordinate_rows=quotient_coordinate_rows,
        chart_universe_id=chart_universe_id,
    )
    return {
        "policy_id": POLICY_ID,
        "semantics": "birational_point_search_chart_order_only",
        "heuristics": list(heuristics),
        "lattice_state_sha256": canonical_hash(state),
        "ordered_chart_count": len(ordered_chart_ids),
        "ordered_chart_ids_sha256": canonical_hash(list(ordered_chart_ids)),
        "state_bound_recomputation_complete": True,
        "valid_only_while_lattice_state_sha256_matches": True,
        "calibration_transfers_to_a_changed_state": False,
    }


def validate_ordering(
    certificate: dict[str, Any],
    *,
    basis_records: Sequence[Any],
    height_gram_rows: Iterable[Iterable[Any]],
    generic_coordinate_rows: Iterable[Iterable[Any]],
    quotient_coordinate_rows: Iterable[Iterable[Any]],
    chart_universe_id: str,
    ordered_chart_ids: Sequence[Any],
) -> None:
    """Fail closed unless an ordering still matches its complete state."""

    if certificate.get("policy_id") != POLICY_ID:
        raise ValueError("chart ordering uses the wrong policy")
    state = lattice_state(
        basis_records=basis_records,
        height_gram_rows=height_gram_rows,
        generic_coordinate_rows=generic_coordinate_rows,
        quotient_coordinate_rows=quotient_coordinate_rows,
        chart_universe_id=chart_universe_id,
    )
    if certificate.get("lattice_state_sha256") != canonical_hash(state):
        raise ValueError(
            "stale chart ordering: recompute representatives, scores, and order "
            "for the changed lattice presentation"
        )
    if certificate.get("ordered_chart_count") != len(ordered_chart_ids):
        raise ValueError("stale chart ordering: chart count changed")
    if certificate.get("ordered_chart_ids_sha256") != canonical_hash(
        list(ordered_chart_ids)
    ):
        raise ValueError("stale chart ordering: ordered chart identities changed")


def policy_document() -> dict[str, Any]:
    """Return the machine-readable semantic and invalidation contract."""

    return {
        "policy_id": POLICY_ID,
        "quartic_role": "pointed_birational_search_chart_of_the_same_elliptic_curve",
        "quartic_is_nontrivial_2_covering": False,
        "quartic_represents_a_selmer_class": False,
        "heuristic_meanings": HEURISTIC_MEANINGS,
        "miss_inferences": MISS_INFERENCES,
        "ordering_invalidation": {
            "invalidating_changes": INVALIDATING_CHANGES,
            "required_action": (
                "discard the cached ordering; recompute chart identities, shortest "
                "representatives, scores, and the complete intended order; bind the "
                "result to the new lattice-state fingerprint before another search"
            ),
            "empirical_effectiveness_transfer": (
                "none: any efficiency or enrichment claim must be revalidated on "
                "blinded controls for the changed presentation"
            ),
        },
        "positive_evidence_rule": (
            "only an exactly mapped rational point followed by the separately required "
            "group-law and independence checks can change a rank lower bound"
        ),
    }


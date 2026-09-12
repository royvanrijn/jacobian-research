#!/usr/bin/env python3
"""Continue the archived ``(72,108)`` Case-1 descent to the full polygons.

The pinned exact replay stops after bracket layer ``-3`` because its thirteen
compatibility equations already give the contradiction.  At that point the
reconstructed Laurent bands stop at ``P_-5`` and ``Q_-4``.  The published
Case-1 polygons continue through ``P_-8`` and ``Q_-12``.

This script reads (and never mutates) the pinned checkpoint, reuses its exact
number-field arithmetic, and continues the same triangular linear recurrence
through a requested terminal layer.  It prints a compact exact summary and
hashes of every reconstructed band and compatibility block.  It is a band
derivation, not a claim that the already-empty coefficient scheme has a
point, and not yet the alternate-chart right-component calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import pickle
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPLAY = (
    ROOT
    / "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)
CHECKPOINT = REPLAY / "case1_checkpoint.pkl"
EXACT_CORE = REPLAY / "exact_core.py"
CHECKPOINT_SHA256 = (
    "2dcf13d924530cdc9a8728e943efdc73"
    "d003ce1c187d5cec273f6f701e0240ba"
)
EXACT_CORE_SHA256 = (
    "3ba2d44e52a8028044dd73a9394449f6"
    "a26c638aebed6cc9f09288e49d77ff82"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_exact_state() -> tuple[Any, dict[int, list[Any]], dict[int, list[Any]], int]:
    """Load the pinned encoded checkpoint into its original exact types."""

    if sha256(CHECKPOINT) != CHECKPOINT_SHA256:
        raise RuntimeError("the pinned Case-1 checkpoint hash changed")
    if sha256(EXACT_CORE) != EXACT_CORE_SHA256:
        raise RuntimeError("the pinned exact_core.py hash changed")

    try:
        from flint import fmpq, fmpq_poly
    except ImportError as error:
        raise RuntimeError("python-flint is required for the Case-1 continuation") from error

    sys.path.insert(0, str(REPLAY))
    exact_core = importlib.import_module("exact_core")

    def decode_field(value: list[tuple[int, int]]) -> Any:
        return exact_core.K(fmpq_poly([fmpq(numerator, denominator) for numerator, denominator in value]))

    def decode_parameter_polynomial(value: list[Any]) -> Any:
        return exact_core.PP(
            {
                tuple(monomial): decode_field(coefficient)
                for monomial, coefficient in value
            }
        )

    def decode_bands(value: dict[Any, list[Any]]) -> dict[int, list[Any]]:
        return {
            int(index): [decode_parameter_polynomial(item) for item in band]
            for index, band in value.items()
        }

    state = pickle.loads(CHECKPOINT.read_bytes())
    if state["done"] != [1, 0, -1, -2, -3] or state["np"] != 6:
        raise RuntimeError("the pinned checkpoint has an unexpected descent state")
    return (
        exact_core,
        decode_bands(state["P"]),
        decode_bands(state["Q"]),
        int(state["np"]),
    )


def add_bracket(exact_core: Any, p_band: list[Any], p_index: int, q_band: list[Any], q_index: int) -> list[Any]:
    """Return ``p_index*P*Q' - q_index*P'*Q`` exactly."""

    return exact_core.tadd(
        exact_core.tscale(
            exact_core.tmul(p_band, exact_core.tder(q_band)),
            exact_core.K(p_index),
        ),
        exact_core.tscale(
            exact_core.tmul(exact_core.tder(p_band), q_band),
            exact_core.K(-q_index),
        ),
    )


def support_p(index: int) -> list[int]:
    return list(range(0, max(0, 8 + index + 1)))


def support_q(index: int) -> list[int]:
    if index == 0:
        return list(range(1, 13))
    return list(range(0, max(0, 12 + index + 1)))


def parameter_term_count(polynomials: list[Any]) -> int:
    return sum(len(polynomial.d) for polynomial in polynomials)


def encode_field(value: Any) -> list[list[int]]:
    return [
        [int(value.p[index].p), int(value.p[index].q)]
        for index in range(len(value.p))
    ]


def encode_parameter_polynomial(value: Any) -> list[Any]:
    return [
        [list(monomial), encode_field(coefficient)]
        for monomial, coefficient in sorted(value.d.items())
    ]


def exact_hash(polynomials: list[Any]) -> str:
    payload = json.dumps(
        [encode_parameter_polynomial(polynomial) for polynomial in polynomials],
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def shape(polynomials: list[Any]) -> dict[str, int | str]:
    monomials = [monomial for polynomial in polynomials for monomial in polynomial.d]
    return {
        "coefficient_slots": len(polynomials),
        "nonzero_coefficients": sum(bool(polynomial.d) for polynomial in polynomials),
        "parameter_terms": parameter_term_count(polynomials),
        "maximum_parameter_degree": max((sum(monomial) for monomial in monomials), default=-1),
        "sha256": exact_hash(polynomials),
    }


def solve_layer(
    exact_core: Any,
    p_bands: dict[int, list[Any]],
    q_bands: dict[int, list[Any]],
    parameter_count: int,
    layer: int,
) -> tuple[list[Any], list[Any], list[Any], int, dict[str, Any]]:
    """Solve one new bracket layer using the archived triangular recurrence."""

    p_index = layer - 2
    q_index = layer - 1
    p_support = support_p(p_index)
    q_support = support_q(q_index)
    column_count = len(p_support) + len(q_support)
    if column_count == 0:
        raise ValueError(f"layer {layer} lies below both published polygons")

    started = time.monotonic()
    known: list[Any] = []
    contributing_pairs: list[tuple[int, int]] = []
    for old_p_index, old_p_band in p_bands.items():
        for old_q_index, old_q_band in q_bands.items():
            if old_p_index + old_q_index - 1 != layer:
                continue
            known = exact_core.tadd(
                known,
                add_bracket(
                    exact_core,
                    old_p_band,
                    old_p_index,
                    old_q_band,
                    old_q_index,
                ),
            )
            contributing_pairs.append((old_p_index, old_q_index))
    assembly_seconds = time.monotonic() - started

    degree_candidates = [len(known) - 1]
    if p_support:
        degree_candidates.append(max(p_support) + 11)
    if q_support:
        degree_candidates.append(7 + max(q_support))
    maximum_degree = max(degree_candidates)
    rows = [
        [exact_core.K(0) for _ in range(column_count)]
        for _ in range(maximum_degree + 1)
    ]
    for column, exponent in enumerate(p_support):
        for q_exponent, top_q_coefficient in enumerate(exact_core.Dvec):
            degree = exponent + q_exponent - 1
            if 0 <= degree <= maximum_degree:
                rows[degree][column] += top_q_coefficient * (
                    p_index * q_exponent - 3 * exponent
                )
    offset = len(p_support)
    for support_index, exponent in enumerate(q_support):
        for p_exponent, top_p_coefficient in enumerate(exact_core.Avec):
            degree = p_exponent + exponent - 1
            if 0 <= degree <= maximum_degree:
                rows[degree][offset + support_index] += top_p_coefficient * (
                    2 * exponent - q_index * p_exponent
                )

    retained_rows: list[list[Any]] = []
    right_hand_side: list[Any] = []
    for degree, row in enumerate(rows):
        known_coefficient = known[degree] if degree < len(known) else exact_core.PP()
        if any(row) or known_coefficient:
            retained_rows.append(row)
            right_hand_side.append(-known_coefficient)

    started = time.monotonic()
    solution, free_columns, compatibility, next_parameter_count = exact_core.linear_solve(
        retained_rows,
        right_hand_side,
        parameter_count,
    )
    solve_seconds = time.monotonic() - started
    p_band = [exact_core.PP() for _ in range(max(p_support, default=-1) + 1)]
    for column, exponent in enumerate(p_support):
        p_band[exponent] = solution[column]
    q_band = [exact_core.PP() for _ in range(max(q_support, default=-1) + 1)]
    for support_index, exponent in enumerate(q_support):
        q_band[exponent] = solution[offset + support_index]

    report = {
        "layer": layer,
        "P_band": p_index if p_support else None,
        "Q_band": q_index if q_support else None,
        "P_support": p_support,
        "Q_support": q_support,
        "contributing_pairs": [list(pair) for pair in contributing_pairs],
        "known_parameter_terms": parameter_term_count(known),
        "row_count": len(retained_rows),
        "column_count": column_count,
        "rank": column_count - len(free_columns),
        "nullity": len(free_columns),
        "parameters_before": parameter_count,
        "parameters_after": next_parameter_count,
        "compatibility_count": len(compatibility),
        "P_shape": shape(p_band) if p_support else None,
        "Q_shape": shape(q_band) if q_support else None,
        "compatibility_shape": shape(compatibility),
        "assembly_seconds": round(assembly_seconds, 3),
        "solve_seconds": round(solve_seconds, 3),
    }
    return p_band, q_band, compatibility, next_parameter_count, report


def continue_descent(
    stop_layer: int,
) -> tuple[dict[str, Any], dict[int, list[Any]], dict[int, list[Any]]]:
    if stop_layer > -4 or stop_layer < -11:
        raise ValueError("stop_layer must lie between -4 and -11")
    exact_core, p_bands, q_bands, parameter_count = load_exact_state()
    reports = []
    for layer in range(-4, stop_layer - 1, -1):
        print(f"START layer {layer}", flush=True)
        p_band, q_band, _, parameter_count, report = solve_layer(
            exact_core,
            p_bands,
            q_bands,
            parameter_count,
            layer,
        )
        if p_band:
            p_bands[layer - 2] = p_band
        if q_band:
            q_bands[layer - 1] = q_band
        reports.append(report)
        print(
            "DONE",
            f"layer={layer}",
            f"rank={report['rank']}/{report['column_count']}",
            f"compatibility={report['compatibility_count']}",
            f"parameters={report['parameters_after']}",
            f"assembly={report['assembly_seconds']}s",
            f"solve={report['solve_seconds']}s",
            flush=True,
        )
    report = {
        "status": "exact continuation of a necessary band recurrence",
        "scope_warning": (
            "the original thirteen compatibility equations already define the "
            "empty scheme; this derives the omitted formal bands but does not "
            "assert existence of a Case-1 point"
        ),
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "exact_core_sha256": EXACT_CORE_SHA256,
        "initial_layer": -4,
        "terminal_layer": stop_layer,
        "initial_parameter_count": 6,
        "final_parameter_count": parameter_count,
        "derived_P_bands": [index for index in sorted(p_bands) if index <= -6],
        "derived_Q_bands": [index for index in sorted(q_bands) if index <= -5],
        "layers": reports,
    }
    return report, p_bands, q_bands


def deterministic_ledger(report: dict[str, Any]) -> dict[str, Any]:
    """Remove wall-clock fields and retain the exact reproducibility ledger."""

    layer_fields = (
        "layer",
        "P_band",
        "Q_band",
        "P_support",
        "Q_support",
        "column_count",
        "rank",
        "nullity",
        "compatibility_count",
        "parameters_before",
        "parameters_after",
        "P_shape",
        "Q_shape",
        "compatibility_shape",
    )
    return {
        "schema_version": 1,
        "status": report["status"],
        "scope_warning": report["scope_warning"],
        "checkpoint_sha256": report["checkpoint_sha256"],
        "exact_core_sha256": report["exact_core_sha256"],
        "initial_layer": report["initial_layer"],
        "terminal_layer": report["terminal_layer"],
        "initial_parameter_count": report["initial_parameter_count"],
        "final_parameter_count": report["final_parameter_count"],
        "derived_P_bands": report["derived_P_bands"],
        "derived_Q_bands": report["derived_Q_bands"],
        "layers": [
            {field: layer[field] for field in layer_fields}
            for layer in report["layers"]
        ],
    }


def encoded_band_checkpoint(
    report: dict[str, Any],
    p_bands: dict[int, list[Any]],
    q_bands: dict[int, list[Any]],
) -> dict[str, Any]:
    """Return a Python-native checkpoint independent of FLINT pickling."""

    def encode_bands(bands: dict[int, list[Any]]) -> dict[int, list[Any]]:
        return {
            index: [encode_parameter_polynomial(item) for item in band]
            for index, band in bands.items()
        }

    return {
        "schema_version": 1,
        "checkpoint_sha256": report["checkpoint_sha256"],
        "exact_core_sha256": report["exact_core_sha256"],
        "terminal_layer": report["terminal_layer"],
        "parameter_count": report["final_parameter_count"],
        "P": encode_bands(p_bands),
        "Q": encode_bands(q_bands),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stop-layer",
        type=int,
        default=-4,
        help="last bracket layer to derive, from -4 through -11",
    )
    parser.add_argument(
        "--ledger-output",
        type=Path,
        help="write the deterministic exact hash ledger as JSON",
    )
    parser.add_argument(
        "--checkpoint-output",
        type=Path,
        help="write the reconstructed bands as a Python-native pickle",
    )
    args = parser.parse_args()
    report, p_bands, q_bands = continue_descent(args.stop_layer)
    if args.ledger_output is not None:
        args.ledger_output.write_text(
            json.dumps(deterministic_ledger(report), indent=2, sort_keys=True)
            + "\n"
        )
    if args.checkpoint_output is not None:
        args.checkpoint_output.write_bytes(
            pickle.dumps(
                encoded_band_checkpoint(report, p_bands, q_bands),
                protocol=5,
            )
        )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

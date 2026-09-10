#!/usr/bin/env python3
"""Fail-closed dispatcher for prospective constructed strict-class signals.

This dispatcher deliberately treats four facts as separate certificates:
an additional constructed strict class, rational solubility of its cover, a
new rational point, and finite-quotient independence of that point.  A
principal-ideal dependency is construction evidence, not any of the latter
three facts.  It reads no historical point, rank, or oracle-labelled input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA = "rank-jump.prospective-constructed-strict-signals.v1"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def evaluate(signal: dict) -> dict:
    """Evaluate one prospective candidate without deriving any new facts."""
    require(signal.get("schema") == SCHEMA, "unsupported signal schema")
    require(signal.get("provenance") == "PROSPECTIVE", "retrospective/oracle inputs are forbidden")
    require(signal.get("historical_point_inputs", 0) == 0, "historical point input is forbidden")
    require(signal.get("generic_strict_dimension") in (0, "UNKNOWN"), "invalid generic strict dimension")
    require(isinstance(signal.get("candidate_id"), str) and signal["candidate_id"], "missing candidate id")

    principal = signal.get("principal_dependency", {})
    strict = signal.get("constructed_strict_class", {})
    cover = signal.get("cover_solubility", {})
    point = signal.get("new_rational_point", {})
    quotient = signal.get("independent_quotient_direction", {})
    require(principal.get("status") in {"CERTIFIED", "UNKNOWN", "ABSENT"}, "invalid principal-dependency status")
    require(strict.get("status") in {"CERTIFIED_NONZERO", "UNKNOWN", "ABSENT"}, "invalid strict-class status")
    require(cover.get("status") in {"CERTIFIED_SOLUBLE", "UNKNOWN", "INSOLUBLE"}, "invalid cover-solubility status")
    require(point.get("status") in {"CERTIFIED_NEW", "UNKNOWN", "ABSENT"}, "invalid rational-point status")
    require(quotient.get("status") in {"CERTIFIED_INDEPENDENT", "UNKNOWN", "DEPENDENT", "ABSENT"}, "invalid quotient status")

    if strict["status"] == "CERTIFIED_NONZERO":
        require(principal["status"] == "CERTIFIED", "a constructed strict class requires a certified principal dependency")
        require(strict.get("certificate"), "missing strict-class certificate")
    if point["status"] == "CERTIFIED_NEW":
        require(cover["status"] == "CERTIFIED_SOLUBLE", "a new rational point requires a certified soluble cover")
        require(point.get("certificate"), "missing rational-point certificate")
    if quotient["status"] == "CERTIFIED_INDEPENDENT":
        require(point["status"] == "CERTIFIED_NEW", "quotient independence requires a certified new rational point")
        require(quotient.get("certificate"), "missing quotient-independence certificate")

    states = {
        "principal_dependency": principal["status"],
        "constructed_strict_class": strict["status"],
        "cover_solubility": cover["status"],
        "new_rational_point": point["status"],
        "independent_quotient_direction": quotient["status"],
    }
    if strict["status"] != "CERTIFIED_NONZERO":
        action = "WAIT_FOR_CONSTRUCTED_STRICT_CLASS"
    elif cover["status"] != "CERTIFIED_SOLUBLE":
        action = "CONSTRUCT_AND_AUDIT_COVER"
    elif point["status"] != "CERTIFIED_NEW":
        action = "SEARCH_OR_PROVE_NEW_RATIONAL_POINT"
    elif quotient["status"] != "CERTIFIED_INDEPENDENT":
        action = "CERTIFY_FINITE_QUOTIENT_INDEPENDENCE"
    else:
        action = "ADMIT_CERTIFIED_QUOTIENT_DIRECTION_TO_FOUNDRY"

    return {
        "schema": "rank-jump.prospective-constructed-strict-seed-gate.v1",
        "status": "PASS_FAIL_CLOSED_GATE",
        "candidate_id": signal["candidate_id"],
        "generic_strict_dimension": signal["generic_strict_dimension"],
        "strict_nonzero_is_outside_generic_strict_subgroup": (
            signal["generic_strict_dimension"] == 0 and strict["status"] == "CERTIFIED_NONZERO"
        ),
        "gates": states,
        "next_action": action,
        "boundary": (
            "A constructed strict class is a prospective construction seed only. "
            "It neither proves cover solubility, supplies a rational point, nor "
            "proves an independent Mordell-Weil quotient direction."
        ),
    }


def build(input_path: Path, output_path: Path) -> None:
    signal = json.loads(input_path.read_text())
    result = evaluate(signal)
    result["bindings"] = {str(input_path): digest(input_path), str(Path(__file__).resolve()): digest(Path(__file__).resolve())}
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x") as handle:
        handle.write(payload)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.input, args.output)

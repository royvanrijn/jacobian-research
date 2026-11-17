#!/usr/bin/env python3
"""Validate the cross-program support-saturation gate ledger."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "verified" / "SUPPORT_SATURATION_PATHS.json"
STATUS = ROOT / "MATH_STATUS.json"

ALLOWED_STAGE_STATUS = {
    "proved",
    "conditional",
    "open",
    "fails_for_current_case",
    "not_applicable",
}

EXPECTED_STAGES = {
    "cubic-keller": {
        "C0": "G0",
        "C1": "G2",
        "C2": "G3",
        "C3": "G4",
    },
    "plane-jc": {
        "P0": "G0",
        "P1": "G1",
        "P2": "G2",
        "P3": "G3",
        "P4": "G4",
    },
    "restricted-weyl": {
        "W0": "G0",
        "W1": "G1",
        "W2": "G0",
        "W3": "G2",
        "W4": "G3",
        "W5": "G4",
    },
}

EXPECTED_STAGE_STATUS = {
    "cubic-keller": {
        "C0": "proved",
        "C1": "open",
        "C2": "open",
        "C3": "conditional",
    },
    "plane-jc": {
        "P0": "open",
        "P1": "open",
        "P2": "open",
        "P3": "open",
        "P4": "conditional",
    },
    "restricted-weyl": {
        "W0": "proved",
        "W1": "fails_for_current_case",
        "W2": "open",
        "W3": "open",
        "W4": "open",
        "W5": "conditional",
    },
}

EXPECTED_NOT_APPLICABLE = {
    "cubic-keller": {"G1"},
    "plane-jc": set(),
    "restricted-weyl": set(),
}

EXPECTED_MODEL_FRONTIERS = {
    "cubic-keller": {"KDSQ6"},
    "plane-jc": set(),
    "restricted-weyl": set(),
}


def main() -> None:
    ledger = json.loads(LEDGER.read_text())
    assert ledger["schema"] == "support-saturation-paths.v1"

    status = json.loads(STATUS.read_text())
    by_id = {entry["id"]: entry for entry in status["entries"]}
    theorem_id = ledger["theorem_status_id"]
    formal_id = ledger["formal_core_status_id"]
    assert by_id[theorem_id]["state"] == "proved"
    assert by_id[formal_id]["state"] == "proved"
    assert formal_id in by_id[theorem_id]["dependencies"]

    canonical = ROOT / ledger["canonical_theorem"]
    assert canonical.is_file()
    canonical_text = canonical.read_text()
    assert LEDGER.name in canonical_text
    common = ledger["common_gate"]
    assert [gate["id"] for gate in common] == [
        "G0",
        "G1",
        "G2",
        "G3",
        "G4",
    ]
    common_ids = {gate["id"] for gate in common}
    assert all(gate["name"] and gate["requirement"] for gate in common)

    programmes = ledger["programmes"]
    assert set(programmes) == set(EXPECTED_STAGES)
    for name, expected in EXPECTED_STAGES.items():
        programme = programmes[name]
        assert programme["claim_status"] == "open"
        source = ROOT / programme["canonical_source"]
        assert source.is_file()
        source_text = source.read_text()
        assert LEDGER.name in source_text
        assert programme["module"]
        assert programme["support_ideal"]
        assert programme["target"]
        assert len(programme["failure_routes"]) >= 3

        stages = programme["stages"]
        stage_ids = [stage["id"] for stage in stages]
        assert len(stage_ids) == len(set(stage_ids))
        assert {stage["id"]: stage["gate"] for stage in stages} == expected
        assert {
            stage["id"]: stage["status"] for stage in stages
        } == EXPECTED_STAGE_STATUS[name]
        assert all(stage["gate"] in common_ids for stage in stages)
        assert all(stage["status"] in ALLOWED_STAGE_STATUS for stage in stages)
        assert all(stage["requirement"] for stage in stages)
        assert all(f"`{stage['id']}`" in source_text for stage in stages)

        omitted = programme.get("not_applicable_gates", [])
        omitted_ids = {item["gate"] for item in omitted}
        assert omitted_ids == EXPECTED_NOT_APPLICABLE[name]
        assert all(item["reason"] for item in omitted)
        covered_ids = {stage["gate"] for stage in stages} | omitted_ids
        assert covered_ids == common_ids

        frontiers = programme.get("model_frontiers", [])
        frontier_ids = {item["status_id"] for item in frontiers}
        assert frontier_ids == EXPECTED_MODEL_FRONTIERS[name]
        for item in frontiers:
            assert by_id[item["status_id"]]["state"] == "proved"
            assert item["scope"] and item["programme_effect"]
            assert item["C1"] and item["C2"]

        conclusion = stages[-1]
        assert conclusion["gate"] == "G4"
        assert conclusion["status"] == "conditional"
        assert any(
            stage["status"] in {"open", "fails_for_current_case"}
            for stage in stages[:-1]
        ), f"{name}: an open programme was accidentally marked complete"

    print("PASS: support-saturation paths use the common G0--G4 gate")
    print("PASS: cubic, plane-JC, and restricted-Weyl outcomes remain open")


if __name__ == "__main__":
    main()

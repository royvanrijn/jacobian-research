#!/usr/bin/env python3
"""Turn a complete ELKIESR28REL Magma transcript into a fail-closed gate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
DEFAULT_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
sys.path.insert(0, str(CAS))

from build_elkies_2026_rank28_relative_descent_magma import (  # noqa: E402
    build_magma,
    load_relative_input,
)
from elkies_residual_selmer_gate import SCHEMA, gate_record  # noqa: E402


PREFIX = "ELKIESR28REL|"


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_protocol(text: str) -> dict[str, Any]:
    records: list[dict[str, str]] = []
    for line in text.splitlines():
        if not line.startswith(PREFIX):
            continue
        record: dict[str, str] = {}
        for field in line[len(PREFIX) :].split("|"):
            if "=" not in field:
                raise ValueError("malformed ELKIESR28REL field")
            key, value = field.split("=", 1)
            if not key or key in record:
                raise ValueError("duplicate or empty ELKIESR28REL key")
            record[key] = value
        records.append(record)
    inputs = [record for record in records if record.get("stage") == "input"]
    completed = [
        record
        for record in records
        if record.get("stage") == "two_selmer" and record.get("status") == "complete"
    ]
    classifications = [record for record in records if "classification" in record]
    if len(inputs) != 1 or len(completed) != 1 or len(classifications) != 1:
        raise ValueError("transcript lacks one complete input/Selmer/classification chain")
    input_record, selmer_record, classification = inputs[0], completed[0], classifications[0]
    if input_record.get("version") != "1" or input_record.get("target_rank") != "32":
        raise ValueError("transcript uses another protocol version or target rank")
    if input_record.get("parameter") != "-9529/5471":
        raise ValueError("transcript belongs to another fibre")
    if input_record.get("generic") != "17" or input_record.get("known_quotient_floor") != "11":
        raise ValueError("transcript used another certified subgroup")
    total = int(selmer_record["total_selmer_dim"])
    residual = int(selmer_record["residual_dim"])
    required = int(selmer_record["required_residual_dim"])
    if total - 17 != residual or residual < 11 or required != 15:
        raise ValueError("transcript contains inconsistent Selmer dimensions")
    gate = gate_record(total_two_selmer_dimension=total)
    if classification.get("classification") != gate["status"]:
        raise ValueError("classification contradicts the exact gate")
    expected_authorization = "true" if gate["expensive_search_authorized"] else "false"
    if classification.get("expensive_search_authorized") != expected_authorization:
        raise ValueError("classification has the wrong search authorization")
    if int(classification["total_selmer_dim"]) != total or int(
        classification["residual_dim"]
    ) != residual:
        raise ValueError("classification repeats different Selmer dimensions")
    return {
        "input": input_record,
        "two_selmer": selmer_record,
        "classification": classification,
        "gate": gate,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", type=Path, default=DEFAULT_CONTROLS)
    parser.add_argument("--program", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = load_relative_input(args.controls)
    expected_program = build_magma(source)
    if args.program.read_text() != expected_program:
        raise SystemExit("the Magma program is not the exact source-pinned generated job")
    parsed = parse_protocol(args.log.read_text())
    if parsed["input"].get("controls_sha256") != source.controls_sha256:
        raise SystemExit("the transcript names another positive-control artifact")
    if parsed["input"].get("generic_sha256") != source.generic_point_sequence_sha256:
        raise SystemExit("the transcript names another generic subgroup")
    if parsed["input"].get("combined_sha256") != source.combined_point_sequence_sha256:
        raise SystemExit("the transcript names another rank-28 positive control")

    document = {
        "schema": SCHEMA,
        "status": parsed["gate"]["status"],
        "parameter": "-9529/5471",
        "global_minimal_model": [str(value) for value in source.model],
        "known_rank_lower_bound": 28,
        "known_additional_directions_beyond_generic_17": 11,
        "directions_still_needed_for_rank_32": 4,
        "positive_control_certificate": {
            "path": str(args.controls.resolve()),
            "sha256": source.controls_sha256,
        },
        "descent_backend": {
            "name": "Magma TwoSelmerGroup",
            "algorithm": "complete 2-Selmer group over Q with Bound=-1",
            "unconditional": True,
            "class_group_completeness_completed": True,
            "all_local_solubility_conditions_completed": True,
            "known_generic_kummer_dimension": 17,
            "proof_boundary": (
                "The completeness claims require the unique protocol completion line "
                "from the exact generated Bound=-1 job. Partial logs are rejected."
            ),
        },
        "backend_result": parsed["two_selmer"],
        "gate": parsed["gate"],
        "provenance": {
            "program": str(args.program.resolve()),
            "program_sha256": file_sha256(args.program),
            "log": str(args.log.resolve()),
            "log_sha256": file_sha256(args.log),
        },
        "stop_rule": (
            "The generated job performs no residual-cover construction before this "
            "gate and contains no rational-point search at any stage."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"ELKIESR28REL|parsed=true|status={document['status']}|"
        f"residual={parsed['gate']['residual_two_selmer_quotient_dimension']}|"
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()

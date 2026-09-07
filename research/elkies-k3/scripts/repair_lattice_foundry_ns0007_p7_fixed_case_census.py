#!/usr/bin/env python3
"""Resolve timed-out NS0007 census cases with checked singleton replays."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import zlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "elkies-k3.lattice-foundry-ns0007-pole0-fixed-case-census-modp.v1"
STATUS_NAMES = {
    0: "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE",
    1: "NONUNIT_IDEAL",
    2: "FINITE_SOLUTION_SET",
    3: "POSITIVE_DIMENSIONAL_SOLUTION_SET",
    4: "TIMEOUT",
    5: "SOLVER_ERROR",
}
TERMINAL_LEADING_CODES = {0, 1}
REPAIRABLE_CODES = {4, 5}


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def decode(row: dict) -> bytearray:
    statuses = bytearray(
        zlib.decompress(base64.b64decode(row["status_vector"]["base64"]))
    )
    if len(statuses) != row["search"]["selected_cases"]:
        raise ArithmeticError("status-vector length mismatch")
    if digest_bytes(bytes(statuses)) != row["status_vector"]["uncompressed_sha256"]:
        raise ArithmeticError("status-vector digest mismatch")
    return statuses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("replays", nargs="+", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    base_path = arguments.base.resolve()
    output = json.loads(base_path.read_text())
    if output.get("schema") != SCHEMA:
        raise ValueError("unexpected base census schema")
    if output["search"].get("solver_mode") != "leading-ideal":
        raise ValueError("base is not a leading-ideal census")
    statuses = decode(output)
    start = output["search"]["selected_start_index"]
    open_indices = {
        start + offset
        for offset, code in enumerate(statuses)
        if code in REPAIRABLE_CODES
    }
    if not open_indices:
        raise ValueError("base census has no timeout or solver-error cases")

    replay_manifest = []
    replayed_indices = set()
    replay_nonunit = []
    for replay_argument in arguments.replays:
        replay_path = replay_argument.resolve()
        replay = json.loads(replay_path.read_text())
        if replay.get("schema") != SCHEMA:
            raise ValueError(f"unexpected replay schema: {replay_path}")
        if replay["input"] != output["input"]:
            raise ArithmeticError(f"replay input mismatch: {replay_path}")
        if replay["prime"] != output["prime"] or replay["lambda"] != output["lambda"]:
            raise ArithmeticError(f"replay field or lambda mismatch: {replay_path}")
        if replay["search"].get("solver_mode") != "leading-ideal":
            raise ValueError(f"replay is not leading-ideal mode: {replay_path}")
        if replay["search"].get("selected_cases") != 1:
            raise ValueError(f"replay is not a singleton: {replay_path}")
        replay_index = replay["search"]["selected_start_index"]
        if replay["search"]["selected_stop_index_exclusive"] != replay_index + 1:
            raise ArithmeticError(f"replay interval mismatch: {replay_path}")
        if replay_index not in open_indices or replay_index in replayed_indices:
            raise ArithmeticError(f"unexpected or duplicate replay index: {replay_index}")
        replay_status = decode(replay)
        if replay_status[0] not in TERMINAL_LEADING_CODES:
            raise ArithmeticError(f"replay did not terminate: {replay_path}")
        statuses[replay_index - start] = replay_status[0]
        replayed_indices.add(replay_index)
        if replay_status[0] == 1:
            replay_nonunit.extend(replay["exceptional_cases"])
        replay_manifest.append(
            {
                "index": replay_index,
                "path": relative(replay_path),
                "sha256": digest(replay_path),
                "resolved_status": STATUS_NAMES[replay_status[0]],
            }
        )
    if replayed_indices != open_indices:
        missing = sorted(open_indices - replayed_indices)
        raise ArithmeticError(f"not all open cases were replayed: {missing}")

    remaining_exceptional = [
        case
        for case in output["exceptional_cases"]
        if case["index"] not in replayed_indices
    ]
    remaining_exceptional.extend(replay_nonunit)
    remaining_exceptional.sort(key=lambda case: case["index"])
    histogram = Counter(statuses)
    if any(code not in TERMINAL_LEADING_CODES for code in statuses):
        raise ArithmeticError("nonterminal status remains after repair")
    compressed = zlib.compress(bytes(statuses), level=9)
    output["status"] = "PASS_BOUNDED_FIXED_CASE_RANGE_LEADING_IDEAL_CENSUS"
    output["search"]["resolved_case_replay_count"] = len(replay_manifest)
    output["accounting"] = {
        "status_histogram": {
            name: histogram[code]
            for code, name in STATUS_NAMES.items()
            if histogram[code]
        },
        "exceptional_case_count": len(remaining_exceptional),
    }
    output["status_vector"] = {
        "encoding": output["status_vector"]["encoding"],
        "base64": base64.b64encode(compressed).decode("ascii"),
        "uncompressed_sha256": digest_bytes(bytes(statuses)),
    }
    output["exceptional_cases"] = remaining_exceptional
    output["repair"] = {
        "base": relative(base_path),
        "base_sha256": digest(base_path),
        "singleton_replays": sorted(replay_manifest, key=lambda row: row["index"]),
    }
    output["reproduce"] = (
        "python3 elkies-k3/scripts/"
        "repair_lattice_foundry_ns0007_p7_fixed_case_census.py --base "
        f"{relative(base_path)} --output {relative(arguments.output.resolve())} "
        + " ".join(relative(path.resolve()) for path in arguments.replays)
    )

    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("repaired NS0007 census artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007REPAIR|"
        f"cases={len(statuses)}|replays={len(replay_manifest)}|"
        f"nonunit={histogram[1]}|status={output['status']}"
    )


if __name__ == "__main__":
    main()

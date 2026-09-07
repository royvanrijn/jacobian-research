#!/usr/bin/env python3
"""Verify and combine contiguous NS0007 fixed-case leading-ideal shards."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import zlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-census-mod7.json"
)
SCHEMA = "elkies-k3.lattice-foundry-ns0007-pole0-fixed-case-census-modp.v1"
UNIT = 0
NONUNIT = 1


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def decode_statuses(row: dict) -> bytes:
    encoded = row["status_vector"]["base64"]
    statuses = zlib.decompress(base64.b64decode(encoded))
    expected = row["search"]["selected_cases"]
    if len(statuses) != expected:
        raise ArithmeticError("shard status-vector length mismatch")
    if digest_bytes(statuses) != row["status_vector"]["uncompressed_sha256"]:
        raise ArithmeticError("shard status-vector digest mismatch")
    return statuses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    loaded = []
    for path in arguments.shards:
        resolved = path.resolve()
        row = json.loads(resolved.read_text())
        if row.get("schema") != SCHEMA:
            raise ValueError(f"unexpected shard schema: {resolved}")
        if row["search"].get("solver_mode") != "leading-ideal":
            raise ValueError(f"shard is not a leading-ideal census: {resolved}")
        if row["search"].get("msolve_polynomial_encoding") != "FULLY_EXPANDED":
            raise ValueError(f"shard does not use expanded certifying syntax: {resolved}")
        if not row["search"].get("certifying_input_syntax"):
            raise ValueError(f"shard input syntax is diagnostic only: {resolved}")
        if row["status"] == "OPEN_TIMEOUT_OR_SOLVER_ERROR_CASES_REMAIN":
            raise ArithmeticError(f"open shard cannot be combined: {resolved}")
        loaded.append((row["search"]["selected_start_index"], resolved, row))
    loaded.sort()

    first = loaded[0][2]
    total_cases = first["search"]["total_cases"]
    invariant = {
        "prime": first["prime"],
        "lambda": first["lambda"],
        "fixed_variables": first["fixed_variables_in_enumeration_order"],
        "input": first["input"],
        "total_cases": total_cases,
    }
    next_index = 0
    status_parts = []
    exceptional = []
    shard_manifest = []
    for start, path, row in loaded:
        current = {
            "prime": row["prime"],
            "lambda": row["lambda"],
            "fixed_variables": row["fixed_variables_in_enumeration_order"],
            "input": row["input"],
            "total_cases": row["search"]["total_cases"],
        }
        if current != invariant:
            raise ArithmeticError(f"shard invariant mismatch: {path}")
        stop = row["search"]["selected_stop_index_exclusive"]
        if start != next_index:
            raise ArithmeticError(
                f"shard gap or overlap: expected start {next_index}, found {start}"
            )
        statuses = decode_statuses(row)
        if stop - start != len(statuses):
            raise ArithmeticError(f"shard interval length mismatch: {path}")
        counted = Counter(statuses)
        declared = row["accounting"]["status_histogram"]
        if counted[UNIT] != declared.get("NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE", 0):
            raise ArithmeticError(f"unit-ideal accounting mismatch: {path}")
        if counted[NONUNIT] != declared.get("NONUNIT_IDEAL", 0):
            raise ArithmeticError(f"nonunit-ideal accounting mismatch: {path}")
        if any(code not in (UNIT, NONUNIT) for code in statuses):
            raise ArithmeticError(f"nonterminal status remains in shard: {path}")
        if len(row["exceptional_cases"]) != counted[NONUNIT]:
            raise ArithmeticError(f"exceptional-case accounting mismatch: {path}")
        for case in row["exceptional_cases"]:
            if not start <= case["index"] < stop:
                raise ArithmeticError(f"exceptional case outside shard: {path}")
        status_parts.append(statuses)
        exceptional.extend(row["exceptional_cases"])
        shard_manifest.append(
            {
                "path": relative(path),
                "sha256": digest(path),
                "start_index": start,
                "stop_index_exclusive": stop,
            }
        )
        next_index = stop

    if next_index != total_cases:
        raise ArithmeticError(
            f"incomplete fixed-case coverage: stopped at {next_index}/{total_cases}"
        )
    statuses = b"".join(status_parts)
    histogram = Counter(statuses)
    compressed = zlib.compress(statuses, level=9)
    nonunit_count = histogram[NONUNIT]
    status = (
        "PASS_EXHAUSTIVE_FIXED_CASE_UNIT_IDEAL_CENSUS"
        if not nonunit_count
        else "PASS_EXHAUSTIVE_FIXED_CASE_LEADING_IDEAL_CENSUS_WITH_NONUNIT_CASES"
    )
    output = {
        "schema": SCHEMA,
        "status": status,
        "input": invariant["input"],
        "prime": invariant["prime"],
        "lambda": invariant["lambda"],
        "fixed_variables_in_enumeration_order": invariant["fixed_variables"],
        "search": {
            "total_cases": total_cases,
            "selected_start_index": 0,
            "selected_stop_index_exclusive": total_cases,
            "selected_cases": total_cases,
            "exhaustive": True,
            "solver_mode": "leading-ideal",
            "verified_contiguous_shard_count": len(loaded),
        },
        "accounting": {
            "status_histogram": {
                "NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE": histogram[UNIT],
                **({"NONUNIT_IDEAL": nonunit_count} if nonunit_count else {}),
            },
            "exceptional_case_count": len(exceptional),
        },
        "status_vector": {
            "encoding": (
                "zlib-compressed bytes in full lexicographic product order; "
                "0=NO_SOLUTION_OVER_ALGEBRAIC_CLOSURE, 1=NONUNIT_IDEAL"
            ),
            "base64": base64.b64encode(compressed).decode("ascii"),
            "uncompressed_sha256": digest_bytes(statuses),
        },
        "exceptional_cases": sorted(exceptional, key=lambda case: case["index"]),
        "shards": shard_manifest,
        "proof_boundary": {
            "proved": (
                "The declared fixed-lambda finite-field chart is partitioned into "
                "all p^6 assignments of the six displayed base-field coordinates. "
                "For every assignment, exact finite-field Groebner-basis computation "
                "classifies the expanded ideal as unit or nonunit."
            ),
            "open": (
                "Any nonunit case requires clean solution decoding, base-field "
                "rationality, exact Kodaira orders, residual squarefreeness, NS "
                "marking, and a characteristic-zero lift."
            ),
        },
        "reproduce": (
            "python3 elkies-k3/scripts/"
            "combine_lattice_foundry_ns0007_p7_fixed_case_census.py "
            + " ".join(relative(path) for _, path, _ in loaded)
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("combined NS0007 fixed-case artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007COMBINE|"
        f"cases={total_cases}|shards={len(loaded)}|"
        f"nonunit={nonunit_count}|status={status}"
    )


if __name__ == "__main__":
    main()

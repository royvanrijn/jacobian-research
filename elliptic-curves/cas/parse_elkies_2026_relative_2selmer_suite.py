#!/usr/bin/env python3
"""Parse completed ELKIESR17REL2 transcripts into an exact suite result."""

from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path
import re
from typing import Any, Iterable, Sequence


INPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-result.v1"
PROTOCOL = "ELKIESR17REL2"
PREFIX = f"{PROTOCOL}|"
BITS_RE = re.compile(r"^\[\s*([01](?:\s*,\s*[01])*)?\s*\]$")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_record(line: str) -> dict[str, str]:
    if not line.startswith(PREFIX):
        raise ValueError("not an ELKIESR17REL2 record")
    record: dict[str, str] = {}
    for field in line[len(PREFIX) :].split("|"):
        if "=" not in field:
            raise ValueError(f"malformed protocol field: {field}")
        key, value = field.split("=", 1)
        if not key or key in record:
            raise ValueError("empty or duplicate protocol key")
        record[key] = value
    return record


def protocol_records(text: str) -> list[dict[str, str]]:
    return [parse_record(line) for line in text.splitlines() if line.startswith(PREFIX)]


def parse_bits(value: str, dimension: int) -> tuple[int, ...]:
    match = BITS_RE.fullmatch(value)
    if match is None:
        raise ValueError(f"malformed GF(2) vector: {value}")
    body = match.group(1)
    bits = () if body is None else tuple(int(item.strip()) for item in body.split(","))
    if len(bits) != dimension:
        raise ValueError(f"expected {dimension} bits, found {len(bits)}")
    return bits


def gf2_basis(rows: Iterable[Sequence[int]]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return basis


def gf2_rank(rows: Iterable[Sequence[int]]) -> int:
    return len(gf2_basis(rows))


def gf2_contains(basis: dict[int, int], row: Sequence[int]) -> bool:
    value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return False
        value ^= basis[pivot]
    return True


def all_bits(dimension: int) -> Iterable[tuple[int, ...]]:
    return itertools.product((0, 1), repeat=dimension)


def one(records: list[dict[str, str]], **filters: str) -> dict[str, str]:
    matches = [
        record
        for record in records
        if all(record.get(key) == value for key, value in filters.items())
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one protocol record for {filters}, found {len(matches)}")
    return matches[0]


def parse_case(text: str, expected: dict[str, Any]) -> dict[str, Any]:
    records = protocol_records(text)
    input_record = one(records, stage="input")
    if (
        input_record.get("version") != "1"
        or input_record.get("case") != expected["case_id"]
        or input_record.get("parameter") != expected["parameter"]
        or input_record.get("role") != expected["role"]
    ):
        raise ValueError("transcript input identity does not match the suite manifest")
    selmer = one(records, stage="two_selmer", status="complete")
    blind_plan = one(records, stage="blind_plan")
    blind_end = one(records, stage="blind_end")
    classification = one(records, stage="classification", status="complete")
    total_dim = int(selmer["total_dim"])
    generic_rank = int(selmer["generic_kummer_rank"])
    residual_dim = int(selmer["residual_dim"])
    if generic_rank != 17 or total_dim - generic_rank != residual_dim:
        raise ValueError("inconsistent generic or residual Selmer dimension")

    generic_rows = [record for record in records if record.get("stage") == "generic_class"]
    quotient_basis_records = [
        record for record in records if record.get("stage") == "quotient_basis"
    ]
    exceptional_records = [
        record for record in records if record.get("stage") == "exceptional_class"
    ]
    cover_records = [record for record in records if record.get("stage") == "blind_cover"]
    if len(generic_rows) != 17 or len(quotient_basis_records) != residual_dim:
        raise ValueError("transcript omitted a generic or quotient basis row")
    expected_exceptional = int(expected["held_out_exceptional_point_count"])
    if len(exceptional_records) != expected_exceptional:
        raise ValueError("transcript omitted a held-out exceptional class")
    target_count = int(blind_plan["target_count"])
    if len(cover_records) != target_count or int(blind_end["target_count"]) != target_count:
        raise ValueError("transcript omitted a planned blind cover")

    quotient_basis = [
        parse_bits(record["quotient_bits"], residual_dim)
        for record in quotient_basis_records
    ]
    if gf2_rank(quotient_basis) != residual_dim:
        raise ValueError("reported quotient basis is not a basis")
    exceptional_rows = [
        parse_bits(record["quotient_bits"], residual_dim)
        for record in exceptional_records
    ]
    exceptional_basis = gf2_basis(exceptional_rows)
    exceptional_rank = len(exceptional_basis)
    if int(classification["exceptional_quotient_rank"]) != exceptional_rank:
        raise ValueError("exceptional quotient rank contradicts its displayed rows")
    unexplained_dim = residual_dim - exceptional_rank
    if int(classification["unexplained_dim"]) != unexplained_dim:
        raise ValueError("unexplained dimension is inconsistent")
    known_realized_count = 1 << exceptional_rank
    unrealized_count = (1 << residual_dim) - known_realized_count
    if (
        int(classification["known_realized_class_count"]) != known_realized_count
        or int(classification["unrealized_class_count"]) != unrealized_count
    ):
        raise ValueError("class counts contradict the quotient subspaces")

    covers = []
    recovered_rows: list[tuple[int, ...]] = []
    for record in cover_records:
        quotient_bits = parse_bits(record["quotient_bits"], residual_dim)
        recovered_bits = None
        if record["search_status"] == "point_found":
            recovered_bits = parse_bits(record["recovered_quotient_bits"], residual_dim)
            recovered_rows.append(recovered_bits)
        elif record["search_status"] != "no_point_within_bound":
            raise ValueError("unknown cover-search status")
        covers.append(
            {
                "index": int(record["index"]),
                "label": record["label"],
                "quotient_bits": list(quotient_bits),
                "known_exceptional_realized": gf2_contains(exceptional_basis, quotient_bits),
                "alpha": record["alpha"],
                "quartic_f": record["quartic_f"],
                "quartic_h": record["quartic_h"],
                "construction_seconds": float(record["construction_seconds"]),
                "search_seconds": float(record["search_seconds"]),
                "search_status": record["search_status"],
                "cover_point": record.get("cover_point"),
                "elliptic_point": record.get("elliptic_point"),
                "recovered_quotient_bits": (
                    None if recovered_bits is None else list(recovered_bits)
                ),
            }
        )
    recovered_rank = gf2_rank(recovered_rows)
    if (
        recovered_rank != int(blind_end["recovered_quotient_rank"])
        or recovered_rank != int(classification["blind_recovered_rank"])
    ):
        raise ValueError("blind recovered rank contradicts its displayed rows")
    intersection_rank = (
        recovered_rank
        + exceptional_rank
        - gf2_rank([*recovered_rows, *exceptional_rows])
    )

    unrealized_classes: list[list[int]] | None = None
    if residual_dim <= 16:
        unrealized_classes = [
            list(bits)
            for bits in all_bits(residual_dim)
            if not gf2_contains(exceptional_basis, bits)
        ]
    enumerate_all = blind_plan["enumerate_all"] == "true"
    if blind_plan["enumerate_all"] not in ("true", "false"):
        raise ValueError("malformed enumerate_all flag")
    if enumerate_all and target_count != (1 << residual_dim) - 1:
        raise ValueError("an exhaustive blind plan has the wrong target count")

    return {
        "case_id": expected["case_id"],
        "role": expected["role"],
        "parameter": expected["parameter"],
        "global_minimal_model": expected["global_minimal_model"],
        "certified_rank_lower_bound": expected["certified_rank_lower_bound"],
        "nagao_record": expected.get("nagao_record"),
        "two_selmer": {
            "dimension": total_dim,
            "two_torsion_dimension": int(selmer["two_torsion_dim"]),
            "seconds": float(selmer["seconds"]),
            "factor_base_size": int(selmer["factor_base_size"]),
            "unconditional_bound_parameter": -1,
        },
        "specialized_generic_subgroup": {
            "point_count": 17,
            "kummer_rank": generic_rank,
            "selmer_coordinates": [record["selmer_bits"] for record in generic_rows],
        },
        "relative_quotient": {
            "dimension": residual_dim,
            "basis": [
                {
                    "selmer_bits": record["selmer_bits"],
                    "quotient_bits": list(parse_bits(record["quotient_bits"], residual_dim)),
                    "alpha": record["alpha"],
                }
                for record in quotient_basis_records
            ],
        },
        "held_out_exceptional_points": {
            "count": expected_exceptional,
            "quotient_rank": exceptional_rank,
            "quotient_classes": [list(row) for row in exceptional_rows],
            "known_realized_class_count_including_zero": known_realized_count,
        },
        "unrealized_quotient": {
            "dimension_beyond_known_exceptional_span": unexplained_dim,
            "class_count": unrealized_count,
            "classes": unrealized_classes,
            "classes_omitted_reason": (
                None if unrealized_classes is not None else "residual dimension exceeds 16"
            ),
        },
        "blind_cover_benchmark": {
            "search_bound": int(input_record["search_bound"]),
            "enumerated_every_nonzero_quotient_class": enumerate_all,
            "target_count": target_count,
            "total_seconds": float(blind_end["seconds"]),
            "recovered_class_count": int(blind_end["recovered_class_count"]),
            "recovered_quotient_rank": recovered_rank,
            "recovered_known_exceptional_direction_rank": intersection_rank,
            "covers": covers,
        },
        "status": "PASS_COMPLETE_RELATIVE_2SELMER_AND_BLIND_COVER_BENCHMARK",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit("unexpected relative 2-Selmer input manifest")
    results = []
    logs = []
    for expected in manifest["cases"]:
        program = Path(expected["program"])
        if file_sha256(program) != expected["program_sha256"]:
            raise SystemExit(f"generated program changed: {program}")
        log = args.log_dir / f"{expected['case_id']}.log"
        result = parse_case(log.read_text(), expected)
        results.append(result)
        logs.append({"path": str(log), "sha256": file_sha256(log)})
    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_COMPLETE_RELATIVE_2SELMER_SUITE",
        "input_manifest": {"path": str(args.manifest), "sha256": file_sha256(args.manifest)},
        "logs": logs,
        "cases": results,
        "case_count": len(results),
        "claim_boundary": (
            "Selmer dimensions and quotient classes are complete unconditional descent "
            "computations. Cover searches are bounded; a miss is not a proof of no rational point."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"{PROTOCOL}|stage=parse_complete|cases={len(results)}|output={args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate the pinned transcription of Banks' Table C.1."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "arithmetic/banks_degree_5_10_candidates.json"
SCHEMA = "banks-strongly-intersective-candidates/v1"
EXPECTED_SHAPES = {
    5: {(2, 3)},
    6: {(2, 2, 2)},
    7: {(2, 2, 3), (2, 5), (3, 4)},
    8: {(2, 2, 2, 2), (2, 2, 4), (2, 3, 3)},
    9: {(2, 2, 2, 3), (2, 2, 5), (2, 3, 4), (2, 7), (4, 5)},
    10: {
        (2, 2, 2, 2, 2),
        (2, 2, 2, 4),
        (2, 2, 3, 3),
        (2, 2, 6),
        (2, 3, 5),
        (2, 4, 4),
        (3, 3, 4),
        (3, 7),
        (4, 6),
    },
}
EXPECTED_GROUP_COUNTS = {
    5: 1,
    6: 1,
    7: 7,
    8: 9,
    9: 33,
    10: 112,
}


def main() -> None:
    data = json.loads(DATA.read_text())
    assert data["schema"] == SCHEMA
    assert "Necessary" in data["scope"]
    assert "not an if-and-only-if" in data["scope"]
    source = data["source"]
    assert source["year"] == 2025
    assert re.fullmatch(r"[0-9a-f]{64}", source["thesis_pdf_sha256"])
    assert re.fullmatch(r"[0-9a-f]{40}", source["code_main_commit"])
    assert source["retrieved"] == "2026-07-27"

    degrees = {entry["degree"]: entry for entry in data["degrees"]}
    assert set(degrees) == set(range(5, 11))
    total_groups = 0
    total_examples = 0
    for degree, expected_shapes in EXPECTED_SHAPES.items():
        rows = degrees[degree]["shapes"]
        actual_shapes = {
            tuple(row["factor_degrees"])
            for row in rows
        }
        assert actual_shapes == expected_shapes
        group_count = 0
        for row in rows:
            shape = row["factor_degrees"]
            assert shape == sorted(shape)
            assert all(part >= 2 for part in shape)
            assert sum(shape) == degree
            group_ids = [group["id"] for group in row["groups"]]
            assert group_ids
            assert len(group_ids) == len(set(group_ids))
            assert all(
                re.fullmatch(r"[0-9]+\.(?:[0-9]+|[a-z]+)", group_id)
                for group_id in group_ids
            )
            assert all(
                isinstance(group["example_in_table"], bool)
                for group in row["groups"]
            )
            group_count += len(group_ids)
            total_examples += sum(
                group["example_in_table"] for group in row["groups"]
            )
        assert group_count == EXPECTED_GROUP_COUNTS[degree]
        total_groups += group_count
        print(
            f"PASS degree {degree}: {len(rows)} shapes, "
            f"{group_count} group rows"
        )
    print(
        f"PASS pinned Banks transcription: {total_groups} candidate rows, "
        f"{total_examples} rows with a displayed thesis example"
    )
    print("PASS every row remains labelled necessary, not sufficient")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit the committed HC4PPG1--9 artifact chain without symbolic replay."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
EXPECTED = {
    "hc4_projective_polar_atlas.json": (
        "hc4-projective-polar-atlas-v5",
        "350bc81b4ba7ac21289d7548f6d46de6526887c4e98a0813596cf20a454b240b",
    ),
    "hc4_quintic_infinity_rees_strata.json": (
        "hc4-quintic-infinity-rees-strata-v1",
        "51ddbf2b7c0c2b9b3f2cd7c1a8dcb4bf0fe97a3e3bce306eeee031cd5c92b99d",
    ),
    "hc4_rank3_vertex_colength.json": (
        "hc4-rank3-vertex-colength-v1",
        "c610f57af67061d0b4eb9523cb018569a7e8220a51dbd2350b71eb7007bfe473",
    ),
    "hc4_codim3_gradient_strata.json": (
        "hc4-codim3-gradient-strata-v1",
        "8759875cf431d18f35321631984d9120c72a2335dcae31d107fa191ae539e5a3",
    ),
    "hc4_binary_root_partition_segre.json": (
        "hc4-binary-root-partition-segre-v1",
        "09f6a57c735b2751d0f890b8cd216822001bae875fd9a5156e2a27550f8e71ad",
    ),
}
HELPER = ROOT / "jcsearch" / "projective_gradient_segre.py"
EXPECTED_HELPER_SHA256 = (
    "53801ca63b1df72d9664336810d7130a9293f1413ffc37dbfcaaa54e6cce45ad"
)


def load(filename: str) -> dict[str, object]:
    expected_format, expected_hash = EXPECTED[filename]
    path = GENERATED / filename
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual_hash == expected_hash, (filename, actual_hash, expected_hash)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["format"] == expected_format, filename
    return payload


def main() -> None:
    helper_hash = hashlib.sha256(HELPER.read_bytes()).hexdigest()
    assert helper_hash == EXPECTED_HELPER_SHA256, (
        helper_hash,
        EXPECTED_HELPER_SHA256,
    )
    atlas = load("hc4_projective_polar_atlas.json")
    strata = load("hc4_quintic_infinity_rees_strata.json")
    vertex = load("hc4_rank3_vertex_colength.json")
    codim3 = load("hc4_codim3_gradient_strata.json")
    binary = load("hc4_binary_root_partition_segre.json")

    remaining = atlas["quintic_coverage_summary"][
        "remaining_numerical_signatures_after_vertex_colength"
    ]
    assert remaining == {
        "affine_degree_2": 318,
        "affine_degree_3": 306,
        "total": 624,
    }
    assert "Lower layers can change" in strata["scope"]
    assert vertex["atlas_intersection"]["affine_degree_2"]["rows_after"] == 318
    assert vertex["atlas_intersection"]["affine_degree_3"]["rows_after"] == 306
    assert "No unconditional codimension-three atlas row" in codim3["scope"]
    assert binary["unconditional_rows_excluded"] == 0
    assert "remains open" in binary["scope"]

    print(
        "PASS: committed HC4PPG1--9 atlas, Rees, vertex, codimension-three, "
        "and binary-root ledgers plus their shared helper match exact hashes; "
        "no SymPy import, Gröbner calculation, or artifact rewrite"
    )
    print(
        "SCOPE: 624 rows are necessary numerical configurations, not existence "
        "results; lower-layer torsion and exceptional codimension-three packets "
        "remain explicit"
    )


if __name__ == "__main__":
    main()

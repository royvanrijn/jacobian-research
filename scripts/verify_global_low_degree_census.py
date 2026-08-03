#!/usr/bin/env python3
"""Recompute and verify the pinned low-degree census ledgers."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.low_degree_pipeline import (  # noqa: E402
    MANIFEST_FILENAME,
    STAGE_FILENAMES,
    build_low_degree_census,
)
from jcsearch.low_degree_census import sha256_json  # noqa: E402


ARTIFACT_ROOT = ROOT / "artifacts/generated-results"


def main() -> None:
    expected = build_low_degree_census(
        progress=lambda message: print(message, flush=True)
    )
    manifest_path = ARTIFACT_ROOT / MANIFEST_FILENAME
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "global-low-degree-census.manifest.v1"
    assert tuple(manifest["stage_sha256"]) == STAGE_FILENAMES

    for filename in STAGE_FILENAMES:
        path = ARTIFACT_ROOT / filename
        content = path.read_bytes()
        actual_hash = hashlib.sha256(content).hexdigest()
        assert actual_hash == manifest["stage_sha256"][filename], filename
        assert json.loads(content) == expected[filename], filename

    profiles = expected[STAGE_FILENAMES[0]]
    supports = expected[STAGE_FILENAMES[1]]
    exact = expected[STAGE_FILENAMES[6]]
    boundary = expected[STAGE_FILENAMES[7]]
    assert profiles["profile_count"] == 74
    assert len(profiles["profile_rank_gates"]) == 74
    assert profiles["profile_rank_gates"][-1] == {
        "profile": [7, 6, 3],
        "rank_above_degree": {
            "0": 3,
            "1": 3,
            "2": 3,
            "3": 2,
            "4": 2,
            "5": 2,
            "6": 1,
            "7": 0,
        },
    }
    assert supports["determinant_balanced_supports_by_size"] == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 30,
        "5": 85,
        "6": 1694,
    }
    assert supports["determinant_balanced_orbits_by_size"] == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 15,
        "5": 47,
        "6": 851,
    }
    buckets = expected[STAGE_FILENAMES[2]]
    assert all(row["buckets"] for row in buckets["representatives"])
    assert all(
        row["bucket_sha256"] == sha256_json(row["buckets"])
        for row in buckets["representatives"]
    )
    assert exact["representative_count"] == 913
    assert exact["unit_ideal_count"] == 913
    assert exact["surviving_support_ids"] == []
    assert exact["sparse_frontier_theorem"]["lower_bound_nonlinear_support"] == 7
    assert exact["sparse_frontier_theorem"]["attainment_proved"] is False
    assert exact["completely_eliminated_profiles"] == [
        [1, 1, 1],
        [2, 1, 1],
        [2, 2, 1],
        [2, 2, 2],
    ]
    assert not boundary["survives_projective_boundary_analysis"]
    print(
        "PASS global low-degree census: 74 profiles; sparse supports "
        "30/85/1694; 913 exact orbit ideals are units; support lower bound seven"
    )


if __name__ == "__main__":
    main()

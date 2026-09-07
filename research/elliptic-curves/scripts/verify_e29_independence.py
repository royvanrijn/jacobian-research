#!/usr/bin/env python3
"""Replay the pinned exact rank-at-least-29 calibration certificate."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.record_rank29 import verify_rank29_manifest  # noqa: E402


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "elkies_klagsbrun_e29_independence_v1.json"
)


def main() -> None:
    artifact_bytes = ARTIFACT.read_bytes()
    manifest = json.loads(artifact_bytes)
    point_data = json.loads(
        (PROGRAM_ROOT / "data" / "elkies_klagsbrun_e29_points.json").read_text(
            encoding="utf-8"
        )
    )
    assert hashlib.sha256(artifact_bytes).hexdigest() == (
        point_data["independence_certificate"]["artifact_sha256"]
    )
    verify_rank29_manifest(manifest)
    certificate = manifest["independence_certificate"]
    assert certificate["relation_prime"] == 2
    assert len(certificate["rows"]) == 29
    assert manifest["target_status"].endswith("no thirtieth point was found")
    print(
        "PASS E29 independence: 29 published points, exact full-rank mod-2 "
        "reduction certificate; no exact-rank or rank-30 claim"
    )


if __name__ == "__main__":
    main()

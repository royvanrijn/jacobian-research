#!/usr/bin/env python3
"""Replay the pinned exact Kihara rank-at-least-14 certificate."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.kihara import verify_kihara_rank14_manifest  # noqa: E402


ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "kihara_rank14_t2_v1.json"
)
FAMILY_METADATA = PROGRAM_ROOT / "families" / "kihara_rank14.json"


def main() -> None:
    artifact_bytes = ARTIFACT.read_bytes()
    manifest = json.loads(artifact_bytes)
    family = json.loads(FAMILY_METADATA.read_text(encoding="utf-8"))
    assert hashlib.sha256(artifact_bytes).hexdigest() == (
        family["certificate_specialization"]["artifact_sha256"]
    )
    verify_kihara_rank14_manifest(manifest)
    certificate = manifest["independence_certificate"]
    assert certificate["relation_prime"] == 5
    assert len(certificate["rows"]) == 14
    assert manifest["target_status"].startswith("baseline only")
    print(
        "PASS Kihara rank-14 replay: fourteen exact points, full-rank mod-5 "
        "finite-reduction certificate; baseline only, no rank-30 claim"
    )


if __name__ == "__main__":
    main()

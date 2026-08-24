#!/usr/bin/env python3
"""Recompute and byte-check the pinned compressed curve-302 certificate."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_icarm_curve302_rank31 import build_certificate  # noqa: E402

PINNED = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve302_rank31_v1.json.gz"
)
EXPECTED_COMPRESSED_SHA256 = "fc50b4b9ec5fe1dd1fe31aa299f13d8bc3476d43f3ed98e2ade5a4fc8972aa04"
EXPECTED_JSON_SHA256 = "3be0d6fe82c58e0f9284df5d9340332944a1d906508ea986d4abe00357036991"


def main() -> None:
    compressed = PINNED.read_bytes()
    if hashlib.sha256(compressed).hexdigest() != EXPECTED_COMPRESSED_SHA256:
        raise SystemExit("compressed curve-302 certificate hash mismatch")
    rendered = gzip.decompress(compressed)
    if hashlib.sha256(rendered).hexdigest() != EXPECTED_JSON_SHA256:
        raise SystemExit("decompressed curve-302 certificate hash mismatch")

    pinned = json.loads(rendered)
    computed = build_certificate()
    # Python's patch release is provenance, not mathematics.  Preserve the
    # pinned producer version before doing an otherwise byte-for-byte replay.
    computed["generation"]["python"] = pinned["generation"]["python"]
    expected = (json.dumps(computed, indent=2, sort_keys=True) + "\n").encode()
    if expected != rendered:
        raise SystemExit("recomputed curve-302 certificate differs from pinned JSON")
    print(
        "R31ICARM|stage=pinned|compressed=PASS|json=PASS|recompute=PASS|"
        "rank_lower_bound=31|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

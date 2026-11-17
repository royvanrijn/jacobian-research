#!/usr/bin/env python3
"""Recompute and compare the pinned multiplicative Hasse count artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from count_multiplicative_hasse_parameters import enumerate_parameters


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "multiplicative_hasse_parameters_1000000.json"
)
BOUND = 1_000_000


def main() -> int:
    expected = json.dumps(
        enumerate_parameters(BOUND),
        indent=2,
        sort_keys=True,
    ) + "\n"
    actual = ARTIFACT.read_text(encoding="utf-8")
    if actual != expected:
        print(
            "FAIL: pinned multiplicative Hasse artifact is stale; "
            "run the documented refresh command"
        )
        return 1
    digest = hashlib.sha256(actual.encode("utf-8")).hexdigest()
    print(
        "PASS: pinned multiplicative Hasse artifact matches the "
        f"dependency-free enumeration through {BOUND}"
    )
    print(f"SHA-256: {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

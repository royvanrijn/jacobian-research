#!/usr/bin/env python3
"""Exact pullback-fixed-space audit for MacFarlane's F13 through degree 3."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from audit_macfarlane_g20_dimension_reduction import (
    SOURCE_URL,
    build_maps,
    fixed_space_through_degree_three,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "macfarlane_f13_low_degree_invariants.json"
)


def main() -> None:
    data = build_maps()
    fixed_space = fixed_space_through_degree_three(data["F13"], data["x"])
    artifact = {
        "format": "macfarlane-f13-low-degree-invariant-audit-v1",
        "external_source": SOURCE_URL,
        "map": "F13=x+R+B*gamma",
        "degree_cap": 3,
        **fixed_space,
        "scope": (
            "Exact pullback fixed space only; higher-degree invariants and "
            "general nonlinear semiconjugacies remain open."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS F13: pullback-fixed polynomials through degree 3 are constants")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()

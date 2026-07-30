#!/usr/bin/env python3
"""Verify that the rotated cut-14 conormal extension splits through q4."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree42_ritt_rotated_conormal_flags import (  # noqa: E402
    source_ideal_cache,
)
from research_degree42_ritt_rotated_tensor_extension import (  # noqa: E402
    matrix_cache,
    rotated_tensor_audit,
)


WORD = (2, 3, 7)
SOURCE_CACHE = source_ideal_cache(WORD)
MATRIX_CACHE = matrix_cache(WORD, 4, 1)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_cut14_tensor_split_q4.json"
)


def main() -> None:
    result = rotated_tensor_audit(WORD, base_order=4, normal_order=1)
    assert result["dimensions"] == {
        "kernel": 9,
        "total": 13,
        "spectator": 4,
    }
    assert result["splits_as_R_module"] is True
    assert result["one_section"] is not None
    assert result["certificate"]["coboundary_rank"] == 32
    assert result["certificate"]["augmented_rank"] == 32
    assert result["certificate"]["obstruction_detected"] is False
    assert result["zero_action_variables"] == [f"n{i}" for i in range(7)]

    output = {
        "schema": "degree42-ritt-cut14-tensor-split-q4.v1",
        "status": "exact finite base-adic module splitting",
        "word": WORD,
        "thick_composite_omission": 14,
        "base_order": 4,
        "normal_order": 1,
        "source_ideal_cache": str(SOURCE_CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": hashlib.sha256(
            SOURCE_CACHE.read_bytes()
        ).hexdigest(),
        "matrix_cache": str(MATRIX_CACHE.relative_to(ROOT)),
        "matrix_cache_sha256": hashlib.sha256(
            MATRIX_CACHE.read_bytes()
        ).hexdigest(),
        "calculation": {
            "dimensions": result["dimensions"],
            "splits_as_R_module": result["splits_as_R_module"],
            "one_section": result["one_section"],
            "zero_action_variables": result["zero_action_variables"],
            "coboundary_rank": result["certificate"]["coboundary_rank"],
            "augmented_rank": result["certificate"]["augmented_rank"],
            "obstruction_detected": result["certificate"][
                "obstruction_detected"
            ],
        },
        "consequence": (
            "The tensor-presented cut-14 first-conormal extension splits "
            "over B/(tau,zeta)^4. Its reductions therefore split at every "
            "base order at most four."
        ),
        "theorem_boundary": (
            "A finite-order section does not prove that compatible sections "
            "exist at every order or that the completed extension splits."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_cut14_tensor_split_q4.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: cut-14 dimensions are 9 -> 13 -> 4 at base order four")
    print("PASS: the exact module extension splits through base order four")
    print("PASS: cocycle and coboundary ranks are both 32")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

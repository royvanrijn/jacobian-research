#!/usr/bin/env python3
"""Verify formal non-splitting of the degree-42 ideal-module extension."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_cellular_extension import CACHE  # noqa: E402
from research_degree42_tensor_extension import audit  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_tensor_extension.json"
)


def main() -> None:
    result = audit(base_order=2, normal_order=2)
    assert result["dimensions"] == {
        "kernel": 8,
        "total": 12,
        "spectator": 4,
    }
    assert result["splits_as_R_module"] is False
    certificate = result["certificate"]
    assert certificate["coboundary_rank"] == 27
    assert certificate["augmented_rank"] == 28
    assert certificate["obstruction_functional"] == {
        "n0[4,0]": 25,
        "n0[5,1]": 10,
        "n1[4,0]": -20,
        "n1[4,1]": 4,
        "n1[5,0]": -10,
    }
    assert certificate["obstruction_value"] == "240"

    output = {
        "schema": "degree42-ritt-tensor-extension.v1",
        "status": "exact finite tensor obstruction",
        "completed_local_ring": "Q[[n0,...,n6,tau,zeta]]",
        "module_sequence": (
            "0 -> I_boundary/I_6 -> K/I_6 -> K/I_boundary -> 0"
        ),
        "tensor_quotient": (
            "R / ((tau,zeta)^2 + (n0,...,n6)^2)"
        ),
        "source_ideal_cache": str(CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": hashlib.sha256(
            CACHE.read_bytes()
        ).hexdigest(),
        "calculation": result,
        "formal_consequence": (
            "the completed R-module sequence is non-split, because any "
            "R-linear splitting would remain a splitting after tensoring "
            "with the displayed finite quotient"
        ),
        "derived_boundary": (
            "this identifies the ideal-module extension, not yet the "
            "transitivity class of the full derived cotangent complexes"
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_tensor_extension.py"
        ),
    }
    ARTIFACT.write_text(
        json.dumps(
            output,
            indent=2,
            default=lambda value: (
                int(value) if value.is_Integer else str(value)
            ),
        )
        + "\n"
    )
    print("PASS: tensor quotient dimensions are 8 -> 12 -> 4")
    print("PASS: all nine coordinate actions commute with the projection")
    print("PASS: no R-linear section exists")
    print("PASS: obstruction rank rises from 27 to 28")
    print("PASS: the primitive obstruction functional evaluates to 240")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

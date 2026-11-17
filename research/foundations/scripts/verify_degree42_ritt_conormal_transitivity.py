#!/usr/bin/env python3
"""Verify the degree-42 cotangent-transitivity conormal obstruction."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_cellular_extension import CACHE  # noqa: E402
from research_degree42_tensor_extension import conormal_audit  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_conormal_transitivity.json"
)


def main() -> None:
    result = conormal_audit(base_order=2)
    assert result["normal_order"] == 1
    assert result["dimensions"] == {
        "kernel": 4,
        "total": 6,
        "spectator": 2,
    }
    assert result["splits_as_R_module"] is False
    assert result["zero_action_variables"] == [
        "n0",
        "n1",
        "n2",
        "n3",
        "n4",
        "n5",
        "n6",
    ]
    certificate = result["certificate"]
    assert all(
        not any(entry for row in matrix for entry in row)
        for variable, matrix in certificate["adapted_coupling"].items()
        if variable.startswith("n")
    )
    assert certificate["coboundary_rank"] == 5
    assert certificate["augmented_rank"] == 6
    assert certificate["obstruction_functional"] == {
        "zeta[0,0]": 2,
        "zeta[1,0]": 5,
    }
    assert certificate["obstruction_value"] == "2"

    output = {
        "schema": "degree42-ritt-conormal-transitivity.v1",
        "status": "exact cotangent Postnikov obstruction",
        "completed_base_ring": "B=Q[[tau,zeta]]=completed R/K",
        "conormal_projection": (
            "H_1(L_(B/A_6))=K/(I_6+K^2) -> "
            "H_1(L_(B/A_boundary))=K/(I_boundary+K^2)"
        ),
        "tensor_quotient": "B/(tau,zeta)^2",
        "source_ideal_cache": str(CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": hashlib.sha256(
            CACHE.read_bytes()
        ).hexdigest(),
        "calculation": result,
        "derived_consequence": (
            "the completed conormal projection has no B-linear section; "
            "therefore the relative cotangent transitivity triangle for "
            "A_6 -> A_boundary -> B is non-split and its connecting "
            "morphism is nonzero"
        ),
        "derived_boundary": (
            "this detects the first-Postnikov shadow of the connecting "
            "morphism; it does not compute the higher cotangent homology "
            "or prove the cellular homotopy-limit comparison"
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_conormal_transitivity.py"
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
    print("PASS: conormal quotient dimensions are 4 -> 6 -> 2")
    print("PASS: all seven normal actions vanish")
    print("PASS: the conormal projection has no base-linear section")
    print("PASS: obstruction rank rises from 5 to 6")
    print("PASS: the primitive obstruction functional evaluates to 2")
    print("PASS: the cotangent transitivity connecting morphism is nonzero")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the degree-42 first-Postnikov conormal exact sequence."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_cellular_extension import CACHE  # noqa: E402
from research_degree42_tensor_extension import (  # noqa: E402
    postnikov_overlap_audit,
)


CONORMAL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_conormal_transitivity.json"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_postnikov_overlap.json"
)


def main() -> None:
    conormal = json.loads(CONORMAL_ARTIFACT.read_text())
    assert conormal["schema"] == "degree42-ritt-conormal-transitivity.v1"
    assert conormal["calculation"]["dimensions"] == {
        "kernel": 4,
        "total": 6,
        "spectator": 2,
    }
    result = postnikov_overlap_audit()
    assert result["sector_source_mod_base_square"] == {
        "artin_rees_cutoff": 4,
        "cutoff_intersection_generators": 497,
        "cutoff_remainder_generators": 0,
        "denominator_quotient_length": 22,
        "numerator_quotient_length": 16,
        "dimension": 6,
    }
    assert result["quadratic_overlap_mod_base_square"] == {
        "numerator_generators": 70,
        "artin_rees_cutoff": 5,
        "cutoff_intersection_generators": 1286,
        "cutoff_remainder_generators": 0,
        "denominator_quotient_length": 38,
        "numerator_quotient_length": 38,
        "dimension": 0,
    }

    output = {
        "schema": "degree42-ritt-postnikov-overlap.v1",
        "status": "exact first-Postnikov conormal sequence",
        "completed_base_ring": "B=Q[[tau,zeta]]=completed R/K",
        "sector_module": (
            "S=H_1(B tensor^L_(A_boundary) "
            "L_(A_boundary/A_6))=I_boundary/(I_6+K*I_boundary)"
        ),
        "quadratic_overlap": (
            "(I_boundary intersect (I_6+K^2))/"
            "(I_6+K*I_boundary)"
        ),
        "source_ideal_cache": str(CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": hashlib.sha256(
            CACHE.read_bytes()
        ).hexdigest(),
        "conormal_artifact": str(CONORMAL_ARTIFACT.relative_to(ROOT)),
        "conormal_artifact_sha256": hashlib.sha256(
            CONORMAL_ARTIFACT.read_bytes()
        ).hexdigest(),
        "calculation": result,
        "formal_consequence": (
            "Nakayama's lemma applied to the zero base-square quotient "
            "of the quadratic overlap proves, after completion, "
            "I_boundary intersect (I_6+K^2) = "
            "I_6+K*I_boundary. Hence 0 -> S -> "
            "H_1(L_(B/A_6)) -> H_1(L_(B/A_boundary)) -> 0 "
            "is exact."
        ),
        "base_change_consequence": (
            "S tensor_B B/(tau,zeta)^2 has dimension 6, while its image "
            "in the conormal quotient has dimension 4; the two-dimensional "
            "kernel is ordinary non-flat base-change Tor, not an image of "
            "higher cotangent homology in the completed H_1 sequence"
        ),
        "derived_boundary": (
            "the individual H_i for i>=2 are not computed; only their "
            "possible image in the first transitivity homology sequence "
            "is proved to vanish"
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_postnikov_overlap.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: the sector conormal base-square quotient has dimension 6")
    print("PASS: the completed quadratic overlap vanishes by Nakayama")
    print("PASS: the first-Postnikov conormal sequence is short exact")
    print("PASS: the two-dimensional truncated kernel is base-change Tor")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

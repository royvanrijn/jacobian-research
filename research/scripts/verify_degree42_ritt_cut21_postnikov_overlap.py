#!/usr/bin/env python3
"""Verify completed first-Postnikov descent on the cut-21 half-braid."""

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
from research_degree42_ritt_rotated_postnikov_overlap import (  # noqa: E402
    rotated_postnikov_overlap_audit,
)


WORD = (3, 2, 7)
SOURCE_CACHE = source_ideal_cache(WORD)
JET_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_rotated_conormal_jet_327.json"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_cut21_postnikov_overlap.json"
)


def main() -> None:
    jet = json.loads(JET_ARTIFACT.read_text())
    assert jet["schema"] == "degree42-ritt-rotated-conormal-jet.v1"
    assert tuple(jet["word"]) == WORD
    assert jet["thick_composite_omission"] == 21
    assert jet["thin_prime_omission"] == 2
    assert jet["jet_lengths_q2"] == [5, 4, 4, 3]

    result = rotated_postnikov_overlap_audit(WORD)
    assert result["thin_boundary_nakayama"] == {
        "boundary_remainder_generators": 0,
        "tautological_reverse_remainder_generators": 0,
        "completed_equality": True,
    }
    assert result["sector_source_mod_base_square"] == {
        "artin_rees_cutoff": 4,
        "cutoff_intersection_generators": 495,
        "cutoff_remainder_generators": 0,
        "denominator_quotient_length": 19,
        "numerator_quotient_length": 16,
        "dimension": 3,
    }
    assert result["quadratic_overlap_mod_base_square"] == {
        "numerator_generators": 61,
        "artin_rees_cutoff": 5,
        "cutoff_intersection_generators": 1278,
        "cutoff_remainder_generators": 0,
        "denominator_quotient_length": 34,
        "numerator_quotient_length": 34,
        "dimension": 0,
    }

    output = {
        "schema": "degree42-ritt-cut21-postnikov-overlap.v1",
        "status": "exact completed first-Postnikov half-braid comparison",
        "word": WORD,
        "opposite_word": tuple(reversed(WORD)),
        "thick_composite_omission": 21,
        "thin_prime_omission": 2,
        "completed_base_ring": "B=Q[[tau,zeta]]=completed R/K",
        "source_ideal_cache": str(SOURCE_CACHE.relative_to(ROOT)),
        "source_ideal_cache_sha256": hashlib.sha256(
            SOURCE_CACHE.read_bytes()
        ).hexdigest(),
        "jet_artifact": str(JET_ARTIFACT.relative_to(ROOT)),
        "jet_artifact_sha256": hashlib.sha256(
            JET_ARTIFACT.read_bytes()
        ).hexdigest(),
        "calculation": result,
        "thin_boundary_consequence": (
            "I_boundary is contained in I_thin+m*I_boundary and I_thin is "
            "contained in I_boundary. Nakayama therefore proves equality "
            "of their completed ideals."
        ),
        "postnikov_consequence": (
            "The explicit cutoff-5 containment and equal finite lengths "
            "prove that the quadratic overlap has zero base-square "
            "quotient. Nakayama then gives "
            "I_boundary intersect (I_thick+K^2) = "
            "I_thick+K*I_boundary after completion, so the associated "
            "first-Postnikov conormal sequence is short exact."
        ),
        "base_change_consequence": (
            "The sector source modulo (tau,zeta)^2 has dimension 3, while "
            "the q2 quotient-length difference shows a one-dimensional "
            "image after finite base change. The two-dimensional loss is "
            "non-flat base-change Tor, not completed quadratic overlap."
        ),
        "derived_boundary": (
            "Together with the cut-6 and cut-14 certificates, this proves "
            "the completed comparison through H_1 on every degree-42 "
            "half-braid. It does not prove non-split extension transport, "
            "braid restriction coherence, or individual H_i for i>=2."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_cut21_postnikov_overlap.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: thin path and boundary agree after completion")
    print("PASS: the completed cut-21 quadratic overlap vanishes")
    print("PASS: the cut-21 first-Postnikov sequence is short exact")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

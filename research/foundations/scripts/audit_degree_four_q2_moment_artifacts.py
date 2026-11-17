#!/usr/bin/env python3
"""Audit committed quartic q2 normal-jet ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/research_degree_four_q2_augmented_nullcone.py":
        "6f49ca37e8e266a6e10e96b5dcb97ee98e0aa2b232eb5a6429e08d04134e4ec7",
    "scripts/research_degree_four_q2_cubic_decomposition.py":
        "1638531492d8ea905c2d4a6efe3e0de8d8cf970078881752e7f334ffbbc26f82",
    "artifacts/generated-results/degree_four_q2_augmented_nullcone_local.json":
        "9363d71f0d35f87b05c2ac370c98a0902f3bc9ce85402b0259261493efb8328b",
    "artifacts/generated-results/degree_four_q2_cubic_decomposition.json":
        "9c812573ba3dd0f014ffc4515a147605297045042ab0487fec9271f52c14d6ec",
}


def load(name: str) -> dict[str, object]:
    return json.loads(
        (ROOT / "artifacts" / "generated-results" / name).read_text(
            encoding="utf-8"
        )
    )


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (
            relative_path,
            actual_hash,
            expected_hash,
        )

    local = load("degree_four_q2_augmented_nullcone_local.json")
    assert local["format"] == "degree-four-q2-augmented-nullcone-local-v1"
    assert local["prime"] == 32003
    assert local["maximum_normal_jet_degree"] == 4
    assert local["formal_series_composition"] == "native"
    assert len(local["allowed_weight_coordinates"]) == 9
    assert len(local["forbidden_weight_coordinates"]) == 12
    assert local["linear_pivot_indices"] == [1, 2, 6, 11]
    assert len(local["free_normal_indices"]) == 8
    bases = local["local_standard_bases"]
    assert bases["2"] == {
        "status": "completed",
        "dimension": 6,
        "vector_space_dimension": -1,
        "standard_basis_size": 4,
    }
    assert bases["3"] == {
        "status": "completed",
        "dimension": 4,
        "vector_space_dimension": -1,
        "standard_basis_size": 85,
    }
    assert bases["4"] == {"status": "timeout", "timeout_seconds": 300}
    assert "no global q2-augmented nullcone equality" in local["scope_warning"]

    decomposition = load("degree_four_q2_cubic_decomposition.json")
    assert decomposition["format"] == "degree-four-q2-cubic-decomposition-v1"
    assert decomposition["prime"] == 32003
    assert decomposition["linear_pivot_indices"] == [1, 2, 6, 11]
    cubic = decomposition["cubic_decomposition"]
    assert cubic["cubic_ideal"] == {
        "standard_basis_size": 85,
        "dimension": 4,
        "multiplicity": 9,
        "x7_cubed_is_zero": True,
    }
    assert cubic["x6_zero_radical"]["containment_of_cubic_ideal"] is True
    assert cubic["x6_zero_radical"]["factor_count_over_finite_field"] == 2
    off_axis = cubic["off_x6_saturation"]
    assert off_axis["dimension"] == 3
    assert off_axis["multiplicity"] == 9
    assert off_axis["generic_fiber"]["degree"] == 9
    assert off_axis["generic_fiber"]["primitive_polynomial_factor_count"] == 1
    quartic = decomposition["dominant_component_quartic_restrictions"]
    assert quartic["common_quartic_radical_in_full_coordinates"] == [
        "x2",
        "x3",
        "x5",
        "x6",
        "x7",
    ]
    assert quartic["linear_sheet"]["ideal_contained_in_candidate_radical"] is True
    assert quartic["quadratic_sheet"]["ideal_contained_in_candidate_radical"] is True
    assert decomposition["off_axis_quartic_status"] == {
        "status": "timeout",
        "timeout_seconds": 240,
        "calculation": (
            "standard basis of the quartic moment ideal plus the cubic "
            "x6-saturation over F_32003(x1,x4,x6)"
        ),
        "inference": "none",
    }
    assert "does not establish global integrality" in decomposition["scope_warning"]

    local_source = (
        ROOT / "scripts/research_degree_four_q2_augmented_nullcone.py"
    ).read_text(encoding="utf-8")
    assert "except subprocess.TimeoutExpired:" in local_source
    assert 'return {"status": "timeout", "timeout_seconds": timeout}' in local_source
    decomposition_source = (
        ROOT / "scripts/research_degree_four_q2_cubic_decomposition.py"
    ).read_text(encoding="utf-8")
    assert '"off_axis_quartic_status": {' in decomposition_source
    assert '"inference": "none"' in decomposition_source

    print(
        "PASS: quartic q2 normal-jet sources and committed ledgers match "
        "pinned hashes, dimensions, radical flags, and timeout boundaries"
    )
    print(
        "BOUNDARY: all support calculations are at F_32003 and one normal "
        "slice; no formal isolation, F2=0 boundary, characteristic-zero "
        "lifting, or global integrality follows"
    )
    print(
        "PROVENANCE: the 240-second off-axis timeout is a retained historical "
        "record serialized by the cubic checker, not a branch rerun by it; "
        "its recorded inference is none"
    )
    print("NO RECOMPUTATION: no SymPy or Singular calculation was run")


if __name__ == "__main__":
    main()

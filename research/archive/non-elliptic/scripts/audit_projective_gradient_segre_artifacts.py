#!/usr/bin/env python3
"""Audit committed projective-gradient Segre ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/verify_projective_gradient_segre_machinery.py":
        "6ea559e6f7a7140439fce0d77bf4f99c757446794ed9f0360b7cb5997b4fef03",
    "scripts/verify_projective_gradient_normal_slices.py":
        "fbcb8ef2d4740c8880fa060d131ba085fa2927fac07cbfe128f7c9abe4784b45",
    "scripts/verify_projective_gradient_singular_slices.py":
        "72236c7d413fca60bff3324de1123ccc25a2fb704fa107c0627295716a37c199",
    "jcsearch/projective_gradient_segre.py":
        "53801ca63b1df72d9664336810d7130a9293f1413ffc37dbfcaaa54e6cce45ad",
    "scripts/verify_projective_polar_calibrations.m2":
        "128b72db2ffca47acc8c6e60bc0aae466494dd6711a3424409f98e941e868e3f",
    "scripts/verify_projective_gradient_segre_families.m2":
        "5131c0024240cc3a42511fc1b89a00f2d659a8caee5e8367dcdb88510226b89e",
    "scripts/verify_projective_gradient_normal_slices.m2":
        "f7cc174b68e8a304f622fa7c2e17a1502c86b79270cc0f4758159fa67f8f62ac",
    "scripts/verify_projective_gradient_singular_slices.m2":
        "4174aeaded3a0d7ac177680cc338dff7ba6241d4edadc11432bcfec4d6262a13",
    "artifacts/generated-results/projective_gradient_segre_registry.json":
        "1678eac19cc8e59a123ec84836f8f2a89f3b697a29c241e3b26e6987180fd00f",
    "artifacts/generated-results/projective_gradient_normal_slices.json":
        "5853c8fa609879663b31f680591a5e612ab944b1637902de5dcd115c9400837b",
    "artifacts/generated-results/projective_gradient_singular_slices.json":
        "c6971874b5359e4aed11a8918328804f9ffdd6e67811f49c9ff79b2a8c5d7b72",
}


def load(name: str) -> dict[str, object]:
    return json.loads(
        (ROOT / "artifacts" / "generated-results" / name).read_text(
            encoding="utf-8"
        )
    )


def audit_registry() -> None:
    registry = load("projective_gradient_segre_registry.json")
    assert registry["format"] == "projective-gradient-segre-registry-v3"
    assert registry["conventions"] == {
        "actual_affine_compactification": "[X0^m:F1^h:...:Fn^h]",
        "formula": "g_i=m^i-sum_{k=1}^i binom(i,k)m^(i-k)sigma_k",
        "full_polar_map_is_separate": True,
        "top_degree": (
            "g_n equals the affine generic degree when the affine map is dominant"
        ),
    }
    assert registry["all_dimension_checks"] == {
        "dimensions": [1, 12],
        "integrability_reconstruction": "exact Euler reconstruction",
        "map_degrees": [1, 7],
        "transform_pairs_checked": 84,
    }

    smooth = registry["smooth_essential_normal_slice_theorem"]
    assert smooth["artifact"] == (
        "artifacts/generated-results/projective_gradient_normal_slices.json"
    )
    assert smooth["checker"] == "scripts/verify_projective_gradient_normal_slices.py"
    assert smooth["independent_checker"] == (
        "scripts/verify_projective_gradient_normal_slices.m2"
    )
    singular = registry["singular_essential_normal_slice_theorem"]
    assert singular["artifact"] == (
        "artifacts/generated-results/projective_gradient_singular_slices.json"
    )
    assert singular["checker"] == (
        "scripts/verify_projective_gradient_singular_slices.py"
    )
    assert singular["independent_checker"] == (
        "scripts/verify_projective_gradient_singular_slices.m2"
    )

    complete = registry["complete_calibrations"]
    assert complete["triangular_r2"]["actual_affine_compactification"][
        "projective_degrees"
    ] == [1, 2, 2, 2, 1]
    assert complete["triangular_r2"]["full_polar_comparison"][
        "projective_degrees"
    ] == [1, 2, 4, 4, 2]
    assert complete["triangular_r3"]["actual_affine_compactification"][
        "projective_degrees"
    ] == [1, 3, 3, 3, 1]
    assert complete["triangular_r3"]["full_polar_comparison"][
        "projective_degrees"
    ] == [1, 3, 6, 6, 3]
    cotangent = registry["cotangent_calibrations"]
    assert cotangent["plane_triangular_r2"]["cotangent_gradient_lift"][
        "projective_degrees"
    ] == [1, 2, 3, 2, 1]
    assert cotangent["plane_triangular_r3"]["cotangent_gradient_lift"][
        "projective_degrees"
    ] == [1, 3, 5, 3, 1]
    stabilized = registry["quadratic_stabilization_calibrations"]
    assert stabilized["triangular_r2"]["after_one_quadratic_variable"][
        "projective_degrees"
    ] == [1, 2, 2, 2, 2, 1]
    assert stabilized["triangular_r3"]["after_one_quadratic_variable"][
        "projective_degrees"
    ] == [1, 3, 3, 3, 3, 1]

    controls = registry["top_degree_controls_and_open_records"]
    assert set(controls) == {
        "homogeneous_cotangent_hn_38",
        "meng_yang_doubled_hc6",
        "meng_yang_schur_hc5",
        "nonhomogeneous_cotangent_hn_40",
        "plane_quartic_packet_cotangent_target",
        "rank_reduced_cotangent_hn_44",
    }
    for name in (
        "homogeneous_cotangent_hn_38",
        "nonhomogeneous_cotangent_hn_40",
        "rank_reduced_cotangent_hn_44",
    ):
        assert controls[name]["affine_degree"] is None
        assert "not yet computed" in controls[name]["status"]
    assert controls["plane_quartic_packet_cotangent_target"]["status"] == (
        "conditional top-degree target; no explicit packet"
    )
    assert "top generic degree determines only the weighted aggregate" in (
        registry["scope"]
    )
    assert "not silently promoted to Segre classes" in registry["scope"]


def audit_smooth_slices() -> None:
    smooth = load("projective_gradient_normal_slices.json")
    assert smooth["format"] == "projective-gradient-normal-slices-v1"
    assert smooth["regression_range"] == {
        "all_essential_ranks": True,
        "ambient_dimensions": [2, 10],
        "map_degrees": [2, 7],
        "records_checked": 270,
    }
    records = smooth["records"]
    actual = {
        (row["ambient_dimension"], row["map_degree"], row["essential_rank"])
        for row in records
    }
    expected = {
        (dimension, degree, rank)
        for dimension in range(2, 11)
        for degree in range(2, 8)
        for rank in range(1, dimension)
    }
    assert len(records) == len(actual) == 270
    assert actual == expected
    assert smooth["specializations"]["HC4PPG7"] == {
        "jacobian_length": 64,
        "non_socle_affine_degree_lower_bound": 6,
        "parameters": {"m": 4, "n": 4, "r": 3},
        "truncated_active_length": 256,
    }
    assert smooth["specializations"]["HC4PPG8"] == {
        "jacobian_length": 16,
        "parameters": {"m": 4, "n": 4, "r": 2},
        "truncated_active_length": 64,
        "unit_penultimate_sigma3": 16,
    }
    assert "controls the first Segre multiplicity" in smooth["scope"]
    assert "does not determine later Segre degrees" in smooth["scope"]


def audit_singular_slices() -> None:
    singular = load("projective_gradient_singular_slices.json")
    assert singular["format"] == "projective-gradient-singular-slices-v1"
    assert singular["regression_range"] == {
        "all_essential_ranks_and_singular_dimensions": True,
        "ambient_dimensions": [3, 10],
        "map_degrees": [2, 7],
        "records_checked": 2160,
    }
    records = singular["records"]
    actual = {
        (
            row["ambient_dimension"],
            row["map_degree"],
            row["essential_rank"],
            row["singular_locus_dimension"],
            row["transverse_jacobian_length"],
        )
        for row in records
    }
    expected = {
        (dimension, degree, rank, singular_dimension, length)
        for dimension in range(3, 11)
        for degree in range(2, 8)
        for rank in range(2, dimension)
        for singular_dimension in range(rank - 1)
        for length in (1, 2, 4)
    }
    assert len(records) == len(actual) == 2160
    assert actual == expected
    profiles = singular["binary_quintic_calibration"]["profiles"]
    assert profiles == {
        "h4=0": {
            "generic_rank": 2,
            "torsion_orders": [],
            "truncated_length": 8,
        },
        "h4=x*y^3": {
            "generic_rank": 0,
            "torsion_orders": [2, 1],
            "truncated_length": 3,
        },
        "h4=y^4": {
            "generic_rank": 0,
            "torsion_orders": [1, 1],
            "truncated_length": 2,
        },
    }
    assert "profile (rho,a_j)" in singular["scope"]
    assert "No universal numerical Segre vector" in singular["scope"]


def audit_fail_closed_writes() -> None:
    for relative_path in (
        "scripts/verify_projective_gradient_segre_machinery.py",
        "scripts/verify_projective_gradient_normal_slices.py",
        "scripts/verify_projective_gradient_singular_slices.py",
    ):
        source = (ROOT / relative_path).read_text(encoding="utf-8")
        assert '"--write"' in source
        assert "if args.write:" in source
        assert "assert OUTPUT.read_text() == serialized" in source
        assert "is stale; regenerate with --write" in source


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (relative_path, actual_hash, expected_hash)
    audit_registry()
    audit_smooth_slices()
    audit_singular_slices()
    audit_fail_closed_writes()
    print(
        "PASS: projective-gradient Segre sources, Macaulay2 calibrations, and "
        "three committed ledgers match pinned hashes and exact coverage"
    )
    print(
        "BOUNDARY: complete vectors, top-degree-only controls, uncomputed "
        "families, generic first-Segre laws, and singular DVR profiles remain "
        "separate evidence types"
    )
    print(
        "WRITE SAFETY: ordinary checker runs compare committed bytes; artifact "
        "replacement now requires explicit --write"
    )
    print("NO RECOMPUTATION: no SymPy or Macaulay2 calculation was run")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit the committed HC4 Fitting-denominator ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py":
        "a57067701721192cdc0b538ebdf850131988ee23e30a9bfd87fb0de7b97788da",
    "artifacts/generated-results/hc4_fitting_denominator_extraction.json":
        "d2ddab3389244192e91a76958e7ea3bf1f48704afa433a89ba73b5294a4a8531",
    "artifacts/generated-results/hc4_fourth_power_support.json":
        "576e3c0f3f0bd686316caf766e554a471866791a33d00e8f5a97762c44508946",
}


def read_json(relative_path: str) -> dict[str, object]:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (
            relative_path,
            actual_hash,
            expected_hash,
        )

    fitting = read_json(
        "artifacts/generated-results/hc4_fitting_denominator_extraction.json"
    )
    assert fitting["base_ring"] == "Q[mu,nu]"
    assert fitting["canonical_module"] == {
        "quadratic_generator_count": 114,
        "quadratic_target_rank": 120,
        "degree_three_source_rank": 1710,
        "degree_three_target_rank": 680,
        "cube_generator_count": 15,
        "definition": (
            "T=(im(Phi)+im(D))/im(Phi), where Phi(e_j tensor s_i)="
            "s_i*f_j and D(e_i)=s_i^3"
        ),
    }
    scans = fitting["finite_field_scans"]
    assert scans["primes"] == [11, 13, 17, 19]
    assert scans["reconstruction"] == {
        "radial": ["1/5", "1/10"],
        "additional": ["-5/3", "-1/6"],
    }
    assert "interpret the 114-equation module on D(nu)" in scans[
        "chart_boundary_warning"
    ]
    mixed = fitting["exact_specializations"]["mixed"]
    assert mixed["nonzero_cube_indices"] == [2, 9, 11]
    assert mixed["all_fourth_powers_zero"] is True
    assert mixed["reduced_fiber_is_origin"] is True
    assert mixed["lower_face_status"].startswith("not run")
    assert fitting["timeouts_seconds"] == {
        "integral_even_block_annihilator": 900,
        "function_field_even_block_lift": 900,
    }
    assert "exact zeroth Fitting ideal and associated-prime equality remain open" in (
        fitting["interpretation"]
    )

    fourth = read_json(
        "artifacts/generated-results/hc4_fourth_power_support.json"
    )
    assert fourth["schema"] == "hc4_fourth_power_support/v1"
    assert fourth["status"] == (
        "exact finite-field scans; not a characteristic-zero support theorem"
    )
    assert fourth["chart"] == "nu != 0"
    symbolic = fourth["symbolic_annihilator_attempt"]
    assert symbolic["result"] == "timeout_before_standard_basis"
    assert symbolic["mathematical_conclusion"] == "none"
    expected_radial = {
        "7": [3, 5],
        "11": [9, 10],
        "13": [8, 4],
    }
    for prime_text, radial in expected_radial.items():
        prime = int(prime_text)
        row = fourth["scans"][prime_text]
        assert row["parameter_points_on_D_nu"] == prime * (prime - 1)
        assert row["certified_empty_points"] == prime * (prime - 1) - 1
        assert row["radial_reduction"] == radial
        assert len(row["exceptional_points"]) == 1
        assert row["exceptional_points"][0]["parameters"] == radial
        assert row["mixed_nilpotence_point_certified_empty"] is True
    limitations = " ".join(fourth["limitations"])
    assert "proper extensions" in limitations
    assert "characteristic-zero" in limitations
    assert "valid only on D(nu)" in limitations

    print(
        "PASS: committed HC4 fitting-denominator and fourth-power ledgers, "
        "their source, and their fail-closed scope match exact hashes"
    )
    print(
        "SCOPE: cube torsion detects a nilpotence jump, not a reduced Schur "
        "component; finite-field support is only on D(nu), and the integral "
        "Fitting ideal and associated primes remain open"
    )
    print("NO RECOMPUTATION: no SymPy import, Singular run, scan, or rewrite")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit committed completed-moment ledgers without recomputing them."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/research_completed_moment_algebra.py":
        "eedc52d34a799c9ae3e164283bc8ff111b57c38d8ff2f509ab708962afc55431",
    "scripts/verify_completed_moment_diagonal_fields.py":
        "94c619f21f923257aaa74fa058f38daa3bf1a84c62c2df18cd884ac9282e9b9e",
    "scripts/verify_completed_moment_single_phase_fields.py":
        "7300a756e630836a79e94a6b87e3661077f2d1c18ba959c0ef0dda4a3607b244",
    "scripts/verify_degree_four_diagonal_moment_field.py":
        "5c948d9ef4189375c761acddf9b58e675d953a525fdd6ad42c93542ad4fb6643",
    "artifacts/generated-results/completed_moment_algebra_bounded_tests.json":
        "590fe262178bc4e8f11f3b633be9649ae82050afab962768e6e7946a1b15aa7c",
    "artifacts/generated-results/automatic_missing_invariants_d3_d6.json":
        "ab13e38ec344194a60fda4ecebfdfe7d6b9c78f8d0be6b7cad1791e9fdb2d824",
    "artifacts/generated-results/completed_moment_diagonal_fields.json":
        "17fdcbd88f261c7aff209e86207d1a3d9170fd5099db592e66cb004da063ff10",
    "artifacts/generated-results/completed_moment_single_phase_fields.json":
        "3626114820cf156c16327827ec9d8659c216964bd74401c2ea1ef8cc647b8f1e",
}
LEGACY_SOURCE_REVISION = "c91498bbca857b568a9961e776416fcba8de6713"
LEGACY_SOURCE_SHA256 = (
    "4ddc8bcaf025ef7992f6ac9a7c6e32ef6974372ac8475bf633f9b6a35a30575c"
)


def load(relative_path: str) -> dict[str, object]:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (
            relative_path,
            actual_hash,
            expected_hash,
        )

    legacy = load(
        "artifacts/generated-results/completed_moment_algebra_bounded_tests.json"
    )
    assert legacy["status"] == (
        "bounded modular nonrelation tests; positive intersections would "
        "require characteristic-zero reconstruction"
    )
    assert legacy["prime"] == 1_000_003
    assert legacy["max_weight"] == 10
    assert legacy["extra_samples"] == 3
    assert set(legacy) == {"status", "prime", "max_weight", "extra_samples", "degrees"}
    assert set(legacy["degrees"]) == {"3", "4", "5"}
    for degree, expected_count in {"3": 57, "4": 83, "5": 101}.items():
        row = legacy["degrees"][degree]
        assert row["moment_jacobian_rank"] == row["quotient_dimension"]
        assert len(row["relation_tests"]) == expected_count
        assert all(
            test["relation_intersection"] == 0
            for test in row["relation_tests"]
        )

    automatic = load(
        "artifacts/generated-results/automatic_missing_invariants_d3_d6.json"
    )
    assert automatic["prime"] == 1_000_003
    assert automatic["max_weight"] == 10
    assert automatic["invariant_cutoff"] == 6
    assert automatic["extra_samples"] == 3
    assert automatic["relation_tests_skipped"] is True
    assert len(automatic["power_witness_scan"]) == 12
    assert len(automatic["all_degree_ladder_regression"]) == 32
    assert set(automatic["degrees"]) == {"3", "4", "5", "6"}
    expected_dimensions = {"3": 13, "4": 22, "5": 33, "6": 46}
    for degree, quotient_dimension in expected_dimensions.items():
        row = automatic["degrees"][degree]
        assert row["quotient_dimension"] == quotient_dimension
        assert row["moment_jacobian_rank"] == quotient_dimension
        assert not row["relation_tests"]
        scan = row["automatic_missing_invariant_scan"]
        assert scan["first_missing_invariant_degree"] == 2
        assert scan["first_missing_dimension"] == int(degree) - 1

    diagonal = load(
        "artifacts/generated-results/completed_moment_diagonal_fields.json"
    )
    assert diagonal["format"] == "completed-moment-diagonal-fields-v1"
    assert set(diagonal["degrees"]) == {"3", "4", "5"}
    for degree, expected_length in {"3": 24, "4": 120, "5": 720}.items():
        row = diagonal["degrees"][degree]
        assert row["finite_field_parameter_quotient_length"] == expected_length
        assert row["characteristic_zero_parameter_quotient_length"] == expected_length
        assert row["fiber_length"] == 2
        assert row["mutual_ideal_reductions_zero"] is True
        assert row["scope_warning"].startswith("slice theorem only")

    single_phase = load(
        "artifacts/generated-results/completed_moment_single_phase_fields.json"
    )
    assert single_phase["format"] == "completed-moment-single-phase-fields-v1"
    assert set(single_phase["degrees"]) == {"3", "5"}
    for degree, phase_count, origin_length in (("3", 3, 54), ("5", 5, 1934)):
        row = single_phase["degrees"][degree]
        assert row["certified_phase_count"] == phase_count
        assert len(row["phases"]) == phase_count
        assert "no determination" in row["scope_warning"]
        for phase in row["phases"]:
            assert phase["special_first_d_plus_3_origin_quotient_length"] == origin_length
            assert phase["characteristic_zero_fiber_length"] == 2
        if degree == "5":
            odd_values = {
                phase["phase"]: phase["first_apolar_odd_invariant_at_lift_b_1_c_z"]
                for phase in row["phases"]
            }
            assert odd_values == {
                1: "-273686400/7",
                2: "-273686400/7",
                3: "0",
                4: "0",
                5: "0",
            }

    current_source = (ROOT / "scripts/research_completed_moment_algebra.py").read_text(
        encoding="utf-8"
    )
    assert "--invariant-cutoff" in current_source
    assert "--ladder-beta-check" in current_source

    print(
        "PASS: committed completed-moment artifacts, direct sources, exact "
        "slice lengths, and stored bounded-search flags match pinned hashes"
    )
    print(
        "BOUNDARY: fixed-field conclusions are slice theorems; Hilbert tests "
        "are necessary only; skipped or positive relation tests prove nothing"
    )
    print(
        "PROVENANCE: completed_moment_algebra_bounded_tests.json is the legacy "
        f"{LEGACY_SOURCE_REVISION[:12]} run (source sha256 {LEGACY_SOURCE_SHA256}); "
        "the current source emits the extended automatic-ledger schema"
    )
    print("NO RECOMPUTATION: no SymPy invariant scan or Singular process was run")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit committed degree-four moment-field ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/verify_two_pair_counterexample_missing_invariant.py":
        "68953889d25c254a94ae726864f20b13487e058788d60ab72d9164e6d2e6c568",
    "scripts/verify_degree_four_tau_even_parameters.py":
        "412c76e9cbddc301922161e44de9407fc1e08646ce4f14f4781bbc960a66eb84",
    "scripts/research_degree_four_moment_field.py":
        "e45a377ef487d3c9e2ccc5b512a7dd8e5ff15d2aa8c0fef4bac1c549575e842c",
    "scripts/verify_degree_four_diagonal_moment_field.py":
        "5c948d9ef4189375c761acddf9b58e675d953a525fdd6ad42c93542ad4fb6643",
    "scripts/verify_degree_four_single_phase_moment_fields.py":
        "b36fe321d4d212dae3b4cc42cf61d9af6b6723ae4d1b2830e1571b7228c5bd69",
    "scripts/research_degree_four_phase_one_chart.py":
        "14b5e68046a6d41911cc49b07b786e0c538e2b403ea4bfe31b4e7d6d70b2dbd7",
    "artifacts/generated-results/two_pair_counterexample_missing_invariant.json":
        "236a927d099e9ac14dad9f36d57d38fd1f8a1e4a11cf0843ce4759474f1d5c34",
    "artifacts/generated-results/degree_four_tau_even_parameters.json":
        "8677dd7457d292d8c0103648397ffdf2803a77914557157823ec4da67e4b4c7c",
    "artifacts/generated-results/degree_four_moment_field_bounded_relations.json":
        "85f61583487da90a5886eac6ce6845ca8c917232f42d6ab75658eacdbdbe6506",
    "artifacts/generated-results/degree_four_diagonal_moment_field.json":
        "fb3854ebb8c7758c94388bae9bc56c967b4fa38e0dc8d7e0cf77ab2c8d872df9",
    "artifacts/generated-results/degree_four_single_phase_moment_fields.json":
        "b6cf6f2fe296f5fb7974c8d3fb5f535952b9c25b3bda040ddafd2c4c8619aadf",
    "artifacts/generated-results/degree_four_phase_one_chart_modular.json":
        "4b9acaa336b8d70b8b1e7581f73b09457566b7df286a2efc363ac06fdb9a74f6",
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

    base = load("two_pair_counterexample_missing_invariant.json")
    assert base["quadratic_invariants"]["values_at_F"] == {
        "q_0": "0",
        "q_2": "-864",
        "q_4": "2016",
        "q_6": "0",
        "q_8": "0",
    }
    assert base["moment_jacobian_certificate"]["rank"] == 22
    adjoint = base["apolar_adjoint"]
    assert adjoint["first_odd_invariant"] == "c_234 in degree 3"
    assert adjoint["conductor"] == "zero"
    assert adjoint["degree_of_d2_first_six_moment_parameter_map"] == 10
    assert adjoint["d2_mu7_not_in_first_six_parameter_ring_certificate"][
        "evaluation_rank_with_mu7"
    ] == 15
    assert adjoint["degree_of_d2_full_moment_field"] == 2

    parameters = load("degree_four_tau_even_parameters.json")
    assert parameters["format"] == "degree-four-tau-even-parameters-v1"
    assert parameters["prime"] == 1_000_003
    assert parameters["parameter_count"] == 22
    assert parameters["jacobian_rank"] == 22
    assert parameters["minor_determinant_mod_prime"] == 531404
    assert parameters["first_22_moment_jacobian"]["rank"] == 22
    assert parameters["combined_even_parameter_and_moment_cotangent_rank"] == 22
    assert "no assertion" in parameters["status"]

    bounded = load("degree_four_moment_field_bounded_relations.json")
    assert bounded["format"] == "degree-four-moment-field-bounded-relations-v1"
    assert bounded["prime"] == 1_000_003
    assert bounded["max_weight"] == 16
    assert bounded["extra_samples"] == 3
    assert bounded["targets"] == "odd-square"
    rows = bounded["results"]["c_234^2"]
    assert [row["weight"] for row in rows] == list(range(6, 17))
    assert all(row["rank"] == row["columns"] and row["nullity"] == 0 for row in rows)
    assert "rank defect would require" in bounded["status"]

    diagonal = load("degree_four_diagonal_moment_field.json")
    assert diagonal["format"] == "degree-four-diagonal-moment-field-v1"
    assert diagonal["first_five_parameter_quotient_length"] == 120
    assert diagonal["fiber_length"] == 2
    assert diagonal["mutual_ideal_reductions_zero"] is True
    assert diagonal["generic_full_moment_degree_on_slice"] == 2
    assert diagonal["scope_warning"].startswith("slice theorem only")

    single_phase = load("degree_four_single_phase_moment_fields.json")
    assert single_phase["format"] == "degree-four-single-phase-moment-fields-v1"
    assert single_phase["slice_count"] == 10
    assert len(single_phase["slices"]) == 10
    assert single_phase["fixed_locus_coordinate_slice_count"] == 6
    moving = single_phase["genuinely_apolar_moving_coordinate_slices"]
    assert len(moving) == 4
    assert {(row["phase"], row["positive_direction_index"], row["negative_direction_index"]) for row in moving} == {
        (1, 0, 1),
        (1, 1, 0),
        (2, 0, 1),
        (2, 1, 0),
    }
    assert "no determination" in single_phase["scope_warning"]

    phase_one = load("degree_four_phase_one_chart_modular.json")
    assert phase_one["format"] == "degree-four-phase-one-chart-modular-v3"
    assert phase_one["prime"] == 101
    assert phase_one["fiber_length"] == 4
    assert phase_one["fiber_reduced"] is True
    assert phase_one["characteristic_zero_orbit_identity_verified"] is True
    assert phase_one["characteristic_zero_reconstructed_points_reduced_and_isolated"] is True
    assert phase_one["characteristic_zero_odd_cubic_values"] == [1728, -1728]
    assert "completeness is proved only over F_101" in phase_one["scope_warning"]
    assert "generic degree remain open" in phase_one["scope_warning"]

    phase_source = (ROOT / "scripts/research_degree_four_phase_one_chart.py").read_text(
        encoding="utf-8"
    )
    assert '"--certify-example"' in phase_source
    assert 'action="store_true"' in phase_source
    assert "assert completed.returncode == 0" in phase_source
    assert "OUTPUT.write_text" in phase_source

    print(
        "PASS: degree-four moment-field sources and committed ledgers match "
        "pinned hashes and their exact stored certificate flags"
    )
    print(
        "BOUNDARY: diagonal/fixed slices are not full-quotient degree tests; "
        "weight-16 nonrelations are support-bounded; F_101 completeness does "
        "not exclude extra characteristic-zero branches"
    )
    print(
        "FAIL-CLOSED: the phase-one pinned certificate asserts solver success "
        "before writing its artifact; timeouts and rank defects prove nothing"
    )
    print("NO RECOMPUTATION: no SymPy, Singular, or msolve calculation was run")


if __name__ == "__main__":
    main()

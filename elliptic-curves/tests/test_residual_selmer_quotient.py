from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "residual_selmer_quotient.py"
spec = importlib.util.spec_from_file_location("residual_selmer_quotient", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)

RELATIVE_MATRIX_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "cas"
    / "build_mw29_relative_2selmer_matrix.py"
)
relative_spec = importlib.util.spec_from_file_location(
    "build_mw29_relative_2selmer_matrix", RELATIVE_MATRIX_SCRIPT
)
relative_module = importlib.util.module_from_spec(relative_spec)
sys.modules[relative_spec.name] = relative_module
assert relative_spec.loader is not None
relative_spec.loader.exec_module(relative_module)

WITNESS_BOUND_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "cas"
    / "audit_mw29_relative_selmer_witness_bound.py"
)
witness_spec = importlib.util.spec_from_file_location(
    "audit_mw29_relative_selmer_witness_bound", WITNESS_BOUND_SCRIPT
)
witness_module = importlib.util.module_from_spec(witness_spec)
sys.modules[witness_spec.name] = witness_module
assert witness_spec.loader is not None
witness_spec.loader.exec_module(witness_module)

PAIRING_SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "audit_residual_cassels_tate.py"
LOCAL_FILTER_SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "filter_bnf_free_local_selmer.py"
pairing_spec = importlib.util.spec_from_file_location("audit_residual_cassels_tate", PAIRING_SCRIPT)
pairing_module = importlib.util.module_from_spec(pairing_spec)
sys.modules[pairing_spec.name] = pairing_module
assert pairing_spec.loader is not None
pairing_spec.loader.exec_module(pairing_module)


class ResidualSelmerQuotientTests(unittest.TestCase):
    def test_dimension_only_obstruction_witnesses_can_close_relative_bound(self):
        result = witness_module.audit(
            {
                "schema": witness_module.INPUT_SCHEMA,
                "case_id": "toy",
                "global_ambient_dimension_upper_bound": 33,
                "known_mw_dimension": 29,
                "condition_blocks": [
                    {"place": "norm", "width": 1},
                    {"place": "2", "width": 2},
                    {"place": "5", "width": 1},
                ],
                "witnesses": [
                    {"label": f"a{index}", "generator": f"alpha{index}",
                     "condition_syndrome": [int(column == index) for column in range(4)]}
                    for index in range(4)
                ],
                "certification": {
                    "method": "certified class-quotient upper bound plus exact local maps",
                    "hypothesis": None,
                    "global_dimension_upper_bound_certified": True,
                    "witnesses_are_exact_global_squareclasses": True,
                    "witnesses_lie_in_global_ambient_certified": True,
                    "condition_syndromes_certified": True,
                    "condition_blocks_are_necessary_selmer_conditions": True,
                    "known_mw_in_condition_kernel_certified": True,
                },
            },
            maximum_cut_size=3,
        )
        self.assertEqual(result["certified_condition_image_rank"], 4)
        self.assertEqual(result["relative_selmer_dimension_upper_bound"], 0)
        self.assertEqual(result["status"], "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO")
        self.assertEqual(result["greedy_place_order"][0]["place"], "2")
        self.assertEqual(result["minimum_closing_place_cut"]["size"], 3)

    def test_dimension_witness_bound_uses_certified_even_parity(self):
        result = witness_module.audit(
            {
                "schema": witness_module.INPUT_SCHEMA,
                "global_ambient_dimension_upper_bound": 31,
                "known_mw_dimension": 29,
                "residual_selmer_dimension_parity": 0,
                "condition_blocks": [{"place": "2", "width": 1}],
                "witnesses": [
                    {"label": "a", "generator": "alpha", "condition_syndrome": [1]}
                ],
                "certification": {
                    "method": "exact dimension bound, local map, and parity",
                    "global_dimension_upper_bound_certified": True,
                    "witnesses_are_exact_global_squareclasses": True,
                    "witnesses_lie_in_global_ambient_certified": True,
                    "condition_syndromes_certified": True,
                    "condition_blocks_are_necessary_selmer_conditions": True,
                    "known_mw_in_condition_kernel_certified": True,
                    "residual_dimension_parity_certified": True,
                },
            }
        )
        self.assertEqual(result["relative_selmer_dimension_upper_bound_raw"], 1)
        self.assertEqual(result["relative_selmer_dimension_upper_bound"], 0)
        self.assertEqual(result["minimum_closing_place_cut"]["places"], ["2"])

    def test_dimension_only_witness_bound_is_fail_closed(self):
        manifest = {
            "schema": witness_module.INPUT_SCHEMA,
            "global_ambient_dimension_upper_bound": 30,
            "known_mw_dimension": 29,
            "condition_blocks": [{"place": "2", "width": 1}],
            "witnesses": [
                {"label": "a", "generator": "alpha", "condition_syndrome": [1]}
            ],
            "certification": {
                "method": "uncertified relation plateau",
                "global_dimension_upper_bound_certified": False,
                "witnesses_are_exact_global_squareclasses": True,
                "witnesses_lie_in_global_ambient_certified": True,
                "condition_syndromes_certified": True,
                "condition_blocks_are_necessary_selmer_conditions": True,
                "known_mw_in_condition_kernel_certified": True,
            },
        }
        result = witness_module.audit(manifest)
        self.assertEqual(result["status"], "INCOMPLETE_WITNESS_BOUND")
        self.assertIsNone(result["relative_selmer_dimension_upper_bound"])
        self.assertIsNone(
            result["greedy_place_order"][0]["raw_relative_selmer_dimension_upper_bound"]
        )

    def test_binary_linear_algebra_rejects_nonbinary_entries(self):
        with self.assertRaisesRegex(module.F2Error, "non-binary"):
            module.f2_row_basis([[2]], 1)

    def test_nullspace_basis_is_exact(self):
        rows = [[1, 1, 0, 1], [0, 1, 1, 0]]
        nullspace = module.f2_nullspace_basis(rows, 4)
        self.assertEqual(len(nullspace), 2)
        self.assertEqual(module.f2_rank_rows(nullspace, 4), 2)
        for vector in nullspace:
            self.assertTrue(all(module.f2_dot(row, vector) == 0 for row in rows))

    def test_relative_matrix_quotients_mw_before_local_conditions(self):
        result = module.build_relative_local_condition_matrix(
            ambient_dimension=5,
            known_mw_rows=[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]],
            places=[
                {
                    "place": "2",
                    "allowed_subspace_basis": [
                        [1, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0],
                        [0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 1],
                    ],
                },
                {
                    "place": "3",
                    "allowed_subspace_basis": [
                        [1, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0],
                        [0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 1],
                    ],
                },
                {
                    "place": "5",
                    "allowed_subspace_basis": [
                        [1, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0],
                        [0, 0, 1, 0, 0],
                        [0, 0, 0, 1, 0],
                    ],
                },
                {
                    "place": "7",
                    "allowed_subspace_basis": [
                        [1, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0],
                        [0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 1],
                    ],
                },
            ],
            maximum_cut_size=4,
        )
        self.assertEqual(result["known_mw_kummer_dimension"], 2)
        self.assertEqual(result["mw_quotient_ambient_dimension"], 3)
        self.assertEqual(result["full_relative_local_condition_matrix_rank"], 3)
        self.assertEqual(result["unexplained_selmer_excess_kernel_dimension"], 0)
        self.assertEqual(
            [row["place"] for row in result["greedy_place_order"]],
            ["2", "3", "5", "7"],
        )
        self.assertEqual(
            result["minimum_annihilating_place_cut"],
            {"size": 3, "places": ["2", "3", "5"], "minimality_proved": True},
        )
        delete = {
            row["deleted_place"]: row for row in result["delete_one_place_ranks"]
        }
        self.assertEqual(delete["2"]["rank_drop"], 0)
        self.assertEqual(delete["3"]["rank_drop"], 1)

    def test_relative_matrix_rejects_known_point_outside_local_image(self):
        with self.assertRaisesRegex(module.F2Error, "violates place 2"):
            module.build_relative_local_condition_matrix(
                ambient_dimension=2,
                known_mw_rows=[[1, 0]],
                places=[{"place": "2", "allowed_subspace_basis": [[0, 1]]}],
            )

    def test_relative_proof_gate_is_fail_closed(self):
        manifest = {
            "schema": relative_module.INPUT_SCHEMA,
            "case_id": "toy",
            "ambient_norm_square_dimension": 2,
            "known_mw_target_rank": 1,
            "known_mw_rows": [{"label": "P", "row": [1, 0]}],
            "places": [
                {"place": "2", "allowed_subspace_basis": [[1, 0]]}
            ],
            "certification": {
                "method": "incomplete relation collection",
                "global_ambient_upper_envelope_certified": False,
                "global_ambient_exact": True,
                "norm_condition_incorporated": True,
                "known_mw_kummer_coordinates_certified": True,
                "supplied_local_conditions_certified": True,
                "supplied_subspaces_are_necessary_selmer_conditions": True,
                "all_required_local_conditions_complete": True,
                "residual_dimension_parity_certified": False,
            },
        }
        incomplete = relative_module.build_certificate(manifest, maximum_cut_size=2)
        self.assertEqual(incomplete["status"], "INCOMPLETE_RELATIVE_2SELMER_MATRIX")
        self.assertTrue(
            incomplete["relative_local_matrix"]["all_residual_candidates_annihilated"]
        )
        manifest["certification"]["global_ambient_upper_envelope_certified"] = True
        complete = relative_module.build_certificate(manifest, maximum_cut_size=2)
        self.assertEqual(
            complete["status"], "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
        )

    def test_certified_even_parity_closes_a_one_dimensional_upper_bound(self):
        result = relative_module.build_certificate(
            {
                "schema": relative_module.INPUT_SCHEMA,
                "case_id": "parity-toy",
                "ambient_norm_square_dimension": 3,
                "known_mw_target_rank": 2,
                "residual_selmer_dimension_parity": 0,
                "known_mw_rows": [
                    {"label": "P1", "row": [1, 0, 0]},
                    {"label": "P2", "row": [0, 1, 0]},
                ],
                "places": [],
                "certification": {
                    "method": "exact global upper bound and parity theorem",
                    "global_ambient_upper_envelope_certified": True,
                    "global_ambient_exact": False,
                    "norm_condition_incorporated": True,
                    "known_mw_kummer_coordinates_certified": True,
                    "supplied_local_conditions_certified": True,
                    "supplied_subspaces_are_necessary_selmer_conditions": True,
                    "all_required_local_conditions_complete": False,
                    "residual_dimension_parity_certified": True,
                },
            },
            maximum_cut_size=0,
        )
        self.assertEqual(
            result["relative_selmer_bound"],
            {
                "raw_upper_bound_from_supplied_local_matrix": 1,
                "certified_residual_dimension_parity": 0,
                "parity_sharpened_upper_bound": 0,
                "global_presentation_exact": False,
                "all_required_local_conditions_complete": False,
                "residual_kernel_exact": False,
            },
        )
        self.assertEqual(
            result["status"], "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
        )

    def test_nonzero_kernel_needs_an_exact_global_presentation(self):
        manifest = {
            "schema": relative_module.INPUT_SCHEMA,
            "ambient_norm_square_dimension": 2,
            "known_mw_target_rank": 1,
            "known_mw_rows": [{"label": "P", "row": [1, 0]}],
            "places": [],
            "certification": {
                "method": "certified dimension envelope only",
                "global_ambient_upper_envelope_certified": True,
                "global_ambient_exact": False,
                "norm_condition_incorporated": True,
                "known_mw_kummer_coordinates_certified": True,
                "supplied_local_conditions_certified": True,
                "supplied_subspaces_are_necessary_selmer_conditions": True,
                "all_required_local_conditions_complete": True,
                "residual_dimension_parity_certified": False,
            },
        }
        result = relative_module.build_certificate(manifest, maximum_cut_size=0)
        self.assertEqual(result["status"], "CERTIFIED_RELATIVE_2SELMER_UPPER_BOUND")
        self.assertFalse(result["relative_selmer_bound"]["residual_kernel_exact"])

    def test_sparse_dependency_retains_actual_generators(self):
        relations = module.SparseF2Relations(ideal_dimension=4)
        self.assertIsNone(relations.add(module.PrincipalRelation("r1", "1 + theta", 0b0011)))
        self.assertIsNone(relations.add(module.PrincipalRelation("r2", "2 - theta", 0b0101)))
        dependency = relations.add(module.PrincipalRelation("r3", "theta^2", 0b0110))
        self.assertIsNotNone(dependency)
        assert dependency is not None
        self.assertEqual(dependency.relation_labels, ("r1", "r2", "r3"))
        self.assertEqual(dependency.generator_product, ("1 + theta", "2 - theta", "theta^2"))

    def test_early_quotient_distinguishes_only_new_signature(self):
        mw = module.SquareclassImage("P1", "m_P1", local=0b01, fingerprint=0b101)
        quotient = module.EarlyQuotient(
            local_dimension=2, fingerprint_dimension=3, known_mw_images=(mw,)
        )
        same = quotient.image(module.SquareclassImage("same", "m_same", 0b01, 0b101))
        new = quotient.image(module.SquareclassImage("new", "m_new", 0b11, 0b101))
        self.assertEqual(quotient.known_mw_rank, 1)
        self.assertTrue(same.killed_by_known_mw)
        self.assertEqual(new.residual_signature, 0b10)

    def test_certification_refuses_rank_stabilization_as_a_proof(self):
        certification = module.ClassQuotientCertification("rank stabilization", 0)
        with self.assertRaises(module.F2Error):
            module.certification_record(certification)

    def test_hypothesis_is_explicit_in_certification_status(self):
        certification = module.ClassQuotientCertification(
            "GRH-certified 2-class computation", 0, hypothesis="GRH"
        )
        record = module.certification_record(certification)
        self.assertEqual(record["status"], "CERTIFIED_UNDER_HYPOTHESIS")
        self.assertEqual(record["remaining_mod2_s_class_dimension_upper_bound"], 0)

    def test_manifest_audit_preserves_generators_and_residuals(self):
        output = module.audit_manifest(
            {
                "local_dimension": 2,
                "fingerprint_dimension": 2,
                "known_mw_images": [
                    {"label": "P", "generator": "delta(P)", "local": "0b01", "fingerprint": 2}
                ],
                "candidate_images": [
                    {"label": "c", "generator": "alpha", "local": 3, "fingerprint": 2}
                ],
                "class_quotient_certification": {
                    "method": "analytic class-number bound",
                    "remaining_dimension_upper_bound": 1,
                },
            }
        )
        self.assertEqual(output["known_mw_target_rank"], 1)
        candidate = output["candidate_images"][0]
        self.assertEqual(candidate["generator"], "alpha")
        self.assertEqual(candidate["residual_support"], (1,))
        self.assertEqual(output["candidate_residual_rank"], 1)
        self.assertEqual(output["independent_candidate_labels"], ("c",))

    def test_nondegenerate_certified_pairing_closes_residual_rank(self):
        result = pairing_module.audit(
            {
                "known_mw_rank": 30,
                "residual_selmer_basis_certified": True,
                "pairing_algorithm": "explicit Cassels--Tate algorithm",
                "pairing_computed_before_cover_search": True,
                "residual_basis": [{"label": "c1"}, {"label": "c2"}],
                "cassels_tate_matrix": [[0, 1], [1, 0]],
                "cover_searches": [],
            }
        )
        self.assertEqual(result["classification"], "CERTIFIED_EXACT_KNOWN_RANK_AFTER_PAIRING")
        self.assertEqual(result["radical_dimension"], 0)
        self.assertEqual(result["rank_upper_after_pairing"], 30)

    def test_pairing_audit_rejects_pre_pairing_cover_search(self):
        with self.assertRaises(pairing_module.PairingError):
            pairing_module.audit(
                {
                    "known_mw_rank": 20,
                    "residual_basis": [{"label": "c"}],
                    "cassels_tate_matrix": [[0]],
                    "cover_searches": ["C"],
                    "pairing_computed_before_cover_search": False,
                }
            )

    def test_local_filter_requires_full_supplied_local_kummer_space(self):
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            images_path = directory_path / "images.json"
            local_map_path = directory_path / "local-map.json"
            output_path = directory_path / "output.json"
            images_path.write_text(
                json.dumps(
                    {
                        "local_dimension": 2,
                        "fingerprint_dimension": 1,
                        "known_mw_images": [
                            {"label": "P", "generator": "delta(P)", "local": 1, "fingerprint": 0}
                        ],
                        "candidate_images": [
                            {"label": "locally-good", "generator": "a", "local": 1, "fingerprint": 0},
                            {"label": "locally-bad", "generator": "b", "local": 2, "fingerprint": 1},
                        ],
                        "class_quotient_certification": {
                            "method": "analytic class-number bound",
                            "remaining_dimension_upper_bound": 0,
                        },
                    }
                )
            )
            local_map_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-local-kummer-map.v1",
                        "local_dimension": 2,
                        "method": "independent local descent",
                        "allowed_local_images": [1],
                    }
                )
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LOCAL_FILTER_SCRIPT),
                    "--images",
                    str(images_path),
                    "--local-kummer-map",
                    str(local_map_path),
                    "--output",
                    str(output_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("local_survivors=1", result.stdout)
            output = json.loads(output_path.read_text())
            self.assertEqual(output["local_survivor_count"], 1)
            self.assertEqual(output["locally_rejected_candidates"][0]["label"], "locally-bad")
            self.assertEqual(
                output["post_local_quotient"]["candidate_images"][0]["label"],
                "locally-good",
            )

    def test_cover_local_filter_rejects_only_certified_obstructions(self):
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            images_path = directory_path / "images.json"
            cover_audit_path = directory_path / "cover-audit.json"
            output_path = directory_path / "output.json"
            images_path.write_text(
                json.dumps(
                    {
                        "local_dimension": 1,
                        "fingerprint_dimension": 0,
                        "known_mw_images": [],
                        "candidate_images": [
                            {"label": "good", "generator": "a", "local": 0, "fingerprint": 0},
                            {"label": "bad", "generator": "b", "local": 1, "fingerprint": 0},
                            {"label": "unknown", "generator": "c", "local": 1, "fingerprint": 0},
                        ],
                        "class_quotient_certification": {
                            "method": "analytic class-number bound",
                            "remaining_dimension_upper_bound": 0,
                        },
                    }
                )
            )
            cover_audit_path.write_text(
                json.dumps(
                    {
                        "protocol": "BNFFREECOVERLOCAL-v1",
                        "covers": [
                            {
                                "label": "good",
                                "finite_places": [
                                    {"classification": "PROVED_QP_POINT_BY_SMOOTH_FP_LIFT"}
                                ],
                            },
                            {
                                "label": "bad",
                                "finite_places": [
                                    {"classification": "PROVED_NO_QP_POINT_BY_EMPTY_FP_REDUCTION"}
                                ],
                            },
                            {
                                "label": "unknown",
                                "finite_places": [
                                    {"classification": "INCONCLUSIVE_SINGULAR_LIFT_PRECISION"}
                                ],
                            },
                        ],
                    }
                )
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LOCAL_FILTER_SCRIPT),
                    "--images",
                    str(images_path),
                    "--cover-local-audit",
                    str(cover_audit_path),
                    "--output",
                    str(output_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("finite_local_survivors=1", result.stdout)
            output = json.loads(output_path.read_text())
            self.assertEqual(output["finite_local_survivor_count"], 1)
            self.assertEqual(output["locally_obstructed_candidates"][0]["label"], "bad")
            self.assertEqual(output["locally_inconclusive_candidates"][0]["label"], "unknown")

    def test_coverage_filter_uses_only_certified_local_places(self):
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            images_path = directory_path / "images.json"
            coverage_path = directory_path / "coverage.json"
            output_path = directory_path / "output.json"
            images_path.write_text(
                json.dumps(
                    {
                        "local_dimension": 2,
                        "fingerprint_dimension": 0,
                        "known_mw_images": [
                            {"label": "P", "generator": "delta(P)", "local": 0, "fingerprint": 0}
                        ],
                        "candidate_images": [
                            {"label": "good", "generator": "a", "local": 0, "fingerprint": 0},
                            {"label": "bad-at-17", "generator": "b", "local": 1, "fingerprint": 0},
                        ],
                        "class_quotient_certification": {
                            "method": "analytic class-number bound",
                            "remaining_dimension_upper_bound": 0,
                        },
                    }
                )
            )
            coverage_path.write_text(
                json.dumps(
                    {
                        "protocol": "BNFFREELOCALCOVERAGE-v1",
                        "signature_local_dimension": 2,
                        "known_mw_local_images": [{"label": "P", "local": "0x0"}],
                        "odd_places": [
                            {
                                "rational_prime": 17,
                                "coordinate_indices": [0],
                                "classification": "CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE",
                            }
                        ],
                        "real_place": {"classification": "UNRESOLVED_REAL_KUMMER_IMAGE_COVERAGE"},
                    }
                )
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LOCAL_FILTER_SCRIPT),
                    "--images",
                    str(images_path),
                    "--local-coverage-audit",
                    str(coverage_path),
                    "--output",
                    str(output_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("coverage_local_survivors=1", result.stdout)
            output = json.loads(output_path.read_text())
            self.assertEqual(output["certified_covered_places"], ["p=17"])
            self.assertEqual(output["locally_rejected_candidates"][0]["label"], "bad-at-17")


if __name__ == "__main__":
    unittest.main()

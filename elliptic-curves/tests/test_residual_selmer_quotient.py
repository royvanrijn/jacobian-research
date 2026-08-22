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

PAIRING_SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "audit_residual_cassels_tate.py"
LOCAL_FILTER_SCRIPT = Path(__file__).resolve().parents[1] / "cas" / "filter_bnf_free_local_selmer.py"
pairing_spec = importlib.util.spec_from_file_location("audit_residual_cassels_tate", PAIRING_SCRIPT)
pairing_module = importlib.util.module_from_spec(pairing_spec)
sys.modules[pairing_spec.name] = pairing_module
assert pairing_spec.loader is not None
pairing_spec.loader.exec_module(pairing_module)


class ResidualSelmerQuotientTests(unittest.TestCase):
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

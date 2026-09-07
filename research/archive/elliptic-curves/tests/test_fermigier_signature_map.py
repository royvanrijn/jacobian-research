from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SIGNATURE_SCRIPT = CAS / "run_fermigier_rank20_auxiliary_fingerprints.py"
QUOTIENT_SCRIPT = CAS / "residual_selmer_quotient.py"
RELATION_POOL_SCRIPT = CAS / "analyze_curve273_relation_pool.py"
CURVE273_SIGNATURE_SCRIPT = CAS / "analyze_curve273_kummer_fingerprint.py"
EVALUATE_SIGNATURE_SCRIPT = CAS / "evaluate_bnf_free_signature_map.py"
EXTRACT_SQUARECLASSES_SCRIPT = CAS / "extract_bnf_free_squareclasses.py"
NORM_FILTER_SCRIPT = CAS / "filter_bnf_free_norm_condition.py"
TWO_COVER_SCRIPT = CAS / "build_bnf_free_two_covers.py"
TWO_COVER_LOCAL_AUDIT_SCRIPT = CAS / "audit_bnf_free_two_cover_reduction.py"
LOCAL_COVERAGE_AUDIT_SCRIPT = CAS / "audit_bnf_free_local_kummer_coverage.py"
CLASS_QUOTIENT_AUDIT_SCRIPT = CAS / "audit_bnf_free_s_class_quotient.py"
CANONICAL_RELATION_AUGMENT_SCRIPT = CAS / "augment_bnf_free_canonical_principal_relations.py"
FERMIGIER_RELATION_SCRIPT = CAS / "run_fermigier_rank20_fixedfb_quadratic_specialq.py"
FERMIGIER_LINEAR_RELATION_SCRIPT = CAS / "run_fermigier_rank20_fixedfb_specialq.py"
FERMIGIER_MINKOWSKI_RELATION_SCRIPT = CAS / "run_fermigier_rank20_minkowski_specialq.py"


class FermigierSignatureMapTests(unittest.TestCase):
    def test_squareclass_extractor_keeps_an_individually_s_supported_relation(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required to multiply cubic-field generators")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            ledger_path = directory_path / "relations.json"
            candidates_path = directory_path / "candidates.json"
            ledger_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
                        "defining_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "S_columns": [0],
                        "generators": [
                            {"power_basis": ["1", "1", "0"]},
                            {"power_basis": ["2", "0", "0"]},
                        ],
                        "closed_relations": [
                            {
                                "kind": "single-large-prime-cycle",
                                "fb_parity_mask_hex": "0x1",
                                "generator_indices": [0, 1],
                            },
                            {
                                "kind": "single-large-prime-cycle",
                                "fb_parity_mask_hex": "0x2",
                                "generator_indices": [0],
                            },
                        ],
                    }
                )
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(EXTRACT_SQUARECLASSES_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(candidates_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            self.assertIn("candidates=1", result.stdout)
            candidates = json.loads(candidates_path.read_text())
            self.assertEqual(candidates["skipped_non_S_supported_relations"], 1)
            self.assertEqual(candidates["candidates"][0]["generator_coefficients"], ["2", "2", "0"])

    def test_squareclass_extractor_cancels_non_s_columns_across_relations(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required to multiply cubic-field generators")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            ledger_path = directory_path / "relations.json"
            candidates_path = directory_path / "candidates.json"
            ledger_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
                        "defining_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "S_columns": [0],
                        "generators": [
                            {"power_basis": ["1", "1", "0"]},
                            {"power_basis": ["2", "0", "0"]},
                        ],
                        "closed_relations": [
                            {
                                "kind": "single-large-prime-cycle",
                                "fb_parity_mask_hex": "0x2",
                                "generator_indices": [0],
                            },
                            {
                                "kind": "double-large-prime-cycle",
                                "fb_parity_mask_hex": "0x3",
                                "generator_indices": [1],
                            },
                        ],
                    }
                )
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(EXTRACT_SQUARECLASSES_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(candidates_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            candidates = json.loads(candidates_path.read_text())
            self.assertEqual(candidates["candidate_generation"], "non_S_projection_kernel")
            self.assertEqual(candidates["s_supported_kernel_dimension"], 1)
            self.assertEqual(candidates["candidates"][0]["source_relation_indices"], [0, 1])
            self.assertEqual(candidates["candidates"][0]["factor_base_support"], [0])
            self.assertEqual(candidates["candidates"][0]["generator_coefficients"], ["2", "2", "0"])

    def test_norm_filter_keeps_only_rational_square_norms(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for cubic-field norms")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            candidates_path = directory_path / "candidates.json"
            output_path = directory_path / "norm-filtered.json"
            candidates_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-squareclass-candidates.v1",
                        "field_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "candidates": [
                            {"label": "theta", "generator_coefficients": ["0", "1", "0"]},
                            {"label": "one", "generator_coefficients": ["1", "0", "0"]},
                            {"label": "two", "generator_coefficients": ["2", "0", "0"]},
                        ],
                    }
                )
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(NORM_FILTER_SCRIPT),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(output_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            self.assertIn("candidates=1", result.stdout)
            output = json.loads(output_path.read_text())
            self.assertEqual(output["candidates"][0]["label"], "theta")
            self.assertEqual(output["norm_rejected"][0]["label"], "two")
            self.assertEqual(output["globally_square_rejected"][0]["label"], "one")

    def test_norm_kernel_finds_square_norm_products_not_visible_individually(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for exact cubic norms")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            candidates_path = directory_path / "candidates.json"
            output_path = directory_path / "norm-kernel.json"
            candidates_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-squareclass-candidates.v1",
                        "field_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "candidates": [
                            {"label": "two", "generator_coefficients": ["2", "0", "0"]},
                            {"label": "two-theta", "generator_coefficients": ["0", "2", "0"]},
                        ],
                    }
                )
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(NORM_FILTER_SCRIPT),
                    "--candidates",
                    str(candidates_path),
                    "--generate-norm-kernel",
                    "--selmer-rational-primes",
                    "2",
                    "--output",
                    str(output_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            self.assertIn("norm_kernel_dimension=1", result.stdout)
            output = json.loads(output_path.read_text())
            self.assertEqual(output["norm_kernel_dimension"], 1)
            self.assertEqual(output["candidates"][0]["generator_coefficients"], ["0", "4", "0"])
            self.assertEqual(output["candidates"][0]["source_candidate_labels"], ["two", "two-theta"])

    def test_two_cover_builder_materializes_homogeneous_quadrics(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required to write exact cover equations")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            candidates_path = directory_path / "norm-filtered.json"
            covers_path = directory_path / "covers.json"
            candidates_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1",
                        "field_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "candidates": [
                            {
                                "label": "unit",
                                "generator_coefficients": ["1", "0", "0"],
                                "norm": "1",
                            }
                        ],
                    }
                )
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(TWO_COVER_SCRIPT),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(covers_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            self.assertIn("covers=1", result.stdout)
            output = json.loads(covers_path.read_text())
            quadrics = output["covers"][0]["quadrics"]
            self.assertEqual(quadrics["theta_squared_coefficient"], "v^2 + 2*u*w + w^2")
            self.assertEqual(quadrics["theta_coefficient_plus_z_squared"], "2*u*v + 2*v*w + w^2 + z^2")

    def test_two_cover_reduction_audit_certifies_a_smooth_local_lift(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for finite-field cover reduction")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            candidates_path = directory_path / "norm-filtered.json"
            covers_path = directory_path / "covers.json"
            audit_path = directory_path / "local-audit.json"
            candidates_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1",
                        "field_polynomial_ascending": ["-1", "-1", "0", "1"],
                        "candidates": [
                            {
                                "label": "unit",
                                "generator_coefficients": ["1", "0", "0"],
                                "norm": "1",
                            }
                        ],
                    }
                )
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(TWO_COVER_SCRIPT),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(covers_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(TWO_COVER_LOCAL_AUDIT_SCRIPT),
                    "--covers",
                    str(covers_path),
                    "--primes",
                    "2,3,5",
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            self.assertIn("PROVED_QP_POINT_BY_SMOOTH_FP_LIFT", result.stdout)
            self.assertIn("PROVED_QP_POINT_BY_SINGULAR_HENSEL_LIFT", result.stdout)
            audit = json.loads(audit_path.read_text())
            classifications = {
                item["rational_prime"]: item["classification"]
                for item in audit["covers"][0]["finite_places"]
            }
            self.assertEqual(
                classifications[2], "PROVED_QP_POINT_BY_SINGULAR_HENSEL_LIFT"
            )
            self.assertEqual(classifications[3], "PROVED_QP_POINT_BY_SMOOTH_FP_LIFT")
            self.assertEqual(classifications[5], "PROVED_QP_POINT_BY_SMOOTH_FP_LIFT")

    def test_odd_local_coverage_audit_does_not_overclaim_two_adic_coverage(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for local reduction data")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            signature_path = directory_path / "signature.json"
            coverage_path = directory_path / "coverage.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(SIGNATURE_SCRIPT),
                    "--prime-bound",
                    "59",
                    "--output",
                    str(signature_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(LOCAL_COVERAGE_AUDIT_SCRIPT),
                    "--signature-map",
                    str(signature_path),
                    "--output",
                    str(coverage_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("certified_odd=5", result.stdout)
            coverage = json.loads(coverage_path.read_text())
            self.assertEqual(
                coverage["two_adic_place"]["classification"],
                "UNRESOLVED_TWO_ADIC_LOCAL_KUMMER_IMAGE_COVERAGE",
            )
            by_prime = {item["rational_prime"]: item for item in coverage["odd_places"]}
            self.assertEqual(
                by_prime[17]["classification"],
                "CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE",
            )
            self.assertEqual(
                coverage["real_place"]["classification"],
                "CERTIFIED_FULL_REAL_KUMMER_IMAGE_COVERAGE",
            )

    def test_signature_evaluator_uses_exact_real_signs_for_large_cancellation(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for algebraic-real sign comparisons")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            signature_path = directory_path / "signature.json"
            candidates_path = directory_path / "candidates.json"
            output_path = directory_path / "images.json"
            root_scale = 10**8
            signature_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-signature-map.v1",
                        "defining_polynomial_ascending": ["1", str(-(root_scale**2)), "0", "1"],
                        "local_coordinates": [
                            {"kind": "real_sign", "embedding_index": index}
                            for index in range(3)
                        ],
                        "fingerprint_coordinates": [],
                        "fingerprint_dimension": 0,
                        "known_mw_images": [],
                        "class_quotient_certification": {"method": "none"},
                    }
                )
            )
            # (theta-root_scale)^2 has a very small positive value at the
            # positive root near root_scale, despite coefficients of size 1e16.
            candidates_path.write_text(
                json.dumps(
                    {
                        "schema": "elliptic-curves.bnf-free-squareclass-candidates.v1",
                        "candidates": [
                            {
                                "label": "square",
                                "generator_coefficients": [
                                    str(root_scale**2),
                                    str(-2 * root_scale),
                                    "1",
                                ],
                            }
                        ],
                    }
                )
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(EVALUATE_SIGNATURE_SCRIPT),
                    "--signature-map",
                    str(signature_path),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(output_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            output = json.loads(output_path.read_text())
            self.assertEqual(output["candidate_images"][0]["local"], "0x0")

    def test_known_kummer_map_is_faithful_and_consumable(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the local-squareclass computation")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            signature_path = directory_path / "signature.json"
            audit_path = directory_path / "audit.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(SIGNATURE_SCRIPT),
                    "--prime-bound",
                    "59",
                    "--output",
                    str(signature_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("faithful_known_kummer_fingerprint", result.stdout)

            signature = json.loads(signature_path.read_text())
            self.assertEqual(
                signature["schema"],
                "elliptic-curves.bnf-free-signature-map.v1",
            )
            self.assertEqual(len(signature["known_mw_images"]), 20)
            self.assertEqual(signature["known_mw_target_rank"], 20)
            self.assertEqual(signature["selected_auxiliary_primes"], [11, 19, 23, 29, 59])
            self.assertEqual(signature["local_dimension"], 51)
            self.assertEqual(signature["fingerprint_dimension"], 24)

            subprocess.run(
                [
                    sys.executable,
                    str(QUOTIENT_SCRIPT),
                    "--input",
                    str(signature_path),
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            audit = json.loads(audit_path.read_text())
            self.assertEqual(audit["known_mw_target_rank"], 20)
            self.assertEqual(audit["candidate_images"], [])

    def test_relation_ledger_labels_an_insufficient_factor_base_uncertified(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the ideal-relation computation")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            ledger_path = directory_path / "relations.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(RELATION_POOL_SCRIPT),
                    "--glob",
                    str(directory_path / "no-logs-*.log"),
                    "--include-full-ideal-chain",
                    "--write-principal-relations",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn(
                "UNCERTIFIED_FACTOR_BASE_DOES_NOT_REACH_MINKOWSKI_BOUND",
                result.stdout,
            )
            ledger = json.loads(ledger_path.read_text())
            self.assertEqual(len(ledger["relations"]), 11)
            self.assertTrue(
                all(row["generator_power_basis"] for row in ledger["relations"])
            )
            self.assertIn("generators", ledger)
            self.assertIn("closed_relations", ledger)
            self.assertFalse(
                ledger["factor_base_completion"]["materialized_complete_factor_base"]
            )
            audit_path = directory_path / "class-audit.json"
            audit_result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CLASS_QUOTIENT_AUDIT_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn(
                "classification=UNCERTIFIED_INCOMPLETE_FACTOR_BASE",
                audit_result.stdout,
            )

    def test_fermigier_quadratic_collector_writes_exact_generator_ledger(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            checkpoint = directory_path / "checkpoint.json"
            ledger_path = directory_path / "relations.json"
            signature_path = directory_path / "signature.json"
            candidates_path = directory_path / "candidates.json"
            images_path = directory_path / "images.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_RELATION_SCRIPT),
                    "--factor-base-bound",
                    "100",
                    "--special-q-min",
                    "3",
                    "--special-q-max",
                    "100",
                    "--max-special-q",
                    "3",
                    "--pairs-per-q",
                    "80",
                    "--coeff-bound",
                    "12",
                    "--checkpoint",
                    str(checkpoint),
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            ledger = json.loads(ledger_path.read_text())
            self.assertEqual(
                ledger["schema"],
                "elliptic-curves.bnf-free-principal-relation-ledger.v1",
            )
            self.assertEqual(ledger["factor_base_bound"], 100)
            self.assertTrue(ledger["factor_base_completion"])
            self.assertTrue(
                all("rational_prime" in row for row in ledger["factor_base"])
            )
            self.assertGreater(len(ledger["generators"]), 0)
            for relation in ledger["closed_relations"]:
                self.assertTrue(
                    all(
                        len(ledger["generators"][index]["power_basis"]) == 3
                        for index in relation["generator_indices"]
                    )
                )
            candidates_path.write_text(
                json.dumps(
                    [
                        {
                            "label": "collector-generator-0",
                            "generator_coefficients": ledger["generators"][0]["power_basis"],
                        }
                    ]
                )
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(SIGNATURE_SCRIPT),
                    "--prime-bound",
                    "59",
                    "--output",
                    str(signature_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(EVALUATE_SIGNATURE_SCRIPT),
                    "--signature-map",
                    str(signature_path),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(images_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            images = json.loads(images_path.read_text())
            self.assertEqual(images["candidate_images"][0]["label"], "collector-generator-0")

            class_audit_path = directory_path / "class-audit.json"
            audit_result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CLASS_QUOTIENT_AUDIT_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(class_audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("classification=UNCERTIFIED_FACTOR_BASE", audit_result.stdout)
            class_audit = json.loads(class_audit_path.read_text())
            self.assertTrue(class_audit["principal_relations_verified"])
            self.assertIsNone(class_audit["k_s_2_dimension_upper_bound"])

            augmented_path = directory_path / "relations-with-canonical.json"
            augmentation = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CANONICAL_RELATION_AUGMENT_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(augmented_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("status=EXACT_PRINCIPAL_ROWS", augmentation.stdout)
            augmented = json.loads(augmented_path.read_text())
            canonical = augmented["canonical_principal_relations"]
            self.assertTrue(canonical["principal_generators_stored"])
            self.assertEqual(
                canonical["completed_relation_count"],
                len({row["rational_prime"] for row in augmented["factor_base"]}),
            )
            canonical_rows = [
                row
                for row in augmented["closed_relations"]
                if row.get("source") == "canonical_rational_prime_principal"
            ]
            self.assertEqual(len(canonical_rows), canonical["completed_relation_count"])

            idempotent_path = directory_path / "relations-with-canonical-again.json"
            repeat = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CANONICAL_RELATION_AUGMENT_SCRIPT),
                    "--relation-ledger",
                    str(augmented_path),
                    "--output",
                    str(idempotent_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("|added=0|", repeat.stdout)
            idempotent = json.loads(idempotent_path.read_text())
            self.assertEqual(idempotent["generators"], augmented["generators"])
            self.assertEqual(
                idempotent["closed_relations"], augmented["closed_relations"]
            )

            augmented_audit_path = directory_path / "canonical-class-audit.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(CLASS_QUOTIENT_AUDIT_SCRIPT),
                    "--relation-ledger",
                    str(augmented_path),
                    "--output",
                    str(augmented_audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            augmented_audit = json.loads(augmented_audit_path.read_text())
            self.assertTrue(augmented_audit["principal_relations_verified"])
            self.assertEqual(
                augmented_audit["canonical_rational_principal_relation_audit"]["status"],
                "COMPLETE_AND_VERIFIED",
            )
            self.assertEqual(
                augmented_audit["relation_source_rank_analysis"][
                    "canonical_rational_principal_relation_rank"
                ],
                canonical["completed_relation_count"],
            )
            self.assertLessEqual(
                augmented_audit["factor_base_quotient_dimension"],
                class_audit["factor_base_quotient_dimension"],
            )

    def test_fermigier_double_large_prime_closures_retain_generator_products(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            ledger_path = Path(directory) / "relations.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_LINEAR_RELATION_SCRIPT),
                    "--factor-base-bound",
                    "5000",
                    "--special-q-min",
                    "3",
                    "--special-q-max",
                    "3",
                    "--max-special-q",
                    "1",
                    "--seed-specials",
                    "5689:5096",
                    "--b-bound",
                    "5000",
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=90,
                check=True,
            )
            ledger = json.loads(ledger_path.read_text())
            closed = ledger["closed_relations"]
            self.assertGreater(len(closed), 0)
            self.assertTrue(all(row["generator_indices"] for row in closed))
            self.assertGreaterEqual(
                max(len(row["generator_indices"]) for row in closed),
                2,
            )

    def test_fermigier_minkowski_ideal_collector_writes_auditable_ledger(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            ledger_path = directory_path / "relations.json"
            audit_path = directory_path / "class-audit.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_MINKOWSKI_RELATION_SCRIPT),
                    "--factor-base-bound",
                    "100",
                    "--seed-specials",
                    "5689:5096",
                    "--max-special-q",
                    "1",
                    "--lattice-combination-bound",
                    "1",
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            ledger = json.loads(ledger_path.read_text())
            self.assertGreater(len(ledger["generators"]), 0)
            self.assertTrue(
                all("shape_twist" in row for row in ledger["generators"])
            )
            audit = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CLASS_QUOTIENT_AUDIT_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("classification=UNCERTIFIED_FACTOR_BASE", audit.stdout)

    def test_minkowski_collector_supports_curve273_complex_embedding(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            ledger_path = Path(directory) / "curve273-relations.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_MINKOWSKI_RELATION_SCRIPT),
                    "--curve273",
                    "--factor-base-bound",
                    "100",
                    "--special-q-min",
                    "101",
                    "--special-q-max",
                    "101",
                    "--max-special-q",
                    "1",
                    "--lattice-combination-bound",
                    "1",
                    "--shape-twists=0:0",
                    "--trial-prime-bound",
                    "100",
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("archimedean_signature=1,1", result.stdout)
            ledger = json.loads(ledger_path.read_text())
            self.assertEqual(ledger["curve_preset"], "icarm-273")
            self.assertEqual(ledger["selmer_rational_primes"][0], 2)
            self.assertIn(3, ledger["selmer_rational_primes"])
            self.assertIn("collection_early_quotient", ledger)
            self.assertLessEqual(
                ledger["collection_early_quotient"]["dimension_after_canonical_rows_and_S"],
                len(ledger["factor_base"]),
            )

            canonical_path = Path(directory) / "curve273-relations-canonical.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(CANONICAL_RELATION_AUGMENT_SCRIPT),
                    "--relation-ledger",
                    str(ledger_path),
                    "--output",
                    str(canonical_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            canonical = json.loads(canonical_path.read_text())
            self.assertEqual(
                canonical["canonical_principal_relations"]["completed_relation_count"],
                len({row["rational_prime"] for row in canonical["factor_base"]}),
            )

    def test_minkowski_collector_supports_double_special_cycle_merging(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            ledger_path = Path(directory) / "pair-cycle-relations.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_MINKOWSKI_RELATION_SCRIPT),
                    "--curve273",
                    "--factor-base-bound",
                    "100",
                    "--special-q-min",
                    "101",
                    "--special-q-max",
                    "300",
                    "--max-special-q",
                    "3",
                    "--special-ideal-mode",
                    "cycle-pairs",
                    "--pair-cycle-length",
                    "3",
                    "--large-prime-merge-mode",
                    "spanning-forest",
                    "--lattice-combination-bound",
                    "1",
                    "--shape-twists=0:0",
                    "--trial-prime-bound",
                    "100",
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("mode=cycle-pairs", result.stdout)
            ledger = json.loads(ledger_path.read_text())
            self.assertEqual(ledger["special_ideal_mode"], "cycle-pairs")
            self.assertEqual(ledger["large_prime_merge_mode"], "spanning-forest")
            self.assertTrue(
                all(
                    len(generator["source_special_ideals"]) == 2
                    for generator in ledger["generators"]
                )
            )

    def test_minkowski_collector_supports_degree_two_special_ideals(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the cubic-field collector")

        with tempfile.TemporaryDirectory() as directory:
            ledger_path = Path(directory) / "degree-two-relations.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(FERMIGIER_MINKOWSKI_RELATION_SCRIPT),
                    "--curve273",
                    "--factor-base-bound",
                    "100",
                    "--special-residue-degree",
                    "2",
                    "--special-q-min",
                    "101",
                    "--special-q-max",
                    "300",
                    "--max-special-q",
                    "1",
                    "--lattice-combination-bound",
                    "1",
                    "--shape-twists=0:0",
                    "--trial-prime-bound",
                    "100",
                    "--relation-ledger",
                    str(ledger_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("count=1", result.stdout)
            ledger = json.loads(ledger_path.read_text())
            self.assertEqual(ledger["special_residue_degree"], 2)
            self.assertTrue(
                all(
                    "[" in generator["source_special_ideals"][0][1]
                    for generator in ledger["generators"]
                )
            )

    def test_curve273_known_kummer_map_is_faithful_and_consumable(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the local-squareclass computation")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            signature_path = directory_path / "signature.json"
            audit_path = directory_path / "audit.json"
            result = subprocess.run(
                [
                    sage,
                    "-python",
                    str(CURVE273_SIGNATURE_SCRIPT),
                    "--prime-bound",
                    "5000",
                    "--output",
                    str(signature_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            self.assertIn("FAITHFUL_KNOWN_KUMMER_FINGERPRINT", result.stdout)
            signature = json.loads(signature_path.read_text())
            self.assertEqual(len(signature["known_mw_images"]), 30)
            self.assertEqual(signature["known_mw_target_rank"], 30)
            self.assertEqual(signature["local_dimension"], 59)
            self.assertEqual(signature["fingerprint_dimension"], 54)

            subprocess.run(
                [
                    sys.executable,
                    str(QUOTIENT_SCRIPT),
                    "--input",
                    str(signature_path),
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            audit = json.loads(audit_path.read_text())
            self.assertEqual(audit["known_mw_target_rank"], 30)

    def test_curve273_evaluator_reproduces_a_known_kummer_image(self):
        sage = shutil.which("sage")
        if sage is None:
            self.skipTest("Sage is required for the local-squareclass computation")

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            signature_path = directory_path / "signature.json"
            candidates_path = directory_path / "candidates.json"
            images_path = directory_path / "images.json"
            audit_path = directory_path / "audit.json"
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(CURVE273_SIGNATURE_SCRIPT),
                    "--prime-bound",
                    "5000",
                    "--output",
                    str(signature_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            signature = json.loads(signature_path.read_text())
            known = signature["known_mw_images"][0]
            candidates_path.write_text(
                json.dumps(
                    [
                        {
                            "label": "replayed-P1",
                            "generator_coefficients": known["generator_coefficients"],
                        },
                        {
                            "label": "rational-2",
                            "generator_coefficients": ["2", "0", "0"],
                        },
                    ]
                )
            )
            subprocess.run(
                [
                    sage,
                    "-python",
                    str(EVALUATE_SIGNATURE_SCRIPT),
                    "--signature-map",
                    str(signature_path),
                    "--candidates",
                    str(candidates_path),
                    "--output",
                    str(images_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
            images = json.loads(images_path.read_text())
            self.assertEqual(len(images["candidate_images"]), 2)
            self.assertGreaterEqual(images["local_dimension"], signature["local_dimension"])
            self.assertEqual(images["candidate_images"][0]["local"], known["local"])
            self.assertEqual(images["candidate_images"][0]["fingerprint"], known["fingerprint"])
            subprocess.run(
                [
                    sys.executable,
                    str(QUOTIENT_SCRIPT),
                    "--input",
                    str(images_path),
                    "--output",
                    str(audit_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            audit = json.loads(audit_path.read_text())
            self.assertTrue(audit["candidate_images"][0]["killed_by_known_mw_in_this_target"])


if __name__ == "__main__":
    unittest.main()

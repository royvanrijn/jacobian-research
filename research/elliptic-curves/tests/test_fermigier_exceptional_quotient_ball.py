#!/usr/bin/env python3
"""Focused replay checks for the signed Fermigier exceptional quotient ball."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
EC = ROOT / "elliptic-curves"
sys.path[:0] = [str(EC), str(CAS)]

from classify_fermigier_exceptional_quotient_ball import (  # noqa: E402
    affine_slice_coefficients,
    polynomial_sha256,
    result_digest,
    sha256_file,
    sha256_lines,
    signed_weight_two_vectors,
)


SCRIPT = CAS / "classify_fermigier_exceptional_quotient_ball.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elliptic_fermigier_exceptional_quotient_ball.json"
)
EXPECTED_ARTIFACT_SHA256 = (
    "fd70a03a49c27322d6b74fcc2dad30576e761b8fbe46f911fa302a5c004192b2"
)
EXPECTED_RESULT_SHA256 = (
    "0d731a8aac70e7238dc20f4e392d0c0d7ee59cea4b5de3588d506110ce05af01"
)


class FermigierExceptionalQuotientBallTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = json.loads(ARTIFACT.read_text())

    def test_signed_vector_populations_are_complete(self) -> None:
        self.assertEqual(len(signed_weight_two_vectors(10)), 200)
        self.assertEqual(len(signed_weight_two_vectors(8)), 128)
        self.assertEqual(
            Counter(sum(value != 0 for value in vector) for vector in signed_weight_two_vectors(10)),
            {1: 20, 2: 180},
        )
        self.assertEqual(
            Counter(sum(value != 0 for value in vector) for vector in signed_weight_two_vectors(8)),
            {1: 16, 2: 112},
        )

    def test_pinned_sources_and_canonical_aliases(self) -> None:
        data = self.artifact
        self.assertEqual(sha256_file(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(data["sources"]["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(
            data["anchors"]["E22"]["canonical_parameter"], {"u": "19754/39"}
        )
        self.assertEqual(
            data["anchors"]["E22"]["aliases"], {"literal_shift_T": "39508/39"}
        )
        self.assertEqual(
            data["anchors"]["rank20"]["canonical_parameter"], {"u": "28917/20"}
        )
        self.assertEqual(
            data["anchors"]["rank20"]["aliases"], {"literal_shift_T": "28917/10"}
        )
        for source_name, hash_name in (
            ("rank22_artifact", "rank22_artifact_sha256"),
            ("rank20_artifact", "rank20_artifact_sha256"),
        ):
            self.assertEqual(
                sha256_file(ROOT / data["sources"][source_name]),
                data["sources"][hash_name],
            )

    def test_all_group_and_generic_coset_checks_are_recorded(self) -> None:
        balls = self.artifact["direction_balls"]
        self.assertEqual(balls["E22"]["signed_direction_count"], 200)
        self.assertEqual(balls["rank20"]["signed_direction_count"], 128)
        for ball in balls.values():
            records = ball["records"]
            self.assertEqual(len(records), ball["signed_direction_count"])
            self.assertEqual(
                len({record["direction_id"] for record in records}), len(records)
            )
            self.assertTrue(ball["all_exact_group_relations_verified"])
            self.assertTrue(ball["all_mod5_separated_from_generic_span"])
            self.assertTrue(
                all(record["exact_pointed_quartic_round_trip"] for record in records)
            )
            self.assertTrue(
                all(record["exact_canonical_group_relation"] for record in records)
            )
            self.assertTrue(
                all(record["mod5_outside_generic_span"] for record in records)
            )

    def test_complete_affine_manifest_and_genus(self) -> None:
        audit = self.artifact["affine_transport"]
        records = audit["records"]
        self.assertEqual(audit["pair_count"], 25_600)
        self.assertEqual(len(records), 25_600)
        self.assertEqual(len(audit["low_genus_candidates"]), 0)
        self.assertEqual(
            audit["factor_squareclass_histogram"],
            [
                {
                    "count": 25_600,
                    "degree": 6,
                    "factor_signature": [[6, 1]],
                    "squareclass_genus": 2,
                    "squareclass_kernel_degree": 6,
                }
            ],
        )
        lines = (
            f"{record['left_direction_id']}|{record['right_direction_id']}|"
            f"{record['degree']}|{record['factor_signature']}|"
            f"{record['irreducible_mod_prime_witness']}|"
            f"{record['squareclass_kernel_degree']}|"
            f"{record['primitive_polynomial_sha256']}"
            for record in records
        )
        self.assertEqual(sha256_lines(lines), audit["manifest_sha256"])
        self.assertTrue(
            all(
                record["degree"] == 6
                and record["factor_signature"] == [[6, 1]]
                and record["squareclass_kernel_degree"] == 6
                and record["squareclass_genus"] == 2
                for record in records
            )
        )

    def test_three_polynomials_reconstruct_from_endpoints(self) -> None:
        data = self.artifact
        left = {
            record["direction_id"]: record["quartic_point"]["x"]
            for record in data["direction_balls"]["E22"]["records"]
        }
        right = {
            record["direction_id"]: record["quartic_point"]["x"]
            for record in data["direction_balls"]["rank20"]["records"]
        }
        records = data["affine_transport"]["records"]
        for record in (records[0], records[len(records) // 2], records[-1]):
            coefficients = affine_slice_coefficients(
                left[record["left_direction_id"]],
                right[record["right_direction_id"]],
            )
            self.assertEqual(len(coefficients) - 1, 6)
            self.assertEqual(
                polynomial_sha256(coefficients),
                record["primitive_polynomial_sha256"],
            )

    def test_result_digest_and_outcome(self) -> None:
        self.assertEqual(self.artifact["result_sha256"], EXPECTED_RESULT_SHA256)
        self.assertEqual(result_digest(self.artifact), EXPECTED_RESULT_SHA256)
        self.assertEqual(
            self.artifact["outcome"],
            {
                "affine_low_genus_candidates": 0,
                "cross_anchor_affine_interpolants": 25_600,
                "new_base_changes": 0,
                "new_sections": 0,
                "new_specializations": 0,
                "signed_E22_exceptional_directions": 200,
                "signed_rank20_exceptional_directions": 128,
                "target_met": False,
            },
        )


if __name__ == "__main__":
    unittest.main()

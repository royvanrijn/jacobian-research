#!/usr/bin/env python3
"""Focused checks for the published-direction fiber-product screen."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))

from fermigier_mestre import FermigierMestreFamily  # noqa: E402


SCRIPT = CAS / "search_fermigier_published_pair_fiber_products.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_published_pair_fiber_products.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "9385499d9a3cf05b5a04ecc6133b8233cb12a7e057026e38a13ff33ef9fde21c"
)
EXPECTED_ARTIFACT_SHA256 = (
    "f6c65ee07f10a915c2073986bc48c2f48f2a7ad55697dbce062e6198da33f849"
)
EXPECTED_SOURCE_LABELS = {"P6", *{f"P{index}" for index in range(13, 23)}}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierPairFiberProductTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_finalized_source(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        source = self.data["source"]
        self.assertEqual(
            source["published_accidental_preimage_sha256"],
            "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15",
        )
        self.assertEqual(set(source["published_accidental_labels"]), EXPECTED_SOURCE_LABELS)

    def test_slice_and_pair_population_is_exact(self) -> None:
        population = self.data["slice_population"]
        self.assertEqual(population["source_direction_count"], 11)
        self.assertEqual(population["signed_slice_count"], 22)
        slices = population["slices"]
        self.assertEqual({row["source_label"] for row in slices}, EXPECTED_SOURCE_LABELS)
        self.assertEqual(
            {(row["source_label"], row["slope"]) for row in slices},
            {(label, slope) for label in EXPECTED_SOURCE_LABELS for slope in (-1, 1)},
        )
        record_t = Fraction(self.data["source"]["record_parameter_t"])
        for row in slices:
            x_value = Fraction(row["slope"]) * record_t + Fraction(row["intercept"])
            self.assertEqual(x_value, Fraction(row["source_quartic_x"]))
            self.assertEqual(
                Fraction(row["source_quartic_z"]) ** 2,
                FermigierMestreFamily.quartic_value(record_t, x_value),
            )
        pairs = self.data["pair_searches"]
        self.assertEqual(len(pairs), 220)
        self.assertEqual(len({row["pair_id"] for row in pairs}), 220)
        self.assertTrue(all(row["left_source_label"] != row["right_source_label"] for row in pairs))
        self.assertTrue(all(row["factor_polynomials_coprime"] for row in pairs))
        self.assertTrue(all(row["pilot"]["product_degree"] == 8 for row in pairs))

    def test_complete_pilot_has_no_product_points_or_escalations(self) -> None:
        pairs = self.data["pair_searches"]
        self.assertTrue(all(row["pilot"]["search"]["status"] == "completed" for row in pairs))
        self.assertTrue(all(row["pilot"]["search"]["signed_point_count"] == 0 for row in pairs))
        self.assertTrue(all(row["pilot"]["incidences"] == [] for row in pairs))
        self.assertTrue(all(row["pilot"]["qualifying_new_parameter_count"] == 0 for row in pairs))
        self.assertTrue(all(row["escalation"] is None for row in pairs))
        outcome = self.data["outcome"]
        self.assertEqual(outcome["pilot_pairs_attempted"], 220)
        self.assertEqual(outcome["pilot_pairs_completed"], 220)
        self.assertEqual(outcome["pilot_pairs_timed_out_or_errored"], 0)
        self.assertEqual(outcome["productive_pilot_pairs"], 0)
        self.assertEqual(outcome["escalation_pairs_attempted"], 0)
        self.assertFalse(outcome["stopped_as_computationally_disproportionate"])

    def test_decontamination_and_negative_target_scope(self) -> None:
        prior = self.data["prior_decontamination"]
        self.assertEqual(prior["unique_prior_parameter_count"], 590)
        self.assertEqual(
            prior["prior_parameter_sha256"],
            "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291",
        )
        self.assertEqual(self.data["candidates"], [])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["genuinely_new_double_forced_fibers"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertEqual(
            outcome["exact_pair_result_sha256"],
            "80413701447b6468a826fa2185528da74057faa02619bdf02f237e7efb8b1b8b",
        )
        self.assertFalse(self.data["target"]["hit"])
        self.assertTrue(self.data["parameters"]["no_retries"])


if __name__ == "__main__":
    unittest.main()

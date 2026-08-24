#!/usr/bin/env python3
"""Focused checks for the terminal H=50000 pair-product extension."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_fermigier_published_pair_fiber_products_h50000.py"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_published_pair_fiber_products_h50000.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "2d5523e16444f08406e23e596c57198cc9a805cbcaf316f197dc47ca664c53c1"
)
EXPECTED_ARTIFACT_SHA256 = (
    "cdca612edcce38150f77785f1fa24cd0208b38e19be6be6a872703c7625f1e52"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierPairFiberProductH50000Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_honest_execution_provenance(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        provenance = self.data["execution_provenance"]
        self.assertEqual(
            provenance["bounded_execution_script_sha256"],
            "61e2435221d523a9c7b3af28b3aad8be8477c57b44dd9709b7f86ddf709030d7",
        )
        self.assertEqual(provenance["stable_replay_script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertTrue(provenance["metadata_normalization_only"])
        self.assertEqual(provenance["bounded_searches_rerun_during_normalization"], 0)

    def test_stable_exact_dependencies_replace_whole_file_guards(self) -> None:
        source = self.data["source"]
        self.assertEqual(
            source["published_accidental_preimage_sha256"],
            "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15",
        )
        self.assertEqual(
            source["H5000_exact_pair_result_sha256"],
            "80413701447b6468a826fa2185528da74057faa02619bdf02f237e7efb8b1b8b",
        )
        prior = self.data["prior_decontamination"]
        self.assertEqual(prior["unique_prior_parameter_count"], 590)
        self.assertEqual(
            prior["prior_parameter_sha256"],
            "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291",
        )
        self.assertEqual(source["record_parameter_t"], "39508/39")
        self.assertEqual(source["record_parameter_projective_height"], 39_508)

    def test_all_pairs_recover_the_record_calibration(self) -> None:
        pairs = self.data["pair_searches"]
        self.assertEqual(len(pairs), 220)
        self.assertEqual(len({row["pair_id"] for row in pairs}), 220)
        self.assertTrue(all(row["search"]["search"]["status"] == "completed" for row in pairs))
        self.assertTrue(all(row["record_T0_positive_calibration_count"] == 1 for row in pairs))
        for row in pairs:
            calibrations = [
                incidence
                for incidence in row["search"]["incidences"]
                if incidence["classification"] == "record-fiber-excluded"
            ]
            self.assertEqual(len(calibrations), 1)
            calibration = calibrations[0]
            self.assertEqual(calibration["canonical_parameter_t"], "39508/39")
            self.assertTrue(calibration["left_factor_is_square"])
            self.assertTrue(calibration["right_factor_is_square"])
            self.assertEqual(len(calibration["exact_forced_quartic_points"]), 2)
            self.assertTrue(
                all(
                    point["exact_membership_checked"]
                    for point in calibration["exact_forced_quartic_points"]
                )
            )
        outcome = self.data["outcome"]
        self.assertEqual(outcome["pairs_completed"], 220)
        self.assertEqual(outcome["pairs_timed_out_or_errored"], 0)
        self.assertEqual(outcome["record_T0_calibrated_pairs"], 220)

    def test_three_extra_product_points_fail_individual_square_gate(self) -> None:
        extras = []
        signed_counts = Counter()
        for row in self.data["pair_searches"]:
            signed_counts[row["search"]["search"]["signed_point_count"]] += 1
            for incidence in row["search"]["incidences"]:
                if incidence["classification"] == "product-square-only":
                    extras.append(
                        (
                            row["pair_id"],
                            incidence["signed_parameter_t"],
                            incidence["left_factor_is_square"],
                            incidence["right_factor_is_square"],
                        )
                    )
        self.assertEqual(signed_counts, Counter({2: 217, 4: 3}))
        self.assertEqual(
            extras,
            [
                ("p6_m1__p13_p1", "-48363/26", False, False),
                ("p6_p1__p13_m1", "23317/6", False, False),
                ("p13_p1__p15_m1", "-42058/25", False, False),
            ],
        )
        outcome = self.data["outcome"]
        self.assertEqual(
            outcome["incidence_classification_counts"],
            {"product-square-only": 3, "record-fiber-excluded": 220},
        )
        self.assertEqual(
            outcome["exact_H50000_pair_result_sha256"],
            "dea8b716c5aec56817a172afd6e894e7748aaddc482a2d29c0a3360abe55bf4b",
        )

    def test_terminal_negative_target_scope(self) -> None:
        self.assertEqual(self.data["candidates"], [])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["genuinely_new_double_forced_fibers"], 0)
        self.assertEqual(outcome["completed_conductors"], 0)
        self.assertEqual(outcome["rank_triage_count"], 0)
        self.assertIsNone(outcome["maximum_stable_numerical_rank"])
        self.assertFalse(self.data["target"]["hit"])
        self.assertTrue(self.data["parameters"]["terminal_bound"])
        self.assertTrue(self.data["parameters"]["no_retries"])


if __name__ == "__main__":
    unittest.main()

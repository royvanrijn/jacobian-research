#!/usr/bin/env python3
"""Focused replay checks for the closed thirteen-family rational scan."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/search_mestre_rank13_multifamily_rational.py"
SCANNER = ROOT / "elliptic-curves/cas/scan_mestre_rank13_multifamily.cpp"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elliptic_mestre_rank13_multifamily_rational.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "cc16c8aeaa2319eb009be5ea6e6b58f041ba646e147b6786ac7ace5e11895335"
)
EXPECTED_SCANNER_SHA256 = (
    "e4363bd5bc03995eec9edd14710a08868073e3af88ced4874125ab7e3a30ff9d"
)
EXPECTED_ARTIFACT_SHA256 = (
    "0f664e937b9983bd7fa1cfb80269b5c734faddbf6d02dbe4dfca0e3b573ac41f"
)

INCLUDED = [
    [0, 12, 33, 142, 150, 169],
    [0, 12, 50, 93, 114, 131],
    [0, 21, 95, 100, 121, 155],
    [0, 23, 89, 124, 147, 181],
    [0, 23, 93, 128, 133, 175],
    [0, 25, 83, 124, 149, 183],
    [0, 26, 53, 70, 88, 117],
    [0, 32, 65, 97, 108, 148],
    [0, 40, 55, 100, 108, 151],
    [0, 5, 110, 111, 115, 133],
    [0, 7, 54, 127, 148, 166],
    [0, 7, 93, 154, 161, 191],
    [0, 8, 60, 93, 108, 125],
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank13MultifamilyRationalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.records = cls.data["selected_records"]

    def test_pinned_files_result_and_clean_execution(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(SCANNER), EXPECTED_SCANNER_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        provenance = self.data["provenance"]
        self.assertEqual(provenance["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(provenance["scanner_sha256"], EXPECTED_SCANNER_SHA256)
        self.assertEqual(provenance["owned_processes_remaining"], 0)
        self.assertEqual(provenance["same_stage_retries"], 0)
        self.assertTrue(provenance["temporary_scanner_binary_and_manifests_removed"])
        self.assertEqual(
            self.data["result_sha256"],
            "cd191a414ed35235169349e8b3a38929da2f494aa775a54ba7a6b5f52580576d",
        )

    def test_exact_family_ingestion_and_exclusions(self) -> None:
        scope = self.data["scope"]
        self.assertEqual(scope["included_families"], INCLUDED)
        self.assertEqual(
            scope["included_family_anchor_manifest_sha256"],
            "8909d91ee1d6caa9ea42cf7d5dd23e31f928bc225cbd14c1ba3af0bc2999037c",
        )
        self.assertEqual(
            scope["excluded_rank14_panel_families"],
            [
                [0, 17, 142, 145, 162, 200],
                [0, 25, 57, 104, 116, 148],
                [0, 7, 121, 128, 183, 194],
            ],
        )
        prior = [row["roots"] for row in scope["excluded_prior_rational_scan_families"]]
        self.assertEqual(
            prior,
            [
                [0, 4, 30, 31, 39, 46],
                [0, 6, 49, 73, 82, 96],
                [0, 17, 142, 145, 162, 200],
                [0, 25, 57, 104, 116, 148],
                [0, 7, 121, 128, 183, 194],
            ],
        )
        self.assertEqual(
            scope["excluded_family_manifest_sha256"],
            "4f5a5d7945e91b13dcc110020ceecc11d0371d4d5922868e46311e352c7f378d",
        )
        self.assertEqual(scope["fixed_panel_parameters_excluded"], list("12345678"))
        self.assertTrue(scope["included_and_prior_rational_families_disjoint"])
        self.assertIn("no T", scope["fixed_fiber_policy"])
        self.assertEqual(len(self.data["frozen_calibrations"]), 13)
        self.assertTrue(
            all(row["source_stable_numerical_rank"] == 13 for row in self.data["frozen_calibrations"])
        )
        self.assertTrue(
            all(row["numerical_rank_is_not_an_independence_certificate"] for row in self.data["frozen_calibrations"])
        )
        for family in self.data["families"]:
            self.assertTrue(all(value == 0 for value in family["A_coefficients_ascending"][1::2]))
            self.assertTrue(all(value == 0 for value in family["B_coefficients_ascending"][1::2]))
            self.assertEqual(len(family["content_free_discriminant_coefficients_ascending"]), 21)

    def test_closed_population_fresh_bands_features_and_selection(self) -> None:
        modular = self.data["modular_scan"]
        self.assertTrue(modular["bands_disjoint"])
        self.assertTrue(modular["fresh_relative_to_all_prior_rational_scanners_through_prime_809"])
        self.assertEqual(modular["discovery_primes"], [811, 821, 823, 827, 829, 839, 853, 857])
        self.assertEqual(modular["held_primes"], [859, 863, 877, 881, 883, 887, 907, 911])
        box = modular["common_box"]
        self.assertEqual(
            (
                box["primitive_positive_rationals"], box["fixed_panel_excluded"],
                box["evaluated_per_family"], box["evaluated_family_parameter_pairs"],
            ),
            (637_913, 8, 637_905, 8_292_765),
        )
        self.assertEqual(
            box["ordered_parameter_manifest_sha256"],
            "9acfbbd212a6d5fd6da96c0e65e4a28bfe280e7cb5629ee4a09a05face7fb235",
        )
        self.assertEqual(len(modular["family_scans"]), 13)
        self.assertTrue(all(row["keep"] == 512 for row in modular["family_scans"]))
        self.assertTrue(all(row["calibration"]["exact_python_replay"] for row in modular["family_scans"]))

        features = self.data["exact_discriminant_feature_screen"]
        self.assertEqual(features["content_free_homogeneous_degree"], 20)
        self.assertEqual(features["trial_division_prime_bound"], 997)
        self.assertEqual(
            [row["admissible_feature_pool_count"] for row in features["pool_audits"].values()],
            [512] * 13,
        )
        self.assertTrue(
            all(row["exact_singular_rejections"] == 0 for row in features["pool_audits"].values())
        )

        selection = self.data["conductor_selection"]
        self.assertTrue(selection["discovery_population_closed_before_held_scores"])
        self.assertTrue(selection["equal_per_family_rule"])
        self.assertFalse(selection["selection_uses_conductor_data"])
        self.assertFalse(selection["selection_uses_point_or_numerical_rank_data"])
        self.assertEqual(selection["selected_population"], 87)
        self.assertEqual(
            selection["selected_population_sha256"],
            "6a392e6bd2a693ed2972f6a269cf7a8fc657ecea53380e0ba29cfab13eabae75",
        )
        self.assertEqual(len(self.records), 87)

    def test_conductor_first_boundary_and_exact_subtarget_counts(self) -> None:
        conductor = self.data["conductor_first_screen"]
        self.assertTrue(conductor["population_closed_before_any_point_or_rank_call"])
        self.assertTrue(conductor["all_selected_received_exact_conductor_attempt"])
        self.assertEqual(
            (
                conductor["completed"], conductor["timeouts"], conductor["errors"],
                conductor["strict_subtarget"], conductor["plausible_below_190"],
            ),
            (80, 7, 0, 55, 61),
        )
        self.assertEqual(
            Counter(row["conductor_phase"]["status"] for row in self.records),
            Counter(
                {
                    "completed exact PARI minimal-model/conductor computation": 80,
                    "timeout": 7,
                }
            ),
        )
        completed = [
            row for row in self.records
            if row["conductor_phase"]["status"].startswith("completed")
        ]
        leader = min(completed, key=lambda row: float(row["conductor_phase"]["log_conductor"]))
        self.assertEqual((leader["family_label"], leader["parameter"]), (
            "r0_23_93_128_133_175", "3535/52"
        ))
        self.assertAlmostEqual(float(leader["conductor_phase"]["log_conductor"]), 101.0114286349402)

    def test_fixed_point_tiers_stop_at_rank15_without_certificate(self) -> None:
        protocol = self.data["point_search_protocol"]
        self.assertEqual(
            [(row["name"], row["attempted"]) for row in protocol["stages"]],
            [("H5000", 80), ("H50000", 26), ("H250000", 13), ("H1000000", 13)],
        )
        self.assertEqual(protocol["completed_stage_calls"], 131)
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 15)
        self.assertEqual(protocol["finite_reduction_attempts"], [])
        self.assertEqual(protocol["finite_reduction_trigger_stable_rank"], 16)
        self.assertTrue(protocol["stop_rule_fired"])
        self.assertEqual(protocol["broadening_calls_after_fixed_protocol"], 0)

        expected = {
            "H5000": Counter({11: 64, 12: 9, 13: 5, 14: 2}),
            "H50000": Counter({11: 12, 12: 6, 13: 4, 14: 4}),
            "H250000": Counter({14: 7, 12: 3, 11: 2, 13: 1}),
            "H1000000": Counter({14: 6, 12: 4, 13: 1, 15: 1}),
        }
        for stage_name, distribution in expected.items():
            observed = Counter(
                row["point_stages"][stage_name]["stable_numerical_rank"]
                for row in self.records
                if row.get("point_stages", {}).get(stage_name, {}).get("status") == "completed"
            )
            self.assertEqual(observed, distribution)

        leaders = []
        for row in self.records:
            stage = row.get("point_stages", {}).get("H1000000", {})
            if stage.get("stable_numerical_rank") == 15:
                leaders.append((row, stage))
        self.assertEqual(len(leaders), 1)
        leader, stage = leaders[0]
        self.assertEqual((leader["family_label"], leader["parameter"]), (
            "r0_7_93_154_161_191", "2731/36"
        ))
        self.assertEqual(
            leader["conductor_phase"]["conductor"],
            "35180263184668233005022967592240992011410114878725318553939716560010",
        )
        self.assertEqual(
            stage["pool_point_sha256"],
            "e1e093d78eacd51713208287158d5bf703df0ea015cadab6a6196444b8355ceb",
        )
        self.assertTrue(stage["numerical_rank_is_not_an_independence_certificate"])
        self.assertEqual(stage["finite_reduction_attempt"]["status"], "not triggered")
        self.assertEqual(self.data["target"]["hits"], [])


if __name__ == "__main__":
    unittest.main()

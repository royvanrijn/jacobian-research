#!/usr/bin/env python3
"""Focused checks for the closed, disjoint T=490/9 annulus."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/search_mestre_rank15_annulus.py"
SCANNER = ROOT / "elliptic-curves/cas/scan_mestre_rank15_annulus.cpp"
ARTIFACT = (
    ROOT / "artifacts/generated-results/elliptic_mestre_rank15_annulus.json"
)
ANCHOR_CERTIFICATE = (
    ROOT / "artifacts/generated-results/elliptic_mestre_rank15_490_9.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "64958776e7e642f74157b77f5accd270a5871936a860b39012f6dec8023c5b9a"
)
EXPECTED_SCANNER_SHA256 = (
    "545bee442bedf6750f51e36d727e1202a7e3477ac1de9a515f551930e4ba4479"
)
EXPECTED_ARTIFACT_SHA256 = (
    "a75a052a6e47b46d7dcf6e5da53fda96b9b31e6669225884689d2a474b22ecf7"
)
EXPECTED_ANCHOR_SHA256 = (
    "50b2b9c8bd24bcb5533534446af6404f3a9a761b5f33e0e28e04dc572227f950"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank15AnnulusTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.records = cls.data["selected_records"]

    def test_pinned_files_and_clean_execution(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(SCANNER), EXPECTED_SCANNER_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(sha256(ANCHOR_CERTIFICATE), EXPECTED_ANCHOR_SHA256)
        provenance = self.data["provenance"]
        self.assertEqual(provenance["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(provenance["scanner_sha256"], EXPECTED_SCANNER_SHA256)
        self.assertEqual(provenance["owned_processes_remaining"], 0)
        self.assertEqual(provenance["same_stage_retries"], 0)
        self.assertTrue(provenance["temporary_scanner_binary_removed"])
        self.assertEqual(
            self.data["result_sha256"],
            "81ba931ade7a6796566d76a76fa2ee0c9c9ebb70fb9059820683d8e2e5a2c17c",
        )

    def test_exact_disjoint_population_manifest(self) -> None:
        scope = self.data["scope"]
        self.assertEqual(scope["family_roots"], [0, 7, 121, 128, 183, 194])
        self.assertEqual(scope["anchor_parameter"], "490/9")
        self.assertTrue(scope["anchor_excluded"])
        population = scope["population"]
        self.assertEqual(population["denominator_interval"], [4097, 16000])
        self.assertEqual(population["raw_population"], 285_696)
        self.assertEqual(population["nonprimitive_rejections"], 112_027)
        self.assertEqual(population["peer_grid_divisor_rejections"], 5_572)
        self.assertEqual(population["peer_farey_rejections"], 0)
        self.assertEqual(population["evaluated_unique_primitive_parameters"], 168_097)
        self.assertEqual(
            population["ordered_parameter_manifest_sha256"],
            "6c2fd81b95bb0d7c0531bbccaeff6bcb7fecd1353420febab32492dd04aef4f2",
        )
        self.assertEqual(population["peer_farey_manifest_count"], 3_976)
        self.assertEqual(population["peer_grid_divisor_count_within_annulus"], 383)
        self.assertIn("denominator<=1000", population["prior_rectangle_exclusion"])

    def test_fresh_local_scan_and_exact_feature_boundary(self) -> None:
        scan = self.data["fresh_local_scan"]
        self.assertTrue(scan["bands_disjoint"])
        self.assertTrue(scan["bands_disjoint_from_companion_search_through_prime_577"])
        self.assertTrue(set(scan["discovery_primes"]).isdisjoint(scan["held_primes"]))
        self.assertGreater(min(scan["discovery_primes"]), 577)
        self.assertEqual(scan["discovery_table_digest"], "17844147381485427251")
        self.assertEqual(scan["held_table_digest"], "3565458856032527775")
        self.assertTrue(scan["Python_table_replay_matches"])
        self.assertEqual(scan["survivor_count"], 8_192)
        self.assertEqual(
            scan["survivor_population_sha256"],
            "338c59a444303ba17ae489e5196f1b57942d8e74a54ba47805b1d21e6c5de860",
        )
        features = self.data["exact_discriminant_features"]
        self.assertEqual(features["survivor_count"], 8_192)
        self.assertEqual(features["singular_rejections"], 0)
        self.assertEqual(features["admissible_count"], 8_192)
        self.assertEqual(features["trial_division_prime_bound"], 997)
        self.assertEqual(
            features["feature_population_sha256"],
            "ab481a4c0b657f33f25dda2ddde0d805b88724d567c450078dbee8edf4ebe2ce",
        )

    def test_leakage_controlled_conductor_first_selection(self) -> None:
        selection = self.data["conductor_selection"]
        self.assertTrue(selection["discovery_survivors_closed_before_held_scores"])
        self.assertTrue(selection["exact_discriminant_features_use_no_conductor_or_rank_data"])
        self.assertFalse(selection["selection_uses_conductor"])
        self.assertFalse(selection["selection_uses_point_or_rank_data"])
        self.assertEqual(selection["selected_population"], 86)
        self.assertEqual(
            selection["selected_population_sha256"],
            "dcd570ed35d047cc92da21bd3c89d58239a969d4f3e3f239ca4297b30b44deff",
        )
        self.assertEqual(len(self.records), 86)
        conductor = self.data["conductor_first_screen"]
        self.assertTrue(conductor["population_closed_before_any_conductor_or_point_call"])
        self.assertEqual(
            (conductor["completed"], conductor["timeouts"], conductor["errors"]),
            (45, 41, 0),
        )
        self.assertEqual(conductor["subtarget"], 0)
        self.assertEqual(
            Counter(record["conductor_phase"]["status"] for record in self.records),
            Counter({
                "completed exact PARI minimal-model/conductor computation": 45,
                "timeout": 41,
            }),
        )
        completed = [
            record for record in self.records
            if record["conductor_phase"]["status"].startswith("completed")
        ]
        self.assertGreater(
            min(float(record["conductor_phase"]["log_conductor"]) for record in completed),
            182.72,
        )

    def test_all_fixed_point_tiers_remain_rank11(self) -> None:
        protocol = self.data["point_search_protocol"]
        self.assertEqual(
            [(row["name"], row["attempted"]) for row in protocol["stages"]],
            [("H5000", 45), ("H50000", 24), ("H250000", 6), ("H1000000", 2)],
        )
        self.assertEqual(protocol["completed_stage_calls"], 77)
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 11)
        self.assertEqual(protocol["finite_reduction_trigger_stable_rank"], 15)
        self.assertEqual(protocol["finite_reduction_certificates"], [])
        self.assertEqual(protocol["same_height_retries"], 0)
        self.assertEqual(protocol["broadening_calls_after_fixed_protocol"], 0)
        for stage_name, expected_count in (
            ("H5000", 45), ("H50000", 24),
            ("H250000", 6), ("H1000000", 2),
        ):
            stages = [
                record["point_stages"][stage_name]
                for record in self.records
                if stage_name in record.get("point_stages", {})
            ]
            self.assertEqual(len(stages), expected_count)
            self.assertTrue(all(stage["status"] == "completed" for stage in stages))
            self.assertEqual(Counter(stage["stable_numerical_rank"] for stage in stages), Counter({11: expected_count}))
        self.assertEqual(self.data["target"]["hits"], [])
        self.assertEqual(
            self.data["status"], "completed fixed annulus; stopped without broadening"
        )


if __name__ == "__main__":
    unittest.main()

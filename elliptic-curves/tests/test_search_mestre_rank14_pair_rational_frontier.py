#!/usr/bin/env python3
"""Focused checks for the closed two-family rational frontier."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/search_mestre_rank14_pair_rational_frontier.py"
SCANNER = ROOT / "elliptic-curves/cas/scan_mestre_rank14_pair.cpp"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elliptic_mestre_rank14_pair_rational_frontier.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "2f6251c67e2eb3cee2eca37d7e866913e9d5de73d30e3bfcb253641454d40d5f"
)
EXPECTED_SCANNER_SHA256 = (
    "85769ce09a991f974f271f4e3913dbc34186a0ca3ba8ac38f84cad5a53b08330"
)
EXPECTED_ARTIFACT_SHA256 = (
    "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank14PairRationalFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.records = cls.data["selected_records"]

    def test_pinned_files_and_clean_execution(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(SCANNER), EXPECTED_SCANNER_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        provenance = self.data["provenance"]
        self.assertEqual(provenance["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(provenance["scanner_sha256"], EXPECTED_SCANNER_SHA256)
        self.assertEqual(provenance["owned_processes_remaining"], 0)
        self.assertEqual(provenance["same_stage_retries"], 0)
        self.assertTrue(provenance["temporary_scanner_binary_removed"])
        self.assertEqual(
            self.data["result_sha256"],
            "d33c0cf0a5e2364bd49e18363e0a1f3ca51512fdf60be7052dd05f1cbfa9d610",
        )

    def test_exact_scope_symmetry_and_calibrations(self) -> None:
        scope = self.data["scope"]
        self.assertEqual(
            scope["included_families"],
            [[0, 17, 142, 145, 162, 200], [0, 7, 121, 128, 183, 194]],
        )
        self.assertEqual(
            scope["explicitly_excluded_family"], [0, 25, 57, 104, 116, 148]
        )
        self.assertIn("T and -T are identical", scope["T_sign_quotient"])
        self.assertEqual(scope["prior_panel_parameters_excluded"], list("12345678"))
        self.assertEqual(
            [family["calibration_parameter"] for family in self.data["families"]],
            ["7", "1"],
        )
        self.assertTrue(
            all("A(-T)=A(T)" in family["exact_symmetry"] for family in self.data["families"])
        )
        calibrations = self.data["frozen_calibrations"]
        self.assertEqual(
            [row["certified_algebraic_rank_lower_bound"] for row in calibrations],
            [14, 14],
        )
        self.assertEqual(
            [
                row["finite_reduction_certificate"]["combined_exact_rank_over_F3"]
                for row in calibrations
            ],
            [14, 14],
        )
        self.assertTrue(all(row["excluded_from_every_scanner_stratum"] for row in calibrations))

    def test_closed_fresh_prime_populations_and_exact_features(self) -> None:
        modular = self.data["modular_scan"]
        self.assertTrue(modular["bands_disjoint"])
        self.assertTrue(modular["fresh_relative_to_prior_broad_scanners_through_prime_199"])
        self.assertTrue(set(modular["discovery_primes"]).isdisjoint(modular["held_primes"]))
        self.assertEqual(
            [box["evaluated"] for box in modular["global_boxes"]],
            [18_244_811, 18_244_811],
        )
        expected_retained = [
            "41ac917ac623d8c70e188242d7cb039613f1470b5d4cde55c3349d407cf5cc21",
            "8bb914deaf352ee0e68f2cd09c8c369a90b6c3c813969f1b151533ec70acaad6",
            "20c574263b310797246a62dc3d7aa13c7480c0448d4229467fefa59fbba8b1b3",
            "1e71d3b6f1e7898bc83814785e77fe995c33098e1ff8108bc6b87656895efbd8",
            "c0378c6a4a406381762455ea088333051bafcdf9cce47d86e7ef4297ce232bfb",
            "d1fb923862c18247969ba933dca3bedd8fd550f37edbb712996cea46091134e4",
        ]
        self.assertEqual(
            [row["retained_candidate_sha256"] for row in modular["strata"]],
            expected_retained,
        )
        features = self.data["exact_discriminant_feature_screen"]
        self.assertEqual(features["content_free_homogeneous_degree"], 20)
        self.assertEqual(features["trial_division_prime_bound"], 997)
        self.assertEqual(
            [row["admissible_feature_pool_count"] for row in features["pool_audits"].values()],
            [6066, 6125],
        )
        self.assertTrue(
            all(row["exact_singular_rejections"] == 0 for row in features["pool_audits"].values())
        )

    def test_leakage_controlled_selection_and_conductor_first_boundary(self) -> None:
        selection = self.data["conductor_selection"]
        self.assertTrue(selection["discovery_survivors_closed_before_held_scores"])
        self.assertTrue(selection["held_scores_rank_only_discovery_survivors"])
        self.assertFalse(selection["selection_uses_conductor"])
        self.assertFalse(selection["selection_uses_point_or_rank_data"])
        self.assertEqual(selection["selected_population"], 170)
        self.assertEqual(selection["selected_per_family"], {
            "r0_17_142_145_162_200": 83,
            "r0_7_121_128_183_194": 87,
        })
        self.assertEqual(
            selection["selected_population_sha256"],
            "98bf682019d51e53117639016d03bc3aaaea2626381509958f0efdec431e88df",
        )
        self.assertEqual(len(self.records), 170)
        conductor = self.data["conductor_first_screen"]
        self.assertTrue(conductor["population_closed_before_any_point_or_rank_call"])
        self.assertEqual(
            (conductor["completed"], conductor["timeouts"], conductor["errors"]),
            (144, 26, 0),
        )
        self.assertEqual(conductor["subtarget"], 97)
        self.assertEqual(
            Counter(row["conductor_phase"]["status"] for row in self.records),
            Counter({
                "completed exact PARI minimal-model/conductor computation": 144,
                "timeout": 26,
            }),
        )

    def test_fixed_point_tiers_stop_at_exact_rank15_leader(self) -> None:
        protocol = self.data["point_search_protocol"]
        self.assertEqual(
            [(row["name"], row["attempted"]) for row in protocol["stages"]],
            [("H5000", 144), ("H50000", 32), ("H250000", 8), ("H1000000", 2)],
        )
        self.assertEqual(protocol["completed_stage_calls"], 186)
        self.assertEqual(protocol["same_height_retries"], 0)
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 15)
        self.assertEqual(protocol["finite_reduction_attempts"], [])
        self.assertTrue(protocol["stop_rule_fired"])
        self.assertEqual(protocol["broadening_calls_after_fixed_protocol"], 0)

        distributions = {}
        for stage in ("H5000", "H50000", "H250000", "H1000000"):
            distributions[stage] = Counter(
                row["point_stages"][stage]["stable_numerical_rank"]
                for row in self.records
                if stage in row.get("point_stages", {})
            )
        self.assertEqual(distributions["H5000"], Counter({11: 125, 12: 10, 13: 8, 14: 1}))
        self.assertEqual(distributions["H50000"], Counter({11: 13, 12: 8, 13: 5, 14: 5, 15: 1}))
        self.assertEqual(distributions["H250000"], Counter({14: 6, 13: 1, 15: 1}))
        self.assertEqual(distributions["H1000000"], Counter({14: 1, 15: 1}))

        leader = next(
            row
            for row in self.records
            if row["family_index"] == 1
            and row["numerator"] == 490
            and row["denominator"] == 9
        )
        self.assertEqual(
            [leader["point_stages"][stage]["stable_numerical_rank"] for stage in (
                "H5000", "H50000", "H250000", "H1000000"
            )],
            [12, 15, 15, 15],
        )
        self.assertEqual(
            leader["point_stages"]["H1000000"]["pool_point_sha256"],
            "2d920def92ff94fb878808f93fa39bee417e8e14b4cfd7667f946c5f28f59e85",
        )
        self.assertEqual(
            leader["conductor_phase"]["conductor"],
            "1468617013201344525723305189723172651230013579564220710",
        )
        self.assertLess(float(leader["conductor_phase"]["log_conductor"]), 182.72)
        self.assertEqual(self.data["target"]["hits"], [])


if __name__ == "__main__":
    unittest.main()

import hashlib
import json
from pathlib import Path
import unittest

import elkies_rank28


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"


class RecordHeightLatticeComparisonTests(unittest.TestCase):
    def test_rank28_public_points_are_exact(self):
        self.assertEqual(len(elkies_rank28.POINTS), 28)
        self.assertTrue(all(map(elkies_rank28.on_curve, elkies_rank28.POINTS)))

    def test_high_precision_height_artifact(self):
        path = ARTIFACTS / "record_height_lattices_28_29_273_302_v1.json"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "aa89dea9eb7cf547633a67b522bd9bef67868b0fe686942853953f0258a1472f",
        )
        payload = json.loads(path.read_text())
        self.assertEqual(payload["decimal_precision_digits"], 100)
        self.assertEqual(
            [record["rank_lower_bound"] for record in payload["curves"]],
            [28, 29, 30, 31],
        )
        for record in payload["curves"]:
            rank = record["rank_lower_bound"]
            self.assertEqual(len(record["height_gram"]), rank)
            self.assertTrue(all(len(row) == rank for row in record["height_gram"]))
            self.assertEqual(len(record["lll_transform_columns"]), rank)
            self.assertEqual(len(record["lll_reduced_gram"]), rank)

    def test_candidate_core_artifact_and_common_shell_profile(self):
        path = ARTIFACTS / "record_rank17_core_candidates_v1.json"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "0e1f32f73bb80033aeb1e2bee55402e685f0bb9a9e4e84c53e3a8026791c2d55",
        )
        payload = json.loads(path.read_text())
        self.assertEqual(payload["target_lattice"]["unoriented_minimal_shell_size"], 1311)
        self.assertEqual(
            [record["candidate_core_lines_at_bound"] for record in payload["curves"]],
            [787, 485, 971, 864],
        )
        profiles = [
            list(map(float, record["approximate_1311_shell_quantiles_over_lambda"]))
            for record in payload["curves"]
        ]
        # Exclude the unusually low first vector.  All other declared nearest-
        # rank quantiles agree across the four curves to within 0.10 lambda.
        for column in range(1, 7):
            values = [profile[column] for profile in profiles]
            self.assertLess(max(values) - min(values), 0.10)

    def test_curve302_bounded_mestre_result(self):
        path = ARTIFACTS / "icarm_construction_fingerprints_v2.json"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "ced63ed67c61bb23484039237259127ffd0864426ae41429cd005e6989bfdc4a",
        )
        payload = json.loads(path.read_text())
        records = {
            record["curve_id"]: record
            for record in payload["six_root_mestre_recognition"]["targets"]
        }
        self.assertEqual(records[302]["families_tested"], 2330)
        self.assertEqual(records[302]["exact_j_matches"], [])
        self.assertEqual(records[302]["survivors_with_rational_square_parameter"], [])

    def test_rank17_fingerprint_fails_known_negative_control(self):
        path = ARTIFACTS / "record_rank17_fingerprint_calibration_v1.json"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "6294fe759b7c03eb908afea12546c896d7137566d884c38caa13541de42ee763",
        )
        payload = json.loads(path.read_text())
        self.assertEqual(payload["discrimination_result"], "FAILS_NEGATIVE_CONTROL")
        self.assertLess(float(payload["control_rms_over_record_mean_rms"]), 1.20)
        truth = payload["negative_control"]["exact_true_generic_subgroup"]
        self.assertEqual(truth["true_generic_rank"], 12)
        self.assertEqual(truth["unique_relation_in_quartic_point_order"], [1] * 12 + [0])
        self.assertTrue(truth["exact_group_law_replay"])
        self.assertEqual(truth["intersection_rank_with_forced_rank17_candidate"], 9)
        for fit in payload["exact_r17_fits"]:
            shell = fit["mapped_full_minimal_shell"]
            self.assertEqual(shell["unoriented_lines"], 1311)
            self.assertGreater(
                float(shell["normalized_rms_deviation_from_4"]),
                2.0,
            )
        pairwise = payload["direct_pairwise_numerical_core_fits"]
        self.assertEqual(pairwise["discrimination_result"], "FAILS_NEGATIVE_CONTROL")
        self.assertLess(
            float(pairwise["curve245_to_302_rms_over_record_source_mean_rms"]),
            1.20,
        )

        fingerprints = {
            record["label"]: record
            for record in payload["out_of_sample_integrality_fingerprint"]
        }
        control_odds = float(
            fingerprints["curve245-negative-control"]["core_over_outside_odds_ratio"]
        )
        for label in ("rank28", "rank29", "curve273", "curve302"):
            self.assertGreater(
                float(fingerprints[label]["core_over_outside_odds_ratio"]),
                control_odds,
            )

    def test_shell_aware_search_also_fails_negative_control(self):
        path = ARTIFACTS / "record_rank17_shell_embedding_search_v1.json"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            "34d9d287d37ee07a7caf862f0db53785a3f197164b3862fb8bb24fe1152f48b3",
        )
        payload = json.loads(path.read_text())
        self.assertEqual(
            payload["calibration"]["discrimination_result"],
            "FAILS_NEGATIVE_CONTROL",
        )
        self.assertLess(
            float(payload["calibration"]["control_cv_over_rank29_cv"]),
            1.20,
        )
        self.assertEqual(
            {record["label"] for record in payload["results"]},
            {"rank28", "rank29", "curve273", "curve302", "curve245-negative-control"},
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "elkies-k3/scripts/build_mw17_jump_v2_zero_gain_rescue.py"
RUNNER = ROOT / "elliptic-curves/cas/run_mw17_jump_v2_zero_gain_rescue.sage"
CAMPAIGN = ROOT / "artifacts/generated-results/elkies-k3-mw17-jump-v2-campaign-v1.json"
PROTOCOL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-mw17-jump-v2-zero-gain-rescue-arm-v1.json"
)


class MW17JumpV2ZeroGainRescueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.campaign = json.loads(CAMPAIGN.read_text())
        cls.protocol = json.loads(PROTOCOL.read_text())

    def test_builder_replays_byte_for_byte(self) -> None:
        module = SourceFileLoader("mw17_jump_v2_rescue_builder_test", str(BUILDER)).load_module()
        self.assertEqual(module.build(), self.protocol)

    def test_assignment_preserves_every_frozen_candidate_in_order(self) -> None:
        assignments = self.protocol["assignments"]
        self.assertEqual(len(assignments), 2_239)
        self.assertEqual(
            [row["sample_id"] for row in assignments],
            [row["sample_id"] for row in self.campaign["rows"]],
        )
        self.assertEqual(
            [row["campaign_index"] for row in assignments], list(range(2_239))
        )
        self.assertTrue(
            self.protocol["source_campaign"]["rows_and_order_reused_without_change"]
        )

    def test_one_of_eight_hash_arm_is_outcome_blind(self) -> None:
        rule = self.protocol["assignment"]
        self.assertEqual(rule["namespace"], "mw17-jump-v2-zero-gain-rescue-v1")
        self.assertEqual(rule["denominator"], 8)
        self.assertEqual(rule["rescue_bucket"], 0)
        self.assertEqual(rule["prospective_treatment_probability"], "1/8")
        self.assertFalse(rule["uses_detector_outcomes"])
        for row in self.protocol["assignments"]:
            expected_hash = sha256(
                f"{rule['namespace']}\0{row['sample_id']}".encode()
            ).hexdigest()
            self.assertEqual(row["assignment_sha256"], expected_hash)
            self.assertEqual(row["assignment_bucket_mod_8"], int(expected_hash, 16) % 8)
            self.assertEqual(
                row["assigned_to_rescue_arm"], row["assignment_bucket_mod_8"] == 0
            )
        treated = [
            row for row in self.protocol["assignments"]
            if row["assigned_to_rescue_arm"]
        ]
        self.assertEqual(rule["assigned_candidate_count"], len(treated))
        self.assertEqual(rule["assigned_candidate_count"], 264)
        self.assertEqual(
            rule["assigned_counts_by_family"],
            {
                "074d9": 109,
                "07ca9": 42,
                "08234": 21,
                "08f72": 36,
                "0e80b": 30,
                "11952": 26,
            },
        )
        self.assertEqual(
            rule["assigned_counts_by_family"],
            dict(sorted(Counter(row["family"] for row in treated).items())),
        )

    def test_rescue_is_disjoint_generic_then_budget_preserving_adaptive(self) -> None:
        detector = self.protocol["rescue_detector"]
        self.assertEqual(detector["generic_rescue_class_ranks_one_based"], [44, 344])
        self.assertEqual(detector["generic_rescue_batch_count"], 7)
        self.assertEqual(detector["generic_rescue_charts_per_batch"], 43)
        self.assertEqual(detector["additional_budget_chart_count"], 301)
        self.assertEqual(detector["maximum_total_chart_count_including_base"], 344)
        self.assertIn("existing adaptive quotient policy", detector["switch_rule"])
        self.assertEqual(detector["ranking_field"], "actual_certified_quotient_rank_gain")

    def test_runner_enforces_clean_zero_and_no_descent_prerequisite(self) -> None:
        source = RUNNER.read_text()
        self.assertIn("generic_rows[BASE_INITIAL_CHARTS:TOTAL_CHART_CAP]", source)
        self.assertIn("range(0, RESCUE_CHARTS, RESCUE_BATCH_SIZE)", source)
        self.assertIn('"NOT_REQUIRED_FOR_BOUNDED_RESCUE"', source)
        self.assertIn("clean_zero_eligible", source)
        policy = self.protocol["proof_and_budget_policy"]
        self.assertFalse(policy["complete_descent_required_to_run_this_bounded_rescue"])
        self.assertFalse(policy["small_field_100_row_laboratory_changed"])


if __name__ == "__main__":
    unittest.main()

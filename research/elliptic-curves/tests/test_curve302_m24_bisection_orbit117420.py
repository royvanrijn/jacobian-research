from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_m24_bisection_orbit117420_v1.json"


class Curve302M24BisectionOrbit117420Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_exact_orbit_and_completed_box(self) -> None:
        data = self.data
        centre = data["exact_half_lattice_degree_two_reconstruction"]
        search = data["completed_search"]
        self.assertEqual(data["status"], "PASS_CERTIFIED_RANK_AT_LEAST_25")
        self.assertEqual(centre["orbit_mask"], 117420)
        self.assertEqual(centre["minimum_generic_MW17_norm"], 10)
        self.assertEqual(centre["M24_parity_mask"], 95909)
        self.assertEqual(len(centre["representative"]), 24)
        self.assertEqual(search["charts"], 1)
        self.assertEqual(search["height"], 125000)
        self.assertEqual(search["search_status"], "bounded_search_complete")
        self.assertEqual(search["mod_2_certified_rank_lower_bound"], 25)
        self.assertEqual(search["mod_3_5_ranks"], {"3": 25, "5": 25})

    def test_execution_blindness_and_historical_exclusion(self) -> None:
        data = self.data
        exclusion = data["why_sealed_17_to_24_policy_excluded_it"]
        audit = data["execution_oracle_audit"]
        self.assertFalse(exclusion["initial_M17_deep_sample_contains_parity"])
        self.assertFalse(exclusion["initial_M17_selected_centres_contain_parity"])
        self.assertEqual(exclusion["padded_parity_shift_right_22"], 0)
        self.assertTrue(audit["rule_is_retrospectively_calibrated"])
        self.assertTrue(audit["execution_is_target_blind"])
        rejected = " ".join(audit["excluded_from_execution"]).lower()
        self.assertIn("residual-visibility diagnostic", rejected)
        for path in (*audit["geometry_reads"], *audit["worker_reads"], *audit["replay_reads"]):
            self.assertNotIn("residual_visibility", path)
            self.assertNotIn("public-span", path)


if __name__ == "__main__":
    unittest.main()

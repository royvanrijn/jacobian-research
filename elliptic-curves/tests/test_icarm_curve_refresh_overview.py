from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SWEEP = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
OVERVIEW = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve_refresh_475_573_overview_v1.json"
)
QUOTIENTS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class IcarmCurveRefreshOverviewTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.sweep = json.loads(SWEEP.read_text())
        cls.overview = json.loads(OVERVIEW.read_text())
        cls.quotients = json.loads(QUOTIENTS.read_text())

    def test_pinned_artifact_hashes(self) -> None:
        self.assertEqual(
            sha256(SWEEP),
            "77a3c051111e7ead5ee2a6f88df4a975c2f5bdb87be1bfe4d88b195f293da50c",
        )
        self.assertEqual(
            sha256(OVERVIEW),
            "1db137c4c006f774ad653b41b8c04ecc7b332d1104905dcb3f5eb23732904e3c",
        )
        self.assertEqual(
            sha256(QUOTIENTS),
            "e0c0e62c4c357a0dffd4c55dcae0ec2b4993a47fcea0a0de9b055b78a6580081",
        )

    def test_complete_573_curve_sweep(self) -> None:
        outcome = self.sweep["outcome"]
        self.assertEqual(outcome["pinned_curve_count"], 573)
        self.assertEqual(outcome["distinct_target_j_invariant_count"], 561)
        self.assertEqual(outcome["class_preimage_decision_count"], 3438)
        self.assertEqual(outcome["class_rational_hit_count"], 86)
        self.assertEqual(outcome["class_miss_count"], 3352)
        self.assertEqual(
            outcome["native_chart_twist_counts"], {"QQ_ISOMORPHIC_UNTWISTED": 479}
        )

    def test_appended_curve_overview_is_exact_and_complete(self) -> None:
        summary = self.overview["summary"]
        self.assertEqual(summary["new_curve_id_interval"], [475, 573])
        self.assertEqual(summary["new_curve_count"], 99)
        self.assertEqual(summary["new_norm12_atlas_hit_count"], 17)
        self.assertEqual(summary["new_norm12_atlas_miss_count"], 82)
        self.assertTrue(summary["all_displayed_points_checked_on_curve"])
        self.assertTrue(summary["all_stored_discriminants_checked"])
        self.assertEqual(summary["priority_independence_not_closed_curve_ids"], [])
        self.assertEqual(summary["priority_curves_requiring_mod3_fallback"], [542])
        self.assertEqual(summary["priority_non_atlas_ids"], [542, 548])
        self.assertEqual(summary["curves_missing_public_conductor"], [537, 543, 545, 568])

    def test_atlas_hit_specialization_audit(self) -> None:
        summary = self.quotients["summary"]
        self.assertEqual(summary["new_atlas_hit_count"], 17)
        self.assertEqual(summary["quotient_curve_count"], 16)
        self.assertEqual(summary["priority_curve_count"], 11)
        self.assertTrue(summary["all_generic_specializations_have_rank_17"])
        self.assertTrue(
            summary["all_quotient_generic_subgroups_primitive_in_displayed_subgroups"]
        )
        self.assertEqual(summary["noninclusive_displayed_subgroup_curve_ids"], [499])
        self.assertEqual(
            summary["displayed_exceptional_quotients"],
            {
                "531": "Z^11",
                "534": "Z^11",
                "535": "Z^11",
                "536": "Z^11",
                "537": "Z^10",
                "540": "Z^8",
                "541": "Z^8",
                "543": "Z^12",
                "544": "Z^11",
                "545": "Z^11",
                "546": "Z^8",
                "498": "Z^6",
                "539": "Z^6",
                "538": "Z^5",
                "478": "Z^4",
                "532": "Z^3",
            },
        )
        obstruction = self.quotients["noninclusive_fibres"][0]
        self.assertEqual(obstruction["curve_id"], 499)
        self.assertFalse(
            obstruction["displayed_subgroup_contains_specialized_generic_subgroup"]
        )
        self.assertEqual(
            obstruction["overgroup_generated_by_displayed_and_generic_modulo_displayed"],
            "Z/3Z",
        )


if __name__ == "__main__":
    unittest.main()

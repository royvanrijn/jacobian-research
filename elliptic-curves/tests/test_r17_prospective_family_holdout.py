from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-prospective-ordinary-family-holdout-v1.json"
)
SCRIPT = (
    ROOT
    / "elkies-k3/scripts/build_r17_norm12_prospective_family_holdout.py"
)


class ProspectiveFamilyHoldoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_commitment_is_unopened_and_label_blind(self) -> None:
        self.assertEqual(
            self.data["status"],
            "FROZEN_UNOPENED_PROSPECTIVE_WHOLE_FAMILY_HOLDOUT",
        )
        commitment = self.data["commitment"]
        self.assertTrue(commitment["frozen_before_any_search_outcome"])
        self.assertFalse(commitment["selection_used_public_rank_or_hit_labels"])
        self.assertFalse(commitment["historical_search_denominator_known"])
        for row in self.data["rows"]:
            self.assertIsNone(row["features"])
            self.assertIsNone(row["search_outcome"])
            self.assertEqual(row["search_status"], "NOT_OPENED")
            self.assertEqual(row["selection_lane"], "ordinary_counter_hash_draw_no_score")

    def test_balanced_new_height_shell_and_unique_fibres(self) -> None:
        rows = self.data["rows"]
        self.assertEqual(len(rows), 1536)
        self.assertEqual(Counter(row["family"] for row in rows), {
            family: 256 for family in self.data["pgl2_family_split"]["family_order_by_frozen_hash"]
        })
        for row in rows:
            numerator, denominator = row["projective_pair"]
            self.assertGreater(denominator, 0)
            self.assertEqual(gcd(abs(numerator), denominator), 1)
            self.assertEqual(row["projective_height"], max(abs(numerator), denominator))
            self.assertGreaterEqual(row["projective_height"], 30_001)
            self.assertLessEqual(row["projective_height"], 60_000)
        self.assertEqual(len({row["sample_id"] for row in rows}), len(rows))
        self.assertEqual(len({row["j_invariant_sha256"] for row in rows}), len(rows))

    def test_outer_split_holds_out_whole_pgl2_families(self) -> None:
        split = self.data["pgl2_family_split"]
        self.assertEqual(len(split["locked_holdout_families"]), 2)
        self.assertEqual(len(split["prospective_development_families"]), 4)
        family_splits = defaultdict(set)
        for row in self.data["rows"]:
            family_splits[row["family"]].add(row["outer_split"])
        self.assertTrue(all(len(values) == 1 for values in family_splits.values()))
        self.assertEqual(
            {
                family
                for family, values in family_splits.items()
                if values == {"locked_family_holdout"}
            },
            set(split["locked_holdout_families"]),
        )
        member_sets = [set(members) for members in split["member_charts_by_representative"].values()]
        self.assertEqual(sum(map(len, member_sets)), 43)
        self.assertEqual(len(set().union(*member_sets)), 43)

    def test_protocol_has_no_adaptive_allocation_or_negative_theorem(self) -> None:
        protocol = self.data["frozen_search_protocol"]
        self.assertEqual(
            protocol["allocation_rule"],
            "identical stages, bounds, and time limit for every scheduled row",
        )
        point_stage = next(
            row
            for row in protocol["per_parameter_stages"]
            if row["stage"] == "uniform_bounded_point_search"
        )
        self.assertFalse(point_stage["adaptive_depth"])
        self.assertFalse(point_stage["presieve"])
        self.assertIn(
            "not a rank-17",
            " ".join(protocol["no_promotion_rules"]),
        )

    def test_generation_hash_is_current(self) -> None:
        self.assertEqual(
            self.data["generation"]["script_sha256"],
            sha256(SCRIPT.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

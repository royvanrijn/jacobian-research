from __future__ import annotations

from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from pointed_regression_sources import historical_digest


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "elkies-k3/scripts/build_mw17_jump_v2_campaign.py"
CAMPAIGN = ROOT / "artifacts/generated-results/elkies-k3-mw17-jump-v2-campaign-v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class MW17JumpV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.campaign = json.loads(CAMPAIGN.read_text())

    def test_historical_builder_replays_byte_for_byte(self):
        module = SourceFileLoader("mw17_jump_v2_campaign_test", str(BUILDER)).load_module()
        with patch.object(module, "digest", side_effect=historical_digest):
            self.assertEqual(module.build(), self.campaign)

    def test_inventory_and_priority_order(self):
        rows = self.campaign["rows"]
        self.assertEqual(len(rows), 2_239)
        self.assertEqual([row["campaign_index"] for row in rows], list(range(2_239)))
        self.assertEqual(len({row["sample_id"] for row in rows}), 2_239)
        self.assertEqual(
            Counter(row["family"] for row in rows),
            Counter({"07ca9": 256, "08234": 256, "08f72": 256, "11952": 256, "0e80b": 256, "074d9": 959}),
        )
        tranches = []
        for row in rows:
            if not tranches or row["priority_tranche"] != tranches[-1]:
                tranches.append(row["priority_tranche"])
        self.assertEqual(tranches, list(range(9)))
        self.assertEqual({row["family"] for row in rows[:256]}, {"07ca9"})
        self.assertEqual({row["family"] for row in rows[256:512]}, {"08234"})
        self.assertEqual({row["frame_class"] for row in rows[512:1024]}, {"alternate-Q80"})

    def test_h10000_is_exactly_requested_191_row_tranche(self):
        rows = [row for row in self.campaign["rows"] if row["source_population"].startswith("bounded_box")]
        self.assertEqual(len(rows), 191)
        self.assertEqual([row["source_index"] for row in rows], list(range(191)))
        source_fact = self.campaign["selection"]["h10000_source_fact"]
        self.assertEqual(source_fact["scanned_parameter_count"], 121_589_944)
        self.assertEqual(source_fact["stored_ranked_finalist_count"], 1_000)
        self.assertIn("No distinct 191-row source object", source_fact["warning"])

    def test_crt_selection_is_balanced_per_anchor(self):
        expected = {"full": 52, "C": 51, "D": 51, "E": 51, "F": 51}
        for anchor in (356, 385):
            rows = [
                row for row in self.campaign["rows"]
                if row["source_population"] == f"crt_anchor_{anchor}_balanced256"
            ]
            self.assertEqual(len(rows), 256)
            self.assertEqual(Counter(row["selection_diagnostic"]["lane"] for row in rows), Counter(expected))
            self.assertTrue(all(not row["selection_diagnostic"]["used_as_detector_filter"] for row in rows))

    def test_sources_are_hash_pinned_and_initial_gain_is_not_a_filter(self):
        for name, expected in self.campaign["immutability"]["source_file_sha256"].items():
            self.assertEqual(digest(ROOT / name), expected, name)
        detector = self.campaign["detector"]
        self.assertEqual(detector["ranking_field"], "actual_certified_quotient_rank_gain")
        self.assertEqual(detector["rank_lower_bound_formula"], "17 + actual_certified_quotient_rank_gain")
        self.assertIn("hyperellratpoints runs", detector["quartic_preprocessing"])
        self.assertIn("preprocessing, not the bounded point search", detector["quartic_preprocessing"])
        self.assertIn("never a candidate-selection or leaderboard filter", detector["initial_gain_policy"])
        self.assertEqual(detector["global_termination"], "write stop sentinel immediately after any certified gain at least 15")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SCRIPT = ROOT / "elliptic-curves/cas/run_curve385_height_compression_pilot.sage"
SOURCE = ART / "curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz"
PROTOCOL = ART / "curve385_height_compression_pilot_protocol_v1.json"
RESULT = ART / "curve385_height_compression_pilot_blind_v1.json"


class Curve385HeightCompressionPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol = json.loads(PROTOCOL.read_text())
        cls.result = json.loads(RESULT.read_text())
        with gzip.open(SOURCE, "rt") as handle:
            cls.source = json.load(handle)

    def test_protocol_is_bound_to_code_and_frozen_m29(self) -> None:
        key = "elliptic-curves/cas/run_curve385_height_compression_pilot.sage"
        self.assertEqual(
            self.protocol["input_hashes"][key],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(self.protocol["source_state"]["basis_rank"], 29)
        self.assertEqual(
            self.protocol["status"], "BUILT_OUTCOME_FREE_BOUNDED_PILOT_PROTOCOL"
        )

    def test_builder_uses_fresh_stable_full_lattice_holes(self) -> None:
        generation = self.protocol["candidate_generation"]
        self.assertEqual(
            generation["universe"],
            "all 2^29 parity classes of the frozen current M29 basis",
        )
        self.assertEqual(generation["seed_count"], 256)
        self.assertEqual(generation["stable_reduced_pool_size"], 32)
        self.assertTrue(
            all(
                row["representative_stable_between_scales"]
                for row in generation["stable_reduced_pool"]
            )
        )
        old_keys = set(self.source["searched_base_point_keys"])
        new_keys = {row["base_point_key"] for row in generation["stable_reduced_pool"]}
        self.assertFalse(old_keys & new_keys)
        self.assertEqual(len(new_keys), 32)

    def test_selected_order_is_presearch_and_diverse(self) -> None:
        selection = self.protocol["selection"]
        selected = selection["selected_charts"]
        self.assertEqual(selection["search_chart_count"], 16)
        self.assertEqual(len(selected), 16)
        self.assertIsNone(selected[0]["minimum_torus_distance_to_earlier_center"])
        self.assertTrue(
            all(
                row["minimum_torus_distance_to_earlier_center"] > 0
                for row in selected[1:]
            )
        )
        encoded = json.dumps(self.protocol, sort_keys=True)
        self.assertNotIn("finite_curve_points", encoded)
        self.assertNotIn("hyperellratpoints", encoded)

    def test_result_is_bound_to_committed_order(self) -> None:
        protocol_key = (
            "artifacts/generated-results/elliptic-curves/"
            "curve385_height_compression_pilot_protocol_v1.json"
        )
        self.assertEqual(
            self.result["input_hashes"][protocol_key],
            hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.result["protocol_selected_order_sha256"],
            self.protocol["selection"]["selected_order_sha256"],
        )
        self.assertEqual(len(self.result["cover_records"]), 16)

    def test_bounded_pilot_has_no_points_or_group_growth(self) -> None:
        self.assertEqual(
            self.result["status"], "COMPLETE_BOUNDED_NO_GROUP_GROWTH"
        )
        self.assertEqual(self.result["basis_rank_before"], 29)
        self.assertEqual(self.result["basis_rank_after"], 29)
        self.assertEqual(self.result["finite_point_occurrence_count"], 0)
        self.assertEqual(self.result["distinct_returned_point_count"], 0)
        self.assertEqual(self.result["timeout_count"], 0)
        self.assertEqual(self.result["pari_failure_count"], 0)
        self.assertEqual(
            self.result["classification"]["status"],
            "PASS_BASIS_EQUALS_DISCOVERED_GROUP",
        )
        self.assertEqual(self.result["classification"]["events"], [])
        comparison = self.result["historical_old_point_hit_comparison"]
        self.assertEqual(
            comparison["prior_natural_weight_one_two"],
            {"chart_count": 3116, "chart_with_finite_point_count": 665},
        )
        self.assertEqual(
            comparison["height_compression_pilot"],
            {"chart_count": 16, "chart_with_finite_point_count": 0},
        )
        self.assertAlmostEqual(
            comparison["descriptive_fixed_margin_fisher_lower_tail"]["decimal"],
            0.021729819748749414,
        )


if __name__ == "__main__":
    unittest.main()

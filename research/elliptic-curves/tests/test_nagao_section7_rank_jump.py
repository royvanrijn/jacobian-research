#!/usr/bin/env python3
"""Focused checks for the Nagao section-7 quotient fingerprint and replay."""

from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
FINGERPRINT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "nagao_section7_rank_jump_fingerprint_v1.json"
)
REPLAY = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "nagao_section7_rank_jump_replay_v1.json"
)


class NagaoSection7RankJumpTests(unittest.TestCase):
    def test_pinned_quotient_fingerprint(self) -> None:
        row = json.loads(FINGERPRINT.read_text())["fingerprints"][0]
        quotient = row["quotient_structure"]
        self.assertEqual(quotient["free_quotient_rank_lower_bound"], 8)
        self.assertEqual(
            quotient["smith_invariant_factors"], [1] + [2] * 11
        )
        self.assertEqual(
            quotient["tensor_dimensions_over_f_ell"], {"2": 19, "3": 8, "5": 8}
        )
        self.assertEqual(
            quotient["specialized_generic_saturation_index_in_displayed_subgroup"],
            2048,
        )
        visibility = row["degree_visibility"][0]
        self.assertEqual(visibility["visible_free_quotient_span_dimension"], 8)
        self.assertEqual(
            visibility["visible_tensor_quotient_span_dimension_over_f2"], 19
        )
        self.assertEqual(visibility["ambient_candidate_point_count"], 224)
        self.assertEqual(visibility["distinct_quotient_classes_over_f2"], 170)

    def test_complete_replay_is_censored(self) -> None:
        document = json.loads(REPLAY.read_text())
        self.assertEqual(document["population"]["primitive_parameter_count"], 18_244_819)
        self.assertEqual(document["population"]["negative_label_count"], 0)
        scores = document["anchors"][0]["scores"]
        self.assertEqual(scores["training"]["rank_position_one_based"], 9_041_935)
        self.assertEqual(scores["validation"]["rank_position_one_based"], 755_065)


if __name__ == "__main__":
    unittest.main()

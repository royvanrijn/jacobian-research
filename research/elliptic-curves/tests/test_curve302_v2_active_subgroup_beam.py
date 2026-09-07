#!/usr/bin/env python3
"""Regression boundary checks for the target-free V2 branching policy."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
CAS = ROOT / "elliptic-curves/cas"


class Curve302V2ActiveSubgroupBeamTest(unittest.TestCase):
    def setUp(self) -> None:
        self.beam_path = ART / "curve302_v2_active_subgroup_beam_v1.json"
        self.evaluation_path = ART / "curve302_v2_active_subgroup_beam_retrospective_v1.json"
        self.beam = json.loads(self.beam_path.read_text())
        self.evaluation = json.loads(self.evaluation_path.read_text())

    def test_target_free_selector_boundary_and_bounded_beam(self) -> None:
        self.assertEqual(self.beam["status"], "PASS_GEOMETRY_ONLY_BASIS_INVARIANT_BEAM_DESIGN")
        self.assertEqual(self.beam["point_searches_run"], 0)
        self.assertLessEqual(len(self.beam["beam"]), 4)
        branches = self.beam["all_branch_scores_in_selection_order"]
        self.assertEqual({entry["branch_mask"] for entry in branches}, set(range(1, 8)))
        self.assertEqual({entry["increment_rank"] for entry in self.beam["beam"]}, {1, 2, 3})
        self.assertIn("exact minimum", self.beam["score_definition"]["cosets"])
        self.assertIn("multiplicity", " ".join(self.beam["score_definition"]["lexicographic_priority"]))
        self.assertIn("reduced-chart", " ".join(self.beam["score_definition"]["lexicographic_priority"]))
        forbidden = ("exceptional_subgroup_landscape", "residual_visibility", "residual-strict", "curve302_visibility")
        for path in self.beam["prospective_boundary"]["selection_read_paths"]:
            self.assertFalse(any(token in path for token in forbidden), path)

    def test_random_unimodular_invariance_is_complete(self) -> None:
        invariance = self.beam["basis_invariance"]
        checks = invariance["checks"]
        self.assertEqual(invariance["random_unimodular_rebases"], 15)
        self.assertEqual(len(checks), 15)
        self.assertEqual({row["status"] for row in checks}, {"PASS_EXACT_CVP_AND_CHART_KEY_INVARIANCE"})
        self.assertEqual({row["branch_mask"] for row in checks}, set(range(1, 8)))
        self.assertTrue(all(abs(row["determinant"]) == 1 for row in checks))
        self.assertIn("exact congruent integer form", self.beam["score_definition"]["finite_metric_transport"])

    def test_retrospective_evaluation_is_separate_and_metric_bound(self) -> None:
        self.assertEqual(self.evaluation["status"], "PASS_SEPARATE_RETROSPECTIVE_EVALUATION")
        self.assertEqual(self.evaluation["point_searches_run"], 0)
        self.assertEqual(self.evaluation["fixed_diagnostic_reference"]["state_count"], 1 << 14)
        self.assertIn("did not read", self.evaluation["boundary"]["selector"])
        self.assertEqual(len(self.evaluation["branches_in_frozen_selector_order"]), 7)
        self.assertEqual(
            {row["branch_mask"] for row in self.evaluation["selected_beam_retrospective_rows"]},
            {row["branch_mask"] for row in self.beam["beam"]},
        )
        for row in self.evaluation["branches_in_frozen_selector_order"]:
            self.assertIn("frozen_16384_lower_tail_percentile_by_fixed_direction", row)
            self.assertGreaterEqual(row["fixed_directions_improved_vs_v2_M24"], 0)

    def test_source_hashes_bind_the_two_separated_programmes(self) -> None:
        selector = CAS / "design_curve302_v2_active_subgroup_beam.sage"
        evaluator = CAS / "evaluate_curve302_v2_active_subgroup_beam.sage"
        self.assertEqual(self.beam["inputs"][str(selector.relative_to(ROOT))], hashlib.sha256(selector.read_bytes()).hexdigest())
        self.assertEqual(self.evaluation["inputs"][str(evaluator.relative_to(ROOT))], hashlib.sha256(evaluator.read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()

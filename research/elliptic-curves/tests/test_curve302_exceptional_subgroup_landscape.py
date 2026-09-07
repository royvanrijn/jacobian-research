#!/usr/bin/env python3
"""Regression checks for the finite 2^14 curve-302 order diagnostic."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
CAS = ROOT / "elliptic-curves/cas"


class Curve302ExceptionalSubgroupLandscapeTest(unittest.TestCase):
    def setUp(self) -> None:
        self.landscape_path = ART / "curve302_exceptional_subgroup_landscape_v1.json"
        self.report_path = ART / "curve302_exceptional_subgroup_landscape_report_v1.json"
        self.landscape = json.loads(self.landscape_path.read_text())
        self.report = json.loads(self.report_path.read_text())

    def test_complete_cube_and_persistent_atlas(self) -> None:
        self.assertEqual(self.landscape["status"], "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE")
        names = self.landscape["fixed_basis"]["direction_ids"]
        states = self.landscape["subset_states"]
        self.assertEqual(len(names), 14)
        self.assertEqual(len(states), 1 << len(names))
        for mask, row in enumerate(states):
            self.assertEqual(row["state_mask"], mask)
            for direction in range(len(names)):
                if (mask >> direction) & 1:
                    self.assertIsNone(row["stage_local_numerators"][direction])
                    self.assertIsNone(row["retained_numerators"][direction])
                else:
                    self.assertIsNotNone(row["stage_local_numerators"][direction])
                    self.assertIsNotNone(row["retained_numerators"][direction])
        for mask, row in enumerate(states):
            for direction in range(len(names)):
                if (mask >> direction) & 1:
                    continue
                for added in range(len(names)):
                    if added == direction or (mask >> added) & 1:
                        continue
                    child = states[mask | (1 << added)]
                    self.assertLessEqual(child["retained_numerators"][direction], row["retained_numerators"][direction])

    def test_report_binds_core_and_preserves_provenance_boundary(self) -> None:
        self.assertEqual(self.report["status"], "PASS_COMPLETE_FINITE_ORDER_DIAGNOSTIC")
        actual = hashlib.sha256(self.landscape_path.read_bytes()).hexdigest()
        self.assertEqual(self.report["inputs"][str(self.landscape_path.relative_to(ROOT))], actual)
        history = self.report["historical_order_provenance"]
        self.assertEqual(len(history["M24_to_M31_attested_tail"]), 7)
        self.assertEqual(history["all_5040_M17_to_M24_compatible_then_historical_tail"]["minimum_bottleneck_numerator"], history["tail_actual"]["maximum_numerator"])
        audit = self.report["raw_vs_cumulative_monotonicity"]
        self.assertGreater(audit["raw_stage_local_increases"], 0)
        self.assertEqual(audit["retained_increases"], 0)
        self.assertIn("not a prospective selector", self.report["scope"])

    def test_core_source_hash_and_static_artifacts_exist(self) -> None:
        source = CAS / "audit_curve302_exceptional_subgroup_landscape.sage"
        self.assertEqual(
            self.landscape["inputs"][str(source.relative_to(ROOT))],
            hashlib.sha256(source.read_bytes()).hexdigest(),
        )
        for suffix in ("curve302_exceptional_subgroup_landscape_v1.tsv", "curve302_exceptional_subgroup_landscape_v1.svg"):
            path = ART / suffix
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 100)


if __name__ == "__main__":
    unittest.main()

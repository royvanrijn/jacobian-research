from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_residual_visibility_geometry_v1.json"


class Curve302ResidualVisibilityGeometryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_panel_and_half_lattice_bijection(self) -> None:
        data = self.data
        self.assertEqual(data["status"], "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC")
        self.assertEqual(data["degree_two_half_lattice_correspondence"]["cardinality"], 2**17)
        self.assertEqual(data["degree_two_half_lattice_correspondence"]["rational_bisection_orbits"], 40917)
        self.assertEqual(len(data["directions"]), 14)
        self.assertEqual(
            [row["cohort"] for row in data["directions"]],
            ["recovered_local"] * 4 + ["recovered_strict"] * 3 + ["residual_strict"] * 7,
        )

    def test_declared_vetting_boundary(self) -> None:
        data = self.data
        schedule = data["protocol"]["parity_schedule"]
        self.assertEqual(schedule["M17_over_2M17_nonzero_classes"], 2**17 - 1)
        self.assertEqual(schedule["cross_precision_CVP_shortlist_per_target"], 512)
        self.assertEqual(schedule["exact_pointed_quartic_rows_per_target"], 32)
        for row in data["directions"]:
            self.assertEqual(row["all_nonzero_M17_parities_babai_scanned"], 2**17 - 1)
            self.assertEqual(len(row["CVP_rows"]), 512)
            self.assertEqual(len(row["exact_pointed_quartic_rows"]), 32)
            finite = row["minimum_finite_reduced_coordinate_height_in_vetted_set"]
            if finite is not None:
                self.assertGreater(int(finite), 0)
            self.assertGreaterEqual(row["projective_infinity_rows_in_vetted_set"], 0)


if __name__ == "__main__":
    unittest.main()

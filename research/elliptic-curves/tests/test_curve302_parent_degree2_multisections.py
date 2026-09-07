from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_degree2_multisection_lattice_v1.json"
)
ORBITS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_degree2_multisection_orbits_v1.tsv"
)
PILOT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_cheapest_lattice_bisection_v1.json"
)
PILOT_SCRIPT = ROOT / "elliptic-curves/cas/construct_curve302_parent_cheapest_lattice_bisection.sage"


class Curve302ParentDegreeTwoMultisectionsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_complete_translation_quotient(self) -> None:
        data = self.data
        self.assertEqual(data["status"], "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT")
        self.assertTrue(data["enumeration"]["complete"])
        self.assertEqual(data["normalization"]["finite_quotient"], "M/2M")
        self.assertEqual(data["normalization"]["total_translation_orbits"], 2**17)
        self.assertEqual(
            sum(data["enumeration"]["minimum_norm_histogram"].values()), 2**17
        )
        self.assertEqual(
            data["enumeration"]["minimum_norm_histogram"],
            {"0": 1, "4": 1218, "6": 24875, "8": 63922, "10": 40917, "12": 139},
        )

    def test_low_genus_filters_and_export(self) -> None:
        data = self.data
        rational = data["rational_bisections"]
        genus_one = data["genus_one_bisection_candidates"]
        self.assertEqual(rational["translation_orbits"], 40917)
        self.assertEqual(rational["section_nonnegative_threshold"], 10)
        self.assertEqual(genus_one["translation_orbits"], 64061)
        self.assertEqual(genus_one["section_nonnegative_threshold"], 8)
        self.assertEqual(
            len(ORBITS.read_text().splitlines()),
            1 + rational["translation_orbits"] + genus_one["translation_orbits"],
        )
        self.assertEqual(
            hashlib.sha256(ORBITS.read_bytes()).hexdigest(), data["orbits_tsv_sha256"]
        )

    def test_blind_cheapest_equation_control(self) -> None:
        pilot = json.loads(PILOT.read_text())
        self.assertEqual(
            pilot["status"], "PASS_EXACT_NONSPLIT_CHEAPEST_LATTICE_BISECTION"
        )
        self.assertEqual(pilot["selection"]["orbit_mask"], 8044)
        self.assertEqual(pilot["selection"]["minimum_norm"], 10)
        self.assertFalse(pilot["selection"]["selection_uses_t0_or_exceptional_points"])
        self.assertEqual(pilot["RR_construction"]["interpolation_rank"], 19)
        self.assertEqual(pilot["RR_construction"]["generic_residual_factor_degrees"], [2])
        self.assertEqual(pilot["zero_fibre"]["factor_degrees"], [2])
        self.assertFalse(pilot["zero_fibre"]["split_over_Q"])
        self.assertEqual(pilot["zero_fibre"]["quotient_attachment"], "VACUOUS_NONSPLIT")
        self.assertEqual(
            pilot["inputs"]["elliptic-curves/cas/construct_curve302_parent_cheapest_lattice_bisection.sage"],
            hashlib.sha256(PILOT_SCRIPT.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

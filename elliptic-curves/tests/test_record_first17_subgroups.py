import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_first17_subgroups_v1.json"
)


class RecordFirst17SubgroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        cls.curves = {item["label"]: item for item in cls.payload["curves"]}

    def test_exact_lattice_dimensions(self) -> None:
        expected = {"curve273": (30, 10, 13), "curve302": (31, 9, 14)}
        for label, (rank, intersection, quotient) in expected.items():
            record = self.curves[label]
            exact = record["exact_coordinate_lattice"]
            self.assertEqual(record["displayed_rank"], rank)
            self.assertEqual(exact["relative_saturation_index_in_displayed_subgroup"], 1)
            self.assertEqual(exact["intersection_rank"], intersection)
            self.assertEqual(exact["remaining_quotient_rank"], quotient)
            self.assertEqual(len(exact["remaining_quotient_classes"]), quotient)

    def test_local_codes_have_the_full_declared_dimensions(self) -> None:
        for label, quotient in (("curve273", 13), ("curve302", 14)):
            record = self.curves[label]
            kummer = record["finite_kummer_code"]
            self.assertEqual(kummer["g17_image_dimension"], 17)
            self.assertEqual(kummer["quotient_image_dimension"], quotient)
            self.assertEqual(len(kummer["remaining_point_quotient_signatures"]), quotient)
            component = record["bad_fibre_component_code"]
            self.assertEqual(component["g17_component_cokernel_order"], 1)
            self.assertEqual(component["pair_sum_replay"], "PASS")

    def test_height_grams_and_theta_profiles_are_complete(self) -> None:
        for record in self.curves.values():
            heights = record["canonical_height_and_theta"]
            gram = heights["g17_gram"]
            self.assertEqual((len(gram), {len(row) for row in gram}), (17, {17}))
            for name in ("g17", "candidate_core"):
                profile = heights[name]
                self.assertEqual(len(profile["lll_gram"]), 17)
                self.assertGreaterEqual(profile["unoriented_lines_through_5_lambda"], 1312)
                self.assertIn("1312", profile["normalized_shell_height_by_nearest_rank"])


if __name__ == "__main__":
    unittest.main()

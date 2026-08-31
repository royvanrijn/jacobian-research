from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)


def rank_mod_two(columns: list[list[int]]) -> int:
    if not columns:
        return 0
    rows = [list(row) for row in zip(*columns)]
    rank = 0
    width = len(columns)
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index][column]:
                rows[index] = [left ^ right for left, right in zip(rows[index], rows[rank])]
        rank += 1
    return rank


class ElkiesBisectionSpecializationControlsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_complete_summary(self) -> None:
        self.assertEqual(self.data["atlas_size"], 39120)
        self.assertEqual(self.data["total_square_tests"], 195600)
        summary = self.data["summary"]
        self.assertEqual(summary["split_counts"], [6, 3, 2, 1, 25])
        self.assertEqual(summary["split_class_span_dimensions"], [5, 3, 2, 1, 4])
        self.assertEqual(summary["finite_quotient_escape_counts"], [0, 0, 0, 0, 0])
        for fibre in self.data["fibres"]:
            rank = fibre["rank_result"]
            self.assertTrue(rank["exact_generated_subgroup_rank_determined"])
            self.assertEqual(
                rank["exact_generated_subgroup_rank_after_adjoining"],
                rank["existing_unconditional_rank_lower_bound"],
            )

    def test_every_hit_replays_on_the_minimal_curve(self) -> None:
        for fibre in self.data["fibres"]:
            model = tuple(Fraction(value) for value in fibre["minimal_model"])
            for hit in fibre["hits"]:
                root = Fraction(hit["canonical_positive_square_root"])
                self.assertEqual(root**2, Fraction(hit["q_value"]))
                for key in ("positive_minimal_point", "negative_minimal_point"):
                    point = tuple(Fraction(value) for value in hit[key])
                    self.assertTrue(is_on_weierstrass_curve(model, point))

    def test_reported_class_spans(self) -> None:
        for fibre in self.data["fibres"]:
            public_dimension = fibre["public_complement"]["dimension"]
            classes = [
                hit["finite_quotient_class_modulo_generic_17"][
                    "coordinates_over_f2"
                ][:public_dimension]
                for hit in fibre["hits"]
            ]
            self.assertEqual(
                rank_mod_two(classes),
                fibre["split_class_span"]["dimension_in_known_public_complement"],
            )


if __name__ == "__main__":
    unittest.main()

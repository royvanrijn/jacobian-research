from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest

from search_nagao_u135_skew_height import (
    CHART_CENTER_COUNT,
    explicit_subset_records,
    projective_naive_height,
    select_chart_centers,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]


class U135SkewHeightSearchTests(unittest.TestCase):
    def test_projective_height(self) -> None:
        self.assertEqual(projective_naive_height((Q(-17, 23), Q(0))), 23)
        self.assertEqual(projective_naive_height((Q(29, 3), Q(0))), 29)

    def test_chart_centres_are_unique_and_largest_first(self) -> None:
        points = (
            (Q(1, 2), Q(1)),
            (Q(7, 3), Q(1)),
            (Q(-11, 5), Q(1)),
            (Q(7, 3), Q(-1)),
        )
        selected = select_chart_centers(points, count=2)
        self.assertEqual([point[0] for point in selected], [Q(-11, 5), Q(7, 3)])

    def test_chart_budget_matches_u42(self) -> None:
        self.assertEqual(CHART_CENTER_COUNT * 2, 76)

    def test_explicit_subset_uses_one_based_indices(self) -> None:
        pool = ((Q(1), Q(2)), (Q(3), Q(4)), (Q(5), Q(6)))
        records = explicit_subset_records(pool, (3, 1))
        self.assertEqual(records[0]["pool_index_one_based"], 3)
        self.assertEqual(records[0]["jacobian_x"], "5")
        self.assertEqual(records[1]["jacobian_y"], "2")

    def test_generated_artifact_is_self_consistent_when_present(self) -> None:
        artifact_path = (
            ROOT
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_u135_skew_height.json"
        )
        if not artifact_path.exists():
            self.skipTest("the long bounded-search artifact is not present")
        data = json.loads(artifact_path.read_text())
        script_path = (
            ROOT
            / "elliptic-curves"
            / "cas"
            / "search_nagao_u135_skew_height.py"
        )
        self.assertEqual(
            data["script_sha256"], hashlib.sha256(script_path.read_bytes()).hexdigest()
        )
        self.assertEqual(data["candidate"]["parameter_u"], "135/2")
        self.assertEqual(data["declared_budget"]["chart_count"], 76)
        self.assertEqual(len(data["skew_search"]["boxes"]), 10)
        self.assertEqual(data["stable_pool_numerical_rank"], 17)
        self.assertFalse(data["stable_rank_at_least_18_observed"])
        points = data["outside_uniform_checkpoint"]["points"]
        self.assertEqual(len(points), 74)
        self.assertTrue(
            all(
                point["exact_relation_replayed_with_fraction_group_law"]
                for point in points
            )
        )


if __name__ == "__main__":
    unittest.main()

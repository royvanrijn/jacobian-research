from __future__ import annotations

from pathlib import Path
import sys
import unittest


ECSEARCH = Path(__file__).resolve().parents[1] / "ecsearch"
sys.path.insert(0, str(ECSEARCH))

from fermigier_baseline_evaluation import evaluate_ranker, validate_config  # noqa: E402


class FermigierBaselineEvaluationTests(unittest.TestCase):
    def test_higher_score_ranking_and_missing_positive_are_explicit(self) -> None:
        ranker = {
            "id": "toy",
            "feature_source": "toy",
            "field": "score",
            "direction": "higher",
            "role": "test",
            "initial_population_count": 10,
            "materialized_population_count": 3,
            "selection": "toy",
            "leakage_boundary": "toy",
        }
        candidates = [
            self._candidate(1.0, None, 12, 0),
            self._candidate(3.0, "positive-a", 15, 4),
            self._candidate(2.0, None, 13, 2),
        ]
        result = evaluate_ranker(
            ranker, candidates, ["positive-a", "positive-b"], [1, 2, 10]
        )
        self.assertEqual(result["positive_positions"]["positive-a"]["position_one_based"], 1)
        self.assertEqual(result["positive_coverage"]["missing_positive_ids"], ["positive-b"])
        self.assertEqual(result["budget_metrics"][0]["recall_of_all_admitted_positives"], 0.5)
        self.assertEqual(
            result["censored_outcome_overlap"]["strata"][
                "legacy_rank_at_least_13_uncertified"
            ]["count"],
            1,
        )

    def test_lower_score_direction(self) -> None:
        ranker = {
            "id": "toy",
            "feature_source": "toy",
            "field": "rank",
            "direction": "lower",
            "role": "test",
            "initial_population_count": 2,
            "materialized_population_count": 2,
            "selection": "toy",
            "leakage_boundary": "toy",
        }
        candidates = [
            self._candidate(2, None, None, None),
            self._candidate(1, "positive-a", None, None),
        ]
        result = evaluate_ranker(ranker, candidates, ["positive-a", "positive-b"], [1])
        self.assertEqual(result["positive_positions"]["positive-a"]["position_one_based"], 1)

    def test_config_rejects_non_retrospective_role(self) -> None:
        with self.assertRaisesRegex(ValueError, "retrospective"):
            validate_config(
                {
                    "schema": "elliptic-curves.fermigier-baseline-rankers.v1",
                    "evaluation_role": "held_out",
                }
            )

    @staticmethod
    def _candidate(value, positive_id, legacy_rank, quartic_points):
        return {
            "value": value,
            "positive_id": positive_id,
            "projective_pair_T": [1, 1],
            "legacy_rank": legacy_rank,
            "quartic_point_count": quartic_points,
            "numerical_rank": None,
        }


if __name__ == "__main__":
    unittest.main()

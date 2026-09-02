from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest


ECSEARCH = Path(__file__).resolve().parents[1] / "ecsearch"
SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(ECSEARCH))
sys.path.insert(0, str(SCRIPTS))

from r17_bisection_ranker import (  # noqa: E402
    FEATURE_NAMES,
    feature_vector,
    fit_weighted_contrast,
    ranking_metrics,
    score,
    semantic_label_sha256,
)
from evaluate_r17_bisection_ranker_prospective import (  # noqa: E402
    average_precision,
    roc_auc,
)


def synthetic_row(value: float) -> dict:
    return {
        "height": 10,
        "features": {
            "level1_nagao": {
                "standardized_block_signals": [value, value, value],
                "worst_block_signal": value,
                "mean_block_signal": value,
                "bad_prime_count": 0,
            },
            "level2_conductor_proxy": {
                "quality_proxy": value,
                "known_prime_log_discriminant_saving": value,
                "log_discriminant_after_known_scaling": 10 + value,
            },
            "level2_quotient_code": {
                "rarity": value,
                "local_E_mod_2_and_mod_3_dimensions": [[0, 0]] * 5,
            },
            "level2_cover_diversity": {
                "distinct_character_patterns": 1 + value,
                "character_pattern_entropy": value,
                "branch_zero_count": 0,
            },
        },
    }


class R17BisectionRankerTests(unittest.TestCase):
    def test_feature_vector_has_frozen_shape(self) -> None:
        self.assertEqual(set(feature_vector(synthetic_row(1))), set(FEATURE_NAMES))

    def test_contrast_scores_positive_direction_higher(self) -> None:
        rows = [synthetic_row(value) for value in (0, 1, 4, 5)]
        model = fit_weighted_contrast(rows, [0, 0, 1, 1], [1, 1, 1, 1])
        self.assertGreater(score(model, rows[-1]), score(model, rows[0]))
        altered = copy.deepcopy(rows[-1])
        altered["features"]["level2_conductor_proxy"]["quality_proxy"] = 9
        self.assertGreater(score(model, altered), score(model, rows[-1]))

    def test_ranking_metrics(self) -> None:
        metrics = ranking_metrics([4, 3, 2, 1], [1, 0, 1, 0], [1, 2])
        self.assertEqual(metrics["budgets"]["1"]["positive_count"], 1)
        self.assertEqual(metrics["budgets"]["2"]["recall"], 0.5)

    def test_semantic_label_hash_ignores_timings(self) -> None:
        first = {"parameter": "1/2", "outcomes": {"finite_quotient_gain_lower_bound": 1, "cpu_seconds": 2.0}}
        second = {"parameter": "1/2", "outcomes": {"finite_quotient_gain_lower_bound": 1, "cpu_seconds": 9.0}}
        self.assertEqual(semantic_label_sha256([first]), semantic_label_sha256([second]))

    def test_prospective_ranking_metrics(self) -> None:
        self.assertEqual(average_precision([4, 3, 2, 1], [1, 1, 0, 0]), 1.0)
        self.assertEqual(roc_auc([4, 3, 2, 1], [1, 1, 0, 0]), 1.0)
        self.assertEqual(roc_auc([1, 1], [1, 0]), 0.5)


if __name__ == "__main__":
    unittest.main()

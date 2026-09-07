from __future__ import annotations

import sys
from pathlib import Path
import unittest


ECSEARCH = Path(__file__).resolve().parents[1] / "ecsearch"
sys.path.insert(0, str(ECSEARCH))

from r17_training_data import (  # noqa: E402
    EMBARGOED_PARAMETERS,
    conductor_proxy_features,
    deterministic_sample,
    development_lane_memberships,
    normalize_quadratic,
    normalized_parameter,
    split_quadratic_indices,
    split_bucket,
)


class R17TrainingDataTests(unittest.TestCase):
    def test_normalization_and_embargoed_sampling(self) -> None:
        self.assertEqual(normalized_parameter(4, -6), (-2, 3))
        first = deterministic_sample(
            count=500, height=100, seed=12345, excluded=EMBARGOED_PARAMETERS
        )
        second = deterministic_sample(
            count=500, height=100, seed=12345, excluded=EMBARGOED_PARAMETERS
        )
        self.assertEqual(first, second)
        self.assertEqual(len(first), len(set(first)))
        self.assertFalse(EMBARGOED_PARAMETERS.intersection(first))
        self.assertTrue(all(b > 0 and __import__("math").gcd(abs(a), b) == 1 for a, b in first))

    def test_split_is_stable(self) -> None:
        self.assertEqual(split_bucket((17, 29), "test"), split_bucket((17, 29), "test"))
        self.assertIn(split_bucket((17, 29), "test"), {"train", "validation", "internal_test"})

    def test_quadratic_normalization(self) -> None:
        self.assertEqual(
            normalize_quadratic(["1/2", "3/4", "-5/6"]),
            ((-6, -9, 10), -12),
        )

    def test_exact_split_scan_distinguishes_ramification(self) -> None:
        quadratics = [
            ((1, 0, 0), 1),
            ((1, 0, 0), -1),
            ((-1, 0, 1), 2),
            ((1, 0, 1), 1),
        ]
        split, ramified = split_quadratic_indices((3, 1), quadratics)
        self.assertEqual(split, [0, 2])
        self.assertEqual(ramified, [])
        split, ramified = split_quadratic_indices((1, 1), quadratics)
        self.assertEqual(split, [0])
        self.assertEqual(ramified, [2])

    def test_conductor_feature_is_explicitly_a_proxy(self) -> None:
        features = conductor_proxy_features((1, 1), [1] + [0] * 8, [1] + [0] * 12, [2, 3, 5, 7])
        self.assertIn("not an exact conductor", features["boundary"])
        self.assertGreater(features["log_discriminant_after_known_scaling"], 0)

    def test_lane_quotas_preserve_each_split(self) -> None:
        splits = ["train"] * 14 + ["validation"] * 3 + ["internal_test"] * 3
        records = []
        for index, split in enumerate(splits):
            records.append(
                {
                    "parameter": f"{index}/101",
                    "height": 101,
                    "split": split,
                    "features": {
                        "level1_nagao": {
                            "worst_block_signal": float(index),
                            "mean_block_signal": float(index),
                        },
                        "level2_conductor_proxy": {"quality_proxy": float(index)},
                        "level2_quotient_code": {"rarity": float(index)},
                        "level2_cover_diversity": {
                            "distinct_character_patterns": index,
                            "character_pattern_entropy": float(index),
                        },
                    },
                }
            )
        lanes = development_lane_memberships(records, 10, "test")
        for indices in lanes.values():
            counts = {name: 0 for name in {"train", "validation", "internal_test"}}
            for index in indices:
                counts[records[index]["split"]] += 1
            self.assertEqual(counts, {"train": 7, "validation": 2, "internal_test": 1})


if __name__ == "__main__":
    unittest.main()

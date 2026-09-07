#!/usr/bin/env python3
"""Tests for the disjoint-prime q12/orbit5867 Nagao rerank."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "rerank_h92_q12o5867_rootless_nagao_holdout.py"
SPEC = importlib.util.spec_from_file_location("q12o5867_holdout", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
HOLDOUT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = HOLDOUT
SPEC.loader.exec_module(HOLDOUT)


class HeldoutNagaoRerankTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = HOLDOUT.load_family_model()

    def test_small_disjoint_rerank_is_a_complete_permutation(self) -> None:
        blocks, rejected = HOLDOUT.build_residue_tables(self.model, ((199, 211),))
        self.assertEqual(rejected, ())
        source = [
            {
                "parameter": "1/1",
                "projective_pair": [1, 1],
                "projective_height": 1,
                "total_score_units_1e12": 30,
                "block_score_units_1e12": [10, 10, 10],
            },
            {
                "parameter": "0/1",
                "projective_pair": [0, 1],
                "projective_height": 1,
                "total_score_units_1e12": 20,
                "block_score_units_1e12": [8, 7, 5],
            },
            {
                "parameter": "infinity",
                "projective_pair": [1, 0],
                "projective_height": 1,
                "total_score_units_1e12": 10,
                "block_score_units_1e12": [4, 3, 3],
            },
        ]
        records, holdout_order, robust_order = HOLDOUT.rerank_population(
            source, blocks[0]
        )
        expected = {"1/1", "0/1", "infinity"}
        self.assertEqual(set(holdout_order), expected)
        self.assertEqual(set(robust_order), expected)
        self.assertEqual(
            sorted(record["holdout_rank"] for record in records), [1, 2, 3]
        )
        self.assertTrue(
            all(
                record["holdout_good_prime_count"]
                + record["holdout_bad_reduction_prime_count"]
                == 2
                for record in records
            )
        )

    def test_overlap_summary_preserves_holdout_order(self) -> None:
        summary = HOLDOUT.overlap_summary(
            ["a", "b", "c", "d"], ["c", "b", "d", "a"], (2, 3)
        )
        self.assertEqual(summary[0]["overlap_parameters_in_holdout_order"], ["b"])
        self.assertEqual(summary[1]["overlap_parameters_in_holdout_order"], ["c", "b"])


if __name__ == "__main__":
    unittest.main()

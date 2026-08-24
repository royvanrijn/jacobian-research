from __future__ import annotations

from decimal import Decimal
from pathlib import Path
import shutil
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from staged_record_rescore import (  # noqa: E402
    Candidate,
    assert_population,
    candidate_sort_key,
    parse_keep_counts,
    score_candidates_with_pari,
)
from search_record_residue_class import primitive_rationals_in_class  # noqa: E402


class StagedRecordRescoreTests(unittest.TestCase):
    def test_keep_counts_and_deterministic_ranking(self) -> None:
        self.assertEqual(parse_keep_counts("50,10,10"), (50, 10, 10))
        with self.assertRaises(Exception):
            parse_keep_counts("10,20")
        records = [
            {"score": "2.0", "height": 9, "numerator": 2, "denominator": 9},
            {"score": "3.0", "height": 50, "numerator": 1, "denominator": 50},
            {"score": "2.0", "height": 8, "numerator": 1, "denominator": 8},
        ]
        records.sort(key=candidate_sort_key)
        self.assertEqual([Decimal(record["score"]) for record in records], [3, 2, 2])
        self.assertEqual(records[1]["height"], 8)

    def test_small_population_validation(self) -> None:
        candidates = tuple(
            Candidate(a, b) for a, b in primitive_rationals_in_class(200)
        )
        assert_population(candidates, 200)
        self.assertEqual(len({candidate.identifier for candidate in candidates}), len(candidates))

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is not installed")
    def test_batch_score_order_and_counts(self) -> None:
        candidates = (Candidate(1666, 9), Candidate(1666, 4227))
        records = score_candidates_with_pari(
            candidates,
            37,
            timeout=10.0,
            stack_bytes=32_000_000,
        )
        self.assertEqual({record["t"] for record in records}, {"1666/9", "1666/4227"})
        self.assertTrue(all(record["primes_used"] > 0 for record in records))


if __name__ == "__main__":
    unittest.main()

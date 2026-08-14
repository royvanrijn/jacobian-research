from __future__ import annotations

import sys
from pathlib import Path
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from screen_multiple_root_frontier import (  # noqa: E402
    frontier_key,
    select_candidates,
)


def record(t: str, height_rank: int, score_rank: int) -> dict:
    return {
        "t": t,
        "numerator": int(t.split("/")[0]),
        "denominator": int(t.split("/")[1]),
        "height": max(abs(int(part)) for part in t.split("/")),
        "height_rank": height_rank,
        "score_rank_within_height_pool": score_rank,
    }


class FrontierSelectionTests(unittest.TestCase):
    def test_union_of_prefixes_is_exact_and_stable(self) -> None:
        records = [
            record("5/7", 3, 1),
            record("1/2", 1, 9),
            record("3/4", 2, 8),
            record("9/10", 4, 2),
        ]
        selected = select_candidates(records, height_count=2, score_count=2)
        self.assertEqual(
            [item["t"] for item in selected], ["1/2", "3/4", "5/7", "9/10"]
        )

    def test_negative_prefix_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            select_candidates([], height_count=-1, score_count=0)

    def test_frontier_prefers_threshold_then_odd_parity(self) -> None:
        candidates = []
        for t, logn, root, height_rank in (
            ("1/2", "170", 1, 1),
            ("3/4", "180", -1, 2),
            ("5/6", "150", -1, 3),
            ("7/8", "190", -1, 4),
        ):
            item = record(t, height_rank, height_rank)
            item["pari"] = {"log_conductor": logn, "root_number": root}
            candidates.append(item)
        self.assertEqual(
            [item["t"] for item in sorted(candidates, key=frontier_key)],
            ["5/6", "3/4", "1/2", "7/8"],
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import sys
from pathlib import Path
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from search_nagao_rank13_integer_u import (  # noqa: E402
    NAGAO_CONSTRUCTION,
    nagao_base_change,
    nagao_curve,
    parse_keep_counts,
    parse_stages,
)


class NagaoIntegerSearchTests(unittest.TestCase):
    def test_base_change_and_u_one_model(self) -> None:
        self.assertEqual(str(nagao_base_change(1)), "23549/2")
        curve = nagao_curve(1)
        self.assertEqual(curve.parameter_t, nagao_base_change(1))
        self.assertEqual(len(curve.coefficients), 5)
        self.assertTrue(NAGAO_CONSTRUCTION.is_quartic_family)
        self.assertNotEqual(curve.coefficients[-1], 0)

    def test_zero_parameter_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            nagao_base_change(0)

    def test_stage_parsers(self) -> None:
        self.assertEqual(parse_stages("200,1000,10000"), (200, 1000, 10000))
        self.assertEqual(parse_keep_counts("40,12,12"), (40, 12, 12))
        with self.assertRaises(Exception):
            parse_stages("1000,200")
        with self.assertRaises(Exception):
            parse_keep_counts("12,40")


if __name__ == "__main__":
    unittest.main()

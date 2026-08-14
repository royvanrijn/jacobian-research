from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import shutil
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from compare_score_cutoffs import (  # noqa: E402
    RECORD_CLASS_CURVES,
    last_primes_for_cutoffs,
    parse_cutoffs,
    score_curve_with_pari,
)
from search_record_residue_class import (  # noqa: E402
    build_score_tables,
    score_rational,
)


class ScoreCutoffComparisonTests(unittest.TestCase):
    def test_cutoff_parsing_and_last_primes(self) -> None:
        self.assertEqual(parse_cutoffs("5,10,100"), (5, 10, 100))
        self.assertEqual(last_primes_for_cutoffs((5, 10, 100)), (5, 7, 97))
        with self.assertRaises(Exception):
            parse_cutoffs("10,10")
        with self.assertRaises(ValueError):
            last_primes_for_cutoffs((4,))
        with self.assertRaises(ValueError):
            last_primes_for_cutoffs((100, 10))

    def test_pinned_primary_candidate(self) -> None:
        self.assertEqual(RECORD_CLASS_CURVES[0].parameter, Fraction(1666, 9))
        self.assertEqual(
            [curve.record_scan_position_at_500 for curve in RECORD_CLASS_CURVES],
            list(range(1, 11)),
        )

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is not installed")
    def test_pari_score_matches_python_reference_at_37(self) -> None:
        parameter = Fraction(1666, 9)
        pari = score_curve_with_pari(
            parameter,
            (37,),
            timeout=10.0,
            stack_bytes=32_000_000,
        )[37]
        reference = score_rational(
            parameter.numerator,
            parameter.denominator,
            build_score_tables(37, "fermigier-good"),
        )
        self.assertAlmostEqual(float(pari["score"]), reference["value"], places=12)
        self.assertEqual(pari["primes_used"], reference["primes_used"])
        self.assertEqual(
            pari["skipped_denominator_primes"],
            reference["skipped_denominator_primes"],
        )
        self.assertEqual(
            pari["skipped_bad_primes"], reference["skipped_bad_primes"]
        )


if __name__ == "__main__":
    unittest.main()

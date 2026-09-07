from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))
CAS_ROOT = PROGRAM_ROOT / "cas"
sys.path.insert(0, str(CAS_ROOT))

from ecsearch.fermigier import evaluate_polynomial, fermigier_quartic  # noqa: E402
from ecsearch.fermigier_offset_search import (  # noqa: E402
    denominator_offset_points,
    normalized_quartic_integer_value,
)
from fermigier_mestre import FermigierMestreFamily  # noqa: E402


class FermigierOffsetSearchTests(unittest.TestCase):
    def test_integer_clearing_matches_fraction_model(self) -> None:
        for parameter, sign, numerator, denominator in (
            (Q(91, 5), 1, 17, 9),
            (Q(1155, 2), -1, -31, 16),
            (Q(385, 12), 1, 99, 37),
        ):
            value, _, square_denominator = normalized_quartic_integer_value(
                parameter, sign, numerator, denominator
            )
            shift = 2 * parameter
            x_value = sign * shift + Q(numerator, denominator)
            expected = FermigierMestreFamily.quartic_value(shift, x_value)
            self.assertEqual(Q(value, square_denominator**2), expected)

    def test_known_denominator_point_is_recovered_and_replayed(self) -> None:
        points = denominator_offset_points(
            Q(91, 5), maximum_denominator=5, maximum_abs_numerator=5_000
        )
        model = fermigier_quartic(Q(182, 5))
        self.assertTrue(points)
        for point in points:
            self.assertGreaterEqual(point.offset_denominator, 2)
            self.assertEqual(
                point.raw_y * point.raw_y,
                evaluate_polynomial(model.quartic, point.x),
            )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3

from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from elkies_klagsbrun_rank29 import PUBLISHED_POINTS, to_short_point  # noqa: E402
from search_elkies_klagsbrun_rank30 import point_add, point_negate, poly_evaluate  # noqa: E402
from search_elkies_klagsbrun_rank30_alternate_covers import (  # noqa: E402
    alternate_parameter,
    alternate_to_short_points,
    build_all_covers,
    build_alternate_charts,
    short_alternate_discriminant,
)


Q = Fraction


class ElkiesKlagsbrunAlternateCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.covers = build_all_covers()

    def test_complete_cover_count_and_first_score(self) -> None:
        self.assertEqual(len(self.covers), 4060)
        self.assertEqual(self.covers[0].subset_indices, (5, 25))
        self.assertEqual(self.covers[0].raw_height_score[0], 1276971481641273)

    def test_weight_two_cover_maps_its_two_summands(self) -> None:
        cover = self.covers[0]
        first = to_short_point(PUBLISHED_POINTS[cover.subset_indices[0]])
        second = to_short_point(PUBLISHED_POINTS[cover.subset_indices[1]])
        parameter = alternate_parameter(cover.short_point, first)
        self.assertEqual(parameter, alternate_parameter(cover.short_point, second))
        discriminant = short_alternate_discriminant(cover.short_point)
        square_root = first[0] - second[0]
        self.assertEqual(poly_evaluate(discriminant, parameter), square_root**2)
        self.assertEqual(
            set(alternate_to_short_points(cover.short_point, parameter, square_root)),
            {first, second},
        )

    def test_weight_two_sum_is_exact(self) -> None:
        cover = self.covers[0]
        expected = point_add(PUBLISHED_POINTS[5], PUBLISHED_POINTS[25])
        self.assertEqual(cover.general_point, expected)
        self.assertEqual(
            point_add(expected, point_negate(PUBLISHED_POINTS[25])),
            PUBLISHED_POINTS[5],
        )

    def test_pinned_chart_count(self) -> None:
        charts = build_alternate_charts(
            cover_count=4,
            offset_count=3,
            cross_ratio_count=4,
            offset_height=10,
            cross_ratio_height=10,
        )
        self.assertEqual(len(charts), 28)


if __name__ == "__main__":
    unittest.main()

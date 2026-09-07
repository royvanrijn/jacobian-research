from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from alternate_quartic_covers import (  # noqa: E402
    alternate_cover,
    mobius_image,
    mobius_preimage,
    point_on_short_curve,
    short_add,
    short_subset_sum,
    three_point_mobius_matrix,
)
from nagao_1994 import rank13_base_changed_short_jacobian_coefficients  # noqa: E402


Q = Fraction


class AlternateQuarticCoverTests(unittest.TestCase):
    def test_small_exact_round_trips(self) -> None:
        coefficients = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        base = (Q(0), Q(1))
        cover = alternate_cover(coefficients, base)
        self.assertEqual(cover.coefficients, (Q(4), Q(-8), Q(0), Q(0), Q(1)))

        for point in ((Q(1), Q(1)), (Q(3), Q(5)), (Q(0), Q(-1))):
            cover_point = cover.curve_point_to_cover(point)
            self.assertEqual(cover.cover_point_to_curve(cover_point), point)
        with self.assertRaises(ValueError):
            cover.curve_point_to_cover(base)

    def test_two_quartic_signs_give_the_two_quadratic_roots(self) -> None:
        coefficients = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        cover = alternate_cover(coefficients, (Q(0), Q(1)))
        self.assertEqual(cover.cover_point_to_curve((Q(2), Q(-2))), (Q(1), Q(1)))
        self.assertEqual(cover.cover_point_to_curve((Q(2), Q(2))), (Q(3), Q(5)))

    def test_exact_group_law_and_subset_sum(self) -> None:
        coefficients = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        point = (Q(0), Q(1))
        doubled = short_add(coefficients, point, point)
        self.assertEqual(doubled, (Q(1, 4), Q(-7, 8)))
        self.assertEqual(short_subset_sum(coefficients, (point, doubled), (0, 1)), short_add(coefficients, point, doubled))
        self.assertIsNone(short_add(coefficients, point, (point[0], -point[1])))

    def test_three_point_chart_and_inverse(self) -> None:
        matrix = three_point_mobius_matrix(Q(-7, 3), Q(11, 5), Q(19, 2))
        self.assertEqual(mobius_image(matrix, Q(0)), Q(-7, 3))
        self.assertEqual(mobius_image(matrix, Q(1)), Q(11, 5))
        self.assertEqual(Q(matrix[0], matrix[2]), Q(19, 2))
        for parameter in (Q(-5, 7), Q(0), Q(1), Q(13, 11)):
            value = mobius_image(matrix, parameter)
            self.assertIsNotNone(value)
            self.assertEqual(mobius_preimage(matrix, value), parameter)

    def test_u135_basis_points_define_distinct_exact_charts(self) -> None:
        artifact = json.loads(
            (
                ROOT
                / "artifacts/generated-results/elliptic-curves/elliptic_nagao_rank17_frontier_certificate.json"
            ).read_text()
        )
        record = next(
            item for item in artifact["certificates"] if item["parameter_u"] == "135/2"
        )
        coefficients = rank13_base_changed_short_jacobian_coefficients(Q(135, 2))
        basis = tuple(
            (Q(point["jacobian_x"]), Q(point["jacobian_y"]))
            for point in record["saturated_basis"]
        )
        first = alternate_cover(coefficients, basis[0])
        second = alternate_cover(coefficients, basis[1])
        self.assertNotEqual(first.coefficients, second.coefficients)
        for cover, point in ((first, basis[2]), (second, basis[3])):
            self.assertTrue(point_on_short_curve(coefficients, point))
            self.assertEqual(cover.cover_point_to_curve(cover.curve_point_to_cover(point)), point)


if __name__ == "__main__":
    unittest.main()

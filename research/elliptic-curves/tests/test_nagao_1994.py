from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from nagao_1994 import (  # noqa: E402
    RANK13_CONSTRUCTION,
    RANK13_CONDUCTOR_FACTORIZATION,
    RANK13_PUBLISHED_CONDUCTOR,
    RANK13_PUBLISHED_MODEL,
    RANK13_PUBLISHED_POINTS,
    RANK21_CONSTRUCTION,
    RANK21_CONDUCTOR_FACTORIZATION,
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_CONDUCTOR,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_PARAMETER,
    RANK21_PUBLISHED_POINTS,
    even_discriminant_polynomial,
    factorization_product,
    point_on_extended_weierstrass,
    polynomial_content,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    rank13_base_changed_discriminant_numerator,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
    rank13_leading_square,
    rank13_published_quartic_coefficients,
    rank21_short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data  # noqa: E402


class Nagao1994Tests(unittest.TestCase):
    def test_both_root_tuples_give_primitive_quartics(self) -> None:
        self.assertEqual(RANK13_CONSTRUCTION.quartic_condition, 0)
        self.assertEqual(RANK21_CONSTRUCTION.quartic_condition, 0)
        self.assertEqual(RANK13_CONSTRUCTION.quartic_content, 1557504)
        self.assertEqual(RANK13_CONSTRUCTION.quartic_square_scale, 1248)
        self.assertEqual(RANK21_CONSTRUCTION.quartic_content, 508953600)
        self.assertEqual(RANK21_CONSTRUCTION.quartic_square_scale, 22560)

    def test_rank13_printed_quartic_and_sections_are_exact(self) -> None:
        for parameter in (Q(1), Q(2), Q(7, 3), Q(-5, 2)):
            self.assertEqual(
                primitive_quartic_coefficients(RANK13_CONSTRUCTION, parameter),
                rank13_published_quartic_coefficients(parameter),
            )
            points = rank13_known_quartic_points(parameter)
            self.assertEqual(len(points), 13)
            self.assertEqual(len(set(points)), 13)
            images = tuple(
                quartic_point_to_short_jacobian(
                    RANK13_CONSTRUCTION, parameter, point
                )
                for point in points
            )
            self.assertEqual(len(images), 13)

    def test_rank21_twelve_visible_points_map_exactly(self) -> None:
        self.assertEqual(
            RANK21_CONSTRUCTOR_PARAMETER, 2 * RANK21_PUBLISHED_PARAMETER
        )
        points = primitive_visible_points(
            RANK21_CONSTRUCTION, RANK21_CONSTRUCTOR_PARAMETER
        )
        self.assertEqual(len(points), 12)
        for point in points:
            quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, RANK21_CONSTRUCTOR_PARAMETER, point
            )

    def test_rank13_base_change_splits_the_points_at_infinity(self) -> None:
        for parameter_u in (Q(1), Q(2), Q(7, 3), Q(-5, 2)):
            parameter_t = rank13_base_parameter(parameter_u)
            leading = primitive_quartic_coefficients(
                RANK13_CONSTRUCTION, parameter_t
            )[-1]
            self.assertEqual(rank13_leading_square(parameter_u) ** 2, leading)
        with self.assertRaises(ValueError):
            rank13_base_parameter(Q(0))
        with self.assertRaises(ValueError):
            rank13_leading_square(Q(0))

    def test_discriminant_geometry_is_pinned(self) -> None:
        rank13 = RANK13_CONSTRUCTION.primitive_discriminant_polynomial
        rank21 = RANK21_CONSTRUCTION.primitive_discriminant_polynomial
        self.assertEqual(len(rank13) - 1, 20)
        self.assertEqual(len(rank21) - 1, 20)
        self.assertTrue(all(not rank13[index] for index in range(1, 21, 2)))
        self.assertTrue(all(not rank21[index] for index in range(1, 21, 2)))
        self.assertEqual(len(even_discriminant_polynomial(RANK13_CONSTRUCTION)) - 1, 10)
        self.assertEqual(len(even_discriminant_polynomial(RANK21_CONSTRUCTION)) - 1, 10)
        self.assertEqual(polynomial_content(rank13), 11664)
        self.assertEqual(polynomial_content(rank21), 8464)
        base_changed = rank13_base_changed_discriminant_numerator()
        self.assertEqual(len(base_changed) - 1, 40)
        self.assertTrue(all(not base_changed[index] for index in range(1, 41, 2)))

    def test_all_printed_points_satisfy_their_models_exactly(self) -> None:
        self.assertEqual(len(RANK13_PUBLISHED_POINTS), 13)
        self.assertEqual(len(RANK21_PUBLISHED_POINTS), 21)
        self.assertTrue(
            all(
                point_on_extended_weierstrass(RANK13_PUBLISHED_MODEL, point)
                for point in RANK13_PUBLISHED_POINTS
            )
        )
        self.assertTrue(
            all(
                point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, point)
                for point in RANK21_PUBLISHED_POINTS
            )
        )
        x, y = RANK21_PUBLISHED_POINTS[7]
        self.assertFalse(
            point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, (x + 1, y))
        )

    def test_conductor_factorizations_multiply_exactly(self) -> None:
        self.assertEqual(
            factorization_product(RANK13_CONDUCTOR_FACTORIZATION),
            RANK13_PUBLISHED_CONDUCTOR,
        )
        self.assertEqual(
            factorization_product(RANK21_CONDUCTOR_FACTORIZATION),
            RANK21_PUBLISHED_CONDUCTOR,
        )

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_base_changed_rank13_u1_reduces_to_printed_model(self) -> None:
        data = minimal_curve_data(
            rank13_base_changed_short_jacobian_coefficients(Q(1)), timeout=10
        )
        self.assertEqual(data["minimal_model"], RANK13_PUBLISHED_MODEL)
        self.assertEqual(data["conductor"], RANK13_PUBLISHED_CONDUCTOR)
        self.assertTrue(data["log_conductor"].startswith("165.406045732330510373"))

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_rank21_specialization_reduces_to_printed_model(self) -> None:
        data = minimal_curve_data(rank21_short_jacobian_coefficients(), timeout=10)
        self.assertEqual(data["minimal_model"], RANK21_PUBLISHED_MODEL)
        self.assertEqual(data["conductor"], RANK21_PUBLISHED_CONDUCTOR)
        self.assertTrue(data["log_conductor"].startswith("196.679545735892153436"))


if __name__ == "__main__":
    unittest.main()


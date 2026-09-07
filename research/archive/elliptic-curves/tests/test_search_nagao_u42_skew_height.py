from fractions import Fraction
import unittest

from nagao_1994 import (
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_value,
    rank13_base_parameter,
)
from search_nagao_u42_skew_height import (
    SEARCH_BOXES,
    centered_unimodular_matrix,
    companion_section_x_values,
    exact_linear_combination,
    map_chart_point,
    short_add,
    short_multiply,
    transform_binary_quartic,
)


Q = Fraction


class SkewHeightSearchTests(unittest.TestCase):
    def test_boxes_are_a_contiguous_disjoint_denominator_staircase(self) -> None:
        lower = 1
        for search_box in SEARCH_BOXES:
            self.assertEqual(search_box.denominator_lower, lower)
            self.assertGreaterEqual(
                search_box.numerator_bound, search_box.denominator_upper
            )
            lower = search_box.denominator_upper + 1
        self.assertEqual(lower, 128_001)

    def test_companion_x_values_at_u42(self) -> None:
        parameter_t = rank13_base_parameter(Q(42))
        self.assertEqual(
            set(companion_section_x_values(parameter_t)),
            {
                Q(6211, 210),
                Q(1829, 10),
                Q(-355, 6),
                Q(47189, 70),
                Q(-39983, 210),
            },
        )

    def test_centered_charts_are_unimodular_and_map_zero_to_center(self) -> None:
        center = Q(948253, 4298)
        for shift in (0, -1, 3):
            matrix = centered_unimodular_matrix(center, shift)
            a_value, b_value, c_value, d_value = matrix
            self.assertEqual(a_value * d_value - b_value * c_value, 1)
            self.assertEqual(
                map_chart_point((Q(0), Q(1)), matrix)[0], center
            )

    def test_binary_quartic_chart_maps_an_exact_point_back(self) -> None:
        parameter_t = rank13_base_parameter(Q(42))
        coefficients = primitive_quartic_coefficients(
            RANK13_CONSTRUCTION, parameter_t
        )
        center = Q(948253, 4298)
        matrix = centered_unimodular_matrix(center)
        transformed = transform_binary_quartic(coefficients, matrix)
        original_z = Q(162342527487816, 4618201)
        transformed_point = Q(0), original_z * matrix[3] ** 2
        self.assertEqual(
            transformed_point[1] ** 2,
            quartic_value(transformed, transformed_point[0]),
        )
        self.assertEqual(
            map_chart_point(transformed_point, matrix), (center, original_z)
        )

    def test_exact_short_group_arithmetic(self) -> None:
        # P=(3,5) lies on y^2=x^3-2.  These checks exercise signs, infinity,
        # double-and-add, and the linear-combination replay used by the artifact.
        point = (Q(3), Q(5))
        self.assertIsNone(short_add(Q(-2), point, (point[0], -point[1])))
        doubled = short_multiply(Q(-2), point, 2)
        self.assertEqual(short_add(Q(-2), point, point), doubled)
        self.assertEqual(
            exact_linear_combination(Q(-2), (point, doubled), (3, -1)),
            point,
        )


if __name__ == "__main__":
    unittest.main()

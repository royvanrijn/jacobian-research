from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_REPORTED_PARAMETER,
    evaluate_polynomial,
    fermigier_canonical_coefficients,
    fermigier_quartic,
    quartic_point_to_canonical_point,
    thirteen_visible_points,
    thirteenth_visible_point,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (  # noqa: E402
    add_rational_points,
    build_independence_certificate,
    is_on_weierstrass_curve,
    matrix_rank_mod_prime,
    negate_rational_point,
    select_independent_subset,
    verify_independence_certificate,
)


class FermigierPointTests(unittest.TestCase):
    def test_thirteenth_point_formula_and_sign(self) -> None:
        # The construction has total degree at most twelve in (x,s); after the
        # linear x-substitution these thirteen exact values certify the
        # polynomial square identity in s.
        for shift in range(-6, 8):
            if shift == 0:
                continue
            model = fermigier_quartic(shift)
            positive = thirteenth_visible_point(model)
            negative = thirteenth_visible_point(model, y_sign=-1)
            self.assertEqual(positive[0], negative[0])
            self.assertEqual(positive[1], -negative[1])
            self.assertEqual(
                positive[1] ** 2,
                evaluate_polynomial(model.quartic, positive[0]),
            )
        reconstruction = fermigier_quartic(
            FERMIGIER_E22_RECONSTRUCTION_SHIFT
        )
        self.assertEqual(
            thirteenth_visible_point(reconstruction),
            (
                Fraction(-46964, 195),
                Fraction(3170976819397626546496, 164775),
            ),
        )

    def test_covariant_map_for_all_thirteen_points(self) -> None:
        for adapter_parameter in (
            Fraction(-3),
            Fraction(1),
            Fraction(7, 3),
            FERMIGIER_REPORTED_PARAMETER,
        ):
            model = fermigier_quartic(2 * adapter_parameter)
            canonical = fermigier_canonical_coefficients(adapter_parameter)
            mapped = tuple(
                quartic_point_to_canonical_point(model, point)
                for point in thirteen_visible_points(model)
            )
            self.assertEqual(len(set(mapped)), 13)
            self.assertTrue(
                all(is_on_weierstrass_curve(canonical, point) for point in mapped)
            )

    def test_e22_section_specialization_has_twelve_differences(self) -> None:
        specialization = specialize_fermigier_rank_sections(
            FERMIGIER_REPORTED_PARAMETER
        )
        self.assertEqual(len(specialization.quartic_points), 13)
        self.assertEqual(len(specialization.section_differences), 12)
        self.assertTrue(
            all(
                is_on_weierstrass_curve(specialization.canonical_model, point)
                for point in specialization.section_differences
            )
        )


class ReductionCertificateTests(unittest.TestCase):
    def test_matrix_rank(self) -> None:
        self.assertEqual(matrix_rank_mod_prime(((1, 2), (2, 4)), 5), 1)
        self.assertEqual(matrix_rank_mod_prime(((1, 2), (2, 0)), 5), 2)
        with self.assertRaisesRegex(ValueError, "equal length"):
            matrix_rank_mod_prime(((1,), (1, 2)), 5)

    def test_exact_rank_two_certificate_and_subset_selection(self) -> None:
        curve = (0, 0, 0, -25, 25)
        first = (Fraction(5), Fraction(5))
        second = (Fraction(-5), Fraction(5))
        certificate = build_independence_certificate(
            curve,
            (first, second),
            relation_prime=2,
            maximum_reduction_prime=100,
        )
        verify_independence_certificate(curve, (first, second), certificate)

        dependent = add_rational_points(curve, first, second)
        assert dependent is not None
        cloud = (first, second, dependent, negate_rational_point(curve, first))
        indices, subset_certificate = select_independent_subset(
            curve,
            cloud,
            relation_prime=2,
            maximum_reduction_prime=100,
        )
        self.assertEqual(indices, (0, 1))
        verify_independence_certificate(
            curve,
            tuple(cloud[index] for index in indices),
            subset_certificate,
        )


if __name__ == "__main__":
    unittest.main()

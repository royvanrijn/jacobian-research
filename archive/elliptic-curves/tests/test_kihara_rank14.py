from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from kihara_rank14 import (  # noqa: E402
    PUBLISHED_HEIGHT_DETERMINANT_APPROX,
    PUBLISHED_INDEPENDENCE_SPECIALIZATION,
    PUBLISHED_RANK_LOWER_BOUND,
    binary_invariants,
    kihara_specialization,
    known_quartic_points,
    published_extra_abscissae,
    short_jacobian_coefficients,
    specialized_parameters,
    verify_rational_specialization,
)
from verify_kihara_rank14_symbolic import symbolic_verification  # noqa: E402


class KiharaRank14Tests(unittest.TestCase):
    def test_rational_specializations_replay_all_fifteen_points(self) -> None:
        for parameter_t in (Q(1), Q(2), Q(3, 2), Q(2, 3), Q(-1)):
            specialization = kihara_specialization(parameter_t)
            points = known_quartic_points(parameter_t)
            self.assertEqual(len(specialization.product_coefficients), 13)
            self.assertEqual(len(specialization.approximant_coefficients), 7)
            self.assertEqual(len(specialization.quartic_coefficients), 5)
            self.assertEqual(len(points), 15)
            self.assertEqual(len(set(points)), 15)
            verification = verify_rational_specialization(parameter_t)
            self.assertTrue(verification["all_points_exact"])
            self.assertTrue(verification["quartic_discriminant_nonzero"])

    def test_t2_extra_abscissae_pin_the_primary_source_transcription(self) -> None:
        specialization = kihara_specialization(Q(2))
        self.assertEqual(
            published_extra_abscissae(specialization),
            (Q(29183481217024, 421), Q(4755487195136, 269), Q(-517996544)),
        )

    def test_short_jacobian_is_nonsingular_at_t2(self) -> None:
        coefficients = short_jacobian_coefficients(Q(2))
        self.assertEqual(coefficients[:3], (Q(0), Q(0), Q(0)))
        coefficient_a, coefficient_b = coefficients[3:]
        self.assertNotEqual(4 * coefficient_a**3 + 27 * coefficient_b**2, 0)
        invariant_i, invariant_j = binary_invariants(
            kihara_specialization(Q(2)).quartic_coefficients
        )
        self.assertEqual(coefficient_a, -27 * invariant_i)
        self.assertEqual(coefficient_b, -27 * invariant_j)

    def test_printed_base_change_excludes_zero(self) -> None:
        with self.assertRaises(ValueError):
            specialized_parameters(Q(0))

    def test_publication_status_is_recorded_without_upgrading_it(self) -> None:
        self.assertEqual(PUBLISHED_RANK_LOWER_BOUND, 14)
        self.assertEqual(PUBLISHED_INDEPENDENCE_SPECIALIZATION, Q(2))
        self.assertEqual(PUBLISHED_HEIGHT_DETERMINANT_APPROX, "221792776617402574.10")

    def test_generic_point_identities_hold_symbolically(self) -> None:
        result = symbolic_verification()
        self.assertEqual(result["coefficient_domain"], "QQ(t)[x]")
        self.assertEqual(result["product_degree"], 12)
        self.assertEqual(result["remainder_degree"], 4)
        self.assertEqual(result["visible_section_count"], 12)
        self.assertTrue(result["visible_section_identities_exact"])
        self.assertTrue(result["all_fifteen_sections_exact"])
        self.assertTrue(result["generic_quartic_discriminant_nonzero"])
        self.assertTrue(
            all(
                signature["is_square"]
                for signature in result["extra_section_square_signatures"].values()
            )
        )


if __name__ == "__main__":
    unittest.main()

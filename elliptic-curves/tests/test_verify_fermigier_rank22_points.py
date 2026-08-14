from __future__ import annotations

import unittest

from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from verify_fermigier_rank22_points import (
    EXPECTED_CERTIFICATE_PRIMES,
    EXPECTED_TWO_TORSION_CERTIFICATE_PRIME,
    PUBLISHED_POINTS,
    curve_residual,
    exact_independence_certificate,
    exact_strict_conductor_comparison,
    minimal_point_to_short,
)
from verify_fermigier_benchmark import PUBLISHED_CONDUCTOR


class FermigierRank22PointCertificateTests(unittest.TestCase):
    def test_published_points_and_exact_short_transport(self) -> None:
        self.assertEqual(len(PUBLISHED_POINTS), 22)
        self.assertTrue(all(curve_residual(point) == 0 for point in PUBLISHED_POINTS))

        _, _, _, coefficient_a, coefficient_b = FermigierMestreFamily.coefficients(
            NORMALIZED_RECORD_PARAMETER
        )
        for point in PUBLISHED_POINTS:
            x_value, y_value = minimal_point_to_short(point)
            self.assertEqual(
                y_value**2,
                x_value**3 + coefficient_a * x_value + coefficient_b,
            )

    def test_finite_reductions_certify_all_22_points(self) -> None:
        certificate = exact_independence_certificate()
        self.assertEqual(certificate["combined_exact_rank_over_F2"], 22)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 22)
        self.assertEqual(
            tuple(certificate["certificate_primes"]),
            EXPECTED_CERTIFICATE_PRIMES,
        )
        self.assertEqual(
            certificate["two_torsion_certificate_prime"],
            EXPECTED_TWO_TORSION_CERTIFICATE_PRIME,
        )

    def test_exact_conductor_target_miss(self) -> None:
        comparison = exact_strict_conductor_comparison(PUBLISHED_CONDUCTOR)
        self.assertTrue(comparison["integer_power_inequality_holds"])
        self.assertFalse(comparison["meets_strict_log_conductor_target"])
        self.assertEqual(
            comparison["exact_conclusion"],
            "log(conductor) > 4568/25 = 182.72",
        )


if __name__ == "__main__":
    unittest.main()

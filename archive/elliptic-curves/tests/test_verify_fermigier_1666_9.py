from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import shutil
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from verify_fermigier_1666_9 import (  # noqa: E402
    BAD_PRIMES,
    EXPECTED_CONDUCTOR,
    EXPECTED_CONDUCTOR_FACTORIZATION,
    EXPECTED_H_VALUATIONS,
    EXPECTED_MINIMAL_DISCRIMINANT,
    EXPECTED_MINIMAL_MODEL,
    EXPECTED_SPECIALIZED_MODEL,
    build_artifact,
    exact_specialization_data,
    factorization_product,
    integral_weierstrass_invariants,
    kodaira_symbol,
    selected_known_points,
)


class Fermigier1666Over9VerifierTests(unittest.TestCase):
    def test_dependency_free_exact_specialization(self) -> None:
        data = exact_specialization_data()
        self.assertEqual(data["parameter"], "1666/9")
        self.assertEqual(
            data["discriminant_factor"]["exact_valuations"],
            {str(prime): exponent for prime, exponent in EXPECTED_H_VALUATIONS.items()},
        )
        self.assertEqual(
            data["specialized_short_weierstrass_model"][
                "coefficients_a1_a2_a3_a4_a6"
            ],
            [str(value) for value in EXPECTED_SPECIALIZED_MODEL],
        )
        self.assertEqual(data["known_jacobian_points"]["count"], 12)
        self.assertTrue(
            data["known_jacobian_points"]["all_checked_exactly_on_curve"]
        )

    def test_selected_points_satisfy_the_exact_short_equation(self) -> None:
        _, _, _, coefficient_a, coefficient_b = EXPECTED_SPECIALIZED_MODEL
        points = selected_known_points()
        self.assertEqual(len(points), 12)
        self.assertEqual(len(set(points)), 12)
        for x_value, y_value in points:
            self.assertEqual(
                y_value**2,
                x_value**3 + coefficient_a * x_value + coefficient_b,
            )

    def test_pinned_conductor_factorization_multiplies_exactly(self) -> None:
        self.assertEqual(BAD_PRIMES[:8], (2, 3, 5, 7, 13, 17, 37, 43))
        self.assertEqual(len(BAD_PRIMES), 10)
        self.assertEqual(
            factorization_product(EXPECTED_CONDUCTOR_FACTORIZATION),
            EXPECTED_CONDUCTOR,
        )

    def test_minimal_model_discriminant_is_independently_exact(self) -> None:
        invariants = integral_weierstrass_invariants(EXPECTED_MINIMAL_MODEL)
        self.assertEqual(invariants["discriminant"], EXPECTED_MINIMAL_DISCRIMINANT)
        self.assertEqual(kodaira_symbol(6), "I2")
        self.assertEqual(kodaira_symbol(-1), "IV*")
        with self.assertRaises(ValueError):
            kodaira_symbol(0)

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is not installed")
    def test_full_default_replay_has_no_rank_claim(self) -> None:
        artifact = build_artifact(timeout=30.0, stack_bytes=256_000_000)
        self.assertEqual(artifact["curve"]["conductor"], str(EXPECTED_CONDUCTOR))
        self.assertEqual(
            set(artifact["curve"]["local_reduction_at_every_bad_prime"]),
            {str(prime) for prime in BAD_PRIMES},
        )
        self.assertEqual(
            artifact["numerical_height_pairing"]["on_curve_count"], 12
        )
        self.assertEqual(artifact["rank_status"]["status"], "unknown")
        self.assertFalse(artifact["rank_status"]["ellrank_invoked"])
        self.assertFalse(
            artifact["targets"]["small_conductor_target"]["certified_hit"]
        )


if __name__ == "__main__":
    unittest.main()

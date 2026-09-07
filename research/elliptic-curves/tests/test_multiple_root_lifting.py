from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from ek_k3 import fraction_mod, legendre_symbol, valuation  # noqa: E402
from fermigier_mestre import (  # noqa: E402
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
)
from multiple_root_lifting import (  # noqa: E402
    RootBall,
    RootLiftCapExceeded,
    affine_variable_coefficients,
    all_roots_mod_prime_power,
    fixed_divisor_valuation,
    scaled_variable_coefficients,
    verify_prime_power_roots,
)
from search_fermigier_mixed_small_prime_crt_gauss import (  # noqa: E402
    CALIBRATION_PROFILE,
    ICARM_282_U,
    discover_groups,
)
from search_mestre_02557104116148_power_root_crt import (  # noqa: E402
    combine_rows,
    kernel_basis,
    matching_ball,
)


class GenericMultipleRootLiftingTests(unittest.TestCase):
    def test_complete_enumeration_matches_brute_force(self) -> None:
        polynomials = (
            (1, 1),
            (0, 0, 1),
            (3, 0, 1),
            (2, -3, 1),
            (0, 1, 1, 1),
        )
        for coefficients in polynomials:
            for prime in (2, 3, 5):
                for exponent in range(1, 5):
                    with self.subTest(
                        coefficients=coefficients,
                        prime=prime,
                        exponent=exponent,
                    ):
                        modulus = prime**exponent
                        expected = tuple(
                            residue
                            for residue in range(modulus)
                            if sum(
                                coefficient * residue**degree
                                for degree, coefficient in enumerate(coefficients)
                            )
                            % modulus
                            == 0
                        )
                        actual = all_roots_mod_prime_power(
                            coefficients, prime, exponent
                        )
                        self.assertEqual(actual.roots, expected)

    def test_multiple_root_splits_and_then_thins(self) -> None:
        # x^2 = 0 modulo 3^3 has exactly the three roots divisible by 3^2.
        result = all_roots_mod_prime_power((0, 0, 1), 3, 3)
        self.assertEqual(result.level_counts, (1, 3, 3))
        self.assertEqual(result.roots, (0, 9, 18))
        self.assertEqual(
            result.maximal_balls(),
            (RootBall(3, 2, 0),),
        )

    def test_simple_roots_are_also_supported(self) -> None:
        result = all_roots_mod_prime_power((-1, 0, 1), 5, 4)
        self.assertEqual(result.level_counts, (2, 2, 2, 2))
        self.assertEqual(result.roots, (1, 624))

    def test_multiple_root_can_disappear(self) -> None:
        result = all_roots_mod_prime_power((3, 0, 1), 3, 2)
        self.assertEqual(result.level_counts, (1, 0))
        self.assertEqual(result.roots, ())
        self.assertEqual(result.maximal_balls(), ())

    def test_cap_never_returns_a_partial_answer(self) -> None:
        with self.assertRaises(RootLiftCapExceeded) as caught:
            all_roots_mod_prime_power((0, 0, 1), 7, 2, max_roots=6)
        self.assertEqual(caught.exception.reached_exponent, 2)

    def test_verifier_rejects_tampering(self) -> None:
        result = all_roots_mod_prime_power((-1, 0, 1), 5, 2)
        bad = result.__class__(
            prime=result.prime,
            exponent=result.exponent,
            modulus=result.modulus,
            roots=(1, 2),
            level_counts=result.level_counts,
            candidate_digits_checked=result.candidate_digits_checked,
        )
        with self.assertRaises(AssertionError):
            verify_prime_power_roots((-1, 0, 1), bad)


class FermigierMultipleRootTests(unittest.TestCase):
    H = DISCRIMINANT_FACTOR_COEFFICIENTS
    FAMILY = FermigierMestreFamily

    def test_mixed_small_prime_projective_ball_counts(self) -> None:
        groups, _ = discover_groups()
        self.assertEqual(
            {prime: len(balls) for prime, balls in groups.items()},
            {5: 2, 7: 1, 11: 3, 13: 4, 17: 3, 19: 8, 23: 6, 29: 6, 31: 6},
        )

    def test_icarm_282_is_radius_52_in_effective_profile(self) -> None:
        groups, _ = discover_groups()
        numerator, denominator = ICARM_282_U.numerator, ICARM_282_U.denominator
        choices = tuple(
            matching_ball(numerator, denominator, groups[prime])
            for prime in CALIBRATION_PROFILE
        )
        coefficient_a, coefficient_b, modulus = combine_rows(choices)
        first, second = kernel_basis(coefficient_a, coefficient_b, modulus)
        determinant = first[0] * second[1] - first[1] * second[0]
        left = (numerator * second[1] - denominator * second[0]) // determinant
        right = (first[0] * denominator - first[1] * numerator) // determinant
        self.assertEqual(modulus, 101_959)
        self.assertEqual((first, second), ((-220, -27), (77, -454)))
        self.assertEqual((left, right), (-52, 3))

    def test_small_prime_lift_profiles(self) -> None:
        expected = {
            5: (5, 25, 25),
            7: (1, 7, 49),
            11: (3, 33, 363),
            13: (6, 78, 338),
            17: (3, 51, 867),
            19: (7, 133, 836),
            23: (6, 138, 0),
            29: (3, 87, 174),
            31: (6, 186, 124),
            37: (6, 222, 2886),
        }
        for prime, counts in expected.items():
            with self.subTest(prime=prime):
                result = all_roots_mod_prime_power(
                    self.H, prime, 3, max_roots=3_000
                )
                self.assertEqual(result.level_counts, counts)

    def test_automatic_high_powers_compress_to_cheap_balls(self) -> None:
        roots_11 = all_roots_mod_prime_power(
            self.H, 11, 4, max_roots=4_000
        )
        self.assertEqual(
            roots_11.maximal_balls(),
            (
                RootBall(11, 1, 0),
                RootBall(11, 1, 5),
                RootBall(11, 1, 6),
            ),
        )

        roots_17 = all_roots_mod_prime_power(
            self.H, 17, 4, max_roots=15_000
        )
        self.assertEqual(
            roots_17.maximal_balls(),
            (
                RootBall(17, 1, 0),
                RootBall(17, 1, 1),
                RootBall(17, 1, 16),
            ),
        )

        # H(t) is always divisible by 5^2: this is a fixed divisor of the
        # polynomial as a function, although its coefficient content is one.
        roots_5 = all_roots_mod_prime_power(self.H, 5, 2)
        self.assertEqual(roots_5.maximal_balls(), (RootBall(5, 0, 0),))
        self.assertEqual(roots_5.reciprocal_density, Q(1))
        self.assertEqual(fixed_divisor_valuation(self.H, 5), 2)

    def test_record_parameter_local_valuations(self) -> None:
        t = NORMALIZED_RECORD_PARAMETER
        invariants = self.FAMILY.invariants(t)
        expected = {
            7: (18, 4, 6),
            17: (4, 2, 2),
            37: (3, 0, 0),
        }
        for prime, (delta_value, c4_value, c6_value) in expected.items():
            with self.subTest(prime=prime):
                self.assertNotEqual(t.denominator % prime, 0)
                self.assertEqual(
                    valuation(invariants["weierstrass_discriminant"], prime),
                    delta_value,
                )
                self.assertEqual(valuation(invariants["c4"], prime), c4_value)
                self.assertEqual(valuation(invariants["c6"], prime), c6_value)

        # At 37, the record lies in one of two classes modulo 37 that force
        # v_37(H)>=3.  Four further, smaller balls modulo 37^2 complete the
        # full root set.
        roots_37 = all_roots_mod_prime_power(self.H, 37, 3)
        balls_37 = roots_37.maximal_balls()
        self.assertEqual(
            balls_37,
            (
                RootBall(37, 1, 4),
                RootBall(37, 1, 33),
                RootBall(37, 2, 153),
                RootBall(37, 2, 661),
                RootBall(37, 2, 708),
                RootBall(37, 2, 1216),
            ),
        )
        self.assertEqual(fraction_mod(t, 37), 33)
        self.assertEqual(
            fixed_divisor_valuation(
                affine_variable_coefficients(self.H, 33, 37), 37
            ),
            3,
        )

    def test_seven_adic_fixed_divisor_after_scaling_parameter(self) -> None:
        # H(7*s) is divisible by 7^18 for every integer s.  Thus the enormous
        # apparent power at the record parameter costs only T=0 (mod 7).
        scaled = scaled_variable_coefficients(self.H, 7)
        self.assertEqual(fixed_divisor_valuation(scaled, 7), 18)
        self.assertEqual(fraction_mod(NORMALIZED_RECORD_PARAMETER, 7), 0)

    def test_small_crt_conductor_fixture_hits_the_same_cheap_balls(self) -> None:
        # This is the shortest representative of T=2142 (mod 7*17*37).
        t = Q(-119, 2)
        expected = {7: (0, 18), 17: (0, 4), 37: (33, 3)}
        for prime, (residue, minimum_valuation) in expected.items():
            with self.subTest(prime=prime):
                self.assertEqual(fraction_mod(t, prime), residue)
                self.assertGreaterEqual(
                    valuation(self.FAMILY.discriminant_factor(t), prime),
                    minimum_valuation,
                )

    def test_clean_multiplicative_and_additive_classes(self) -> None:
        # These base-prime classes force unexpectedly high powers while c4
        # remains a unit.  All are clean split-multiplicative fibers.
        clean_classes = {
            11: (4, (0, 5, 6)),
            17: (4, (1, 16)),
            19: (4, (5, 14)),
            37: (3, (4, 33)),
        }
        for prime, (forced_exponent, residues) in clean_classes.items():
            for residue in residues:
                with self.subTest(prime=prime, residue=residue):
                    self.assertGreaterEqual(
                        fixed_divisor_valuation(
                            affine_variable_coefficients(
                                self.H, residue, prime
                            ),
                            prime,
                        ),
                        forced_exponent,
                    )
                    local = self.FAMILY.local_data(residue, prime)
                    self.assertTrue(local.split_multiplicative)
                    self.assertNotEqual(
                        fraction_mod(
                            self.FAMILY.invariants(Q(residue))["c4"], prime
                        ),
                        0,
                    )

        # The zero class at 17 also forces a fourth power, but c4 vanishes in
        # an already minimal model, so it is additive rather than cleanly
        # multiplicative.
        additive = self.FAMILY.local_data(0, 17)
        self.assertIsNone(additive.split_multiplicative)
        self.assertEqual(
            fraction_mod(self.FAMILY.invariants(Q(0))["c4"], 17), 0
        )

    def test_seven_class_is_nonminimal_before_scaling(self) -> None:
        # local_data sees c4=0 on T=0 (mod 7), but that chart is nonminimal.
        # Refining to T=7*s (mod 49), s!=0, permits x=7^2*x', y=7^3*y'.
        # The scaled model has unit c4 and split multiplicative reduction.
        for residue_after_division in range(1, 7):
            t = Q(7 * residue_after_division)
            invariants = self.FAMILY.invariants(t)
            self.assertEqual(valuation(invariants["c4"], 7), 4)
            self.assertEqual(valuation(invariants["c6"], 7), 6)
            self.assertGreaterEqual(
                valuation(invariants["weierstrass_discriminant"], 7), 18
            )
            coefficients = self.FAMILY.coefficients(t)
            scaled_a = coefficients[3] / 7**4
            scaled_b = coefficients[4] / 7**6
            character_sum = sum(
                legendre_symbol(
                    x**3
                    + fraction_mod(scaled_a, 7) * x
                    + fraction_mod(scaled_b, 7),
                    7,
                )
                for x in range(7)
            )
            self.assertEqual(-character_sum, 1)

        # T=0 (mod 49) does not acquire unit c4 after the same scaling and is
        # genuinely additive at this first refined residue.
        invariants_at_zero = self.FAMILY.invariants(Q(0))
        self.assertGreater(valuation(invariants_at_zero["c4"], 7), 4)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from math import gcd
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.calibration import (  # noqa: E402
    binary_discriminant_factor,
    integral_model,
    model_discriminant,
)
from ecsearch.crt_lattice import (  # noqa: E402
    crt,
    enumerate_rationals_in_height_box,
    evaluate_polynomial,
    first_rationals_by_height,
    gauss_reduce_congruence_lattice,
    hensel_lift_roots,
)
from ecsearch.local_data import (  # noqa: E402
    calibration_family_local_data,
    short_weierstrass_local_data,
    weierstrass_local_data,
)
from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS,
    FERMIGIER_REPORTED_PARAMETER,
    FERMIGIER_ROOTS,
    evaluate_polynomial as evaluate_fraction_polynomial,
    fermigier_canonical_coefficients,
    fermigier_discriminant_factor,
    fermigier_quartic,
    multiply_polynomials,
    twelve_visible_points,
    weierstrass_discriminant,
    weierstrass_c_invariants,
)


class HenselTests(unittest.TestCase):
    def test_simple_roots(self) -> None:
        roots = hensel_lift_roots((-2, 0, 1), 7, 3)
        self.assertEqual(roots, [108, 235])
        self.assertTrue(all((root * root - 2) % 343 == 0 for root in roots))

    def test_singular_branching_and_death(self) -> None:
        self.assertEqual(hensel_lift_roots((0, 0, 1), 3, 1), [0])
        self.assertEqual(hensel_lift_roots((0, 0, 1), 3, 2), [0, 3, 6])
        self.assertEqual(hensel_lift_roots((0, 0, 1), 3, 3), [0, 9, 18])
        self.assertEqual(hensel_lift_roots((3, 0, 1), 3, 1), [0])
        self.assertEqual(hensel_lift_roots((3, 0, 1), 3, 2), [])


class CrtAndLatticeTests(unittest.TestCase):
    def test_generalized_crt(self) -> None:
        self.assertEqual(crt(((1, 125), (48, 49))), (1126, 6125))
        self.assertEqual(crt(((2, 6), (5, 9))), (14, 18))
        with self.assertRaisesRegex(ValueError, "incompatible"):
            crt(((0, 4), (1, 2)))

    def test_gauss_basis_and_primitive_trap(self) -> None:
        first, second = gauss_reduce_congruence_lattice(1126, 6125)
        self.assertEqual((first, second), ((-49, -49), (38, -87)))
        self.assertEqual(abs(first[0] * second[1] - first[1] * second[0]), 6125)
        # The shortest lattice vector reduces to t=1, which is not in the
        # original projective congruence because 49 shares a factor with M.
        self.assertEqual(gcd(abs(first[0]), abs(first[1])), 49)
        self.assertNotEqual((-1 - 1126 * -1) % 6125, 0)

    def test_gauss_skew_weights_and_weight_validation(self) -> None:
        for weights in ((1, 1000), (1000, 1)):
            first, second = gauss_reduce_congruence_lattice(
                1126, 6125, weights=weights
            )
            norm = lambda vector: (
                weights[0] * vector[0] ** 2
                + weights[1] * vector[1] ** 2
            )
            self.assertLessEqual(norm(first), norm(second))
            self.assertLessEqual(
                2 * abs(
                    weights[0] * first[0] * second[0]
                    + weights[1] * first[1] * second[1]
                ),
                norm(first),
            )
            self.assertEqual(
                abs(first[0] * second[1] - first[1] * second[0]),
                6125,
            )
        for weights in (
            (0, 1),
            (1, -1),
            (0.5, 1),
            (True, 1),
            (1,),
            (1, 1, 1),
        ):
            with self.assertRaisesRegex(ValueError, "two positive integers"):
                gauss_reduce_congruence_lattice(7, 31, weights=weights)

    def test_complete_height_box(self) -> None:
        representatives = enumerate_rationals_in_height_box(
            1126,
            6125,
            numerator_bound=87,
            denominator_bound=87,
        )
        pairs = {
            (item.numerator, item.denominator) for item in representatives
        }
        self.assertEqual(pairs, {(-38, 87), (-87, 38)})
        self.assertEqual((-38) ** 2 - 87**2, -6125)

    def test_height_enumerator_against_direct_scan(self) -> None:
        for modulus in range(2, 25):
            for residue in range(modulus):
                bound = 12
                exact = enumerate_rationals_in_height_box(
                    residue,
                    modulus,
                    numerator_bound=bound,
                    denominator_bound=bound,
                )
                exact_pairs = {
                    (item.numerator, item.denominator) for item in exact
                }
                direct_pairs = {
                    (numerator, denominator)
                    for numerator in range(-bound, bound + 1)
                    for denominator in range(1, bound + 1)
                    if gcd(abs(numerator), denominator) == 1
                    and gcd(denominator, modulus) == 1
                    and (numerator - residue * denominator) % modulus == 0
                }
                self.assertEqual(exact_pairs, direct_pairs)

    def test_asymmetric_height_boxes(self) -> None:
        for numerator_bound, denominator_bound, weights in (
            (3, 19, (1, 1000)),
            (19, 3, (1000, 1)),
            (7, 13, (3, 11)),
        ):
            exact = enumerate_rationals_in_height_box(
                7,
                31,
                numerator_bound=numerator_bound,
                denominator_bound=denominator_bound,
                weights=weights,
            )
            exact_pairs = {
                (item.numerator, item.denominator) for item in exact
            }
            direct_pairs = {
                (numerator, denominator)
                for numerator in range(-numerator_bound, numerator_bound + 1)
                for denominator in range(1, denominator_bound + 1)
                if gcd(abs(numerator), denominator) == 1
                and gcd(denominator, 31) == 1
                and (numerator - 7 * denominator) % 31 == 0
            }
            self.assertEqual(exact_pairs, direct_pairs)

    def test_first_rationals_count_order_and_failure(self) -> None:
        representatives, checked_height = first_rationals_by_height(
            2, 7, count=4, maximum_height=9
        )
        direct = enumerate_rationals_in_height_box(
            2,
            7,
            numerator_bound=checked_height,
            denominator_bound=checked_height,
        )
        self.assertEqual(representatives, direct[:4])
        self.assertLessEqual(checked_height, 9)
        with self.assertRaisesRegex(ValueError, "fewer than"):
            first_rationals_by_height(2, 7, count=100, maximum_height=9)


class CalibrationFamilyTests(unittest.TestCase):
    def test_shortest_calibration_pair(self) -> None:
        representatives, checked_height = first_rationals_by_height(
            3249333373,
            143227016087,
            maximum_height=262144,
        )
        self.assertEqual(checked_height, 131072)
        pair = representatives[0]
        self.assertEqual((pair.numerator, pair.denominator), (-110627, 84367))
        self.assertEqual(
            binary_discriminant_factor(pair.numerator, pair.denominator),
            23**3 * 47**2 * 73**2,
        )

    def test_integral_model_identity(self) -> None:
        numerator, denominator = -110627, 84367
        a4, a6 = integral_model(numerator, denominator)
        x_coordinate = 0
        y_coordinate = numerator * denominator**2
        self.assertEqual(
            y_coordinate**2,
            x_coordinate**3 + a4 * x_coordinate + a6,
        )
        self.assertEqual(
            model_discriminant(numerator, denominator),
            -16 * (4 * a4**3 + 27 * a6**2),
        )
        self.assertEqual(
            evaluate_polynomial((27, 0, -4), numerator * pow(denominator, -1, 23**3), 23**3),
            0,
        )


class LocalTableTests(unittest.TestCase):
    def test_good_trace_matches_direct_point_count(self) -> None:
        local = calibration_family_local_data(2, 5)
        self.assertEqual(local.reduction, "good")
        direct_count = 1 + sum(
            1
            for x_coordinate in range(5)
            for y_coordinate in range(5)
            if (
                y_coordinate**2
                - (x_coordinate**3 - 4 * x_coordinate + 4)
            )
            % 5
            == 0
        )
        self.assertEqual(local.point_count, direct_count)
        self.assertEqual(local.trace, 5 + 1 - direct_count)

    def test_bad_primes_do_not_report_good_trace(self) -> None:
        for prime, root in ((23, 1), (47, 18), (73, 5)):
            local = calibration_family_local_data(root, prime)
            self.assertEqual(local.reduction, "split_multiplicative")
            self.assertIsNone(local.trace)
            self.assertIsNone(local.point_count)
            self.assertEqual(local.local_euler_coefficient, 1)
        unresolved = short_weierstrass_local_data(0, 0, 5)
        self.assertEqual(unresolved.reduction, "unresolved_bad")
        self.assertIsNone(unresolved.local_euler_coefficient)

    def test_general_weierstrass_trace_matches_direct_count(self) -> None:
        coefficients = (1, 2, 1, 3, 4)
        checked = 0
        for prime in (5, 7, 11, 13):
            local = weierstrass_local_data(coefficients, prime)
            if local.reduction != "good":
                continue
            a1, a2, a3, a4, a6 = coefficients
            direct_count = 1 + sum(
                1
                for x_coordinate in range(prime)
                for y_coordinate in range(prime)
                if (
                    y_coordinate**2
                    + a1 * x_coordinate * y_coordinate
                    + a3 * y_coordinate
                    - x_coordinate**3
                    - a2 * x_coordinate**2
                    - a4 * x_coordinate
                    - a6
                )
                % prime
                == 0
            )
            self.assertEqual(local.point_count, direct_count)
            self.assertEqual(local.trace, prime + 1 - direct_count)
            checked += 1
        self.assertGreater(checked, 0)


class FermigierFamilyTests(unittest.TestCase):
    def test_product_square_minus_quartic_identity(self) -> None:
        model = fermigier_quartic(Fraction(7, 3))
        square = multiply_polynomials(model.square_part, model.square_part)
        padded_quartic = model.quartic + (Fraction(0),) * 8
        self.assertEqual(
            square,
            tuple(
                model.product[index] + padded_quartic[index]
                for index in range(13)
            ),
        )

    def test_twelve_visible_points(self) -> None:
        model = fermigier_quartic(FERMIGIER_REPORTED_PARAMETER)
        points = twelve_visible_points(model)
        self.assertEqual(len(points), 2 * len(FERMIGIER_ROOTS))
        self.assertEqual(len(set(points)), 12)
        for x_coordinate, y_coordinate in points:
            self.assertEqual(
                y_coordinate * y_coordinate,
                evaluate_fraction_polynomial(model.quartic, x_coordinate),
            )

    def test_published_and_reconstruction_shifts_are_distinct(self) -> None:
        self.assertEqual(
            FERMIGIER_E22_RECONSTRUCTION_SHIFT,
            2 * FERMIGIER_REPORTED_PARAMETER,
        )
        literal = fermigier_quartic(FERMIGIER_REPORTED_PARAMETER)
        reconstruction = fermigier_quartic(
            FERMIGIER_E22_RECONSTRUCTION_SHIFT
        )
        self.assertNotEqual(literal.quartic, reconstruction.quartic)

    def test_canonical_discriminant_identity(self) -> None:
        # Both sides have degree at most 24.  Equality at 25 distinct values
        # proves the polynomial identity, while keeping the implementation
        # dependency-free.
        for parameter in range(-12, 13):
            coefficients = fermigier_canonical_coefficients(parameter)
            self.assertEqual(
                weierstrass_discriminant(coefficients),
                fermigier_discriminant_factor(parameter),
            )

    def test_literal_and_adapter_discriminant_coordinates(self) -> None:
        self.assertEqual(
            len(FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS), 21
        )
        for degree, (literal, adapter) in enumerate(
            zip(
                FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS,
                FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
                strict=True,
            )
        ):
            self.assertEqual(literal * 2**degree, 16 * adapter)
        expected_lifts = {
            89: [199136, 505833],
            131: [360891, 1887200],
            137: [713225, 1858128],
        }
        for prime, roots in expected_lifts.items():
            self.assertEqual(
                hensel_lift_roots(
                    FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS,
                    prime,
                    3,
                ),
                roots,
            )
            self.assertTrue(
                all(
                    sum(
                        coefficient * root**degree
                        for degree, coefficient in enumerate(
                            FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS
                        )
                    )
                    % prime**3
                    == 0
                    for root in roots
                )
            )

    def test_quartic_to_canonical_invariant_bridge(self) -> None:
        scale = 101232
        for adapter_parameter in range(-7, 8):
            if adapter_parameter == 0:
                continue
            e, d, c, b, a = fermigier_quartic(
                2 * adapter_parameter
            ).quartic
            quartic_i = 12 * a * e - 3 * b * d + c * c
            quartic_j = (
                72 * a * c * e
                + 9 * b * c * d
                - 27 * a * d * d
                - 27 * b * b * e
                - 2 * c**3
            )
            c4, c6 = weierstrass_c_invariants(
                fermigier_canonical_coefficients(adapter_parameter)
            )
            self.assertEqual(
                quartic_i,
                scale**4 * adapter_parameter**4 * c4,
            )
            self.assertEqual(
                quartic_j,
                2 * scale**6 * adapter_parameter**6 * c6,
            )

    def test_clean_simple_split_shaping_roots(self) -> None:
        expected_roots = {
            89: [23, 66],
            131: [7, 124],
            137: [67, 70],
            191: [37, 154],
        }
        derivative = tuple(
            degree * FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS[degree]
            for degree in range(1, len(FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS))
        )
        for prime, roots in expected_roots.items():
            self.assertEqual(
                hensel_lift_roots(
                    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
                    prime,
                    1,
                ),
                roots,
            )
            for root in roots:
                self.assertNotEqual(
                    evaluate_polynomial(derivative, root, prime), 0
                )
                coefficients = fermigier_canonical_coefficients(root)
                integral_coefficients = tuple(
                    int(coefficient) for coefficient in coefficients
                )
                local = weierstrass_local_data(integral_coefficients, prime)
                self.assertEqual(local.reduction, "split_multiplicative")
                self.assertEqual(local.local_euler_coefficient, 1)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for the fixed-cubic-field varying-curve experiment."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from fixed_cubic_field_curve_family import (  # noqa: E402
    bounded_integer_parameters,
    bounded_rational_parameters,
    covering_values,
    cubic_discriminant,
    discriminant_multiplier,
    f2_kernel_masks,
    f2_rank,
    field_multiply,
    fixed_field_cubic_coefficients,
    inverse_theta_coefficients,
)


Q = Fraction
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
)


class FixedCubicFieldAlgebraTests(unittest.TestCase):
    def test_discriminant_identity_and_anchor(self) -> None:
        A, B = -7, 11
        base = (Q(B), Q(A), Q(0), Q(1))
        for u in bounded_rational_parameters(4):
            coefficients = fixed_field_cubic_coefficients(A, B, u)
            self.assertEqual(
                cubic_discriminant(coefficients),
                cubic_discriminant(base) * discriminant_multiplier(A, B, u) ** 2,
            )
        self.assertEqual(fixed_field_cubic_coefficients(A, B, 0), base)

    def test_explicit_inverse_preserves_the_labelled_roots(self) -> None:
        A, B, u = Q(-7), Q(11), Q(2, 3)
        alpha = (Q(0), Q(1), u)
        alpha_squared = field_multiply(alpha, alpha, A, B)
        inverse = inverse_theta_coefficients(A, B, u)
        recovered = tuple(
            inverse[0] * (1 if index == 0 else 0)
            + inverse[1] * alpha[index]
            + inverse[2] * alpha_squared[index]
            for index in range(3)
        )
        self.assertEqual(recovered, (Q(0), Q(1), Q(0)))

    def test_covering_equations_for_an_anchor_kummer_point(self) -> None:
        A, B = Q(-7), Q(11)
        x_value = Q(3)
        beta = (x_value, Q(-1), Q(0))
        first, second, recovered_x_numerator = covering_values(
            beta, (Q(1), Q(0), Q(0)), Q(1), A, B, Q(0)
        )
        self.assertEqual((first, second), (0, 0))
        self.assertEqual(recovered_x_numerator, x_value)

    def test_whole_span_kernel_keeps_surviving_combinations(self) -> None:
        # Neither displayed row is locally zero, but their sum is.  A
        # basis-element-only filter would incorrectly return the zero space.
        rows = [[1], [1]]
        kernel = f2_kernel_masks(rows)
        self.assertEqual(f2_rank(rows), 1)
        self.assertEqual(kernel, [0b11])

    def test_parameter_policies_are_deterministic(self) -> None:
        self.assertEqual(
            bounded_integer_parameters(2),
            tuple(map(Q, (-2, -1, 0, 1, 2))),
        )
        self.assertEqual(
            bounded_rational_parameters(2),
            tuple(map(Q, (-2, -1, Q(-1, 2), 0, Q(1, 2), 1, 2))),
        )


@unittest.skipUnless(ARTIFACT.is_file(), "fixed-cubic-field artifact is absent")
class FixedCubicFieldArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_run_is_class_group_free_and_all_places_complete(self) -> None:
        self.assertEqual(
            self.data["status"],
            "PASS_EXACT_FULL_SPAN_LOCAL_INTERSECTIONS_NO_CLASS_GROUP",
        )
        self.assertFalse(self.data["class_group_computation_performed"])
        self.assertFalse(self.data["point_realization_computation_performed"])
        self.assertTrue(
            all(run["all_local_kummer_images_complete"] for run in self.data["runs"])
        )

    def test_whole_span_dimensions_and_zero_control(self) -> None:
        dimensions = {
            run["parameter_u"]: run["W_u_dimension"] for run in self.data["runs"]
        }
        self.assertEqual(dimensions, {"-2": 13, "-1": 18, "0": 20, "1": 13, "2": 13})
        zero = next(run for run in self.data["runs"] if run["parameter_u"] == "0")
        self.assertEqual(zero["combined_condition_rank"], 0)
        self.assertEqual(
            [row["one_based_anchor_basis_indices"] for row in zero["W_u_basis"]],
            [[index] for index in range(1, 21)],
        )

    def test_new_bad_primes_are_in_the_checked_support(self) -> None:
        for run in self.data["runs"]:
            support = set(run["complete_finite_place_support"])
            newly_bad = set(run["newly_bad_primes_relative_to_anchor"])
            self.assertLessEqual(newly_bad, support)
            local = {row["prime"]: row for row in run["finite_local_conditions"]}
            self.assertTrue(all(local[prime]["bad_reduction"] for prime in newly_bad))

    def test_curves_are_genuinely_different_in_the_bounded_run(self) -> None:
        self.assertTrue(self.data["family"]["pairwise_distinct_j_invariants"])
        self.assertTrue(self.data["family"]["u_zero_is_anchor_short_model"])
        self.assertEqual(self.data["family"]["parameter_count"], 5)


if __name__ == "__main__":
    unittest.main()

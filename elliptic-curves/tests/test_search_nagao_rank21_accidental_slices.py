#!/usr/bin/env python3

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from nagao_1994 import quartic_value  # noqa: E402
from nagao_1994_section7 import section7_primitive_quartic_coefficients  # noqa: E402
from search_nagao_rank21_accidental_slices import (  # noqa: E402
    EXPECTED_RECONSTRUCTION_SHA256,
    PINNED_IDENTITY_HEIGHT_200000_PARAMETERS,
    T0,
    build_priority_search_plans,
    build_slices,
    generic_labels_for_jacobian_image,
    generic_quartic_points,
    homogenized_transform,
    map_transformed_point,
    normalize_slice,
    rational_square_root,
    select_minimum_intercept_priority_slices,
    slice_polynomial,
    validate_known_linear_sections,
)


Q = Fraction
SCRIPT = CAS / "search_nagao_rank21_accidental_slices.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank21_accidental_slices.json"
EXAMPLE_POINT = (
    Q(1512322, 15933),
    Q(351536935858073, 84620163),
)
EXAMPLE_INTERCEPT = Q(-4471, 339)
EXAMPLE_NORMALIZED_QUARTIC = (
    1168595207736118032132129,
    -29899908013820942677608,
    308236804344177057976,
    -1537208192800455072,
    3140209905972624,
)


class NagaoAccidentalSliceTests(unittest.TestCase):
    def test_bivariate_substitution_matches_section7_quartic(self) -> None:
        for parameter, slope, intercept in (
            (Q(3, 2), Q(1), Q(-7, 3)),
            (T0, Q(-1), Q(130932, 491)),
            (Q(-11, 5), Q(4), Q(17, 9)),
        ):
            x_value = slope * parameter + intercept
            self.assertEqual(
                Q(quartic_value(section7_primitive_quartic_coefficients(parameter), x_value)),
                Q(sum(value * parameter**index for index, value in enumerate(slice_polynomial(slope, intercept)))),
            )

    def test_example_m1_slice_is_the_pinned_genus_one_quartic(self) -> None:
        normalized = normalize_slice(slice_polynomial(Q(1), EXAMPLE_INTERCEPT))
        self.assertEqual(normalized.raw_degree, 4)
        self.assertEqual(normalized.normalized_degree, 4)
        self.assertEqual(normalized.genus, 1)
        self.assertEqual(normalized.normalized_coefficients, EXAMPLE_NORMALIZED_QUARTIC)
        root = rational_square_root(normalized.normalized_value(T0))
        self.assertIsNotNone(root)
        self.assertEqual(
            normalized.original_ordinate(T0, root) ** 2,
            EXAMPLE_POINT[1] ** 2,
        )

    def test_known_linear_sections_normalize_to_squares(self) -> None:
        records = validate_known_linear_sections()
        self.assertEqual(len(records), 18)
        self.assertTrue(
            all(record["normalized_degree"] == 0 and record["genus"] == 0 for record in records)
        )

    def test_priority_choice_and_mobius_map_are_exact(self) -> None:
        slices = build_slices((EXAMPLE_POINT,))
        priority = select_minimum_intercept_priority_slices(slices)
        self.assertEqual(len(priority), 1)
        self.assertEqual(priority[0].slope, 1)
        plans = build_priority_search_plans(
            priority,
            identity_height=10,
            chart_height=10,
            chart_shifts=(0,),
            include_identity=False,
        )
        self.assertEqual(len(plans), 1)
        plan = plans[0]
        root = rational_square_root(priority[0].normalized.normalized_value(T0))
        self.assertIsNotNone(root)
        # T0 is the image of local coordinate zero in a centered chart.
        a_value, b_value, c_value, d_value = plan.matrix
        self.assertEqual(Q(b_value, d_value), T0)
        transformed_ordinate = root * d_value ** (plan.homogenized_degree // 2)
        self.assertEqual(
            sum(value * Q(0) ** index for index, value in enumerate(plan.polynomial)),
            transformed_ordinate**2,
        )
        self.assertEqual(
            map_transformed_point(
                (Q(0), transformed_ordinate),
                plan.matrix,
                total_degree=plan.homogenized_degree,
            ),
            (T0, root),
        )

    def test_generic_jacobian_decontamination(self) -> None:
        label, point = generic_quartic_points(Q(5605, 94))[6]
        labels = generic_labels_for_jacobian_image(Q(5605, 94), point)
        self.assertIn(label, labels)

    def test_pinned_pilot_and_generated_artifact(self) -> None:
        self.assertEqual(len(PINNED_IDENTITY_HEIGHT_200000_PARAMETERS), 81)
        self.assertEqual(len(set(PINNED_IDENTITY_HEIGHT_200000_PARAMETERS)), 81)
        self.assertTrue(ARTIFACT.exists())
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["reconstruction"]["reconstructed_point_sha256"],
            EXPECTED_RECONSTRUCTION_SHA256,
        )
        self.assertEqual(data["decontamination_at_T0"]["accidental_point_count"], 16)
        classification = data["complete_slice_classification"]
        self.assertEqual(classification["slice_count"], 272)
        self.assertEqual(
            [(row["raw_degree"], row["genus"], row["count"]) for row in classification["classification_counts"]],
            [(4, 1, 32), (6, 2, 240)],
        )
        auxiliary = data["auxiliary_search"]
        self.assertEqual(auxiliary["positive_rank_slice_count_by_point_cardinality"], 5)
        self.assertTrue(
            all(not record["retried"] for record in auxiliary["run_records"])
        )
        decontamination = data["new_parameter_decontamination"]
        self.assertEqual(decontamination["distinct_new_parameter_count_before_generic_filter"], 81)
        self.assertEqual(decontamination["generic_only_parameter_count"], 71)
        self.assertEqual(decontamination["forced_non_generic_parameter_count"], 10)


if __name__ == "__main__":
    unittest.main()

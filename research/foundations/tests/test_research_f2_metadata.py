"""Retained metadata aliases must not hide changed formal-jet evidence."""

from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plane-jc/cas"))
from f2_formal_homotopy_metadata import (
    REGULAR_GAUGE_VARIABLES,
    boundary_text,
    comparison_record,
)


class FormalMetadataTests(unittest.TestCase):
    def record(self):
        return {
            "field": {"prime": 61, "rho": 19, "y": 19},
            "claim_boundary": boundary_text(31),
            "higher_order_gauge": {
                "name": "seven-pole-coordinate-zero",
                "variables_prescribed_zero_from_order_two": list(REGULAR_GAUGE_VARIABLES),
            },
            "variable_series": {"x": [1, 2]},
            "root_series_digest_sha256": "retained-root-digest",
            "truncated_lambda_one_evaluation": {"nonzero_total": 50},
            "software": {"sympy": "recorded-version"},
        }

    def test_known_metadata_corrections_preserve_original_record(self):
        old = self.record()
        saved = deepcopy(old)
        current = deepcopy(old)
        current["claim_boundary"] = boundary_text(61)
        current["higher_order_gauge"]["name"] = "selected-coordinate-zero"
        self.assertEqual(comparison_record(old), comparison_record(current))
        self.assertEqual(old, saved)

    def test_changed_arithmetic_and_field_are_still_rejected(self):
        old = self.record()
        for key, replacement in (
            ("field", {"prime": 31, "rho": 14, "y": 3}),
            ("variable_series", {"x": [1, 3]}),
            ("root_series_digest_sha256", "changed"),
            ("truncated_lambda_one_evaluation", {"nonzero_total": 0}),
        ):
            with self.subTest(key=key):
                self.assertNotEqual(comparison_record(old), comparison_record({**old, key: replacement}))

    def test_gauge_alias_requires_the_exact_seven_coordinates(self):
        old = self.record()
        old["higher_order_gauge"]["variables_prescribed_zero_from_order_two"].pop()
        self.assertEqual(comparison_record(old)["higher_order_gauge"]["name"], "seven-pole-coordinate-zero")

    def test_unrecognized_boundary_or_gauge_is_not_normalized(self):
        old = self.record()
        old["claim_boundary"] = "an exact modular point"
        old["higher_order_gauge"]["name"] = "another-gauge"
        result = comparison_record(old)
        self.assertEqual(result["claim_boundary"], "an exact modular point")
        self.assertEqual(result["higher_order_gauge"]["name"], "another-gauge")


if __name__ == "__main__":
    unittest.main()

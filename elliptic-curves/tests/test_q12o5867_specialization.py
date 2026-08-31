#!/usr/bin/env python3
"""Exact unit tests for the q12o5867 specialization adapter."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    homogeneous_value,
    load_q12o5867_data,
    normalize_projective_parameter,
    short_certificate_model,
)
from elliptic_candidate_record import (  # noqa: E402
    change_weierstrass_model,
    is_on_weierstrass_curve,
    source_point_to_target,
)


class Q12O5867SpecializationTests(unittest.TestCase):
    def test_projective_normalization_and_homogeneous_charts(self) -> None:
        self.assertEqual(normalize_projective_parameter(-534, 1694), (-267, 847))
        self.assertEqual(normalize_projective_parameter(3, -6), (-1, 2))
        self.assertEqual(normalize_projective_parameter(-9, 0), (1, 0))
        with self.assertRaises(ValueError):
            normalize_projective_parameter(0, 0)
        coefficients = (Fraction(2), Fraction(3), Fraction(5))
        self.assertEqual(homogeneous_value(coefficients, 2, 3, 2), 56)
        self.assertEqual(homogeneous_value(coefficients, 1, 0, 2), 5)

    def test_minimal_general_model_to_integral_short_model(self) -> None:
        model = tuple(Fraction(value) for value in (1, -1, 1, -10, 20))
        point = (Fraction(0), Fraction(4))
        self.assertTrue(is_on_weierstrass_curve(model, point))
        short_model, change = short_certificate_model(model)
        self.assertEqual(short_model[:3], (0, 0, 0))
        self.assertEqual(change_weierstrass_model(model, change), short_model)
        self.assertTrue(
            is_on_weierstrass_curve(
                short_model, source_point_to_target(point, change)
            )
        )

    def test_published_sections_reconstruct_from_quadratic_chords(self) -> None:
        data = load_q12o5867_data(
            ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json",
            ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json",
        )
        self.assertEqual(data.coordinate, "elkies_2026_published_t")
        for a, b in ((-2, 377), (-308, 251), (2456, 135), (-9529, 5471), (1, 0)):
            specialization = evaluate_projective_specialization(data, a, b)
            self.assertEqual(len(specialization.points), 17)
            self.assertTrue(
                all(
                    is_on_weierstrass_curve(specialization.model, point)
                    for point in specialization.points
                )
            )

    @unittest.skipUnless(
        (ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json").exists()
        and (
            ROOT
            / "artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json"
        ).exists(),
        "local exact q12o5867 artifacts are unavailable",
    )
    def test_real_finite_and_infinity_section_identities(self) -> None:
        sys.set_int_max_str_digits(0)
        data = load_q12o5867_data(
            ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json",
            ROOT
            / "artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json",
        )
        for a, b in ((-267, 847), (1, 0)):
            specialization = evaluate_projective_specialization(data, a, b)
            self.assertEqual(len(specialization.points), 17)
            self.assertTrue(
                all(
                    is_on_weierstrass_curve(specialization.model, point)
                    for point in specialization.points
                )
            )


if __name__ == "__main__":
    unittest.main()

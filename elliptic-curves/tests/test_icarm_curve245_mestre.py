#!/usr/bin/env python3
"""Focused exact tests for the recovered Mestre parent of ICARM #245."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    change_weierstrass_model,
    weierstrass_invariants,
)
from icarm_curve245_mestre import (  # noqa: E402
    AFFINE_SCALE,
    ANCHOR_SHORT_TO_PUBLIC_CHANGE,
    CANONICAL_PARAMETER,
    CANONICAL_ROOTS,
    CONSTRUCTION,
    FERMIGIER_U,
    FERMIGIER_V,
    NATIVE_PARAMETER,
    NATIVE_ROOTS,
    PUBLIC_MODEL,
    extra_quartic_point,
    fermigier_roots,
    primitive_short_model,
)


Q = Fraction


class IcarmCurve245MestreTests(unittest.TestCase):
    def test_recovered_fermigier_parameters_and_affine_normalization(self) -> None:
        self.assertEqual(
            tuple(sorted(fermigier_roots(FERMIGIER_U, FERMIGIER_V))),
            NATIVE_ROOTS,
        )
        translated = tuple(
            int(AFFINE_SCALE * (root - min(NATIVE_ROOTS)))
            for root in sorted(NATIVE_ROOTS)
        )
        self.assertEqual(translated, CANONICAL_ROOTS)
        self.assertEqual(AFFINE_SCALE * NATIVE_PARAMETER, CANONICAL_PARAMETER)
        self.assertEqual(CONSTRUCTION.quartic_condition, 0)

    def test_pinned_short_family_coefficients(self) -> None:
        for parameter in (Q(1), Q(7, 3), Q(14), Q(-17, 5), CANONICAL_PARAMETER):
            self.assertEqual(
                primitive_short_model(parameter),
                CONSTRUCTION.primitive_jacobian_coefficients(parameter),
            )

    def test_extra_generic_section(self) -> None:
        for parameter in (Q(1), Q(7, 3), Q(14), CANONICAL_PARAMETER):
            x_value, y_value = extra_quartic_point(parameter)
            coefficients = CONSTRUCTION.primitive_quartic_coefficients(parameter)
            value = Q(0)
            for coefficient in reversed(coefficients):
                value = value * x_value + coefficient
            self.assertEqual(y_value**2, value)

    def test_anchor_is_exactly_the_public_model(self) -> None:
        short = primitive_short_model(CANONICAL_PARAMETER)
        change = WeierstrassChange.from_values(ANCHOR_SHORT_TO_PUBLIC_CHANGE)
        self.assertEqual(change_weierstrass_model(short, change), PUBLIC_MODEL)
        short_invariants = weierstrass_invariants(short)
        public_invariants = weierstrass_invariants(PUBLIC_MODEL)
        self.assertEqual(
            short_invariants["c4"]**3 / short_invariants["discriminant"],
            public_invariants["c4"]**3 / public_invariants["discriminant"],
        )


if __name__ == "__main__":
    unittest.main()

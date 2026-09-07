#!/usr/bin/env python3
"""Pure exact tests for the q12o5867 BNF-free descent bridge."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/build_q12o5867_bnf_free_signature.py"
SPEC = importlib.util.spec_from_file_location("q12_bnf_signature", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Q12O5867BNFFreeSignatureTests(unittest.TestCase):
    def test_generalized_to_monic_cubic_identity(self) -> None:
        model = (1, 2, 3, -4, 5)
        coefficients = MODULE.monic_cubic_coefficients(model)
        self.assertEqual(coefficients, (464, -40, 9, 1))
        point = (Fraction(0), Fraction(5))
        # This point is not required to lie on the example curve; test the
        # coordinate formula separately from the identity below.
        self.assertEqual(MODULE.point_on_monic_cubic(model, point), (0, 52))

    def test_point_transport_on_integral_curve(self) -> None:
        # y^2 + y = x^3 - x has P=(0,0).
        model = (0, 0, 1, -1, 0)
        coefficients = MODULE.monic_cubic_coefficients(model)
        transformed = MODULE.point_on_monic_cubic(model, (0, 0))
        self.assertEqual(coefficients, (16, -16, 0, 1))
        self.assertEqual(transformed, (0, 4))
        self.assertEqual(
            transformed[1] ** 2,
            MODULE.evaluate_cubic(coefficients, transformed[0]),
        )


if __name__ == "__main__":
    unittest.main()

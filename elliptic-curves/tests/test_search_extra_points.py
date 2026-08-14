#!/usr/bin/env python3
"""Focused dependency-free tests for the bounded extra-point search helpers."""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from fermigier_mestre import FermigierMestreFamily  # noqa: E402
from search_extra_points import (  # noqa: E402
    exact_point_record,
    parse_parameters,
    parse_point_vector,
    parse_precisions,
    signless_quartic_points,
)


class ExtraPointHelperTests(unittest.TestCase):
    def test_parameter_and_precision_parsing(self) -> None:
        self.assertEqual(parse_parameters("1666/9,70/223"), ("1666/9", "70/223"))
        self.assertEqual(parse_precisions("96,192"), (96, 192))
        for value in ("", "not-rational"):
            with self.assertRaises(argparse.ArgumentTypeError):
                parse_parameters(value)
        for value in ("96", "192,96", "96,96", "96,not-an-integer"):
            with self.assertRaises(argparse.ArgumentTypeError):
                parse_precisions(value)

    def test_gp_point_parsing_and_signless_deduplication(self) -> None:
        output = "noise [[1/2, 3/4], [1/2, -3/4], [-5, 7], [-5, -7]] tail"
        parsed = parse_point_vector(output)
        self.assertEqual(
            parsed,
            (
                (Fraction(1, 2), Fraction(3, 4)),
                (Fraction(1, 2), Fraction(-3, 4)),
                (Fraction(-5), Fraction(7)),
                (Fraction(-5), Fraction(-7)),
            ),
        )
        self.assertEqual(signless_quartic_points(parsed), (parsed[0], parsed[2]))

    def test_exact_record_checks_both_models(self) -> None:
        parameter = Fraction(1666, 9)
        quartic_point = FermigierMestreFamily.known_quartic_points(parameter)[0]
        jacobian_point = FermigierMestreFamily.quartic_point_to_jacobian(
            parameter, quartic_point
        )
        record = exact_point_record(parameter, quartic_point, jacobian_point)
        self.assertEqual(record["quartic_residual"], "0")
        self.assertEqual(record["jacobian_residual"], "0")
        self.assertEqual(record["quartic_x"], str(quartic_point[0]))
        self.assertEqual(record["jacobian_x"], str(jacobian_point[0]))

        invalid_jacobian = (jacobian_point[0], jacobian_point[1] + 1)
        with self.assertRaises(AssertionError):
            exact_point_record(parameter, quartic_point, invalid_jacobian)


if __name__ == "__main__":
    unittest.main()

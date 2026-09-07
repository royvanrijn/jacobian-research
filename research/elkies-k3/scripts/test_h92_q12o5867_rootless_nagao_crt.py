#!/usr/bin/env python3
"""Tests for the q12/orbit5867 projective CRT/Gauss constructor."""

from __future__ import annotations

import importlib.util
from math import gcd
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "construct_h92_q12o5867_rootless_nagao_crt.py"
SPEC = importlib.util.spec_from_file_location("q12o5867_crt", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CRT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CRT
SPEC.loader.exec_module(CRT)


class ProjectiveCRTNagaoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = CRT.load_family_model()

    def test_mixed_finite_infinity_lattice_constraints(self) -> None:
        # a=2b mod 5 and b=0 mod 7.
        basis = CRT.congruence_basis(2, 5, 7)
        self.assertEqual(abs(basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0]), 35)
        representatives = CRT.short_representatives(2, 5, 7, coefficient_radius=4)
        self.assertTrue(representatives)
        for numerator, denominator, _ in representatives:
            self.assertEqual((numerator - 2 * denominator) % 5, 0)
            self.assertEqual(denominator % 7, 0)
            self.assertEqual(gcd(abs(numerator), denominator), 1)

    def test_small_projective_beam_is_deterministic(self) -> None:
        groups = (
            (
                CRT.ProjectiveChoice(5, 2, 100, -1),
                CRT.ProjectiveChoice(5, 5, 90, 1),
            ),
            (
                CRT.ProjectiveChoice(7, 3, 110, -2),
                CRT.ProjectiveChoice(7, 7, 80, 2),
            ),
        )
        first = CRT.beam_combine_projective(
            groups, beam_width=4, height_weight=0.01, beam_coefficient_radius=3
        )
        second = CRT.beam_combine_projective(
            groups, beam_width=4, height_weight=0.01, beam_coefficient_radius=3
        )
        self.assertEqual(first, second)
        self.assertTrue(any(state.infinity_modulus > 1 for state in first))
        parameters = CRT.enumerate_beam_parameters(
            first, coefficient_radius=4, minimum_height=1
        )
        self.assertEqual(len(parameters), len({parameter[:2] for parameter in parameters}))

    def test_exact_specialization_clears_denominators_and_is_nonsingular(self) -> None:
        result = CRT.exact_specialization(
            self.model,
            {"parameter": "1/1", "projective_pair": [1, 1]},
        )
        source = result["source_short_model"]
        integral = result["denominator_cleared_integral_short_model"]
        self.assertTrue(source["nonsingular"])
        self.assertTrue(integral["nonsingular"])
        self.assertNotEqual(int(integral["discriminant"]), 0)
        self.assertFalse(integral["minimality_claimed"])

    def test_prime_interval_cli_helper_is_inclusive(self) -> None:
        self.assertEqual(CRT.primes_in_interval(199, 211), (199, 211))


if __name__ == "__main__":
    unittest.main()

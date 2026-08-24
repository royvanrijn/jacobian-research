#!/usr/bin/env python3
"""Lightweight replay tests for the pinned Fermigier bidegree-(2,1) pilot."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
SCRIPT = ROOT / "elliptic-curves/cas/analyze_fermigier_bidegree21_pilot.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "elliptic_fermigier_bidegree21_p13_r20e1_pilot.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "44a3f90c416c6d0cba5e63cc7dab9ca13fb4300b9570a162b47b6e27a6c339bf"
)
EXPECTED_ARTIFACT_SHA256 = (
    "423bec6bd9545783da0a550c1abb44bb6ac096c361011976ab3b209028341bae"
)

SPEC = importlib.util.spec_from_file_location("fermigier_bidegree21_pilot", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierBidegree21PilotTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.polynomial, cls.symbols = MODULE.pencil_polynomial(
            Fraction(256714, 39), Fraction(-8545)
        )

    def test_pinned_files_and_result_digest(self) -> None:
        self.assertEqual(sha256_file(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256_file(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        stable = dict(self.data)
        digest = stable.pop("result_sha256")
        stable.pop("generated_at_utc")
        self.assertEqual(
            digest,
            hashlib.sha256(
                json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
        )

    def test_generic_degree_and_discriminant_factor_degree(self) -> None:
        T, c, k = (
            self.symbols["T"],
            self.symbols["c"],
            self.symbols["k"],
        )
        generic = sp.Poly(
            self.polynomial.as_expr(), T, domain=sp.QQ.frac_field(c, k)
        )
        self.assertEqual(generic.degree(), 10)
        self.assertEqual(
            sp.factor(generic.LC()),
            6666883836179368888567109693610000
            * (c - k) ** 2
            * (c + k) ** 2,
        )
        factors = self.data["discriminant_factorization_over_QQ_c_k"]["factors"]
        self.assertEqual(
            sum(item["total_degree"] * item["exponent"] for item in factors),
            72,
        )
        discriminant = self.data["discriminant_factorization_over_QQ_c_k"]
        self.assertEqual(discriminant["nonlinear_nonconstant_factor_count"], 1)
        self.assertEqual(discriminant["nonlinear_total_degree"], 32)
        self.assertEqual(discriminant["nonlinear_exponent"], 1)

    def test_rational_component_squareclasses(self) -> None:
        records = {
            item["component"]: item
            for item in self.data["rational_linear_components"]
        }
        self.assertEqual(records["k-c=0"]["squareclass_kernel_degree"], 8)
        self.assertEqual(records["k+c=0"]["squareclass_kernel_degree"], 8)
        self.assertEqual(
            records["5899690*c+732683*k=0"]["squareclass_kernel_degree"],
            6,
        )
        self.assertTrue(
            all(not item["genus_at_most_one"] for item in records.values())
        )

    def test_cancellation_and_invalid_pole_identities(self) -> None:
        T, c, k = (
            self.symbols["T"],
            self.symbols["c"],
            self.symbols["k"],
        )
        numerator = self.symbols["numerator"]
        denominator = self.symbols["denominator"]
        cancellation = numerator.subs(
            k, -sp.Rational(5899690, 732683) * c
        )
        remainder = sp.rem(
            sp.Poly(cancellation, T, domain=sp.QQ.frac_field(c)),
            sp.Poly(denominator, T, domain=sp.QQ.frac_field(c)),
        )
        self.assertTrue(remainder.is_zero)

        for c_value, anchor in (
            (-sp.Rational(10, 28917), sp.Rational(28917, 10)),
            (-sp.Rational(39, 39508), sp.Rational(39508, 39)),
        ):
            self.assertEqual(denominator.subs({c: c_value, T: anchor}), 0)
            self.assertEqual(numerator.subs({c: c_value, T: anchor}), 0)

    def test_outcome_is_negative_and_bounded(self) -> None:
        self.assertEqual(self.data["scope"]["completed_pair_count"], 1)
        self.assertFalse(self.data["scope"]["all_pairs_classified"])
        self.assertFalse(self.data["outcome"]["target_met"])
        self.assertEqual(self.data["outcome"]["genus_at_most_one_components"], 0)


if __name__ == "__main__":
    unittest.main()

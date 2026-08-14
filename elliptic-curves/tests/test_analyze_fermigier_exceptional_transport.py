#!/usr/bin/env python3
"""Focused tests for exact Fermigier exceptional transport."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/analyze_fermigier_exceptional_transport.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic_fermigier_exceptional_transport.json"
SPEC = importlib.util.spec_from_file_location("fermigier_exceptional_transport", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ExceptionalTransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_exact_quotient_dimensions_and_relation(self) -> None:
        e22 = self.data["exceptional_quotients"]["E22"]
        rank20 = self.data["exceptional_quotients"]["rank20"]
        self.assertEqual(e22["exceptional_quotient_rank_lower_bound"], 10)
        self.assertEqual(rank20["exceptional_quotient_rank_lower_bound"], 8)
        self.assertTrue(e22["P6_exact_relation"]["verified_by_exact_group_law"])
        self.assertEqual(
            e22["independent_exceptional_labels_modulo_generic"],
            [f"P{index}" for index in range(13, 23)],
        )

    def test_affine_and_quadratic_structural_gates(self) -> None:
        self.assertEqual(self.data["transport"]["cross_anchor_pair_count"], 88)
        affine = self.data["transport"]["affine"]
        quadratic = self.data["transport"]["quadratic"]
        self.assertEqual(
            affine["histogram"],
            [{"count": 88, "degree": 6, "factor_signature": [[6, 1]], "genus": 2}],
        )
        self.assertEqual(quadratic["histogram"][0]["count"], 88)
        self.assertEqual(
            quadratic["histogram"][0]["discriminant_factor_signature"],
            [[1, 16], [32, 1]],
        )
        self.assertEqual(
            quadratic["histogram"][0]["rational_collision_parameters"], ["0"]
        )
        self.assertEqual(quadratic["low_genus_candidates"], [])

    def test_complete_mobius_structural_gate(self) -> None:
        mobius = self.data["transport"]["mobius"]
        self.assertEqual(mobius["pair_count"], 88)
        self.assertEqual(len(mobius["records"]), 88)
        self.assertEqual(mobius["unexpected_rational_degenerations"], [])
        self.assertEqual(mobius["low_genus_candidates"], [])
        histogram = mobius["histogram"]
        self.assertEqual(len(histogram), 1)
        self.assertEqual(histogram[0]["count"], 88)
        self.assertEqual(histogram[0]["generic_degree"], 10)
        self.assertEqual(
            histogram[0]["discriminant_factor_signature"],
            [[1, 12], [1, 12], [1, 16], [32, 1]],
        )
        self.assertEqual(histogram[0]["infinity_degree"], 10)
        self.assertEqual(histogram[0]["infinity_factor_signature"], [[10, 1]])

    def test_actual_fiber_products(self) -> None:
        products = self.data["fiber_products"]
        self.assertEqual(products["pair_count"], 3160)
        self.assertEqual(
            products["endpoint_strata"],
            {
                "distinct-at-both-anchors": 2520,
                "shared-E22-endpoint": 280,
                "shared-rank20-endpoint": 360,
            },
        )
        self.assertEqual(
            products["histogram"],
            [
                {
                    "common_branch_gcd_degree": 0,
                    "count": 3160,
                    "fiber_product_genus": 9,
                    "third_quotient_genus": 5,
                    "third_quotient_squarefree_degree": 12,
                }
            ],
        )
        self.assertEqual(products["low_genus_third_quotients"], [])

    def test_result_digest(self) -> None:
        self.assertEqual(self.data["result_sha256"], MODULE.result_digest(self.data))


if __name__ == "__main__":
    unittest.main()

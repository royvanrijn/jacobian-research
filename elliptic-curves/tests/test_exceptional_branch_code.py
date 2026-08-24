#!/usr/bin/env python3
"""Focused tests for the exact exceptional branch-divisor code."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))
from branch_divisor_code import (  # noqa: E402
    SquareCondition,
    analyze_square_conditions,
    conditions_from_json_records,
)

SCRIPT = CAS / "analyze_exceptional_branch_code.py"
ARTIFACT = ROOT / "artifacts/generated-results/elliptic-curves/elliptic_exceptional_branch_code.json"
SPEC = importlib.util.spec_from_file_location("exceptional_branch_code", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class BranchDivisorCodeTest(unittest.TestCase):
    def test_synthetic_cancellation_and_genus(self) -> None:
        parameter = sp.symbols("u")
        first = sp.Poly(parameter * (parameter - 1), parameter, domain=sp.QQ)
        second = sp.Poly(parameter * (parameter - 2), parameter, domain=sp.QQ)
        third = sp.Poly(first.as_expr() * second.as_expr(), parameter, domain=sp.QQ)
        data = analyze_square_conditions(
            [
                SquareCondition.polynomial("a", first),
                SquareCondition.polynomial("b", second),
                SquareCondition.polynomial("a_plus_b", third),
            ]
        )
        self.assertEqual(data["code_dimension"], 2)
        self.assertEqual(
            data["complete_branch_cancellations"],
            [{"combination": ["a", "b", "a_plus_b"], "geometric_branch_count": 0, "quadratic_quotient_genus": None}],
        )
        self.assertIsNone(data["minimum_quadratic_quotient"]["genus"])
        self.assertEqual(data["multiquadratic_cover"]["genus"], 0)

    def test_exact_json_condition_interchange(self) -> None:
        conditions = conditions_from_json_records(
            [
                {
                    "label": "rational",
                    "numerator_coefficients_ascending": ["-1", "0", "1"],
                    "denominator_coefficients_ascending": ["1", "1"],
                }
            ],
            parameter_name="q",
        )
        self.assertEqual(str(conditions[0].numerator.as_expr()), "q**2 - 1")
        self.assertEqual(str(conditions[0].denominator.as_expr()), "q + 1")
        self.assertEqual(
            analyze_square_conditions(conditions)["minimum_quadratic_quotient"]["genus"],
            0,
        )

    def test_geometric_code_automorphisms_expand_factor_degree(self) -> None:
        parameter = sp.symbols("u")
        data = analyze_square_conditions(
            [
                SquareCondition.polynomial(
                    "quadratic_orbit",
                    sp.Poly(parameter**2 + 1, parameter, domain=sp.QQ),
                )
            ]
        )
        self.assertEqual(
            data["incidence_indistinguishable_closed_branch_places"], []
        )
        self.assertEqual(
            data["geometric_incidence_indistinguishable_branch_places"],
            [
                {
                    "incidence_pattern_binary": "1",
                    "closed_branch_places": ["B1"],
                    "geometric_branch_place_count": 2,
                }
            ],
        )
        self.assertEqual(
            data["geometric_incidence_code_automorphisms"][
                "coordinate_permutation_group_order"
            ],
            2,
        )

    def test_bounded_subcutoff_preimages_recheck_at_their_anchors(self) -> None:
        parameter, x = sp.symbols("T X")
        for roots, anchor, preimages in (
            (
                MODULE.ICARM243_ROOTS,
                MODULE.ICARM243_ANCHOR,
                MODULE.ICARM243_BOUNDED_ACCIDENTAL_X,
            ),
            (
                MODULE.ICARM226_ROOTS,
                MODULE.ICARM226_ANCHOR,
                MODULE.ICARM226_BOUNDED_ACCIDENTAL_X,
            ),
        ):
            expression = MODULE.mestre_quartic_expression(roots, parameter, x)
            conditions = MODULE.bounded_observed_conditions(
                [(f"test_{index}", value) for index, value in enumerate(preimages)],
                expression,
                parameter,
                x,
                anchor,
            )
            self.assertEqual(len(conditions), len(preimages))

    def test_pinned_artifact(self) -> None:
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(data["result_sha256"], MODULE.result_digest(data))
        for name, expected_dimension in (
            ("fermigier_E22", 10),
            ("fermigier_rank20_anchor", 8),
            ("fermigier_E22_rank20_combined", 18),
            ("icarm_245", 8),
            ("icarm_275", 8),
        ):
            family = data["families"][name]
            self.assertEqual(family["code_dimension"], expected_dimension)
            self.assertEqual(family["minimum_quadratic_quotient"]["genus"], 2)
            self.assertEqual(family["low_genus_quadratic_quotients"], [])
            self.assertEqual(family["complete_branch_cancellations"], [])
            self.assertEqual(family["shared_branch_pairs"], [])
        for name, expected_dimension in (
            ("icarm_243_observed_preimages_h200000", 16),
            ("icarm_226_observed_preimages_h200000", 19),
        ):
            family = data["supplemental_bounded_input_families"][name]
            self.assertEqual(family["code_dimension"], expected_dimension)
            self.assertEqual(family["minimum_quadratic_quotient"]["genus"], 2)
            self.assertEqual(family["low_genus_quadratic_quotients"], [])
            self.assertEqual(family["complete_branch_cancellations"], [])
            self.assertEqual(family["shared_branch_pairs"], [])


if __name__ == "__main__":
    unittest.main()

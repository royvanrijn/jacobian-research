#!/usr/bin/env python3
"""Regression tests for the partial finite-class representative screen."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from search_nagao_rank21_productive_auxiliary_partial_classes import (  # noqa: E402
    PAIR_REPRESENTATIVES_PER_TARGET_PAIR,
    PARENT_ENGINE_SHA256,
    SINGLE_REPRESENTATIVES_PER_TARGET_POINT,
    vector_priority,
)


SCRIPT = CAS / "search_nagao_rank21_productive_auxiliary_partial_classes.py"
PARENT = CAS / "search_nagao_rank21_productive_auxiliary_local_classes.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_productive_auxiliary_partial_classes.json"
)


class ProductiveAuxiliaryPartialClassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_priority_and_parent_hash(self) -> None:
        self.assertLess(vector_priority((0, 1, 0)), vector_priority((1, 1, 0)))
        self.assertLess(vector_priority((-1, 0)), vector_priority((1, 0)))
        self.assertEqual(hashlib.sha256(PARENT.read_bytes()).hexdigest(), PARENT_ENGINE_SHA256)

    def test_complete_product_and_selection_counts(self) -> None:
        finite = self.data["finite_product"]
        self.assertEqual(finite["group_orders"], [20, 48, 90])
        self.assertEqual(finite["complete_state_count"], 86400)
        selection = self.data["selection"]
        self.assertEqual(
            selection["single_representatives_per_target_finite_point"],
            SINGLE_REPRESENTATIVES_PER_TARGET_POINT,
        )
        self.assertEqual(
            selection["pair_representatives_per_target_finite_pair"],
            PAIR_REPRESENTATIVES_PER_TARGET_PAIR,
        )
        self.assertEqual(
            selection["single_target_finite_point_counts"],
            {"p13": 4, "p37": 4, "p83": 12},
        )
        self.assertEqual(
            selection["pair_target_finite_point_counts"],
            {"p13_p37": 16, "p13_p83": 48, "p37_p83": 48},
        )
        self.assertEqual(selection["unique_selected_vector_count"], 1315)

    def test_exact_decontaminated_frontier(self) -> None:
        exact = self.data["exact_lifts"]
        self.assertEqual(exact["exceptional_vector_count"], 0)
        self.assertEqual(exact["distinct_parameter_x_incidence_count"], 1274)
        self.assertEqual(exact["generic_visible_incidence_count"], 4)
        self.assertEqual(exact["published_fiber_incidence_count"], 1)
        self.assertEqual(exact["decontaminated_distinct_parameter_count"], 1269)
        self.assertEqual(exact["proxy_below_gate_count"], 0)
        self.assertAlmostEqual(
            exact["minimum_log_radical_upper_proxy"], 327.55666183774576
        )
        self.assertEqual(exact["top_32"][0]["parameter"], "44246851/613596")

    def test_artifact_hash_and_scope(self) -> None:
        self.assertFalse(self.data["target_hit"])
        self.assertEqual(self.data["conclusion"]["exact_conductor_call_count"], 0)
        self.assertEqual(
            self.data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()

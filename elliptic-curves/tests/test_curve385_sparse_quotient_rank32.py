#!/usr/bin/env python3
"""Regression tests for the curve-385 sparse quotient rank-32 protocol."""

from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
PROFILE = ART / "curve385_quotient_weight_profile_v1.json"
PROTOCOL = ART / "curve385_sparse_quotient_rank32_protocol_v1.json"
POLICY_SOURCE = CAS / "curve385_sparse_quotient_policy.py"
RUNNER = CAS / "run_curve385_sparse_quotient_rank32_search.sage"

PROFILE_SHA256 = "c321d1b40d9e5fc77ebff64e5d6584feeab5f503b13eadda4f6d524d0e38162a"
PROTOCOL_SHA256 = "2c9150f50f305b8aa3763590cd5e81c4d7e121f9373177827780789ce472834f"
DEFINITION_SHA256 = "5723679da2907e036095f90376cdabde457a4f7ba5bc284ad4a4ca3edea1aa37"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_module():
    spec = importlib.util.spec_from_file_location("curve385_sparse_policy_test", POLICY_SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load sparse quotient policy")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Curve385SparseQuotientRank32Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.profile = json.loads(PROFILE.read_text())
        cls.protocol = json.loads(PROTOCOL.read_text())
        cls.policy = load_module()

    def test_frozen_artifacts(self) -> None:
        self.assertEqual(digest(PROFILE), PROFILE_SHA256)
        self.assertEqual(digest(PROTOCOL), PROTOCOL_SHA256)
        self.assertEqual(self.protocol["protocol_definition_hash"], DEFINITION_SHA256)
        self.assertEqual(
            self.policy.canonical_hash(self.protocol["protocol_definition"]),
            DEFINITION_SHA256,
        )
        runner_source = RUNNER.read_text()
        self.assertIn(f'EXPECTED_PROTOCOL_DEFINITION_HASH = "{DEFINITION_SHA256}"', runner_source)

    def test_exact_weight_profile(self) -> None:
        self.assertEqual(self.profile["status"], "PASS_EXACT_POSTHOC_WEIGHT_PROFILE")
        split = self.profile["primitive_split"]
        self.assertEqual(
            (split["generic_rank"], split["initial_discovered_rank"], split["new_complement_rank"]),
            (17, 20, 3),
        )
        self.assertTrue(split["initial_basis_is_generic_basis_followed_by_complement"])
        cumulative = self.profile["cylinder"]["cumulative_weight_profile"]
        self.assertEqual([row["chart_count"] for row in cumulative], [129, 258, 301])
        self.assertEqual([row["quotient_rank_over_M20"] for row in cumulative], [7, 9, 9])
        events = self.profile["cylinder"]["basis_extension_events_in_search_priority_order"]
        self.assertEqual(len(events), 9)
        self.assertTrue(all(row["minimum_quotient_weight"] <= 2 for row in events))

    def test_basis_dependence_is_explicit(self) -> None:
        audit = self.profile["basis_dependence_audit"]
        self.assertEqual(audit["unordered_F2_basis_count"], 28)
        self.assertEqual(audit["natural_basis_weight_one_rank"], 7)
        self.assertEqual(audit["best_weight_one_rank"], 9)
        self.assertEqual(audit["weight_at_most_two_full_rank_basis_count"], 20)
        self.assertEqual(audit["weight_at_most_two_deficient_basis_count"], 8)

    def test_primary_12_bit_stages_are_small_and_complete(self) -> None:
        plan = self.protocol["stage_plans_by_quotient_bit_count"]["12"]
        self.policy.validate_stage_plan(plan, 12, 43)
        self.assertEqual(plan[0]["id"], "natural-weight-1")
        self.assertEqual(plan[0]["new_chart_count"], 516)
        self.assertEqual(plan[1]["id"], "natural-weight-2")
        self.assertEqual(plan[1]["new_chart_count"], 2_838)
        self.assertEqual(plan[1]["cumulative_chart_count"], 3_354)
        self.assertLess(plan[3]["cumulative_chart_count"], 12_814)
        full = self.profile["rank32_scale_projection"]["staged_chart_counts"][-1]
        self.assertEqual(full["chart_count"], 176_085)

    def test_alternate_bases_are_precommitted_and_disjoint(self) -> None:
        for bit_count in (12, 13, 14):
            plan = self.protocol["stage_plans_by_quotient_bit_count"][str(bit_count)]
            seen = set()
            for row in plan:
                words = set(row["new_physical_words"])
                self.assertFalse(seen.intersection(words))
                seen.update(words)
                self.assertEqual(
                    self.policy.gf2_rank(row["basis_words_in_natural_coordinates"]),
                    bit_count,
                )
            self.assertEqual(plan[2]["basis_id"], "alternate-a")
            self.assertEqual(plan[3]["basis_id"], "alternate-b")

    def test_restart_and_failure_rules_are_fail_closed(self) -> None:
        definition = self.protocol["protocol_definition"]
        self.assertIn("restart at stage 1", definition["restart_rule"])
        self.assertEqual(
            definition["primary_campaign"]["maximum_stage_each_lattice_state"], 2
        )
        self.assertEqual(
            definition["fail_closed"]["completed_sparse_miss"],
            "bounded negative result only",
        )
        self.assertEqual(definition["fail_closed"]["unclassified_point"], "stop UNKNOWN")


if __name__ == "__main__":
    unittest.main()

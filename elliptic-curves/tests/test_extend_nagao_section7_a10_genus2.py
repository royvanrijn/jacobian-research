from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "extend_nagao_section7_a10_genus2.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_a10_genus2_extension.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "extend_nagao_section7_a10_genus2", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class A10GenusTwoExtensionTests(unittest.TestCase):
    def test_exact_inputs_and_seven_declared_plans(self) -> None:
        item, prior, metadata = MODULE.load_extension_inputs()
        self.assertEqual(item.identifier, "a10_sp08")
        self.assertEqual(len(prior), 7_631)
        self.assertEqual(
            metadata["canonical_prior_union_sha256"],
            "dac7e9cd950f6d96d17d95322be66593d672e583222f7f58720cfe642d3ff0b3",
        )
        plans = MODULE.search_plans(item)
        self.assertEqual(len(plans), 7)
        self.assertEqual(plans[0]["id"], "direct_H1000000")
        self.assertIsNone(plans[0]["matrix"])
        self.assertEqual(
            [plan["matrix"] for plan in plans[1:]],
            [
                (-162, 163, -1, 1),
                (1, 163, 0, 1),
                (164, 163, 1, 1),
                (-944, 1049, -9, 10),
                (105, 1049, 1, 10),
                (1154, 1049, 11, 10),
            ],
        )
        for plan in plans[1:]:
            a_value, b_value, c_value, d_value = plan["matrix"]
            self.assertEqual(a_value * d_value - b_value * c_value, 1)
            self.assertEqual(MODULE.Q(b_value, d_value), plan["center"])
            root = MODULE.rational_square_root(
                item.normalized.normalized_value(plan["center"])
            )
            self.assertIsNotNone(root)
            self.assertEqual(
                MODULE.polynomial_value(plan["polynomial"], MODULE.Q(0)),
                (root * d_value**3) ** 2,
            )

    def test_prior_and_generic_decontamination_are_exact(self) -> None:
        item, prior, _ = MODULE.load_extension_inputs()
        prior_root = MODULE.rational_square_root(
            item.normalized.normalized_value(MODULE.Q(163))
        )
        self.assertIsNotNone(prior_root)
        prior_record = MODULE.classify_extension_point(
            item,
            (MODULE.Q(163), prior_root),
            prior,
            plan_id="test",
        )
        self.assertEqual(prior_record["classification"], "prior_parameter_population")
        self.assertFalse(prior_record["accepted_new_parameter"])

        generic_root = MODULE.rational_square_root(
            item.normalized.normalized_value(MODULE.Q(69))
        )
        self.assertIsNotNone(generic_root)
        generic_record = MODULE.classify_extension_point(
            item,
            (MODULE.Q(69), generic_root),
            prior,
            plan_id="test",
        )
        self.assertEqual(
            generic_record["classification"],
            "known_generic_section_intersection",
        )
        self.assertIn("quartic-x:visible-00", generic_record["generic_labels"])

    def test_artifact_records_seven_timeouts_without_negative_claim(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.file_sha256(SCRIPT)
        )
        budget = data["search_budget"]
        self.assertEqual(budget["direct_height"], 1_000_000)
        self.assertEqual(budget["chart_height"], 1_000_000)
        self.assertEqual(budget["declared_auxiliary_call_cap"], 7)
        self.assertEqual(budget["actual_auxiliary_call_count"], 7)
        self.assertTrue(budget["one_call_per_plan_no_retry"])
        runs = data["search_runs"]
        self.assertEqual(runs["status_counts"], {"timeout": 7})
        self.assertEqual(len(runs["records"]), 7)
        self.assertTrue(
            all(
                record["status"] == "timeout"
                and record["mapped_signless_point_count"] == 0
                and record["one_call_no_retry"]
                for record in runs["records"]
            )
        )
        self.assertEqual(
            data["returned_point_population"]["mapped_signless_incidence_count"],
            0,
        )
        self.assertEqual(data["new_candidate_population"]["unique_parameter_count"], 0)
        self.assertEqual(data["exact_conductors"]["attempted"], 0)
        self.assertEqual(data["rank_triage"]["parameter_count"], 0)
        self.assertFalse(data["outcome"]["breakthrough_curve_found"])
        # Seven timeouts make this a completed capped attempt, not a negative
        # enumeration of any of the seven declared height boxes.
        self.assertEqual(data["claim_scope"]["rank_certificate"], False)


if __name__ == "__main__":
    unittest.main()

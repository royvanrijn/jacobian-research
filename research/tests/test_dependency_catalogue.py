"""Check real retained proof routes against the catalogue's promised dependencies."""

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dependency_graph import DependencyGraph
from dependency_sources import missing_source_closure


class RetainedDependencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "MATH_STATUS.json").read_text())
        cls.graph = DependencyGraph(cls.data, ROOT)

    def test_every_principal_scope_has_a_route_and_local_proof_closure(self):
        for topic in self.data["topics"]:
            self.assertTrue(topic["scope_notes"], topic["id"])
            for result in topic["results"]:
                self.assertTrue(self.graph.incoming[result["id"]], result["id"])
                self.assertIn(result["dependency_coverage"], {"partial", "complete"})
        self.assertEqual(missing_source_closure(self.graph), [])

    def test_included_arithmetic_proof_changes_reach_the_stated_consequence(self):
        path = "../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex"
        affected = self.graph.source_impact(path)
        self.assertIn("ARITHMETIC-COMMON-ARITHMETIC-FIBERS", affected)
        self.assertIn("ARITHMETIC-JC-NEGATION", affected)
        changed = deepcopy(self.data)
        for node in changed["graph"]["nodes"]:
            for source in node.get("sources", []):
                if source["path"] == path:
                    source["sha256"] = "0" * 64
        graph = DependencyGraph(changed, ROOT)
        self.assertEqual(graph.evaluate()["ARITHMETIC-COMMON-ARITHMETIC-FIBERS"]["state"], "stale")
        # An independent cubic collision may continue to support the same negated
        # conjecture, so its conclusion must not be unconditionally downgraded here.
        self.assertTrue(any(route["state"] == "stale"
                            for route in graph.evaluate()["ARITHMETIC-JC-NEGATION"]["routes"]))

    def test_known_tang_and_trace_inputs_cannot_be_bypassed(self):
        cases = [
            ("HESS-TANG", "HESSIAN-RICCI-FLAT"),
            ("HESS-LW", "HESSIAN-MOMENT-TRACE"),
            ("HESS-LY", "HESSIAN-MOMENT-TRACE"),
            ("HESS-KN", "HESSIAN-MOMENT-TRACE"),
        ]
        for premise, conclusion in cases:
            with self.subTest(premise=premise):
                scenario = self.graph.simulate({premise: "stale"})
                self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_real_multistep_progress_stays_hypothetical_and_preserves_assurance(self):
        before = deepcopy(self.data)
        # Fine decomposition exposes a retained-input verification gap. A
        # hypothetical proof must carry that assumption too, rather than bypass it.
        overrides = {"GAP-AUTOMATIC-T3-PERFECTION": "proved"}
        for blocker in self.graph.blockers("HESS-LZ"):
            if blocker["type"] == "node":
                overrides[blocker["id"]] = "proved"
        result = self.graph.simulate(overrides)
        for conclusion in ("GOAL-JC2", "GOAL-HC4"):
            self.assertEqual(result["nodes"][conclusion]["state"], "conditional")
            self.assertIn("GAP-AUTOMATIC-T3-PERFECTION",
                          result["nodes"][conclusion]["assumptions"])
        self.assertEqual(self.data, before)
        self.assertFalse(self.graph.evaluate()["GOAL-HC4"]["hypothetical"])

    def test_reverse_hessian_equivalence_propagates_without_extra_tang_context(self):
        result = self.graph.simulate({"GOAL-HC4": "proved"})
        self.assertEqual(result["nodes"]["GOAL-JC2"]["state"], "conditional")
        self.assertEqual(result["nodes"]["GOAL-JC2"]["assumptions"], ["GOAL-HC4"])


if __name__ == "__main__":
    unittest.main()

"""Exercise implication logic and evidence invalidation without proving mathematics."""

from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_graph import DependencyGraph, GraphValidationError


class DependencyGraphTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        self.sources = {}
        for name in ("import", "own", "reduction", "alternative", "second-step"):
            path = self.root / (name + ".md")
            path.write_text("Exact scoped evidence for " + name + ".\n")
            self.sources[name] = {
                "path": name + ".md", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "locator": "Theorem 1 and its hypotheses",
            }
        self.result = {
            "id": "CLAIM", "title": "A scoped principal claim", "state": "proved",
            "scope": "The precise conclusion in the recorded setting.",
            "canonical_source": "own.md", "dependency_coverage": "partial",
            "formal_verification": False, "external_review": False,
            "independent_replay": False,
        }
        self.data = {
            "topics": [{"id": "topic", "results": [self.result]}],
            "graph": {"schema_version": 1, "nodes": [], "implications": []},
        }

    def node(self, node_id, status="proved", source="import", **extra):
        node = {"id": node_id, "title": node_id, "kind": "external_theorem",
                "scope": "Precisely the scoped statement of " + node_id,
                "status": status, **extra}
        if source is not None:
            node["sources"] = [deepcopy(self.sources[source])]
        self.data["graph"]["nodes"].append(node)
        return node

    def rule(self, rule_id, premises, conclusion="CLAIM", status="proved", source="reduction", **extra):
        rule = {"id": rule_id, "premises": premises, "conclusion": conclusion,
                "status": status, **extra}
        if source is not None:
            rule["proof"] = deepcopy(self.sources[source])
        self.data["graph"]["implications"].append(rule)
        return rule

    def graph(self):
        return DependencyGraph(self.data, self.root, self.repository)

    def test_recorded_proved_is_not_automatically_a_graph_proof(self):
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual(record["recorded_state"], "proved")
        self.assertEqual(record["state"], "unmapped")
        self.assertEqual(record["dependency_coverage"], "partial")
        self.assertFalse(record["support"])

    def test_and_rule_needs_all_premises_and_explicit_checked_implication(self):
        self.node("IMPORTED")
        obligation = self.node("GAP", "open", None, kind="open_obligation")
        self.rule("REDUCTION", ["IMPORTED", "GAP"])
        graph = self.graph()
        self.assertEqual(graph.evaluate()["CLAIM"]["state"], "open")
        self.assertEqual({item["id"] for item in graph.blockers("CLAIM")}, {"GAP"})
        obligation["status"] = "proved"
        obligation["sources"] = [deepcopy(self.sources["own"])]
        self.data["graph"]["implications"][0]["status"] = "unverified"
        self.assertEqual(self.graph().evaluate()["CLAIM"]["state"], "open")
        self.data["graph"]["implications"][0]["status"] = "proved"
        self.assertEqual(self.graph().evaluate()["CLAIM"]["state"], "proved")

    def test_alternative_routes_are_or_and_do_not_require_every_route(self):
        self.node("GAP", "open", None)
        self.node("OWN", source="own")
        self.rule("BLOCKED", ["GAP"])
        self.rule("WORKING", ["OWN"], source="alternative")
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual(record["state"], "proved")
        self.assertEqual(record["support"], [{"type": "implication", "id": "WORKING"}])
        self.assertEqual({route["state"] for route in record["routes"]}, {"proved", "blocked"})
        self.assertEqual(record["blockers"], [])

    def test_stale_leaf_evidence_invalidates_the_entire_dependent_chain(self):
        self.node("IMPORTED")
        self.node("GOAL", "open", None, kind="goal")
        self.rule("FIRST", ["IMPORTED"])
        self.rule("SECOND", ["CLAIM"], "GOAL", source="second-step")
        graph = self.graph()
        self.assertEqual(graph.evaluate()["GOAL"]["state"], "proved")
        (self.root / "import.md").write_text("Changed argument.\n")
        evaluated = graph.evaluate()
        self.assertEqual([evaluated[key]["state"] for key in ("IMPORTED", "CLAIM", "GOAL")],
                         ["stale", "stale", "stale"])
        self.assertEqual(graph.impact("IMPORTED"), ["CLAIM", "GOAL"])
        self.assertEqual(graph.source_impact("import.md"), ["IMPORTED", "CLAIM", "GOAL"])
        self.assertEqual({item["id"] for item in graph.blockers("GOAL", evaluated)}, {"IMPORTED"})
        self.assertEqual(graph.stale_evidence(evaluated)[0]["owner"], "IMPORTED")

    def test_stale_implication_proof_blocks_route_and_records_exact_change(self):
        self.node("IMPORTED")
        self.rule("REDUCTION", ["IMPORTED"])
        graph = self.graph()
        (self.root / "reduction.md").write_text("A changed reduction.\n")
        record = graph.evaluate()["CLAIM"]
        self.assertEqual(record["state"], "stale")
        blocker = record["blockers"][0]
        self.assertEqual((blocker["id"], blocker["type"]), ("REDUCTION", "implication"))
        self.assertEqual(blocker["stale_sources"][0]["path"], "reduction.md")
        self.assertEqual(graph.source_impact("reduction.md"), ["CLAIM"])

    def test_proved_declaration_on_intermediate_node_cannot_bypass_its_premises(self):
        self.node("IMPORTED")
        self.node("INTERMEDIATE", source="own")
        self.rule("FIRST", ["IMPORTED"], "INTERMEDIATE")
        self.rule("SECOND", ["INTERMEDIATE"], source="second-step")
        graph = self.graph()
        self.assertEqual(graph.evaluate()["CLAIM"]["state"], "proved")
        (self.root / "import.md").write_text("Changed imported argument.\n")
        evaluation = graph.evaluate()
        self.assertEqual(evaluation["INTERMEDIATE"]["state"], "stale")
        self.assertEqual(evaluation["CLAIM"]["state"], "stale")

    def test_fresh_alternative_survives_and_stale_evidence_remains_visible(self):
        self.node("OLD")
        self.node("FRESH", source="own")
        self.rule("FIRST", ["OLD"])
        self.rule("ALTERNATIVE", ["FRESH"], source="alternative")
        graph = self.graph()
        (self.root / "import.md").unlink()
        evaluation = graph.evaluate()
        self.assertEqual(evaluation["CLAIM"]["state"], "proved")
        self.assertEqual(evaluation["CLAIM"]["support"][0]["id"], "ALTERNATIVE")
        stale = graph.stale_evidence(evaluation)
        self.assertEqual((stale[0]["owner"], stale[0]["reason"]), ("OLD", "missing"))

    def test_false_premise_disables_route_without_negating_its_conclusion(self):
        self.node("FALSE", "falsified")
        self.rule("REDUCTION", ["FALSE"])
        evaluation = self.graph().evaluate()
        self.assertEqual(evaluation["FALSE"]["state"], "falsified")
        self.assertEqual(evaluation["CLAIM"]["state"], "open")
        self.assertEqual(evaluation["CLAIM"]["blockers"][0]["state"], "falsified")

    def test_refuted_route_never_establishes_a_conclusion(self):
        self.node("TRUE")
        self.rule("BROKEN", ["TRUE"], status="refuted")
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual(record["state"], "open")
        self.assertEqual(record["blockers"][0]["state"], "refuted")

    def test_assumptions_propagate_until_an_unconditional_route_exists(self):
        self.node("CONDITIONAL", assumptions=["characteristic-zero hypothesis"])
        self.rule("CONDITIONAL-ROUTE", ["CONDITIONAL"], assumptions=["descent hypothesis"])
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual(record["state"], "conditional")
        self.assertEqual(record["assumptions"], ["characteristic-zero hypothesis", "descent hypothesis"])
        self.node("UNCONDITIONAL", source="own")
        self.rule("UNCONDITIONAL-ROUTE", ["UNCONDITIONAL"], source="alternative")
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual((record["state"], record["assumptions"]), ("proved", []))

    def test_assumptions_on_a_derived_conclusion_are_not_discharged_by_its_route(self):
        self.node("IMPORTED")
        self.result["assumptions"] = ["unverified descent"]
        self.rule("REDUCTION", ["IMPORTED"])
        record = self.graph().evaluate()["CLAIM"]
        self.assertEqual(record["state"], "conditional")
        self.assertEqual(record["assumptions"], ["unverified descent"])

    def test_simulation_is_read_only_and_hypothetical_assumptions_propagate(self):
        self.node("GAP", "open", None)
        self.node("GOAL", "open", None, kind="goal")
        self.rule("REDUCTION", ["GAP"])
        self.rule("NEXT", ["CLAIM"], "GOAL", source="second-step")
        original = deepcopy(self.data)
        graph = self.graph()
        scenario = graph.simulate({"GAP": "proved"})
        self.assertTrue(scenario["hypothetical"])
        for key in ("GAP", "CLAIM", "GOAL"):
            self.assertEqual(scenario["nodes"][key]["state"], "conditional")
            self.assertEqual(scenario["nodes"][key]["assumptions"], ["GAP"])
            self.assertTrue(scenario["nodes"][key]["hypothetical"])
        self.assertEqual(self.data, original)
        self.assertEqual(graph.data, original)
        self.assertEqual(graph.evaluate()["GOAL"]["state"], "open")
        self.assertFalse(self.result["formal_verification"])
        self.assertFalse(self.result["external_review"])

    def test_simulated_refutation_downgrades_dependents_without_false_negations(self):
        self.node("IMPORTED")
        self.rule("REDUCTION", ["IMPORTED"])
        scenario = self.graph().simulate({"IMPORTED": "falsified"})
        self.assertEqual(scenario["nodes"]["CLAIM"]["state"], "open")
        self.assertTrue(scenario["nodes"]["CLAIM"]["hypothetical"])
        self.assertEqual({item["id"] for item in scenario["changes"]}, {"IMPORTED", "CLAIM"})

    def test_explanation_has_exact_recursive_provenance_and_and_diamonds(self):
        self.node("IMPORT")
        self.node("OWN", source="own")
        self.rule("REDUCTION", ["IMPORT", "OWN"])
        graph = self.graph()
        explanation = graph.explain("CLAIM")
        self.assertEqual(explanation["proof"]["locator"], "Theorem 1 and its hypotheses")
        self.assertEqual({node["id"] for node in explanation["premise_provenance"]}, {"IMPORT", "OWN"})
        self.assertEqual(explanation["premise_provenance"][0]["sources"][0]["sha256"],
                         self.sources["import"]["sha256"])
        dot = graph.to_dot()
        self.assertIn('"REDUCTION" [label="REDUCTION (AND)", shape=diamond]', dot)
        self.assertIn('"IMPORT" -> "REDUCTION"', dot)
        self.assertIn('"OWN" -> "REDUCTION"', dot)

    def test_invalid_graphs_reject_dangling_duplicates_self_dependencies_and_cycles(self):
        for failure in ("dangling", "duplicate-node", "duplicate-rule", "duplicate-premise", "self", "cycle"):
            with self.subTest(failure=failure):
                self.data["graph"]["nodes"] = []
                self.data["graph"]["implications"] = []
                self.node("A")
                self.node("B", "open", None)
                if failure == "dangling":
                    self.rule("R", ["MISSING"])
                elif failure == "duplicate-node":
                    self.node("A")
                elif failure == "duplicate-rule":
                    self.rule("R", ["A"])
                    self.rule("R", ["B"])
                elif failure == "duplicate-premise":
                    self.rule("R", ["A", "A"])
                elif failure == "self":
                    self.rule("R", ["CLAIM"])
                else:
                    self.rule("R1", ["A"], "B")
                    self.rule("R2", ["B"], "A", status="unverified", source=None)
                with self.assertRaises(GraphValidationError):
                    self.graph()

    def test_established_leaves_and_rules_require_pinned_exact_evidence(self):
        for failure in ("missing-source", "missing-locator", "bad-hash", "missing-proof", "empty-premises"):
            with self.subTest(failure=failure):
                self.data["graph"]["nodes"] = []
                self.data["graph"]["implications"] = []
                node = self.node("A")
                if failure == "missing-source":
                    node["sources"] = []
                elif failure == "missing-locator":
                    node["sources"][0].pop("locator")
                elif failure == "bad-hash":
                    node["sources"][0]["sha256"] = "a receipt is not an exact pin"
                elif failure == "missing-proof":
                    self.rule("R", ["A"], source=None)
                else:
                    self.rule("R", [])
                with self.assertRaises(GraphValidationError):
                    self.graph()

    def test_graph_evidence_cannot_escape_repository_on_either_platform(self):
        node = self.node("A")
        for path in ("../../outside.md", "/tmp/outside.md", "C:/outside.md", "C:outside.md"):
            with self.subTest(path=path):
                node["sources"][0]["path"] = path
                with self.assertRaises(GraphValidationError):
                    self.graph()

    def test_result_node_ids_are_reserved_and_unknown_scenarios_are_rejected(self):
        self.node("CLAIM")
        with self.assertRaisesRegex(GraphValidationError, "duplicate node"):
            self.graph()
        self.data["graph"]["nodes"] = []
        graph = self.graph()
        with self.assertRaisesRegex(GraphValidationError, "unknown simulation node"):
            graph.simulate({"UNKNOWN": "proved"})
        with self.assertRaisesRegex(GraphValidationError, "invalid simulation state"):
            graph.simulate({"CLAIM": "published"})


if __name__ == "__main__":
    unittest.main()

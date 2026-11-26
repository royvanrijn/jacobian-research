"""Exercise implication logic and evidence invalidation without proving mathematics."""

from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_graph import DependencyGraph, GraphValidationError
from dependency_views import _blocker_text


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


    def equivalence(self, members=("A", "B"), contexts=None):
        for member in members:
            if not any(node["id"] == member for node in self.data["graph"]["nodes"]):
                self.node(member, "open", None, kind="goal")
        directions = []
        contexts = contexts or {}
        for index, premise in enumerate(members):
            conclusion = members[(index + 1) % len(members)]
            directions.append({
                "id": premise + "-TO-" + conclusion,
                "premises": [premise, *contexts.get(premise, [])],
                "conclusion": conclusion, "status": "proved",
                "proof": deepcopy(self.sources["reduction"]),
            })
        group = {"id": "EQUIVALENT-GOALS", "title": "Exact directional equivalence",
                 "members": list(members), "directions": directions}
        self.data["graph"]["equivalences"] = [group]
        return group

    def test_equivalent_open_goals_cannot_prove_each_other(self):
        self.equivalence()
        graph = self.graph()
        evaluation = graph.evaluate()
        self.assertEqual([evaluation[key]["state"] for key in ("A", "B")], ["open", "open"])
        self.assertIn("EQUIVALENT-GOALS", {blocker["id"] for blocker in graph.blockers("A")})
        self.assertEqual(graph.component_of["A"], graph.component_of["B"])
        self.assertEqual(graph.impact("A"), ["B"])
        quotient = graph.quotient()
        self.assertIn({"id": "EQUIVALENT-GOALS", "members": ["A", "B"], "equivalence": True},
                      quotient["components"])
        self.assertIn("independent grounding required", graph.to_dot())

    def test_navigation_preserves_first_blocker_witness_across_routes(self):
        self.node("GAP", "open", None)
        self.equivalence(contexts={"A": ["GAP"]})
        self.node("LEFT", "open", None)
        self.node("RIGHT", "open", None)
        self.rule("TO-LEFT", ["A"], "LEFT", assumptions=["extra condition"])
        self.rule("TO-RIGHT", ["B"], "RIGHT")
        self.rule("JOIN", ["LEFT", "RIGHT"])
        graph = self.graph()
        full = graph.evaluate()
        compact = graph.evaluate(distinct_blockers=True)
        self.assertLess(len(compact["CLAIM"]["blockers"]), len(full["CLAIM"]["blockers"]))
        for node_id, record in full.items():
            first_witnesses = {}
            for blocker in record["blockers"]:
                identity = (blocker["type"], blocker["id"], blocker["state"])
                first_witnesses.setdefault(identity, blocker)
            self.assertEqual(compact[node_id]["blockers"], list(first_witnesses.values()))
            self.assertEqual(_blocker_text(compact[node_id]), _blocker_text(record))

        overrides = {"GAP": "open"}
        scenario = graph.simulate(overrides)["nodes"]
        full_scenario = graph.evaluate(overrides)
        self.assertLess(len(scenario["CLAIM"]["blockers"]), len(full_scenario["CLAIM"]["blockers"]))
        for node_id, record in scenario.items():
            self.assertEqual(record["blockers"], compact[node_id]["blockers"])
            self.assertEqual({**record, "blockers": full_scenario[node_id]["blockers"]}, full_scenario[node_id])

    def test_refuted_equivalent_goals_have_no_open_blockers(self):
        self.node("A", "falsified")
        self.node("B", "falsified")
        self.equivalence()
        self.rule("USES-REFUTED", ["A"])
        evaluation = self.graph().evaluate()
        for member in ("A", "B"):
            self.assertEqual(evaluation[member]["state"], "falsified")
            self.assertEqual(evaluation[member]["blockers"], [])
        self.assertEqual([(blocker["id"], blocker["state"])
                          for blocker in evaluation["CLAIM"]["blockers"]], [("A", "falsified")])

    def test_proving_either_equivalent_member_propagates_without_circular_provenance(self):
        for ground in ("A", "B"):
            with self.subTest(ground=ground):
                self.data["graph"]["nodes"] = []
                self.data["graph"]["implications"] = []
                self.equivalence()
                self.node("SOURCE", source="own")
                self.rule("ZZ-GROUND", ["SOURCE"], ground, source="own")
                graph = self.graph()
                evaluated = graph.evaluate()
                for member in ("A", "B"):
                    self.assertEqual(evaluated[member]["state"], "proved")
                    explanation = graph.explain(member)
                    self.assertTrue(explanation["premise_provenance"])
                    selected = graph.rules[evaluated[member]["support"][0]["id"]]
                    self.assertTrue(all(evaluated[p]["proof_depth"] < evaluated[member]["proof_depth"]
                                        for p in selected["premises"]))

    def test_each_equivalence_direction_retains_its_own_context_premises(self):
        self.node("TANG", "open", None)
        self.equivalence(contexts={"A": ["TANG"]})
        self.node("SOURCE", source="own")
        self.rule("GROUND-A", ["SOURCE"], "A", source="own")
        graph = self.graph()
        self.assertEqual(graph.evaluate()["A"]["state"], "proved")
        self.assertEqual(graph.evaluate()["B"]["state"], "open")
        self.assertIn("TANG", {blocker["id"] for blocker in graph.blockers("B")})
        # The reverse direction needs no Tang input.
        self.data["graph"]["implications"][0]["conclusion"] = "B"
        graph = self.graph()
        self.assertEqual(graph.evaluate()["A"]["state"], "proved")
        self.assertEqual(graph.simulate({"TANG": "falsified"})["nodes"]["A"]["state"], "proved")

    def test_changed_equivalence_seed_cannot_survive_by_circulating(self):
        self.equivalence()
        self.node("SOURCE", source="own")
        self.rule("GROUND", ["SOURCE"], "A", source="own")
        graph = self.graph()
        self.assertEqual(graph.evaluate()["B"]["state"], "proved")
        (self.root / "own.md").write_text("Withdrawn source evidence.\n")
        evaluated = graph.evaluate()
        self.assertEqual([evaluated[key]["state"] for key in ("A", "B")], ["stale", "stale"])
        self.assertIn("B", graph.source_impact("own.md"))

    def test_changed_direction_does_not_invalidate_independently_grounded_member(self):
        group = self.equivalence()
        group["directions"][0]["proof"] = deepcopy(self.sources["alternative"])
        self.node("SOURCE", source="own")
        self.rule("GROUND", ["SOURCE"], "A", source="own")
        graph = self.graph()
        (self.root / "alternative.md").write_text("Changed forward equivalence proof.\n")
        evaluated = graph.evaluate()
        self.assertEqual(evaluated["A"]["state"], "proved")
        self.assertEqual(evaluated["B"]["state"], "stale")
        self.assertIn("A-TO-B", {source["owner"] for source in graph.stale_evidence()})

    def test_hypothetical_equivalence_proof_never_erases_its_assumption(self):
        self.equivalence()
        self.node("GOAL", "open", None)
        self.rule("DOWNSTREAM", ["B"], "GOAL")
        original = deepcopy(self.data)
        graph = self.graph()
        for seed in ("A", "B"):
            simulated = graph.simulate({seed: "proved"})
            for node_id in ("A", "B", "GOAL"):
                record = simulated["nodes"][node_id]
                self.assertEqual(record["state"], "conditional")
                self.assertEqual(record["assumptions"], [seed])
                self.assertTrue(record["hypothetical"])
            graph.explain("GOAL", simulated["nodes"])
        self.assertEqual(self.data, original)

    def test_conditional_equivalence_support_improvements_remain_acyclic(self):
        self.equivalence(("A", "B", "C"))
        self.node("SOURCE-A", source="own", assumptions=["X", "Y"])
        self.node("SOURCE-C", source="alternative", assumptions=["X"])
        self.rule("GROUND-A", ["SOURCE-A"], "A", source="own")
        self.rule("GROUND-C", ["SOURCE-C"], "C", source="alternative")
        graph = self.graph()
        evaluated = graph.evaluate()
        for member in ("A", "B", "C"):
            self.assertEqual(evaluated[member]["assumptions"], ["X"])
            graph.explain(member)

    def test_equivalence_explanation_freezes_exact_conditional_premise_proofs(self):
        self.node("SA", source="own", assumptions=["a"])
        self.node("SZ", source="alternative", assumptions=["z"])
        group = self.equivalence(("A", "B", "C"))
        group["directions"] = [
            {"id": "B-contextZ-to-A", "premises": ["B", "SZ"], "conclusion": "A",
             "status": "proved", "proof": deepcopy(self.sources["reduction"])},
            {"id": "C-to-B", "premises": ["C"], "conclusion": "B",
             "status": "proved", "proof": deepcopy(self.sources["reduction"])},
            {"id": "A-to-C", "premises": ["A"], "conclusion": "C",
             "status": "proved", "proof": deepcopy(self.sources["reduction"])},
        ]
        self.rule("seed-B", ["SZ"], "B", source="alternative")
        self.rule("seed-C", ["SA"], "C", source="own")
        graph = self.graph()
        evaluated = graph.evaluate()
        self.assertEqual(evaluated["A"]["assumptions"], ["z"])
        self.assertEqual(evaluated["B"]["assumptions"], ["a"])
        explanation = graph.explain("A", evaluated)
        frozen_b, frozen_sz = explanation["premise_provenance"]
        self.assertEqual(frozen_b["assumptions"], ["z"])
        self.assertEqual(frozen_b["support"][0]["id"], "seed-B")
        self.assertEqual(frozen_b["preferred_assumptions"], ["a"])
        self.assertEqual(frozen_sz["assumptions"], ["z"])
        selected_route = next(route for route in evaluated["A"]["routes"] if route.get("selected"))
        self.assertEqual(selected_route["assumptions"], ["z"])
        self.assertEqual(selected_route["preferred_premise_assumptions"], ["a", "z"])

        def verify_union(record, ancestors=frozenset()):
            self.assertNotIn(record["id"], ancestors)
            ancestors = ancestors | {record["id"]}
            support = record["support"][0]
            if support["type"] != "implication":
                return
            rule = graph.rules[support["id"]]
            selected_route = next(route for route in record["routes"] if route.get("selected"))
            self.assertEqual(selected_route["id"], support["id"])
            self.assertEqual(selected_route["assumptions"], record["assumptions"])
            required = set(graph.nodes[record["id"]].get("assumptions", []))
            required.update(rule.get("assumptions", []))
            self.assertEqual([p["id"] for p in record["premise_provenance"]], rule["premises"])
            for premise in record["premise_provenance"]:
                required.update(premise["assumptions"])
                verify_union(premise, ancestors)
            self.assertEqual(set(record["assumptions"]), required)
        for member in ("A", "B", "C"):
            verify_union(graph.explain(member, evaluated))
        # The old z-based proof cannot survive fresh loss of its seed.
        (self.root / "alternative.md").write_text("Changed z-based source.\n")
        fresh = graph.evaluate()
        self.assertEqual(fresh["A"]["state"], "stale")
        self.assertIsNone(fresh["A"]["proof_certificate"])

    def test_undeclared_and_cross_component_cycles_still_fail_validation(self):
        self.equivalence()
        self.node("OUTSIDE", "open", None)
        self.rule("OUT", ["A"], "OUTSIDE")
        self.rule("BACK", ["OUTSIDE"], "B")
        with self.assertRaisesRegex(GraphValidationError, "cycle outside equivalence"):
            self.graph()
        self.data["graph"]["implications"] = []
        self.rule("UNDOCUMENTED", ["A"], "B")
        with self.assertRaisesRegex(GraphValidationError, "must be an equivalence direction"):
            self.graph()

    def test_invalid_equivalence_groups_are_rejected(self):
        for failure in ("one-way", "overlap", "multiple-members", "dangling"):
            with self.subTest(failure=failure):
                self.data["graph"]["nodes"] = []
                group = self.equivalence()
                if failure == "one-way":
                    group["directions"].pop()
                elif failure == "overlap":
                    self.data["graph"]["equivalences"].append({**deepcopy(group), "id": "SECOND"})
                elif failure == "multiple-members":
                    group["directions"][0]["premises"] = ["A", "B"]
                else:
                    group["members"][0] = "MISSING"
                with self.assertRaises(GraphValidationError):
                    self.graph()


if __name__ == "__main__":
    unittest.main()

"""Catalogue-to-navigation propagation and read-only CLI behavior."""

from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research


class DependencyIntegrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        (self.root / "topic").mkdir()
        (self.root / "topic/README.md").touch()
        (self.root / "proof.md").write_text("Lemma A. A and B imply C. C implies D.\n")
        self.source = {
            "path": "proof.md",
            "sha256": hashlib.sha256((self.root / "proof.md").read_bytes()).hexdigest(),
            "locator": "Lemma A and its stated implications",
        }
        for name, value in (("ROOT", self.root), ("REPOSITORY", self.repository),
                            ("STATUS", self.root / "MATH_STATUS.json")):
            changed = patch.object(research, name, value)
            changed.start()
            self.addCleanup(changed.stop)
        self.data = {
            "schema_version": 2, "authority": "MATH_STATUS.json",
            "topics": [{
                "id": "topic", "title": "Example topic", "summary": "A scoped result.",
                "status": "B is still open.", "readme": "topic/README.md",
                "scope_notes": ["Only this exact family is considered."], "gap_ids": ["B"],
                "open_questions": [{"question": "Does B hold?", "source": "proof.md",
                                    "node": "B", "kind": "mathematical"}],
                "results": [{
                    "id": "D", "title": "Restricted target", "scope": "The exact restricted D.",
                    "state": "partial", "canonical_source": "proof.md", "checker": None,
                    "artifact_hash": None, "evidence": [], "proof_type": "Written implication",
                    "independent_replay": False, "formal_verification": False,
                    "external_review": False, "dependency_coverage": "partial",
                }],
            }],
            "graph": {"schema_version": 1, "nodes": [
                {"id": "A", "kind": "source", "title": "Known lemma", "scope": "A.",
                 "status": "proved", "sources": [self.source], "topics": ["topic"]},
                {"id": "B", "kind": "gap", "title": "Missing lemma", "scope": "B.",
                 "status": "open", "sources": [self.source], "topics": ["topic"],
                 "gap_kind": "mathematical"},
                {"id": "C", "kind": "claim", "title": "Intermediate", "scope": "C.",
                 "status": "open", "topics": ["topic"]},
            ], "implications": [
                {"id": "A-B-to-C", "premises": ["A", "B"], "conclusion": "C",
                 "status": "proved", "proof": self.source},
                {"id": "C-to-D", "premises": ["C"], "conclusion": "D",
                 "status": "proved", "proof": self.source},
            ]},
        }

    def write_navigation(self, data):
        for path, text in research.navigation_outputs(data).items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)

    def cli(self, *arguments):
        self.root.joinpath("MATH_STATUS.json").write_text(json.dumps(self.data))
        output = StringIO()
        with patch.object(sys, "argv", ["research.py", *arguments]), redirect_stdout(output):
            research.main()
        return output.getvalue()

    def test_recording_new_evidence_updates_all_downstream_navigation(self):
        research.validate(self.data)
        self.write_navigation(self.data)
        research.check_navigation(self.data)
        changed = deepcopy(self.data)
        changed["graph"]["nodes"][1]["status"] = "proved"
        # It is the explicit proof-source update that discharges B, not a CLI simulation.
        graph = research.dependency_graph(changed)
        evaluation = graph.evaluate()
        self.assertEqual(evaluation["C"]["state"], "proved")
        self.assertEqual(evaluation["D"]["state"], "proved")
        self.assertEqual(changed["topics"][0]["results"][0]["state"], "partial")
        self.assertFalse(changed["topics"][0]["results"][0]["formal_verification"])
        with self.assertRaisesRegex(AssertionError, "is stale"):
            research.check_navigation(changed)
        self.write_navigation(changed)
        research.check_navigation(changed)
        self.assertIn("| proved |", (self.root / "STATUS.md").read_text())
        self.assertIn("D", (self.root / "DEPENDENCIES.md").read_text())
        topic_text = (self.root / "topic/README.md").read_text()
        self.assertIn("mathematical / proved", topic_text)
        self.assertIn("Assessment update needed", topic_text)
        self.assertNotIn("Open questions:", topic_text)

    def test_generated_navigation_cannot_be_used_as_a_proof_source(self):
        data = deepcopy(self.data)
        path = self.root / "topic/README.md"
        data["graph"]["nodes"][0]["sources"] = [{
            "path": "topic/README.md", "locator": "Generated summary",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }]
        with self.assertRaisesRegex(AssertionError, "generated navigation"):
            research.validate(data)

    def test_changed_source_invalidates_a_previously_completed_chain(self):
        self.data["graph"]["nodes"][1]["status"] = "proved"
        graph = research.dependency_graph(self.data)
        self.assertEqual(graph.evaluate()["D"]["state"], "proved")
        self.write_navigation(self.data)
        (self.root / "proof.md").write_text("Corrected proof; previous evidence needs review.\n")
        changed = graph.evaluate()
        self.assertEqual(changed["D"]["state"], "stale")
        self.assertEqual(graph.source_impact("proof.md"), ["A", "B", "C", "D"])
        with self.assertRaisesRegex(AssertionError, "is stale"):
            research.check_navigation(self.data)

    def test_cli_simulation_explains_assumptions_and_never_writes_proof_status(self):
        baseline = deepcopy(self.data)
        result = json.loads(self.cli("simulate", "--prove", "B", "--json"))
        self.assertTrue(result["hypothetical"])
        self.assertEqual(result["nodes"]["D"]["state"], "conditional")
        self.assertEqual(result["nodes"]["D"]["assumptions"], ["B"])
        self.assertEqual(self.data, baseline)
        self.assertEqual(json.loads(research.STATUS.read_text()), baseline)
        actual = json.loads(self.cli("explain", "D", "--json"))
        self.assertEqual(actual["state"], "open")
        self.assertIn("B", [blocker["id"] for blocker in actual["blockers"]])

    def test_focused_export_marks_required_premises_outside_the_view(self):
        exported = json.loads(self.cli("graph", "--format", "json", "--focus", "B",
                                       "--depth", "1", "--direction", "descendants"))
        self.assertEqual(set(exported["nodes"]), {"B", "C"})
        self.assertEqual(exported["omitted_premises"], ["A"])
        self.assertFalse(exported["selection_is_complete"])
        text = self.cli("graph", "--focus", "B", "--depth", "1", "--direction", "descendants")
        self.assertIn("omitted required premises", text)

    def test_cli_exports_full_dag_and_reports_structural_impact(self):
        exported = json.loads(self.cli("graph", "--format", "json"))
        self.assertEqual(set(exported["nodes"]), {"A", "B", "C", "D"})
        impact = json.loads(self.cli("impact", "B", "--json"))
        self.assertEqual([node["id"] for node in impact], ["C", "D"])
        gaps = json.loads(self.cli("gaps", "--kind", "mathematical", "--json"))
        self.assertEqual([node["id"] for node in gaps], ["B"])
        self.assertEqual(gaps[0]["kind"], "mathematical")
        self.assertEqual(gaps[0]["node_kind"], "gap")

    def test_endpoint_inventory_covers_the_union_without_conflating_roles_or_states(self):
        self.data["topics"][0]["gap_ids"].append("IMPORTED")
        self.data["topics"].append({"id": "second", "results": [], "gap_ids": ["B"]})
        self.data["graph"]["nodes"] += [
            {"id": "UNLISTED", "kind": "goal", "title": "Open theorem", "scope": "Unlisted goal.",
             "status": "open", "topic": "topic", "next_action": "Check the full coefficient family."},
            {"id": "UNTAGGED", "kind": "goal", "title": "Untagged theorem", "scope": "Untagged goal.",
             "status": "open"},
            {"id": "ALTERNATIVE", "kind": "gap", "title": "Alternative method", "scope": "Optional method.",
             "status": "open", "endpoint_role": "alternative_method", "topics": ["second"]},
            {"id": "IMPORTED", "kind": "source", "title": "Imported theorem", "scope": "Audit the application.",
             "status": "unverified"},
            {"id": "DISCHARGED", "kind": "goal", "title": "Proved theorem", "scope": "Already proved.",
             "status": "open"},
            {"id": "REFUTED", "kind": "goal", "title": "Falsified target", "scope": "Already falsified.",
             "status": "falsified", "sources": [self.source]},
        ]
        self.data["graph"]["implications"] += [
            {"id": "A-to-discharged", "premises": ["A"], "conclusion": "DISCHARGED",
             "status": "proved", "proof": self.source},
            {"id": "alternative-to-discharged", "premises": ["ALTERNATIVE"], "conclusion": "DISCHARGED",
             "status": "proved", "proof": self.source},
        ]

        inventory = json.loads(self.cli("gaps", "--all", "--json"))
        rows = {row["id"]: row for row in inventory}
        self.assertEqual(len(inventory), 7)
        self.assertEqual(set(rows), {"B", "ALTERNATIVE", "DISCHARGED", "IMPORTED", "REFUTED", "UNLISTED", "UNTAGGED"})
        self.assertEqual(rows["B"]["topics"], ["second", "topic"])
        self.assertEqual(rows["B"]["listed_in_topics"], ["second", "topic"])
        self.assertEqual(rows["UNLISTED"]["topics"], ["topic"])
        self.assertEqual(rows["UNLISTED"]["listed_in_topics"], [])
        self.assertEqual(rows["UNLISTED"]["next_action"], "Check the full coefficient family.")
        self.assertEqual(rows["UNTAGGED"]["topics"], [])
        self.assertEqual(rows["DISCHARGED"]["state"], "proved")
        self.assertEqual(rows["DISCHARGED"]["recorded_state"], "open")
        self.assertEqual(rows["DISCHARGED"]["disposition"], "resolved")
        self.assertEqual(rows["REFUTED"]["disposition"], "falsified")
        self.assertEqual(rows["ALTERNATIVE"]["role"], "alternative_method")
        self.assertEqual(rows["IMPORTED"]["role"], "source_audit")
        self.assertEqual(rows["IMPORTED"]["node_kind"], "source")
        self.assertEqual(rows["UNLISTED"]["role"], "active_goal")

        self.assertEqual(json.loads(self.cli("gaps", "--json")), inventory)
        active = json.loads(self.cli("gaps", "--unresolved", "--role", "active_goal", "--json"))
        self.assertEqual([row["id"] for row in active], ["B", "UNLISTED", "UNTAGGED"])
        second = json.loads(self.cli("gaps", "--topic", "second", "--json"))
        self.assertEqual([row["id"] for row in second], ["ALTERNATIVE", "B"])
        text = self.cli("gaps", "--all")
        self.assertIn("proved (recorded open)", text)
        self.assertIn("UNTAGGED [untagged", text)
        self.assertIn("Next action — Check the full coefficient family.", text)

    def test_plain_explain_includes_selected_recursive_source_pins(self):
        self.data["graph"]["nodes"][1]["status"] = "proved"
        text = self.cli("explain", "D")
        self.assertIn("Source: proof.md :: Lemma A", text)
        self.assertIn("sha256:" + self.source["sha256"], text)
        self.assertIn("A: proved", text)


if __name__ == "__main__":
    unittest.main()

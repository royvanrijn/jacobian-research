"""Generated reports must expose unresolved scope and the actual implication logic."""

from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_graph import DependencyGraph
from dependency_views import render_dependency_report, render_scope_report


class DependencyViewsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        proof = self.root / "proof.md"
        proof.write_text("An exact source statement.\n")
        self.source = {"path": "proof.md", "locator": "Theorem 2, restricted to degree two",
                       "sha256": hashlib.sha256(proof.read_bytes()).hexdigest()}
        self.result = {
            "id": "CLAIM", "title": "Restricted construction", "state": "proved",
            "scope": "Only degree two over characteristic zero; no all-degree assertion.",
            "canonical_source": "proof.md", "dependency_coverage": "partial",
            "independent_replay": False, "formal_verification": True, "external_review": False,
            "verification_note": "Formal verification covers only the finite identity.",
        }
        self.data = {"topics": [{
            "id": "topic", "title": "A research topic", "summary": "Degree-bounded constructions.",
            "scope_notes": ["The all-degree implication remains unverified."], "gap_ids": ["GAP"],
            "results": [self.result],
            "open_questions": [{"question": "Does the all-degree extension hold?", "source": "proof.md",
                                "node": "GAP", "kind": "mathematical"}],
        }], "graph": {"schema_version": 1, "nodes": [
            {"id": "SOURCE", "kind": "source", "title": "Retained written proof",
             "scope": "The exact written degree-two argument.", "status": "proved",
             "sources": [deepcopy(self.source)], "topics": ["topic"]},
            {"id": "GAP", "kind": "gap", "title": "Missing extension",
             "scope": "The missing all-degree extension argument.", "status": "open", "topics": ["topic"],
             "gap_kind": "mathematical"},
            {"id": "GOAL", "kind": "goal", "title": "All-degree target",
             "scope": "The complete all-degree statement.", "status": "open", "topics": ["topic"]},
        ], "implications": [
            {"id": "SCOPED-ROUTE", "premises": ["SOURCE"], "conclusion": "CLAIM",
             "status": "proved", "proof": deepcopy(self.source)},
            {"id": "EXTEND", "premises": ["CLAIM", "GAP"], "conclusion": "GOAL",
             "status": "proved", "proof": deepcopy(self.source)},
        ]}}

    def graph(self):
        return DependencyGraph(self.data, self.root, self.repository)

    def test_scope_report_keeps_exact_scope_coverage_and_assurance_separate(self):
        graph = self.graph()
        text = render_scope_report(self.data, graph)
        self.assertIn(self.result["scope"], text)
        self.assertIn("**Graph state:** proved", text)
        self.assertIn("**Dependency coverage:** partial", text)
        self.assertIn("Formal verification: recorded", text)
        self.assertIn("External review: not recorded", text)
        self.assertIn(self.result["verification_note"], text)
        self.assertIn("The all-degree implication remains unverified.", text)
        self.assertIn("**GAP** (mathematical; open): The missing all-degree extension argument.", text)

    def test_diagrams_show_and_rules_or_alternatives_and_hidden_source_count(self):
        self.data["graph"]["implications"].append({
            "id": "ALTERNATIVE", "premises": ["GAP"], "conclusion": "CLAIM",
            "status": "unverified",
        })
        text = render_dependency_report(self.data, self.graph())
        self.assertIn("```mermaid\nflowchart LR", text)
        self.assertIn('AND<br/>SCOPED-ROUTE<br/>+1 pinned source', text)
        self.assertIn("CLAIM<br/>proved<br/>OR routes", text)
        self.assertNotIn('SOURCE<br/>proved', text)
        self.assertIn('GAP<br/>open', text)
        self.assertIn('GOAL<br/>open', text)
        self.assertIn('AND<br/>EXTEND', text)
        self.assertIn("including hidden sources, are required", text)
        self.assertIn("simulate --prove ID", text)

    def test_changed_source_reappears_and_reports_downgrade_without_assurance_upgrade(self):
        graph = self.graph()
        (self.root / "proof.md").write_text("A materially changed proof.\n")
        dependencies = render_dependency_report(self.data, graph)
        scopes = render_scope_report(self.data, graph)
        self.assertIn('SOURCE<br/>stale', dependencies)
        self.assertIn("Changed or missing evidence", dependencies)
        self.assertIn("| SOURCE | [proof.md](proof.md) | changed |", dependencies)
        self.assertIn("**Authored state:** proved. **Graph state:** stale", scopes)
        self.assertIn("Formal verification: recorded", scopes)
        self.assertIn("External review: not recorded", scopes)
        self.assertIn("SOURCE (stale)", scopes)

    def test_rendering_is_deterministic_and_does_not_edit_the_authority(self):
        original = deepcopy(self.data)
        graph = self.graph()
        for renderer in (render_dependency_report, render_scope_report):
            first = renderer(self.data, graph)
            self.assertEqual(first, renderer(self.data, graph))
            self.assertTrue(first.endswith("\n"))
        self.assertEqual(self.data, original)
        self.assertEqual(graph.data, original)

    def test_discharged_question_renders_computed_proved_state_and_authored_assessment(self):
        gap = next(node for node in self.data["graph"]["nodes"] if node["id"] == "GAP")
        gap["status"] = "proved"
        gap["sources"] = [deepcopy(self.source)]
        graph = self.graph()
        scopes = render_scope_report(self.data, graph)
        dependencies = render_dependency_report(self.data, graph)
        self.assertNotIn("Open questions", scopes)
        self.assertIn("Tracked obligations from the catalogue:", scopes)
        self.assertIn("**GAP — mathematical; proved / discharged.** Authored question:", scopes)
        self.assertIn("**Computed obligation progress:** 1 proved / discharged.", scopes)
        self.assertIn("**Authored scope assessment:** The all-degree implication remains unverified.", scopes)
        self.assertIn("| GAP | mathematical | proved / discharged |", dependencies)


if __name__ == "__main__":
    unittest.main()

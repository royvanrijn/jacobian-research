"""Keep bibliographic identity and proof assurance separate."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_graph import DependencyGraph
from dependency_references import bibliography_items, check_references, reference_state


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        self.source = self.root / "proof.tex"
        self.source.write_text(r"\bibitem{Input}A. Author, Exact input." + "\n")
        self.graph = DependencyGraph({"topics": [], "graph": {"schema_version": 1, "nodes": [
            {"id": "A", "title": "Claim", "scope": "The exact proposition.", "kind": "claim",
             "status": "open"}], "implications": []}}, self.root)
        self.data = {"graph": {"external_references": {"schema_version": 1, "papers": [{
            "id": "PAPER", "title": "Exact input", "citation": "A. Author, Exact input.",
            "verification": "citation_only", "local_sources": [], "uses": [{
                "path": "proof.tex", "locator": "bibitem Input", "role": "application_unverified",
                "nodes": ["A"]}],
        }], "internal_citations": []}}}

    def test_citation_mapping_does_not_prove_its_application(self):
        before = deepcopy(self.data)
        state = check_references(self.data, self.graph, ["proof.tex"])
        self.assertTrue(state["complete"])
        self.assertEqual(self.graph.evaluate()["A"]["state"], "open")
        self.assertEqual(self.data, before)

    def test_new_bibliography_item_requires_a_resource_mapping(self):
        self.source.write_text(self.source.read_text() + r"\bibitem{New}A new paper." + "\n")
        with self.assertRaisesRegex(AssertionError, "unmapped bibliography item.*New"):
            check_references(self.data, self.graph, ["proof.tex"])

    def test_bibtex_entries_are_covered_and_macros_are_not_papers(self):
        path = self.root / "references.bib"
        path.write_text('@article{Alpha, title={An input}}\n'
                        '@book(Beta, title={A book})\n'
                        '@string{journal = "Name"}\n'
                        '@comment{ignored, anything}\n')
        self.assertEqual(bibliography_items(self.root, ["references.bib"]),
                         {("references.bib", "Alpha"), ("references.bib", "Beta")})
        with self.assertRaisesRegex(AssertionError, "references.bib :: Alpha"):
            check_references(self.data, self.graph, ["proof.tex", "references.bib"])

    def test_changed_retained_source_is_reported(self):
        cached = self.root / "input.txt"
        cached.write_text("The precise primary source.\n")
        paper = self.data["graph"]["external_references"]["papers"][0]
        paper["local_sources"] = [{"path": "input.txt", "locator": "Theorem 1",
                                  "sha256": hashlib.sha256(cached.read_bytes()).hexdigest()}]
        check_references(self.data, self.graph, ["proof.tex"])
        cached.write_text("Changed input statement.\n")
        state = reference_state(self.data, self.graph, ["proof.tex"])
        self.assertIn("changed external source: input.txt", state["errors"])

    def test_new_markdown_url_requires_identity_and_preserves_doi_parentheses(self):
        note = self.root / "note.md"
        url = "https://doi.org/10.1016/S0019-3577(98)80020-7"
        note.write_text("[Precise input](" + url + ").\n")
        self.graph.data["graph"]["coverage"] = {"documents": [{"path": "note.md", "role": "proof"}]}
        state = reference_state(self.data, self.graph, ["proof.tex", "note.md"])
        self.assertIn("unmapped external URL: note.md :: " + url, state["errors"])
        self.data["graph"]["external_references"]["papers"][0]["urls"] = [url]
        check_references(self.data, self.graph, ["proof.tex", "note.md"])

    def test_external_primary_input_needs_a_cited_resource_mapping(self):
        primary = self.root / "sources/input.pdf"
        primary.parent.mkdir()
        primary.write_bytes(b"Retained primary input bytes")
        self.graph.nodes["A"]["origin"] = "external"
        self.graph.nodes["A"]["sources"] = [{
            "path": "sources/input.pdf", "locator": "Theorem 1",
            "sha256": hashlib.sha256(primary.read_bytes()).hexdigest()}]
        state = reference_state(self.data, self.graph, ["proof.tex"])
        self.assertTrue(any("lacks cited-resource mapping" in error for error in state["errors"]))
        self.data["graph"]["external_references"]["papers"][0]["local_sources"] = self.graph.nodes["A"]["sources"]
        check_references(self.data, self.graph, ["proof.tex"])

    def test_inline_tex_reference_is_not_hidden_by_missing_bibitem(self):
        url = "https://dlmf.nist.gov/1.10.iv"
        self.source.write_text(self.source.read_text() + r"\href{" + url + "}{Rouche theorem}")
        self.graph.data["graph"]["coverage"] = {"documents": [{"path": "proof.tex", "role": "proof"}]}
        state = reference_state(self.data, self.graph, ["proof.tex"])
        self.assertIn("unmapped external URL: proof.tex :: " + url, state["errors"])


if __name__ == "__main__":
    unittest.main()

"""Detect omitted sources/declarations without confusing coverage with proof."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_coverage import (CoverageError, check_coverage, coverage_state,
                                 declarations, inventory_paths)


class CoverageTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        self.path = self.root / "proof.tex"
        self.path.write_text(r"\begin{theorem}\label{thm:one}A.\end{theorem}" + "\n")
        self.graph = SimpleNamespace(root=self.root, repository=self.repository,
                                     nodes={"A": {"status": "open"}})
        self.data = {"graph": {"coverage": {"schema_version": 1, "documents": [{
            "path": "proof.tex", "sha256": hashlib.sha256(self.path.read_bytes()).hexdigest(),
            "role": "proof", "claims": ["A"], "dependencies": [], "gaps": [],
            "explanation": "Exact retained statement, awaiting its proof input.",
            "statements": [{"locator": "thm:one", "nodes": ["A"]}],
        }]}}}

    def test_mapped_open_claim_is_inventory_complete_without_becoming_proved(self):
        before = deepcopy(self.data)
        state = check_coverage(self.data, self.graph, ["proof.tex"])
        self.assertTrue(state["complete"])
        self.assertEqual(state["declarations"], 1)
        self.assertEqual(self.graph.nodes["A"]["status"], "open")
        self.assertEqual(self.data, before)

    def test_new_retained_document_is_not_silently_omitted(self):
        with self.assertRaisesRegex(CoverageError, "unmapped retained document: new.md"):
            check_coverage(self.data, self.graph, ["proof.tex", "new.md"])

    def test_new_statement_and_changed_bytes_require_reconciliation(self):
        self.path.write_text(self.path.read_text() + r"\begin{lemma}\label{lem:new}B.\end{lemma}" + "\n")
        state = coverage_state(self.data, self.graph, ["proof.tex"])
        self.assertTrue(any("changed coverage source" in error for error in state["errors"]))
        self.assertTrue(any("lem:new" in error for error in state["errors"]))
        # Accepting a byte hash alone cannot hide the unregistered statement.
        self.data["graph"]["coverage"]["documents"][0]["sha256"] = hashlib.sha256(self.path.read_bytes()).hexdigest()
        with self.assertRaisesRegex(CoverageError, "lem:new"):
            check_coverage(self.data, self.graph, ["proof.tex"])

    def test_unknown_claim_and_reference_without_reason_are_rejected(self):
        record = self.data["graph"]["coverage"]["documents"][0]
        record["claims"] = ["missing"]
        record["explanation"] = ""
        state = coverage_state(self.data, self.graph, ["proof.tex"])
        self.assertTrue(any("unknown coverage node missing" in error for error in state["errors"]))
        self.assertTrue(any("lacks explanation" in error for error in state["errors"]))

    def test_formal_axiom_and_theorem_are_indexed_with_actual_names(self):
        path = self.root / "Test.lean"
        path.write_text("/- theorem ignored : True := by trivial -/\n"
                        "axiom imported_input : True\n"
                        "theorem exported : True := imported_input\n")
        self.assertEqual([(row["locator"], row["kind"]) for row in declarations(path)],
                         [("imported_input", "axiom"), ("exported", "theorem")])

    def test_manifest_lists_ignored_artifacts_and_requires_the_expected_symlink_and_digest(self):
        subprocess.run(["git", "init", "-q", str(self.repository)], check=True)
        (self.repository / ".gitignore").write_text("/research/bulk.xz\n")
        external = tempfile.TemporaryDirectory()
        self.addCleanup(external.cleanup)
        external_root = Path(external.name)
        target = external_root / "research/bulk.xz"
        target.parent.mkdir()
        target.write_bytes(b"original evidence")
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        (self.repository / "external-artifacts.json").write_text(json.dumps({
            "schema_version": 1, "external_root": str(external_root),
            "files": [{"path": "research/bulk.xz", "sha256": digest}],
        }))
        path = self.root / "bulk.xz"
        try:
            path.symlink_to(target)
        except (NotImplementedError, OSError):
            self.skipTest("symbolic links are unavailable")
        record = deepcopy(self.data["graph"]["coverage"]["documents"][0])
        record.update(path="bulk.xz", sha256=digest, role="evidence", statements=[])
        self.data["graph"]["coverage"]["documents"].append(record)

        self.assertIn("bulk.xz", inventory_paths(self.root, self.repository))
        self.assertTrue(check_coverage(self.data, self.graph)["complete"])
        target.write_bytes(b"changed evidence")
        with self.assertRaisesRegex(CoverageError, "changed external artifact bytes"):
            check_coverage(self.data, self.graph)

        wrong = target.with_name("wrong.xz")
        wrong.write_bytes(b"original evidence")
        path.unlink()
        path.symlink_to(wrong)
        with self.assertRaisesRegex(CoverageError, "coverage path escapes repository"):
            check_coverage(self.data, self.graph)

    def test_inventory_excludes_generated_and_archived_editions_but_keeps_current_sources(self):
        subprocess.run(["git", "init", "-q", str(self.repository)], check=True)
        current = self.repository / "results/paper"
        (current / "published").mkdir(parents=True)
        (current / "main.tex").write_text("Current theorem.\n")
        (current / "published/main.tex").write_text("Archived theorem.\n")
        (self.root / "README.md").write_text("Generated index.\n")
        (self.root / "MATH_STATUS.json").write_text("{}\n")
        (self.root / "moments").mkdir()
        (self.root / "moments/README.md").write_text("Generated guide.\n")
        paths = inventory_paths(self.root, self.repository)
        self.assertIn("../results/paper/main.tex", paths)
        self.assertNotIn("../results/paper/published/main.tex", paths)
        self.assertNotIn("README.md", paths)
        self.assertNotIn("MATH_STATUS.json", paths)
        self.assertNotIn("moments/README.md", paths)


if __name__ == "__main__":
    unittest.main()

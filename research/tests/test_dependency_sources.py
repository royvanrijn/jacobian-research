"""Verify local proof input/import coverage and read-only closure inspection."""

from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_sources import missing_source_closure, source_closure


class DependencySourceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()

    def write(self, path, content):
        target = self.root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")
        return target

    def closure(self, *paths):
        return source_closure(list(paths), self.root, self.repository)

    def test_tex_entrypoint_and_current_directory_nested_comments_and_cycles(self):
        self.write("paper/main.tex", "\\input{sections/proof}\n% \\include{missing}\n")
        self.write("paper/sections/proof.tex", "\\include{helper}\n\\input{../main}\n")
        self.write("paper/helper.tex", "Entry-root helper\n")
        self.write("paper/sections/helper.tex", "Including-directory helper\n")
        paths = {item["path"] for item in self.closure("paper/main.tex")}
        self.assertEqual(paths, {"paper/main.tex", "paper/sections/proof.tex",
                                 "paper/helper.tex", "paper/sections/helper.tex"})

    def test_tex_escaped_percent_literal_input_external_and_macro_paths(self):
        self.write("main.tex", "\\% text \\input leaf\n\\input{\\dynamic}\n\\input{/etc/passwd}\n\\input{../../outside}\n")
        self.write("leaf.tex", "Leaf\n")
        self.assertEqual({item["path"] for item in self.closure("main.tex")},
                         {"main.tex", "leaf.tex"})
        with self.assertRaisesRegex(ValueError, "escapes repository"):
            self.closure("../../outside")

    def test_lean_local_imports_custom_srcdir_configuration_nested_comments(self):
        self.write("formal/lakefile.toml", 'name = "sample"\nsrcDir = "src"\n')
        self.write("formal/lean-toolchain", "leanprover/lean4:v4.24.0\n")
        self.write("formal/lake-manifest.json", '{"packages":[]}\n')
        self.write("formal/src/Entry.lean",
                   "import Local.A Mathlib.Algebra\n/- import Missing\n/- nested -/\n-/\n"
                   'def text := "import Missing"\npublic import Local.B -- comment\n')
        self.write("formal/src/Local/A.lean", "import Local.B\n")
        self.write("formal/src/Local/B.lean", "import Local.A\n")
        rows = self.closure("formal/src/Entry.lean")
        self.assertEqual({item["path"] for item in rows},
                         {"formal/src/Entry.lean", "formal/src/Local/A.lean",
                          "formal/src/Local/B.lean", "formal/lakefile.toml",
                          "formal/lean-toolchain", "formal/lake-manifest.json"})
        self.assertTrue(any("Local Lean import Local.A" in item["locator"] for item in rows))

    def test_nearest_lean_project_wins_and_no_external_package_walk(self):
        self.write("lakefile.toml", 'name = "outer"\n')
        self.write("Local/A.lean", "Outer module\n")
        self.write("nested/lakefile.toml", 'name = "inner"\n')
        self.write("nested/Entry.lean", "import Local.A\n")
        self.assertEqual({item["path"] for item in self.closure("nested/Entry.lean")},
                         {"nested/Entry.lean", "nested/lakefile.toml"})

    def test_established_leaf_requires_pin_closure_and_does_not_mutate(self):
        main = self.write("main.tex", "\\input{proof}\n")
        proof = self.write("proof.tex", "Proof\n")
        source = {"path": "main.tex", "sha256": hashlib.sha256(main.read_bytes()).hexdigest(),
                  "locator": "Theorem 1"}
        nodes = {
            "LEAF": {"kind": "source", "status": "proved", "sources": [source]},
            "GOAL": {"kind": "goal", "status": "proved", "sources": [source]},
            "OPEN": {"kind": "source", "status": "unverified", "sources": [source]},
            "DERIVED": {"kind": "source", "status": "proved", "sources": [source]},
        }
        graph = SimpleNamespace(nodes=nodes, incoming={"DERIVED": ["RULE"]},
                                root=self.root, repository=self.repository)
        baseline = deepcopy(nodes)
        self.assertEqual(missing_source_closure(graph), [{"owner": "LEAF", "path": "proof.tex"}])
        self.assertEqual(nodes, baseline)
        nodes["LEAF"]["sources"].append(
            {"path": "proof.tex", "sha256": hashlib.sha256(proof.read_bytes()).hexdigest(),
             "locator": "The proof"})
        self.assertEqual(missing_source_closure(graph), [])

    def test_changed_input_bytes_do_not_get_rewritten_or_rehashed(self):
        main = self.write("main.tex", "\\input{proof}\n")
        proof = self.write("proof.tex", "Old argument\n")
        before = main.read_bytes()
        proof.write_text("Changed argument\n", encoding="utf-8")
        first = self.closure("main.tex")
        self.assertEqual(first, self.closure("main.tex"))
        self.assertEqual(main.read_bytes(), before)
        self.assertEqual(proof.read_text(), "Changed argument\n")


if __name__ == "__main__":
    unittest.main()


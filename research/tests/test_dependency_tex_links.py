"""Explicit TeX proof links preserve paper scope and source provenance."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_tex_links import tex_proof_link_report, tex_proof_links


class TexProofLinksTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.root.mkdir()
        self.nodes = {}
        self.data = {"graph": {"coverage": {"documents": []}}}

    def write(self, path, text):
        target = self.root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
        self.data["graph"]["coverage"]["documents"].append({"path": path})

    def claim(self, node_id, path, label, locator=None):
        self.nodes[node_id] = {"id": node_id, "kind": "claim",
                               "sources": [{"path": path, "locator": locator or label}]}

    def graph(self):
        return SimpleNamespace(root=self.root, repository=self.repository,
                               nodes=self.nodes, rules={})

    def test_cross_input_proof_refs_have_original_file_and_line(self):
        self.write("paper/main.tex", r"""\documentclass{article}
\input{lemma}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{proof}
Use Lemma~\ref{lem:input}.
\end{proof}
""")
        self.write("paper/lemma.tex", r"""\begin{lemma}\label{lem:input}Input.\end{lemma}
\begin{proof}Direct.\end{proof}
""")
        self.claim("MOM-MAIN", "paper/main.tex", "thm:main")
        self.claim("MOM-INPUT", "paper/lemma.tex", "lem:input")
        baseline = deepcopy(self.data), deepcopy(self.nodes)
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["issues"], [])
        self.assertEqual(len(report["links"]), 1)
        link = report["links"][0]
        self.assertEqual((link["conclusion"], link["dependency"]), ("MOM-MAIN", "MOM-INPUT"))
        self.assertEqual(link["dependency_path"], "paper/lemma.tex")
        self.assertEqual(link["source_path"], "paper/main.tex")
        self.assertIn("line 5", link["locator"])
        self.assertEqual((self.data, self.nodes), baseline)

    def test_forward_reference_and_included_proof_body(self):
        self.write("p/main.tex", r"""\documentclass{article}
\begin{theorem}\label{thm:early}Early.\end{theorem}
\begin{proof}\input{proof}\end{proof}
\input{late}
""")
        self.write("p/proof.tex", r"Use \ref{lem:later}.")
        self.write("p/late.tex", r"\begin{lemma}\label{lem:later}Later.\end{lemma}")
        self.claim("MOM-EARLY", "p/main.tex", "thm:early")
        self.claim("MOM-LATER", "p/late.tex", "lem:later")
        links = tex_proof_links(self.data, self.graph())
        self.assertEqual(len(links), 1)
        self.assertEqual(links[0]["source_path"], "p/proof.tex")
        self.assertEqual(links[0]["dependency"], "MOM-LATER")

    def test_identical_labels_are_separate_in_different_papers(self):
        for paper in ("a", "b"):
            self.write(paper + "/main.tex", r"""\documentclass{article}
\input{lemma}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{proof}Use\ref{lem:input}.\end{proof}
""")
            self.write(paper + "/lemma.tex",
                       r"\begin{lemma}\label{lem:input}Input.\end{lemma}")
            self.claim("MOM-" + paper + "-MAIN", paper + "/main.tex", "thm:main")
            self.claim("MOM-" + paper + "-INPUT", paper + "/lemma.tex", "lem:input")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["issues"], [])
        self.assertEqual({(x["conclusion"], x["dependency"]) for x in report["links"]},
                         {("MOM-a-MAIN", "MOM-a-INPUT"), ("MOM-b-MAIN", "MOM-b-INPUT")})

    def test_equation_and_section_labels_are_not_claim_dependencies(self):
        self.write("main.tex", r"""\documentclass{article}
\section{Start}\label{sec:start}
\begin{lemma}\label{lem:input}Input.\end{lemma}
\begin{equation}a=b\label{eq:identity}\end{equation}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{proof}Use\ref{eq:identity},\eqref{eq:missing},\ref{sec:start},\ref{lem:input}.\end{proof}
""")
        self.claim("MOM-MAIN", "main.tex", "thm:main")
        self.claim("MOM-INPUT", "main.tex", "lem:input")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["issues"], [])
        self.assertEqual([x["reference"] for x in report["links"]], ["lem:input"])

    def test_unknown_and_ambiguous_references_are_reported_without_guesses(self):
        self.write("main.tex", r"""\documentclass{article}
\begin{lemma}\label{lem:duplicate}One.\end{lemma}
\begin{lemma}\label{lem:duplicate}Two.\end{lemma}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{proof}Use\ref{lem:unknown},\ref{lem:duplicate},\ref{lem:unmapped}.\end{proof}
\begin{lemma}\label{lem:unmapped}Not inventoried.\end{lemma}
""")
        self.claim("MOM-MAIN", "main.tex", "thm:main")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["links"], [])
        self.assertEqual({x["kind"] for x in report["issues"]},
                         {"unresolved_reference", "ambiguous_label", "unmapped_theorem"})

    def test_precise_fine_nodes_outrank_aggregate_mentions_and_ties_are_explicit(self):
        self.write("main.tex", r"""\documentclass{article}
\begin{lemma}\label{lem:input}Input.\end{lemma}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{proof}Use\ref{lem:input}.\end{proof}
""")
        self.claim("MOM-MAIN", "main.tex", "thm:main")
        self.claim("BUNDLE", "main.tex", "lem:input",
                   "Principal results thm:main and lem:input")
        self.claim("MOM-INPUT", "main.tex", "lem:input", "lem:input;line2")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["links"][0]["dependency"], "MOM-INPUT")
        self.claim("MOM-OTHER-INPUT", "main.tex", "lem:input", "lem:input;line2")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["links"], [])
        self.assertEqual(report["issues"][0]["kind"], "ambiguous_node")

    def test_explicit_proof_title_resolves_owner_after_other_theorems(self):
        self.write("main.tex", r"""\documentclass{article}
\begin{theorem}\label{thm:main}Target.\end{theorem}
\begin{lemma}\label{lem:input}Input.\end{lemma}
\begin{proof}Direct.\end{proof}
\begin{proof}[Proof of Theorem~\ref{thm:main}]Use\ref{lem:input}.\end{proof}
""")
        self.claim("MOM-MAIN", "main.tex", "thm:main")
        self.claim("MOM-INPUT", "main.tex", "lem:input")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["issues"], [])
        self.assertEqual([(x["conclusion"], x["dependency"]) for x in report["links"]],
                         [("MOM-MAIN", "MOM-INPUT")])

    def test_missing_or_ambiguous_inputs_and_comments_are_visible(self):
        self.write("main.tex", r"""\documentclass{article}
% \input{ignored}
\input{missing}
\begin{verbatim}\input{ignored_again}\end{verbatim}
""")
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual([x["kind"] for x in report["issues"]], ["missing_input"])
        self.assertEqual(report["issues"][0]["reference"], "missing")

    def test_ambiguous_local_input_is_reported_instead_of_choosing_a_file(self):
        self.write("paper/main.tex", r"\documentclass{article}\input{sub/part}")
        self.write("paper/sub/part.tex", r"\input{shared}")
        self.write("paper/shared.tex", "Entry-root interpretation.")
        self.write("paper/sub/shared.tex", "Including-directory interpretation.")
        report = tex_proof_link_report(self.data, self.graph())
        ambiguous = [item for item in report["issues"] if item["kind"] == "ambiguous_input"]
        self.assertEqual(len(ambiguous), 1)
        self.assertEqual(ambiguous[0]["candidates"],
                         ["paper/shared.tex", "paper/sub/shared.tex"])

    def test_deposited_snapshots_and_primary_mirrors_remain_leaves(self):
        text = (r"\documentclass{article}\begin{theorem}\label{thm:old}X.\end{theorem}"
                r"\begin{proof}Use\ref{lem:missing}.\end{proof}")
        self.write("paper/published/main.tex", text)
        self.write("sources/external.tex", text)
        self.assertEqual(tex_proof_link_report(self.data, self.graph()),
                         {"links": [], "issues": []})

    def test_late_named_global_proof_does_not_become_an_input_proposition_proof(self):
        self.write("main.tex", r"""\documentclass{article}
\begin{theorem}\label{thm:global}Global.\end{theorem}
\begin{lemma}\label{lem:case}Case.\end{lemma}
\begin{proof}Direct.\end{proof}
\begin{proposition}\label{prop:local}Local.\end{proposition}
\input{local_proof}
\begin{proof}[Proof of Theorem~\ref{thm:global}]
Use\ref{prop:local}.
\end{proof}
""")
        self.write("local_proof.tex",
                   r"\begin{proof}[Proof of Proposition~\ref{prop:local}]Use\ref{lem:case}.\end{proof}")
        for id, path, label in [
                ("MOM-GLOBAL", "main.tex", "thm:global"),
                ("MOM-CASE", "main.tex", "lem:case"),
                ("MOM-LOCAL", "main.tex", "prop:local")]:
            self.claim(id, path, label)
        report = tex_proof_link_report(self.data, self.graph())
        self.assertEqual(report["issues"], [])
        self.assertEqual({(item["conclusion"], item["dependency"])
                          for item in report["links"]},
                         {("MOM-GLOBAL", "MOM-LOCAL"), ("MOM-LOCAL", "MOM-CASE")})


if __name__ == "__main__":
    unittest.main()

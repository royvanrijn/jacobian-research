"""Regression cases discovered while auditing catalogue completeness controls."""

from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_bounds import bound_comparisons, bounds_records
from dependency_coverage import declarations


def graph_for_bounds(bounds):
    nodes = {str(index): {"id": str(index), "title": str(index), "bounds": [bound]}
             for index, bound in enumerate(bounds)}
    states = {node_id: {"state": "proved", "assumptions": [],
                        "dependency_coverage": "partial"} for node_id in nodes}
    return SimpleNamespace(nodes=nodes, evaluate=lambda: states, impact=lambda node_id: [])


class ToolingReviewTests(unittest.TestCase):
    def test_blockquoted_and_admonition_theorems_are_not_omitted(self):
        forms = [
            "> **Theorem (quoted).** Every admissible x satisfies P.\n",
            "> > **Theorem (nested).** Every admissible x satisfies P.\n",
            "> ### Theorem (heading)\n> Every admissible x satisfies P.\n",
            "> [!IMPORTANT]\n> **Theorem (admonition).** Every admissible x satisfies P.\n",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "proof.md"
            for form in forms:
                with self.subTest(form=form):
                    path.write_text(form)
                    found = declarations(path)
                    self.assertEqual(len(found), 1)
                    self.assertEqual(found[0]["kind"], "theorem")

    def test_expression_does_not_make_an_invalid_numeric_value_comparable(self):
        for value in (True, False, float("nan"), float("inf"), -float("inf"), "12", None):
            with self.subTest(value=value):
                graph = graph_for_bounds([{
                    "quantity": "dimension", "domain": "the same fixed class of maps",
                    "sense": "upper", "value": value, "expression": "a symbolic bound",
                }])
                with self.assertRaises(AssertionError):
                    bounds_records(graph)

    def test_units_and_unspecified_units_form_separate_comparison_domains(self):
        common = {
            "quantity": "length", "domain": "the same mathematical object",
            "parameters": {"family": "fixed"}, "sense": "upper",
        }
        graph = graph_for_bounds([
            {**common, "value": 100, "units": "centimetres"},
            {**common, "value": 2, "units": "metres"},
            {**common, "value": 1},
        ])
        groups = bound_comparisons(graph)
        self.assertEqual(len(groups), 3)
        self.assertTrue(all(len(group["claims"]) == 1 for group in groups))
        self.assertTrue(all(not group["opportunities"] for group in groups))


if __name__ == "__main__":
    unittest.main()

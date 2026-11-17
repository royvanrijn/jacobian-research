"""Compare bounds only when their mathematical scopes agree."""

from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_bounds import bound_comparisons


class BoundsTests(unittest.TestCase):
    def graph(self, rows):
        nodes, evaluation = {}, {}
        for node_id, domain, sense, value, state, assumptions in rows:
            nodes[node_id] = {"title": node_id, "bounds": [{
                "quantity": "collision dimension", "domain": domain, "sense": sense, "value": value}]}
            evaluation[node_id] = {"state": state, "assumptions": assumptions, "dependency_coverage": "partial"}
        return SimpleNamespace(nodes=nodes, evaluate=lambda: evaluation, impact=lambda node_id: [node_id + "-consumer"])

    def test_stronger_bound_reports_consumers_without_changing_status(self):
        graph = self.graph([("old", "cubic-linear", "upper", 90, "proved", []),
                            ("new", "cubic-linear", "upper", 89, "proved", [])])
        compared = bound_comparisons(graph)
        self.assertEqual(compared[0]["best_upper"], 89)
        self.assertEqual(compared[0]["opportunities"][0]["affected"], ["old-consumer"])

    def test_general_and_cubic_classes_are_not_compared(self):
        graph = self.graph([("general", "ordinary polynomial maps", "upper", 3, "proved", []),
                            ("cubic", "cubic-linear maps", "upper", 89, "proved", [])])
        self.assertEqual(len(bound_comparisons(graph)), 2)
        self.assertTrue(all(not row["opportunities"] for row in bound_comparisons(graph)))

    def test_conditional_improvement_cannot_replace_unconditional_bound(self):
        graph = self.graph([("known", "same", "upper", 12, "proved", []),
                            ("hypothesis", "same", "upper", 10, "conditional", ["H"])])
        compared = bound_comparisons(graph)
        self.assertEqual(len(compared), 2)
        self.assertTrue(all(not row["opportunities"] for row in compared))

    def test_inconsistent_bounds_are_exposed(self):
        graph = self.graph([("upper", "same", "upper", 3, "proved", []),
                            ("lower", "same", "lower", 4, "proved", [])])
        self.assertTrue(bound_comparisons(graph)[0]["conflict"])


if __name__ == "__main__":
    unittest.main()

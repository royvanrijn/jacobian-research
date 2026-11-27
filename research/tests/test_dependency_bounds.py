"""Compare bounds only when their mathematical scopes agree."""

from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dependency_bounds import bound_comparisons, bounds_records


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

    def test_proved_exclusion_retires_comparison_while_retaining_bound(self):
        graph = self.graph([("old", "excluded class", "lower", 4, "proved", []),
                            ("new", "excluded class", "lower", 6, "proved", []),
                            ("other", "separate class", "lower", 2, "proved", [])])
        graph.nodes["exclusion"] = {"title": "Exclusion theorem"}
        for node_id in ("old", "new"):
            graph.nodes[node_id]["bounds"][0]["domain_excluded_by"] = "exclusion"
        for state in ("proved", "conditional", "stale"):
            with self.subTest(state=state):
                graph.evaluate()["exclusion"] = {
                    "state": state, "assumptions": [], "dependency_coverage": "partial"}
                self.assertEqual(len(bounds_records(graph)), 3)
                self.assertEqual(len(bound_comparisons(graph)), 1 if state == "proved" else 2)

    def test_recorded_replacement_keeps_the_older_bound_as_a_fallback(self):
        graph = self.graph([("old", "same", "lower", 4, "proved", []),
                            ("replacement", "same", "lower", 6, "proved", [])])
        graph.nodes["old"]["bounds"][0]["superseded_by"] = "replacement"

        comparison = bound_comparisons(graph)[0]
        self.assertEqual(comparison["opportunities"], [])
        self.assertEqual(comparison["recorded_updates"], [
            {"weaker_claim": "old", "sense": "lower", "old_value": 4,
             "stronger_value": 6, "replacement": "replacement"}])
        self.assertEqual(len(bounds_records(graph)), 2)

        graph.evaluate()["replacement"]["state"] = "stale"
        comparison = bound_comparisons(graph)[0]
        self.assertEqual(comparison["best_lower"], 4)
        self.assertEqual(comparison["recorded_updates"], [])

    def test_recorded_replacement_does_not_hide_a_further_improvement(self):
        graph = self.graph([("old", "same", "upper", 10, "conditional", ["H"]),
                            ("replacement", "same", "upper", 9, "conditional", ["H"]),
                            ("new", "same", "upper", 8, "conditional", ["H"])])
        graph.nodes["old"]["bounds"][0]["superseded_by"] = "replacement"

        comparison = bound_comparisons(graph)[0]
        self.assertEqual(comparison["recorded_updates"], [])
        self.assertEqual([row["weaker_claim"] for row in comparison["opportunities"]], ["old", "replacement"])

        graph.evaluate()["new"]["state"] = "stale"
        self.assertEqual(len(bound_comparisons(graph)[0]["recorded_updates"]), 1)
        graph.evaluate()["replacement"]["assumptions"] = ["K"]
        self.assertTrue(all(not row["recorded_updates"] for row in bound_comparisons(graph)))

    def test_recorded_replacement_must_strictly_improve_an_identical_bound(self):
        for target_domain, target_sense, target_value in (
                ("different", "lower", 6), ("same", "upper", 6), ("same", "lower", 4)):
            with self.subTest(domain=target_domain, sense=target_sense, value=target_value):
                graph = self.graph([("old", "same", "lower", 4, "proved", []),
                                    ("replacement", target_domain, target_sense, target_value, "proved", [])])
                graph.nodes["old"]["bounds"][0]["superseded_by"] = "replacement"
                with self.assertRaisesRegex(AssertionError, "strictly improve"):
                    bounds_records(graph)

    def test_recorded_replacements_propagate_through_a_chain(self):
        graph = self.graph([("old", "same", "lower", 4, "proved", []),
                            ("middle", "same", "lower", 6, "proved", []),
                            ("new", "same", "lower", 8, "proved", [])])
        graph.nodes["old"]["bounds"][0]["superseded_by"] = "middle"
        graph.nodes["middle"]["bounds"][0]["superseded_by"] = "new"

        comparison = bound_comparisons(graph)[0]
        self.assertEqual(comparison["opportunities"], [])
        self.assertEqual([row["replacement"] for row in comparison["recorded_updates"]],
                         ["new", "new"])
        self.assertEqual(len(bounds_records(graph)), 3)

        # A changed intermediate proof must expose the unpropagated improvement.
        graph.evaluate()["middle"]["state"] = "stale"
        comparison = bound_comparisons(graph)[0]
        self.assertEqual(comparison["recorded_updates"], [])
        self.assertEqual(comparison["opportunities"][0]["weaker_claim"], "old")

        # An intermediate conditional theorem cannot bridge unconditional bounds.
        graph.evaluate()["middle"].update(state="conditional", assumptions=["H"])
        comparison = next(row for row in bound_comparisons(graph) if row["state"] == "proved")
        self.assertEqual(comparison["recorded_updates"], [])
        self.assertEqual(comparison["opportunities"][0]["weaker_claim"], "old")

    def test_replacement_chain_checks_the_matching_assertion_on_multi_bound_nodes(self):
        graph = self.graph([("old", "same", "lower", 4, "proved", []),
                            ("middle", "same", "lower", 6, "proved", []),
                            ("new", "same", "lower", 8, "proved", [])])
        graph.nodes["old"]["bounds"][0]["superseded_by"] = "middle"
        graph.nodes["middle"]["bounds"].append({
            "quantity": "collision dimension", "domain": "other", "sense": "lower",
            "value": 7, "superseded_by": "new"})
        graph.nodes["new"]["bounds"].append({
            "quantity": "collision dimension", "domain": "other", "sense": "lower", "value": 9})

        comparison = next(row for row in bound_comparisons(graph) if row["domain"] == "same")
        self.assertEqual(comparison["recorded_updates"], [])
        self.assertEqual([row["weaker_claim"] for row in comparison["opportunities"]],
                         ["old", "middle"])

    def test_inconsistent_bounds_are_exposed(self):
        graph = self.graph([("upper", "same", "upper", 3, "proved", []),
                            ("lower", "same", "lower", 4, "proved", [])])
        self.assertTrue(bound_comparisons(graph)[0]["conflict"])


if __name__ == "__main__":
    unittest.main()

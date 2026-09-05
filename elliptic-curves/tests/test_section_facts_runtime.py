"""Shared finite-field facts and theorem gates for direct section invocations."""
from dataclasses import asdict
from fractions import Fraction
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
from research_runtime.finite_fields import family_traces
from research_runtime.regulator import Surface
from research_runtime.section_gate import guard_export
from research_runtime.store import FactStore, FiniteFieldFacts


class SectionFactsTests(unittest.TestCase):
    def test_family_trace_cache_and_exact_scope(self):
        from mod2_reduction_independence import finite_curve_points
        with TemporaryDirectory() as directory:
            facts = FiniteFieldFacts(FactStore(directory))
            record = family_traces([-1], [1], 7, a_degree=0, b_degree=0, facts=facts)
            expected = len(finite_curve_points(6, 1, 7))
            self.assertEqual([row["group_order"] for row in record["fibres"]], [expected]*8)
            self.assertEqual(family_traces([Fraction(-2, 2)], [1], 7, a_degree=0, b_degree=0,
                                         facts=facts, discover=False), record)
            with self.assertRaises(FileNotFoundError):
                family_traces([-1], [2], 7, a_degree=0, b_degree=0, facts=facts, discover=False)
            singular = family_traces([0], [0], 7, a_degree=0, b_degree=0, facts=facts)
            self.assertTrue(all(row["trace"] is None for row in singular["fibres"]))

    def test_rank_zero_gate_precedes_direct_solver_and_allows_reduction_proof(self):
        results = ROOT/"artifacts/generated-results"
        model = json.loads((results/"elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json").read_text())["weierstrass_model"]
        pairs = json.loads((results/"elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json").read_text())["pairs"]
        pair = next(row for row in pairs if row["pair_key"] == "alternate-orbit-19bad:alternate-orbit-083ad")
        surface = Surface(model["A_coefficients_low_to_high"], model["B_coefficients_low_to_high"],
                          pair["product_quartic_coefficients_low_to_high"])
        export = {"exact_surface": asdict(surface), "inputs": {}}
        with self.assertRaisesRegex(SystemExit, "EXCLUDED_BY_THEOREM"):
            guard_export(export, ROOT, limits={"seconds": 1})
        proof = guard_export(export, ROOT, limits={"seconds": 1}, reduction_only=True)
        self.assertEqual(proof["purpose"], "finite_field_proof")
        self.assertEqual(guard_export({}, ROOT)["status"], "UNKNOWN")


if __name__ == "__main__":
    unittest.main()

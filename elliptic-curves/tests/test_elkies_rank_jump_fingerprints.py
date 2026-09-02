#!/usr/bin/env python3
"""Focused tests for quotient-first Elkies 2026 control fingerprints."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/cas/build_elkies_2026_rank_jump_fingerprints.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank_jump_fingerprints_v1.json"
)
SPEC = importlib.util.spec_from_file_location("rank_jump_fingerprints", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RankJumpFingerprintTests(unittest.TestCase):
    def test_smith_tensor_dimensions_include_prime_torsion(self) -> None:
        record = MODULE.smith_quotient_record(
            ambient_rank=5,
            subgroup_rank=3,
            smith_factors=(1, 2, 6),
        )
        self.assertEqual(record["free_quotient_rank_lower_bound"], 2)
        self.assertEqual(record["torsion_invariant_factors"], [2, 6])
        self.assertEqual(record["specialized_generic_saturation_index_in_displayed_subgroup"], 12)
        self.assertEqual(
            record["tensor_dimensions_over_f_ell"], {"2": 4, "3": 3, "5": 2}
        )

    def test_binary_matroid_detects_circuit_and_duplicate(self) -> None:
        labels = ("a", "b", "c", "d")
        vectors = ((1, 0), (0, 1), (1, 1), (1, 0))
        circuits = MODULE._minimal_binary_circuits(labels, vectors)
        self.assertIn(["a", "d"], circuits)
        self.assertIn(["a", "b", "c"], circuits)
        histogram = MODULE._subset_rank_histogram(vectors)
        self.assertEqual(histogram["4"], {"2": 1})

    def test_pinned_control_summary(self) -> None:
        document = json.loads(ARTIFACT.read_text())
        self.assertEqual(
            [
                row["quotient_structure"]["free_quotient_rank_lower_bound"]
                for row in document["fingerprints"]
            ],
            [4, 8, 9, 10, 11],
        )
        self.assertEqual(
            [
                row["degree_visibility"][0][
                    "visible_exceptional_span_dimension_over_f2"
                ]
                for row in document["fingerprints"]
            ],
            [4, 5, 3, 2, 1],
        )
        for row in document["fingerprints"]:
            quotient = row["quotient_structure"]
            gain = quotient["free_quotient_rank_lower_bound"]
            self.assertEqual(
                quotient["tensor_dimensions_over_f_ell"],
                {"2": gain, "3": gain, "5": gain},
            )
            self.assertEqual(
                len(row["quotient_height_geometry"]["successive_minima"]), gain
            )
            self.assertTrue(
                row["quotient_height_geometry"]["enumeration"][
                    "complete_for_successive_minima"
                ]
            )


if __name__ == "__main__":
    unittest.main()

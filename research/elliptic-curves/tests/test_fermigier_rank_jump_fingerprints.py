#!/usr/bin/env python3
"""Focused checks for Fermigier quotient fingerprints and frozen replay."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
FINGERPRINT_SCRIPT = (
    ROOT / "elliptic-curves/cas/build_fermigier_rank_jump_fingerprints.py"
)
FINGERPRINT_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fermigier_rank_jump_fingerprints_v1.json"
)
REPLAY_SCRIPT = ROOT / "elliptic-curves/cas/build_fermigier_rank_jump_replay.py"
REPLAY_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fermigier_rank_jump_replay_v1.json"
)


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


FINGERPRINTS = load_module("fermigier_jump_fingerprints", FINGERPRINT_SCRIPT)
REPLAY = load_module("fermigier_jump_replay", REPLAY_SCRIPT)


class FermigierRankJumpTests(unittest.TestCase):
    def test_complete_graph_direction_matroid(self) -> None:
        vectors = []
        for index in range(3):
            row = [0, 0, 0]
            row[index] = 1
            vectors.extend((row, row))
        for left in range(3):
            for right in range(left + 1, 3):
                row = [0, 0, 0]
                row[left] = row[right] = 1
                vectors.extend((row, row, row, row))
        record = FINGERPRINTS.cycle_matroid_intersections(
            dimension=3, coefficient_vectors=vectors
        )
        self.assertEqual(record["matroid_identification"], "cycle matroid of complete graph K_4")
        self.assertEqual(record["distinct_unoriented_mod2_classes"], 6)
        self.assertEqual(
            record["minimal_dependency_circuit_census_by_size"]["3"][
                "distinct_class_circuits"
            ],
            4,
        )

    def test_pinned_quotient_fingerprints(self) -> None:
        document = json.loads(FINGERPRINT_ARTIFACT.read_text())
        e22, rank20 = document["fingerprints"]
        self.assertEqual(
            e22["quotient_structure"]["free_quotient_rank_lower_bound"], 10
        )
        self.assertEqual(
            e22["quotient_structure"]["tensor_dimensions_over_f_ell"],
            {"2": 22, "3": 11, "5": 10},
        )
        self.assertEqual(
            e22["quotient_structure"][
                "specialized_generic_saturation_index_in_displayed_subgroup"
            ],
            24576,
        )
        self.assertEqual(
            rank20["quotient_structure"]["free_quotient_rank_lower_bound"], 8
        )
        self.assertEqual(
            rank20["quotient_structure"]["tensor_dimensions_over_f_ell"],
            {"2": 8, "3": 8, "5": 8},
        )
        self.assertEqual(
            [
                row["degree_visibility"][0][
                    "visible_free_quotient_span_dimension"
                ]
                for row in document["fingerprints"]
            ],
            [10, 8],
        )

    def test_replay_is_censored_and_complete(self) -> None:
        document = json.loads(REPLAY_ARTIFACT.read_text())
        self.assertEqual(document["population"]["primitive_parameter_count"], 60815684)
        self.assertEqual(document["population"]["negative_label_count"], 0)
        self.assertEqual(document["population"]["certified_positive_count"], 2)
        self.assertEqual(
            [
                row["scores"]["discovery_rank"]["rank_position_one_based"]
                for row in document["anchors"]
            ],
            [2755127, 3070200],
        )
        for metric in REPLAY.METRICS:
            self.assertEqual(
                document["retrieval"][metric]["gain_at_least_8"][
                    "recall_at_budget"
                ]["100000"],
                0.0,
            )


if __name__ == "__main__":
    unittest.main()

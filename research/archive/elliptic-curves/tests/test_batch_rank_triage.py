#!/usr/bin/env python3
"""Dependency-free tests for deterministic batch rank-triage helpers."""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from batch_rank_triage import (  # noqa: E402
    DEFAULT_CANDIDATE_IDS,
    candidate_lookup,
    parse_candidate_ids,
    select_escalations,
    stable_rank,
    unique_new_points,
)
from fermigier_mestre import FermigierMestreFamily  # noqa: E402


class BatchRankTriageTests(unittest.TestCase):
    def test_candidate_id_validation(self) -> None:
        identifiers = parse_candidate_ids(
            "primary-1666-9,low-conductor-644-87,low-conductor-154-103,"
            "low-conductor-847-184,low-conductor-70-223,"
            "low-conductor-1057-218"
        )
        self.assertEqual(
            identifiers,
            (
                "primary-1666-9",
                "low-conductor-644-87",
                "low-conductor-154-103",
                "low-conductor-847-184",
                "low-conductor-70-223",
                "low-conductor-1057-218",
            ),
        )
        for value in ("", "primary-1666-9,primary-1666-9", "unknown"):
            with self.assertRaises(argparse.ArgumentTypeError):
                parse_candidate_ids(value)

        candidate = candidate_lookup()["low-conductor-644-87"]
        self.assertIn("low-conductor-644-87", DEFAULT_CANDIDATE_IDS)
        self.assertEqual(candidate.parameter, Fraction(644, 87))
        self.assertEqual(
            candidate.known_log_conductor,
            "153.964400023010213083969770128430120323448601004558877669472",
        )
        self.assertEqual(candidate.known_global_root_number, -1)

        discovered = candidate_lookup()["low-conductor-154-103"]
        self.assertIn("low-conductor-154-103", DEFAULT_CANDIDATE_IDS)
        self.assertEqual(discovered.parameter, Fraction(154, 103))
        self.assertEqual(
            discovered.known_log_conductor,
            "162.234032455648408902235970522085591159150828298194308496774",
        )
        self.assertEqual(discovered.known_global_root_number, 1)
        self.assertEqual(discovered.provenance, "automatic local-condition discovery")
        self.assertIn("p=7,11,13,17,19", discovered.selection_metadata)

        frontier = candidate_lookup()
        expected = {
            "low-conductor-847-184": (
                Fraction(847, 184),
                "162.852739805513698244891454752320221986022685402317389213062",
            ),
            "low-conductor-1057-218": (
                Fraction(1057, 218),
                "169.606323764358951772312077715131200953848832572805084059891",
            ),
        }
        for identifier, (parameter, log_conductor) in expected.items():
            with self.subTest(identifier=identifier):
                self.assertIn(identifier, DEFAULT_CANDIDATE_IDS)
                self.assertEqual(frontier[identifier].parameter, parameter)
                self.assertEqual(frontier[identifier].known_log_conductor, log_conductor)
                self.assertEqual(frontier[identifier].known_global_root_number, 1)
                self.assertEqual(
                    frontier[identifier].provenance,
                    "height-first multiple-root frontier screen",
                )

    def test_new_point_filter_removes_visible_and_sign_duplicates(self) -> None:
        parameter = Fraction(1666, 9)
        visible = FermigierMestreFamily.known_quartic_points(parameter)[0]
        extra = (Fraction(-5336, 9), Fraction(89189931418, 81))
        points = (visible, (visible[0], -visible[1]), extra, (extra[0], -extra[1]))
        new_points = unique_new_points(parameter, points)
        self.assertEqual(len(new_points), 1)
        self.assertEqual(new_points[0][0], extra)
        self.assertEqual(
            new_points[0][1],
            FermigierMestreFamily.quartic_point_to_jacobian(parameter, extra),
        )

    def test_stable_rank_and_escalation_order(self) -> None:
        height_runs = (
            {"numerical_rank": 16, "subset_indices_one_based": [1, 2]},
            {"numerical_rank": 16, "subset_indices_one_based": [1, 2]},
        )
        self.assertEqual(stable_rank(height_runs), 16)
        inconsistent = (
            height_runs[0],
            {"numerical_rank": 15, "subset_indices_one_based": [1, 2]},
        )
        with self.assertRaises(AssertionError):
            stable_rank(inconsistent)

        results = [
            {
                "candidate_id": "rank-15-many",
                "stages": [
                    {
                        "stable_pool_numerical_rank": 15,
                        "numerical_rank_gain_over_seed": 3,
                        "new_x_values_beyond_visible_sections": 30,
                    }
                ],
            },
            {
                "candidate_id": "rank-16-few",
                "stages": [
                    {
                        "stable_pool_numerical_rank": 16,
                        "numerical_rank_gain_over_seed": 4,
                        "new_x_values_beyond_visible_sections": 4,
                    }
                ],
            },
            {
                "candidate_id": "ineligible",
                "stages": [
                    {
                        "stable_pool_numerical_rank": 12,
                        "numerical_rank_gain_over_seed": 0,
                        "new_x_values_beyond_visible_sections": 2,
                    }
                ],
            },
        ]
        self.assertEqual(
            select_escalations(
                results,
                minimum_rank_gain=3,
                minimum_new_x=20,
                limit=2,
            ),
            ("rank-16-few", "rank-15-many"),
        )


if __name__ == "__main__":
    unittest.main()

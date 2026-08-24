#!/usr/bin/env python3
"""Focused replay checks for the disjoint Fermigier deep tranche artifact."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves" / "cas" / "search_fermigier_record_residue_deep_tranche.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_record_residue_deep_tranche.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "4520daee70358ae147038a109f2276ff1d10705cb3b6137d4171b254e0c5ea06"
)
EXPECTED_ARTIFACT_SHA256 = (
    "fd17b69198950c1c26aca1b8b87cd37ccea894831b6952c00df08d1be251b74b"
)
EXPECTED_SELECTED_SHA256 = (
    "22452d49f569aee1bc2b0cf6a3206445271d82c704f1b997f5f19cd7c03b3c44"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierDeepTrancheTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_population_sources_and_files(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        replay = self.data["population_replay"]
        self.assertEqual(replay["candidates_checked_exactly"], 23_769)
        self.assertEqual(replay["nonsingular_candidates"], 23_769)
        self.assertEqual(
            replay["population_and_valuation_sha256"],
            "b3c6751977ca35febe9d4f1d974183226bd0d77ec34321ed9d59244cffe4f086",
        )
        self.assertEqual(self.data["execution"]["phase"], "complete")
        self.assertTrue(self.data["execution"]["no_retries"])

    def test_selection_is_disjoint_and_held_forward(self) -> None:
        candidates = self.data["candidates"]
        parameters = [Fraction(row["t"]) for row in candidates]
        self.assertEqual(len(candidates), 48)
        self.assertEqual(len(set(parameters)), 48)
        digest = hashlib.sha256(
            "".join(f"{parameter}\n" for parameter in parameters).encode()
        ).hexdigest()
        self.assertEqual(digest, EXPECTED_SELECTED_SHA256)
        self.assertEqual(
            self.data["selection"]["selected_parameter_sha256"],
            EXPECTED_SELECTED_SHA256,
        )
        prior = self.data["prior_exclusions"]
        self.assertEqual(prior["batch_rank_triage_parameters"], 11)
        self.assertEqual(prior["unique_prior_parameters_including_benchmark"], 111)
        self.assertEqual(
            prior["fermigier_accidental_parameter_sha256"],
            "e8410fbcba4491165fd114e86cff11a1eabc1373f74cebda8dbc856bbcf0045f",
        )
        self.assertEqual(self.data["scoring"]["prior_parameter_intersection_count"], 4)
        self.assertEqual(self.data["scoring"]["eligible_population_count"], 23_765)
        self.assertEqual(
            self.data["scoring"]["selection_does_not_use"],
            "the cumulative B=500 score or any point/rank result",
        )
        reasons = [row["selection_reason"] for row in candidates]
        self.assertEqual(
            sum(reason.startswith("global-low-height") for reason in reasons), 8
        )
        self.assertEqual(sum("robust-held-maximin" in reason for reason in reasons), 16)
        self.assertEqual(sum(reason.endswith(":held-one") for reason in reasons), 12)
        self.assertEqual(sum(reason.endswith(":held-two") for reason in reasons), 12)

    def test_conductor_first_gate_is_exactly_replayed(self) -> None:
        candidates = self.data["candidates"]
        status_counts = {
            status: sum(
                row["conductor_probe"]["status"] == status for row in candidates
            )
            for status in ("completed", "timeout", "error")
        }
        self.assertEqual(status_counts, {"completed": 28, "timeout": 20, "error": 0})
        subtarget = {
            row["t"]: row["conductor_probe"]
            for row in candidates
            if row["conductor_probe"].get("below_strict_log_conductor_target")
        }
        self.assertEqual(set(subtarget), {"1925/157", "3206/265"})
        self.assertEqual(subtarget["1925/157"]["root_number"], 1)
        self.assertEqual(subtarget["3206/265"]["root_number"], -1)
        self.assertEqual(
            subtarget["3206/265"]["log_conductor"],
            "168.031754726474094634576523333554069039986353437969946262615",
        )
        point_searched = {row["t"] for row in candidates if "point_search" in row}
        self.assertEqual(point_searched, set(subtarget))
        self.assertEqual(self.data["outcome"]["completed_conductors"], 28)
        self.assertEqual(self.data["outcome"]["subtarget_completed_conductors"], 2)

    def test_staged_rank_evidence_and_negative_target(self) -> None:
        by_t = {row["t"]: row for row in self.data["candidates"]}
        shallow = by_t["1925/157"]["point_search"]["stages"]
        deep = by_t["3206/265"]["point_search"]["stages"]
        self.assertEqual(
            [(row["quartic_height_bound"], row["stable_pool_numerical_rank"]) for row in shallow],
            [(50_000, 12)],
        )
        self.assertEqual(
            [
                (
                    row["quartic_height_bound"],
                    row["stable_pool_numerical_rank"],
                    row["unique_new_jacobian_images"],
                )
                for row in deep
            ],
            [(50_000, 13, 1), (250_000, 14, 3), (1_000_000, 14, 9)],
        )
        for stage in shallow + deep:
            self.assertEqual(
                {run["numerical_rank"] for run in stage["height_matrix_runs"]},
                {stage["stable_pool_numerical_rank"]},
            )
            self.assertNotIn(
                "finite_reduction_attempt",
                by_t["1925/157"]["point_search"],
            )
        self.assertEqual(self.data["outcome"]["maximum_stable_numerical_rank"], 14)
        self.assertFalse(self.data["target"]["hit"])
        self.assertEqual(self.data["outcome"]["certified_target_hits"], [])


if __name__ == "__main__":
    unittest.main()

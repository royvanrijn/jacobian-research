#!/usr/bin/env python3
"""Focused checks for the stopped direct H=1000000 pair-product control."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_fermigier_published_pair_fiber_products_h1000000.py"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_published_pair_fiber_products_h1000000.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "9aeac3f4b2640ebf46819b4349ff4304ee9bd602eca2c0c744dce2e4b4c5a9f0"
)
EXPECTED_ARTIFACT_SHA256 = (
    "706cc2933641f3fb1edfe255a027d16c3bc17c44259798320e92b2c853360b26"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierPairFiberProductH1000000ControlTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_pinned_files_and_stable_dependencies(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(self.data["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(
            self.data["source"]["published_preimage_sha256"],
            "6224da9ce4db3150a197a2cf1d9bc6c1a7d0cc6f01245b3f834945f76775ab15",
        )
        self.assertEqual(
            self.data["source"]["H50000_exact_pair_result_sha256"],
            "dea8b716c5aec56817a172afd6e894e7748aaddc482a2d29c0a3360abe55bf4b",
        )
        prior = self.data["prior_decontamination"]
        self.assertEqual(prior["terminal_prior_parameter_count"], 593)
        self.assertEqual(
            prior["terminal_prior_parameter_sha256"],
            "a4d06e4662d2e30c1a0f8873f91d8d348dae10f2abaffce88dcc0f480cfeede0",
        )

    def test_direct_box_was_stopped_by_declared_disproportionate_gate(self) -> None:
        rows = self.data["pair_searches"]
        self.assertEqual([row["pair_id"] for row in rows], [
            "p6_m1__p13_m1",
            "p6_m1__p13_p1",
        ])
        for row in rows:
            search = row["search"]["search"]
            self.assertEqual(search["height_bound"], 1_000_000)
            self.assertEqual(search["status"], "timeout")
            self.assertEqual(search["timeout_seconds"], 20.0)
            self.assertFalse(search["retried"])
        outcome = self.data["outcome"]
        self.assertEqual(outcome["declared_pair_count"], 220)
        self.assertEqual(outcome["pairs_attempted"], 2)
        self.assertEqual(outcome["pairs_completed"], 0)
        self.assertEqual(outcome["pairs_timed_out"], 2)
        self.assertTrue(outcome["stopped_as_computationally_disproportionate"])
        self.assertEqual(self.data["execution"]["owned_processes_remaining"], 0)

    def test_no_unearned_arithmetic_claim(self) -> None:
        self.assertEqual(self.data["candidates"], [])
        self.assertEqual(self.data["outcome"]["conductor_calls"], 0)
        self.assertEqual(self.data["outcome"]["rank_triage_calls"], 0)
        self.assertFalse(self.data["target"]["hit"])
        self.assertTrue(self.data["parameters"]["no_retries"])


if __name__ == "__main__":
    unittest.main()

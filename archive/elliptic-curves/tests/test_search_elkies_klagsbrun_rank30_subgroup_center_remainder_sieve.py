#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from search_elkies_klagsbrun_rank30_companion_center_sieve import (  # noqa: E402
    exact_nonoverlap_proof,
    load_companion_centers,
)
from search_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve import (  # noqa: E402
    EXPECTED_PARENT_ARTIFACT_SHA256,
    EXPECTED_PARENT_SCRIPT_SHA256,
    PARENT_ARTIFACT,
    PARENT_SCRIPT,
    exact_cross_separation,
    selected_remainder,
)
from search_elkies_klagsbrun_rank30_subgroup_center_sieve import (  # noqa: E402
    exact_prior_center_separation,
)


SCRIPT = CAS / "search_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve.py"
ARTIFACT = (
    GENERATED
    / "elliptic_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve.json"
)
EXPECTED_REMAINDER_SHA256 = (
    "7a25613b3e5bbd9f8fb16f5c2b78a40c88c8c476bd8c9e83303099367f9bfb6d"
)
EXPECTED_PRIMITIVE_COUNT = 663_014_076_508
EXPECTED_SURVIVOR_COUNT = 3_031_077
EXPECTED_SURVIVOR_SHA256 = (
    "98695a55f14368692c218932e7b3e7e61c7930562fb8a4050ad3a7f27754462b"
)


class Rank30SubgroupCenterRemainderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.first, cls.remainder, cls.manifest = selected_remainder()

    def test_predecessor_and_exact_complement_are_pinned(self) -> None:
        self.assertEqual(
            hashlib.sha256(PARENT_SCRIPT.read_bytes()).hexdigest(),
            EXPECTED_PARENT_SCRIPT_SHA256,
        )
        self.assertEqual(
            hashlib.sha256(PARENT_ARTIFACT.read_bytes()).hexdigest(),
            EXPECTED_PARENT_ARTIFACT_SHA256,
        )
        self.assertEqual(len(self.first), 128)
        self.assertEqual(len(self.remainder), 713)
        self.assertEqual(
            self.manifest["remainder_sha256"], EXPECTED_REMAINDER_SHA256
        )
        self.assertEqual(
            {center.x for center in self.first}.intersection(
                center.x for center in self.remainder
            ),
            set(),
        )

    def test_all_nonoverlap_gates_replay_exactly(self) -> None:
        internal = exact_nonoverlap_proof(
            self.remainder, denominator_min=50_001, offset_radius=16_384
        )
        self.assertTrue(internal["all_exact_nonoverlap_checks_passed"])
        old, _ = load_companion_centers()
        prior = exact_prior_center_separation(
            self.remainder,
            old,
            denominator_min=50_001,
            offset_radius=16_384,
        )
        self.assertTrue(prior["all_prior_center_boxes_disjoint"])
        cross = exact_cross_separation(
            self.remainder,
            self.first,
            denominator_min=50_001,
            offset_radius=16_384,
        )
        self.assertTrue(cross["passed"])

    def test_generated_artifact_is_complete_and_negative(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            data["center_population"]["remainder_sha256"],
            EXPECTED_REMAINDER_SHA256,
        )
        result = data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(result["completed_center_count"], 713)
        self.assertEqual(
            result["declared_primitive_candidate_count"], EXPECTED_PRIMITIVE_COUNT
        )
        self.assertEqual(
            result["processed_primitive_candidate_count"], EXPECTED_PRIMITIVE_COUNT
        )
        self.assertEqual(
            result["modular_survivor_count_after_primitivity"],
            EXPECTED_SURVIVOR_COUNT,
        )
        self.assertEqual(
            result["modular_survivor_manifest_sha256"],
            EXPECTED_SURVIVOR_SHA256,
        )
        self.assertEqual(
            result["exact_nonsquare_count_after_sieve"], EXPECTED_SURVIVOR_COUNT
        )
        self.assertEqual(result["exact_square_abscissa_count"], 0)
        self.assertFalse(result["rank30_target_hit"])


if __name__ == "__main__":
    unittest.main()

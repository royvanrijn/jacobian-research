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

from search_elkies_klagsbrun_rank30 import exact_linear_combination  # noqa: E402
from search_elkies_klagsbrun_rank30_companion_center_sieve import (  # noqa: E402
    exact_nonoverlap_proof,
)
from search_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve import (  # noqa: E402
    exact_cross_separation,
)
from search_elkies_klagsbrun_rank30_subgroup_center_sieve import (  # noqa: E402
    generate_subgroup_centers,
)
from search_elkies_klagsbrun_rank30_weight3_center_sieve import (  # noqa: E402
    generate_weight3_centers,
)


SCRIPT = CAS / "search_elkies_klagsbrun_rank30_weight3_center_sieve.py"
ARTIFACT = GENERATED / "elliptic_elkies_klagsbrun_rank30_weight3_center_sieve.json"
EXPECTED_POPULATION_SHA256 = (
    "ad08e841c09b78d1b82b502c8b5b355c5ac58ce49ddee34501459e45589b3139"
)
EXPECTED_SELECTED_SHA256 = (
    "bb914ebbfc84455bc9b4a8fd8a5530f2c01edd43164e64d38c2552fe8d99b3ef"
)
EXPECTED_PRIMITIVE_COUNT = 118_096_311_742
EXPECTED_SURVIVOR_COUNT = 534_043
EXPECTED_SURVIVOR_SHA256 = (
    "f9bf8f87822be9789cc4f202f0e47af905fb8be6f4833a1ff994bd4cb2fdb8ed"
)


class Rank30Weight3CenterSieveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.population, cls.selected, cls.manifest = generate_weight3_centers()

    def test_exact_population_decontamination_and_selection(self) -> None:
        self.assertEqual(self.manifest["raw_count"], 14_616)
        self.assertEqual(self.manifest["prior_overlap_count"], 5)
        self.assertEqual(len(self.population), 14_611)
        self.assertEqual(len(self.selected), 128)
        self.assertEqual(
            self.manifest["decontaminated_population_sha256"],
            EXPECTED_POPULATION_SHA256,
        )
        self.assertEqual(
            self.manifest["selected_sha256"], EXPECTED_SELECTED_SHA256
        )
        self.assertTrue(
            all(
                record["prior_sources"] == ["prior_32_companion"]
                for record in self.manifest["prior_overlap_records"]
            )
        )
        for center in self.selected:
            self.assertEqual(exact_linear_combination(center.relation), center.point)

    def test_selected_boxes_are_disjoint_from_all_lower_weight_boxes(self) -> None:
        internal = exact_nonoverlap_proof(
            self.selected, denominator_min=50_001, offset_radius=16_384
        )
        self.assertTrue(internal["all_exact_nonoverlap_checks_passed"])
        lower, _, _ = generate_subgroup_centers()
        cross = exact_cross_separation(
            self.selected,
            lower,
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
            data["center_population"]["decontaminated_population_sha256"],
            EXPECTED_POPULATION_SHA256,
        )
        self.assertEqual(
            data["center_population"]["selected_sha256"],
            EXPECTED_SELECTED_SHA256,
        )
        result = data["search_result"]
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["wall_cap_reached"])
        self.assertEqual(result["completed_center_count"], 128)
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

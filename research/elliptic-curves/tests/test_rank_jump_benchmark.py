#!/usr/bin/env python3
"""Tests for the rank-jump evaluator and laboratory registry."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC_ROOT / "ecsearch"))

from rank_jump_benchmark import (  # noqa: E402
    GROUP_FIELDS,
    LAB_RESULT_SCHEMA,
    SCHEMA,
    evaluate_lab_manifest,
    evaluate_manifest,
    validate_lab_manifest,
    validate_manifest,
)


FEATURE_NAMES = (
    "family_residual_s0_cumulative_B200",
    "family_residual_s5_cumulative_B200",
    "family_residual_s0_window_101_200",
    "family_residual_s5_window_101_200",
    "conductor_scaled_s0_B200",
    "root_number",
    "predicted_local_conductor",
    "quotient_escape",
)
LAB_MANIFEST = ELLIPTIC_ROOT / "data" / "rank_jump_laboratory_v1.json"


def candidate(
    identifier: str, family: str, jump: int, signal: float, twist: str
) -> dict[str, object]:
    discovery = {
        name: signal if "root_number" not in name else (-1 if signal > 0 else 1)
        for name in FEATURE_NAMES
    }
    held = {name: value * 0.9 for name, value in discovery.items()}
    return {
        "id": identifier,
        "family": family,
        "root_shape": f"shape-{family}",
        "parametrization_component": f"component-{family}",
        "quadratic_twist_class": twist,
        "exceptional_quotient_rank_lower_bound": jump,
        "discovery_features": discovery,
        "held_forward_features": held,
    }


def legacy_manifest() -> dict[str, object]:
    return {
        "schema": SCHEMA,
        "discovery_primes_disjoint_from_held_forward": True,
        "feature_names": list(FEATURE_NAMES),
        "candidates": [
            candidate("a-plus", "a", 10, 5.0, "twist-a"),
            candidate("a-control", "a", 0, -2.0, "twist-a-control"),
            candidate("b-plus", "b", 8, 4.0, "twist-b"),
            candidate("b-control", "b", 0, -1.0, "twist-b-control"),
            candidate("c-plus", "c", 6, 3.0, "twist-c"),
            candidate("c-control", "c", 0, -3.0, "twist-c-control"),
            candidate("d-plus", "d", 10, 6.0, "twist-d"),
            candidate("d-control", "d", 0, -4.0, "twist-d-control"),
        ],
    }


class LegacyRankJumpBenchmarkTests(unittest.TestCase):
    def test_manifest_requires_declared_residual_and_auxiliary_features(self) -> None:
        self.assertEqual(validate_manifest(legacy_manifest()), FEATURE_NAMES)
        incomplete = legacy_manifest()
        incomplete["feature_names"] = ["family_residual_s0_cumulative_B200"]
        with self.assertRaisesRegex(ValueError, "family_residual_s5"):
            validate_manifest(incomplete)

    def test_all_structural_protocols_exclude_the_held_out_group(self) -> None:
        manifest = legacy_manifest()
        result = evaluate_manifest(manifest, (1, 2))
        self.assertEqual(
            [item["group_field"] for item in result["protocols"]],
            list(GROUP_FIELDS),
        )
        for protocol in result["protocols"]:
            field = protocol["group_field"]
            for fold in protocol["folds"]:
                self.assertTrue(
                    fold["leakage_check"]["training_excludes_entire_held_out_group"]
                )
                held = fold["held_out_group"]
                training_ids = set(fold["training_candidate_ids"])
                for row in manifest["candidates"]:
                    if row["id"] in training_ids:
                        self.assertNotEqual(row[field], held)


class RankJumpLaboratoryTests(unittest.TestCase):
    def load_manifest(self) -> dict[str, object]:
        return json.loads(LAB_MANIFEST.read_text(encoding="utf-8"))

    def test_active_registry_and_provenance_validate(self) -> None:
        validate_lab_manifest(self.load_manifest(), ROOT)

    def test_r17_complete_population_metrics(self) -> None:
        result = evaluate_lab_manifest(self.load_manifest(), ROOT)
        self.assertEqual(result["schema"], LAB_RESULT_SCHEMA)
        self.assertEqual(result["summary"]["ranked_family_count"], 3)
        self.assertEqual(result["summary"]["certified_positive_count"], 7)
        self.assertEqual(result["summary"]["ranked_positive_count"], 18)

        r17 = result["families"][0]["ranking_runs"][0]
        self.assertEqual(r17["population_count"], 121_589_944)
        self.assertEqual(
            [row["population_rank"] for row in r17["positive_observations"]],
            [54_624, 593_936, 422_873, 55_387],
        )
        top_half_percent = r17["all_known_positives"]["budget_metrics"][3]
        self.assertEqual(top_half_percent["population_fraction"], "1/200")
        self.assertEqual(top_half_percent["known_positive_recall"], 1.0)
        self.assertEqual(
            r17["all_known_positives"]["worst_known_positive_population_fraction"][
                "exact"
            ],
            "10606/2171249",
        )
        held_out = result["families"][0]["ranking_runs"][1]
        self.assertEqual(held_out["evaluation_role"], "held_out")
        self.assertEqual(held_out["population_count"], 100_000)
        self.assertEqual(
            [row["population_rank"] for row in held_out["positive_observations"]],
            [580, 232, 784, 24_210],
        )

        fermigier = result["families"][1]["ranking_runs"]
        self.assertEqual(len(fermigier), 4)
        self.assertTrue(all(run["population_count"] == 60_815_684 for run in fermigier))
        self.assertEqual(
            [row["population_rank"] for row in fermigier[0]["positive_observations"]],
            [2_755_127, 3_070_200],
        )
        self.assertTrue(
            all(
                run["all_known_positives"]["budget_metrics"][4][
                    "known_positive_recall"
                ]
                == 0.0
                for run in fermigier
            )
        )

        nagao = result["families"][2]["ranking_runs"]
        self.assertEqual(len(nagao), 2)
        self.assertTrue(all(run["population_count"] == 18_244_819 for run in nagao))
        self.assertEqual(
            [run["positive_observations"][0]["population_rank"] for run in nagao],
            [9_041_935, 755_065],
        )

    def test_leakage_and_stale_provenance_fail_closed(self) -> None:
        manifest = self.load_manifest()
        leaked = copy.deepcopy(manifest)
        leaked["families"][0]["ranking_runs"][0]["uses_point_search_features"] = True
        with self.assertRaisesRegex(ValueError, "point-search features are forbidden"):
            validate_lab_manifest(leaked, ROOT)

        stale = copy.deepcopy(manifest)
        stale["families"][0]["ranking_runs"][0]["extractor"]["source"][
            "sha256"
        ] = "0" * 64
        with self.assertRaisesRegex(ValueError, "stale source hash"):
            validate_lab_manifest(stale, ROOT)

        ungated = copy.deepcopy(manifest)
        del ungated["families"][0]["ranking_runs"][1]["extractor"][
            "arithmetic_group_gate"
        ]
        with self.assertRaisesRegex(ValueError, "arithmetic_group_gate"):
            validate_lab_manifest(ungated, ROOT)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Focused checks for the fixed rational neighborhood of the new rank-15 fiber."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from search_mestre_02136217261290_t2_rational_neighborhood import (  # noqa: E402
    ANCHOR,
    ROOTS,
    family_coefficients,
    scope_audit,
)
from search_mestre_root_tuple_scale import (  # noqa: E402
    point_digest,
    point_on_short_curve,
)


SOURCE = CAS / "scan_mestre_02136217261290_t2_neighborhood.cpp"
SCRIPT = CAS / "search_mestre_02136217261290_t2_rational_neighborhood.py"
ARTIFACT = (
    GENERATED / "elliptic_mestre_02136217261290_t2_rational_neighborhood.json"
)
EXPECTED_SOURCE_SHA256 = (
    "b1ffda487851260a302c5d9a28bcf7e529a74f28f25ff1cccd03785144f0a76f"
)
EXPECTED_SCRIPT_SHA256 = (
    "49a1d88ced4ba5eccb18126aabaf2a929df5e2598dab6d6c900efd4df7e9ba0d"
)
EXPECTED_ARTIFACT_SHA256 = (
    "9970b2a947821d33dabc67824294c72ec67f4d6f17ddf90cd4a9adacab3d164b"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre02136217261290RationalNeighborhoodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_files_result_and_anchor_are_pinned(self) -> None:
        self.assertEqual(sha256(SOURCE), EXPECTED_SOURCE_SHA256)
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.data["result_sha256"],
            "eff8abb06aeb5fdf99d645066337d347d6ab761f2b91cec4781e05020e0cd8a5",
        )
        self.assertEqual(tuple(self.data["family"]["roots"]), ROOTS)
        self.assertEqual(
            self.data["anchor"]["log_conductor"],
            "106.931803973405473954774824140429339998614531372730705726131",
        )
        self.assertEqual(ANCHOR, 2)
        self.assertTrue(
            self.data["anchor"]["excluded_from_every_searchable_population"]
        )

    def test_exact_disjoint_raw_populations_and_manifests(self) -> None:
        scope = self.data["scope"]
        self.assertTrue(scope["raw_populations_are_disjoint"])
        self.assertEqual(scope["raw_union_primitive_parameter_count"], 1_932_123)
        self.assertEqual(
            scope["raw_union_manifest_sha256"],
            "12e7a4346242789efd823fa8798ddb95ef6d303959b3eeaeddc968b02f504d8b",
        )
        near = scope_audit("near", 20_000)
        ordinary = scope_audit("ordinary", 2_000)
        self.assertEqual(near["primitive_parameter_count"], 793_645)
        self.assertEqual(
            near["parameter_manifest_sha256"],
            "14fb3a771c23b60757ff69992ec3d749264e0993d6cfd6a68cabbf95bd31019a",
        )
        self.assertEqual(ordinary["primitive_parameter_count"], 1_138_478)
        self.assertEqual(
            ordinary["parameter_manifest_sha256"],
            "01f13f9b8d887a2ba325f0e31141f01219f4b831a02579874ba32700aae57657",
        )

    def test_leakage_boundaries_and_fixed_followup(self) -> None:
        features = self.data["exact_discriminant_feature_screen"]
        self.assertEqual(features["admissible_feature_pool_count"], 20_267)
        self.assertEqual(features["exact_singular_rejections"], 0)
        self.assertEqual(
            features["exact_feature_population_sha256"],
            "142e0f378475aa9e14d86462089c86dba32227652bf1bf4a4179b00514569ef0",
        )
        selection = self.data["conductor_selection"]
        self.assertFalse(selection["selection_uses_conductor"])
        self.assertFalse(selection["selection_uses_point_or_rank_search_data"])
        self.assertTrue(selection["discovery_survivors_closed_before_held_scores"])
        self.assertEqual(selection["selected_population"], 115)
        self.assertEqual(selection["selected_per_raw_stratum"], {"near": 64, "ordinary": 51})
        self.assertEqual(
            selection["selected_population_sha256"],
            "73fef669b05034ba47713aac67c269bc3c86c1b681906b02241142f299498f35",
        )
        conductor = self.data["conductor_first_screen"]
        self.assertTrue(conductor["all_conductor_calls_completed_before_any_point_call"])
        self.assertEqual(conductor["completed"], 87)
        self.assertEqual(conductor["timeouts"], 28)
        self.assertEqual(conductor["errors"], 0)
        self.assertEqual(conductor["subtarget"], 34)
        points = self.data["point_search_protocol"]
        self.assertEqual(
            [stage["attempted"] for stage in points["stages"]], [87, 32, 8, 2]
        )
        self.assertEqual(points["maximum_stable_numerical_rank"], 15)
        self.assertEqual(points["same_height_retries"], 0)
        self.assertEqual(points["broadening_calls_after_fixed_protocol"], 0)
        self.assertEqual(self.data["target"]["hits"], [])

    def test_two_new_exact_rank15_fibers(self) -> None:
        expected = {
            "13/10": {
                "conductor": "9577365148718819631533974079429020799341218098640810",
                "log_conductor": "119.691242260181396476394440177432235539182143832174340838467",
                "point_sha256": "5285b094b82dba105cf9197f706f4a1a71d63e54753fd707ea2aebd6e0bcf15f",
            },
            "8/7": {
                "conductor": "716655325145183252004644963198717811664700168495725110",
                "log_conductor": "124.006434749645891980207193194330140573258004346510313420694",
                "point_sha256": "f07e0b1019046ef1c19e30e801baaacf0c473e312b0b7b3a9dbf8b1d4c17207c",
            },
        }
        records = {
            record["parameter"]: record
            for record in self.data["selected_records"]
            if max(
                (
                    stage.get("stable_numerical_rank", -1)
                    for stage in record.get("point_stages", {}).values()
                ),
                default=-1,
            )
            >= 15
        }
        self.assertEqual(set(records), set(expected))
        for parameter, pinned in expected.items():
            record = records[parameter]
            self.assertEqual(record["conductor_phase"]["conductor"], pinned["conductor"])
            self.assertEqual(
                record["conductor_phase"]["log_conductor"], pinned["log_conductor"]
            )
            self.assertEqual(record["conductor_phase"]["root_number"], -1)
            coefficients = family_coefficients(Fraction(parameter))
            for stage_name in ("H50000", "H250000", "H1000000"):
                stage = record["point_stages"][stage_name]
                self.assertEqual(stage["stable_numerical_rank"], 15)
                exact = stage["finite_reduction_attempt"]
                self.assertEqual(exact["certified_algebraic_rank_lower_bound"], 15)
                self.assertEqual(exact["point_sha256"], pinned["point_sha256"])
                points = tuple(
                    (Fraction(point["x"]), Fraction(point["y"]))
                    for point in stage["numerical_subset"]
                )
                self.assertEqual(point_digest(points), pinned["point_sha256"])
                self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))


if __name__ == "__main__":
    unittest.main()

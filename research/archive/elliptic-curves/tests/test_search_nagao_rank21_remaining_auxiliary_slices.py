from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_rank21_remaining_auxiliary_slices.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_remaining_auxiliary_slices.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_rank21_remaining_auxiliary_slices", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class NagaoRank21RemainingAuxiliarySliceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        MODULE.validate_bivariate_quartic()
        cls.preimages, cls.trials = MODULE.select_reconstruction_convention()
        visible_x = {
            point[0]
            for point in MODULE.primitive_visible_points(
                MODULE.RANK21_CONSTRUCTION, MODULE.T0
            )
        }
        cls.accidental = tuple(
            point for point in cls.preimages if point[0] not in visible_x
        )
        cls.slices = MODULE.build_slices(cls.accidental)

    def test_exact_printed_point_reconstruction_and_slice_exclusion(self) -> None:
        self.assertEqual(len(self.preimages), 21)
        self.assertEqual(
            MODULE.point_sequence_digest(self.preimages),
            MODULE.EXPECTED_PREIMAGE_SHA256,
        )
        self.assertEqual(len(self.accidental), 11)
        self.assertEqual(
            MODULE.point_sequence_digest(self.accidental),
            MODULE.EXPECTED_ACCIDENTAL_SHA256,
        )
        winners = [
            record
            for record in self.trials
            if record["generic_abscissa_matches"] == 10
        ]
        self.assertEqual(
            [(record["base_index_one_based"], record["sign"]) for record in winners],
            [(1, 1)],
        )
        self.assertEqual(len(self.slices), 22)
        productive = [
            item
            for item in self.slices
            if item.source_point[0] == MODULE.PRODUCTIVE_SOURCE_X
            and item.slope == MODULE.PRODUCTIVE_SLOPE
        ]
        self.assertEqual(len(productive), 1)
        self.assertEqual(productive[0].intercept, MODULE.PRODUCTIVE_INTERCEPT)
        self.assertEqual(
            len([item for item in self.slices if item != productive[0]]), 21
        )

    def test_pointed_models_and_declared_coefficient_boxes(self) -> None:
        for item in self.slices:
            auxiliary = MODULE.PointedQuartic.from_slice(item)
            self.assertIsNone(
                auxiliary.forward((MODULE.T0, auxiliary.base_ordinate))
            )
            conjugate = auxiliary.conjugate_base_image()
            self.assertEqual(
                auxiliary.forward((MODULE.T0, -auxiliary.base_ordinate)),
                conjugate,
            )
            self.assertTrue(
                MODULE.point_on_short_curve(
                    auxiliary.weierstrass_coefficients, conjugate
                )
            )
            self.assertIsNone(auxiliary.inverse(conjugate))
        self.assertEqual(
            [
                len(MODULE.coefficient_vectors(dimension, MODULE.TERNARY_ALPHABET))
                for dimension in range(7)
            ],
            [0, 2, 8, 26, 80, 242, 728],
        )
        self.assertFalse(
            MODULE.expansion_is_warranted(
                stable_rank=6, accepted_incidences=728, best_proxy=290.0
            )
        )
        self.assertTrue(
            MODULE.expansion_is_warranted(
                stable_rank=4, accepted_incidences=64, best_proxy=249.99
            )
        )

    def test_generated_artifact_closes_the_declared_panel(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.file_sha256(SCRIPT)
        )
        self.assertEqual(data["scope"]["searched_nonproductive_slice_count"], 21)
        self.assertEqual(data["scope"]["slice_height"], 200_000)
        self.assertEqual(
            data["scope"]["excluded_productive_slice"]["intercept"],
            "57361/139",
        )
        expected_ranks = {
            "a01_sm01": 5,
            "a01_sp01": 1,
            "a02_sm01": 1,
            "a02_sp01": 1,
            "a03_sm01": 1,
            "a03_sp01": 1,
            "a04_sm01": 1,
            "a04_sp01": 4,
            "a05_sm01": 1,
            "a05_sp01": 4,
            "a06_sp01": 1,
            "a07_sm01": 1,
            "a07_sp01": 1,
            "a08_sm01": 6,
            "a08_sp01": 1,
            "a09_sm01": 1,
            "a09_sp01": 2,
            "a10_sm01": 1,
            "a10_sp01": 2,
            "a11_sm01": 1,
            "a11_sp01": 3,
        }
        observed_ranks = {
            record["id"]: record["height_selection"]["stable_numerical_rank"]
            for record in data["slices"]
        }
        self.assertEqual(observed_ranks, expected_ranks)
        self.assertEqual(
            sum(
                record["generation"]["ternary"]["coefficient_vectors"]
                for record in data["slices"]
            ),
            1_200,
        )
        self.assertEqual(data["generation"]["unique_accepted_parameter_count"], 883)
        self.assertEqual(data["generation"]["accepted_incidence_count"], 1_118)
        self.assertEqual(data["generation"]["expanded_slice_ids"], [])
        self.assertEqual(data["decontamination"]["global_H200000_or_T0_parameter_count"], 20)
        self.assertGreater(
            data["proxy_filter"]["minimum_log_radical_upper_proxy"], 290.0
        )
        self.assertLess(
            data["proxy_filter"]["minimum_log_radical_upper_proxy"], 291.0
        )
        self.assertEqual(data["proxy_filter"]["below_gate_count"], 0)
        self.assertTrue(
            data["exact_conductors"]["all_proxy_below_190_attempted"]
        )
        self.assertEqual(data["exact_conductors"]["strict_subtarget_parameters"], [])
        self.assertFalse(data["conclusion"]["rank_signal_at_least_18"])
        self.assertFalse(data["conclusion"]["target_hit"])


if __name__ == "__main__":
    unittest.main()

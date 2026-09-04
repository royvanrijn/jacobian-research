from __future__ import annotations

from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-small-field-class-quotient-cohort-v1.json"
)
BUILDER = ROOT / "elkies-k3/scripts/build_r17_small_field_class_quotient_cohort.sage"
FEATURES = ROOT / "elkies-k3/scripts/run_r17_small_field_class_quotient_features.sage"
FREEZER = ROOT / "elkies-k3/scripts/freeze_r17_small_field_class_quotient_detector_protocol.sage"
DETECTOR = ROOT / "elkies-k3/scripts/run_r17_small_field_class_quotient_detector.sage"
ANALYZER = ROOT / "elkies-k3/scripts/analyze_r17_small_field_class_quotient_experiment.py"

EXPECTED_CANDIDATE_HASH = "47cb093c0cbf4803e2cb6c176a45579f241ecf574a48c7c742dce2f053c4d95e"


class R17SmallFieldClassQuotientLabTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(ARTIFACT.read_text())

    def test_exactly_one_hundred_rank_blind_rows_are_frozen(self):
        self.assertEqual(
            self.data["status"], "FROZEN_RANK_BLIND_PRE_CLASS_GROUP_COHORT"
        )
        self.assertEqual(len(self.data["rows"]), 100)
        commitment = self.data["commitment"]
        self.assertEqual(commitment["candidate_list_sha256"], EXPECTED_CANDIDATE_HASH)
        self.assertFalse(commitment["selection_used_rank_or_point_search_labels"])
        self.assertFalse(commitment["selection_used_exceptional_point_coordinates"])
        self.assertFalse(commitment["selection_used_nagao_or_selmer_scores"])

    def test_small_height_distinct_irreducible_structural_cohort(self):
        rows = self.data["rows"]
        self.assertEqual(len({row["sample_id"] for row in rows}), len(rows))
        self.assertEqual(len({row["j_invariant"] for row in rows}), len(rows))
        for row in rows:
            numerator, denominator = row["projective_pair"]
            self.assertGreater(denominator, 0)
            self.assertEqual(gcd(abs(numerator), denominator), 1)
            self.assertLessEqual(max(abs(numerator), denominator), 24)
            self.assertEqual(row["generic_rank"], 17)
        absolute_discriminants = [
            abs(int(row["cubic_order_discriminant"])) for row in rows
        ]
        self.assertEqual(absolute_discriminants, sorted(absolute_discriminants))
        self.assertEqual(
            (
                min(row["absolute_cubic_order_discriminant_bits"] for row in rows),
                max(row["absolute_cubic_order_discriminant_bits"] for row in rows),
            ),
            (286, 353),
        )

    def test_feature_and_outcome_cells_are_sealed(self):
        for row in self.data["rows"]:
            self.assertIsNone(row["class_quotient_features"])
            self.assertEqual(row["feature_status"], "NOT_OPENED")
            self.assertIsNone(row["detector_outcome"])
            self.assertEqual(
                row["outcome_status"], "SEALED_UNTIL_ALL_FEATURES_FREEZE"
            )
        self.assertEqual(
            self.data["phase_order"],
            [
                "0_freeze_rank_blind_cohort",
                "1_compute_and_freeze_complete_unconditional_class_quotient_features",
                "2_freeze_uniform_blind_half_lattice_protocol_bound_to_phase_1_hash",
                "3_run_detector_without_loading_feature_values",
                "4_join_features_and_certified_outcomes_once",
            ],
        )

    def test_mw16_expansion_is_fail_closed(self):
        gate = self.data["family_expansion_gate"]
        self.assertFalse(gate["included_in_v1"])
        self.assertIn("30 seconds", gate["reason"])
        self.assertIn("complete certified BNF", gate["admission_rule"])

    def test_generator_hash_and_replay_inputs_are_pinned(self):
        self.assertEqual(
            self.data["generation"]["script_sha256"],
            sha256(BUILDER.read_bytes()).hexdigest(),
        )
        for path in (
            ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json",
            ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json",
            ROOT
            / "artifacts/generated-results/elliptic-curves"
            / "icarm_curve398_hidden_a1_mw16_v1.json",
        ):
            self.assertEqual(
                self.data["inputs"][str(path.relative_to(ROOT))],
                sha256(path.read_bytes()).hexdigest(),
            )

    def test_phase_boundary_is_executable_and_fail_closed(self):
        feature_source = FEATURES.read_text()
        freezer_source = FREEZER.read_text()
        detector_source = DETECTOR.read_text()
        self.assertIn("pari.bnfcertify(bnf)", feature_source)
        self.assertIn("localized_2_torsion_quotient", feature_source)
        self.assertIn(
            "localized_s_class_group_2_torsion_residual_v1", feature_source
        )
        self.assertIn("dim_Q", feature_source)
        self.assertNotIn(
            '"proves_generic_subgroup_primitive_and_independent"', feature_source
        )
        self.assertIn(
            '"proves_generic_subgroup_2_primitive_and_independent"',
            feature_source,
        )
        self.assertNotIn("run_quartic_search", feature_source)
        self.assertIn("all unconditional Q_t features must freeze", freezer_source)
        self.assertIn("obsolete class quotient", freezer_source)
        self.assertIn("digest(FEATURES)", detector_source)
        self.assertIn("Deliberately do not json.load(FEATURES)", detector_source)

    def test_localized_2_torsion_fixtures(self):
        completed = subprocess.run(
            [
                "sage",
                "-python",
                str(FEATURES),
                "--self-test-localized-class-group",
            ],
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
        self.assertIn(
            "R17LOCALIZED2TORSIONFIXTURES|Z4=PASS|Z8=PASS|S_QUOTIENT=PASS"
            "|NON_TORSION_REJECTION=PASS",
            completed.stdout,
        )

    def test_kendall_implementation_has_expected_direction(self):
        module = SourceFileLoader("small_field_analysis_test", str(ANALYZER)).load_module()
        increasing, _counts = module.kendall_tau_b([0, 1, 2, 3], [0, 2, 2, 5])
        decreasing, _counts = module.kendall_tau_b([0, 1, 2, 3], [5, 2, 2, 0])
        self.assertGreater(increasing, 0)
        self.assertLess(decreasing, 0)


if __name__ == "__main__":
    unittest.main()

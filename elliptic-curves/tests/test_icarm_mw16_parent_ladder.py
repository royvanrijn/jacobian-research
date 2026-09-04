from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
BLIND_INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
INITIAL = ART / "icarm_mw16_parent_ladder_blind_v1.json"
PRESENTATIONS = ART / "icarm_mw16_parent_presentation_audit_v1.json"
NAGAO = ART / "icarm_mw16_parent_nagao_prefilter_h300_v1.json"
SPECIALIZATIONS = ART / "icarm_mw16_nagao_finalist_specializations_h300_v1.json"
PROSPECTIVE = ART / "icarm_mw16_nagao_finalist_half_lattice_h300_v1.json"
CURVE400 = ART / "icarm_mw16_curve400_adaptive_calibration_v1.json"
CALIBRATION = ART / "icarm_mw16_blind_ladder_calibration_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class IcarmMW16ParentLadderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.blind_input = json.loads(BLIND_INPUT.read_text())
        cls.initial = json.loads(INITIAL.read_text())
        cls.presentations = json.loads(PRESENTATIONS.read_text())
        cls.nagao = json.loads(NAGAO.read_text())
        cls.specializations = json.loads(SPECIALIZATIONS.read_text())
        cls.prospective = json.loads(PROSPECTIVE.read_text())
        cls.curve400 = json.loads(CURVE400.read_text())
        cls.calibration = json.loads(CALIBRATION.read_text())

    def test_blind_fixture_and_initial_curve_responses(self) -> None:
        self.assertEqual(
            self.blind_input["status"],
            "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS",
        )
        self.assertEqual(len(self.blind_input["parents"]), 9)
        boundary = self.blind_input["blindness_boundary"]
        self.assertFalse(boundary["public_point_lists_loaded"])
        self.assertFalse(boundary["public_complement_coordinates_loaded"])
        self.assertFalse(boundary["target_rank_lower_bounds_loaded"])
        self.assertEqual(
            self.initial["status"],
            "PASS_COMPLETE_NINE_PARENT_INITIAL_HALF_LATTICE_LADDER",
        )
        responses = {}
        for row in self.initial["parents"]:
            responses.setdefault(int(row["curve_id"]), set()).add(
                int(row["exact_quotient_rank_recovered"])
            )
            self.assertEqual(
                row["discovered_group_saturation"]["status"],
                "PASS_BASIS_EQUALS_DISCOVERED_GROUP",
            )
            self.assertTrue(
                all(
                    cover["search"]["status"] == "bounded_search_complete"
                    for cover in row["cover_records"]
                )
            )
        self.assertEqual(
            responses,
            {398: {5}, 400: {5}, 401: {10}, 542: {10}, 548: {8}},
        )

    def test_nine_presentations_are_five_fibrations(self) -> None:
        self.assertEqual(
            self.presentations["status"],
            "PASS_EXACT_NINE_PRESENTATIONS_FIVE_FIBRATIONS",
        )
        self.assertEqual(self.presentations["presentation_count"], 9)
        self.assertEqual(self.presentations["exact_fibration_class_count"], 5)
        self.assertEqual(len(self.presentations["clusters"]), 5)
        self.assertEqual(
            sorted(
                len(row["presentation_ids"])
                for row in self.presentations["clusters"]
            ),
            [1, 1, 2, 2, 3],
        )
        self.assertEqual(len(self.presentations["cross_cluster_non_equivalences"]), 10)

    def test_curve400_complete_adaptive_recovery(self) -> None:
        self.assertEqual(
            self.curve400["status"],
            "PASS_COMPLETE_CURVE400_ADAPTIVE_CALIBRATION",
        )
        initial = self.curve400["initial"]
        adaptive = self.curve400["adaptive"]
        self.assertEqual((initial["basis_rank_after"], adaptive["basis_rank_after"]), (21, 28))
        self.assertEqual(self.curve400["exact_quotient_rank_recovered_total"], 12)
        self.assertEqual(len(initial["cover_records"]), 4)
        self.assertEqual(len(adaptive["cover_records"]), 124)
        self.assertEqual(
            Counter(row["search"]["status"] for row in adaptive["cover_records"]),
            Counter({"bounded_search_complete": 124}),
        )
        self.assertEqual(
            Counter(
                row["type"]
                for row in adaptive["discovered_group_saturation"]["events"]
            ),
            Counter({"NEW_Q_INDEPENDENT_DIRECTION": 7}),
        )

    def test_five_curve_heldout_calibration(self) -> None:
        self.assertEqual(
            self.calibration["status"],
            "PASS_EXACT_FIVE_CURVE_BLIND_LADDER_CALIBRATION",
        )
        summary = self.calibration["curve_level_summary"]
        self.assertEqual(summary["demonstrated_jump_directions_total"], 55)
        self.assertEqual(summary["initial_exact_directions_recovered_total"], 38)
        self.assertEqual(summary["best_blind_ladder_exact_directions_recovered_total"], 54)
        self.assertEqual(summary["initial_full_recovery_curve_count"], 2)
        self.assertEqual(summary["best_blind_ladder_full_recovery_curve_count"], 4)
        self.assertEqual(
            summary["only_remaining_demonstrated_direction"],
            {
                "curve_id": 401,
                "count": 1,
                "complete_next_wave_chart_count": 8184,
                "complete_next_wave_run": False,
            },
        )

    def test_prospective_gate_is_wholly_censored(self) -> None:
        self.assertEqual(
            self.nagao["status"],
            "PASS_BOUNDED_NINE_PRESENTATION_NAGAO_PREFILTER",
        )
        self.assertEqual(
            self.nagao["search"][
                "total_finalist_rows_before_exact_curve_deduplication"
            ],
            104,
        )
        self.assertEqual(
            self.specializations["status"],
            "PASS_EXACT_MW16_NAGAO_FINALIST_SPECIALIZATIONS",
        )
        self.assertEqual(self.specializations["successful_specialization_count"], 104)
        self.assertEqual(self.specializations["exact_q_isomorphism_class_count"], 104)
        self.assertEqual(
            self.prospective["status"],
            "PASS_COMPLETE_FROZEN_NAGAO_FINALIST_HALF_LATTICE_GATE",
        )
        self.assertEqual(self.prospective["completed_candidate_count"], 104)
        self.assertEqual(self.prospective["positive_candidate_count"], 0)
        self.assertEqual(self.prospective["failed_closed_candidate_count"], 0)
        self.assertEqual(
            self.prospective["chart_status_counts"],
            {"bounded_search_timeout": 856},
        )
        self.assertTrue(
            self.prospective["arithmetic_size_diagnostic"][
                "all_chart_attempts_timed_out"
            ]
        )
        self.assertEqual(
            self.prospective["selmer_gate"][
                "complete_residual_2_selmer_calls_required"
            ],
            0,
        )

    def test_generated_hash_chains(self) -> None:
        for document in (
            self.blind_input,
            self.initial,
            self.presentations,
            self.nagao,
            self.specializations,
            self.curve400,
            self.calibration,
        ):
            inputs = document.get("inputs", document.get("input_hashes", {}))
            for name, expected in inputs.items():
                path = ROOT / name
                if path.is_file():
                    self.assertEqual(digest(path), expected, name)
        prospective_inputs = self.prospective["inputs"]
        for path in (SPECIALIZATIONS, ROOT / "elliptic-curves/cas/merge_icarm_mw16_nagao_finalist_half_lattice_shards.py"):
            name = str(path.relative_to(ROOT))
            self.assertEqual(prospective_inputs[name], digest(path))


if __name__ == "__main__":
    unittest.main()

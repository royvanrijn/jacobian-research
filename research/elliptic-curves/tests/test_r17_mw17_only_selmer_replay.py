from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import build_r17_mw17_only_selmer_replay as builder  # noqa: E402
import run_r17_mw17_only_selmer_replay as runner  # noqa: E402
from build_elkies_2026_relative_2selmer_suite import magma_points  # noqa: E402


class MW17OnlySelmerReplayTests(unittest.TestCase):
    def test_generated_programs_are_fixture_separated(self) -> None:
        document, programs = builder.build()
        self.assertEqual(document["case_count"], 2)
        self.assertFalse(
            document["operational_gate"][
                "selmer_candidate_gate_operationally_calibrated"
            ]
        )
        cases = builder.load_record_pair_cases()
        by_curve = {
            int(case.case_id.rsplit("-", 1)[1]): case for case in cases
        }
        for record, program in zip(document["cases"], programs.values()):
            self.assertTrue(all(record["program_audit"].values()))
            self.assertEqual(program.count("E!["), 17)
            self.assertNotIn("exceptional", program)
            self.assertNotIn("MW29", program)
            self.assertNotIn("TwoCover(", program)
            self.assertNotIn("Points(", program)
            self.assertTrue(program.rstrip().endswith("quit;"))
            for point in by_curve[record["curve_id"]].exceptional_points:
                declaration = magma_points("probe", [point])
                start = declaration.index("E![")
                row = declaration[start : declaration.index(", 1]", start) + 4]
                self.assertNotIn(row, program)

    @staticmethod
    def complete_transcript(case_id: str, residual_dimension: int = 12) -> str:
        total_dimension = 17 + residual_dimension

        def bits(row: list[int]) -> str:
            return "[ " + ", ".join(map(str, row)) + " ]"

        lines = [
            (
                f"{runner.PROTOCOL}|version=1|stage=input|case={case_id}|"
                "role=prospective-mw17-only-selmer-control|generic_count=17|"
                "known_subgroup=MW17|fixture_access=false|magma=V2.29-1"
            ),
            (
                f"{runner.PROTOCOL}|stage=two_selmer|status=complete|"
                f"total_dim={total_dimension}|two_torsion_dim=0|"
                f"generic_kummer_rank=17|residual_dim={residual_dimension}|"
                "seconds=1|factor_base_size=1"
            ),
        ]
        for index in range(17):
            row = [int(column == index) for column in range(total_dimension)]
            lines.append(
                f"{runner.PROTOCOL}|stage=generic_class|index={index + 1}|"
                f"selmer_bits={bits(row)}"
            )
        for index in range(residual_dimension):
            selmer = [
                int(column == 17 + index) for column in range(total_dimension)
            ]
            quotient = [
                int(column == index) for column in range(residual_dimension)
            ]
            lines.append(
                f"{runner.PROTOCOL}|stage=quotient_basis|index={index + 1}|"
                f"selmer_bits={bits(selmer)}|quotient_bits={bits(quotient)}|alpha=a"
            )
        lines.append(
            f"{runner.PROTOCOL}|stage=blind_freeze|status=complete|"
            f"total_dim={total_dimension}|generic_kummer_rank=17|"
            f"residual_dim={residual_dimension}"
        )
        return "\n".join(lines) + "\n"

    def test_complete_transcript_passes_only_at_blind_freeze(self) -> None:
        expected = {"case_id": "mw17-only-control-356", "curve_id": 356}
        result = runner.parse_complete_transcript(
            self.complete_transcript(expected["case_id"]), expected, 12
        )
        self.assertEqual(result["selmer_modulo_mw17_dimension"], 12)
        self.assertEqual(result["magma_version"], "V2.29-1")
        self.assertTrue(result["blind_control_detected_required_jump"])

        leaked = self.complete_transcript(expected["case_id"]) + (
            f"{runner.PROTOCOL}|stage=post_freeze_fixture|status=complete\n"
        )
        with self.assertRaises(ValueError):
            runner.parse_complete_transcript(leaked, expected, 12)

    def test_committed_run_is_fail_closed(self) -> None:
        result = json.loads(runner.OUTPUT.read_text())
        self.assertEqual(
            result["status"], "INCOMPLETE_BLINDED_REPLAY_NOT_EXECUTED"
        )
        self.assertEqual(result["completed_record_replays"], 0)
        self.assertFalse(result["selmer_candidate_gate_operationally_calibrated"])
        self.assertFalse(result["prospective_sample_stage_authorized"])


if __name__ == "__main__":
    unittest.main()

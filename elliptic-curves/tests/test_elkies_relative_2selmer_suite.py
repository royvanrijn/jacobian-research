from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import build_elkies_2026_relative_2selmer_suite as builder  # noqa: E402
import parse_elkies_2026_relative_2selmer_suite as result_parser  # noqa: E402


class ElkiesRelative2SelmerSuiteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rank21 = builder.load_rank21_case()
        cls.high_rank = builder.load_high_rank_cases()

    def test_all_five_controls_replay(self) -> None:
        cases = [self.rank21, *self.high_rank]
        self.assertEqual([case.certified_rank_lower_bound for case in cases], [21, 25, 26, 27, 28])
        self.assertEqual([len(case.exceptional_points) for case in cases], [4, 8, 9, 10, 11])
        for case in cases:
            self.assertEqual(len(case.generic_points), 17)
            self.assertEqual(
                case.certified_rank_lower_bound,
                17 + len(case.exceptional_points),
            )

    def test_public_points_are_declared_after_blind_search(self) -> None:
        program = builder.build_magma(
            self.rank21, search_bound=1000, enumerate_class_limit=4095
        )
        blind_end = program.index("stage=blind_end")
        exceptional_declaration = program.index("exceptional :=")
        self.assertLess(blind_end, exceptional_declaration)
        self.assertIn("Hints := generic_hints", program)
        self.assertIn("AtoS(mu(P))", program)
        self.assertIn("SelmerGroup(", program)
        self.assertIn("MultiplicationByMMap(E, 2)", program)
        self.assertIn("Bound := -1", program)
        self.assertIn("Raw := true", program)
        self.assertIn("TwoCover(alpha : E := E)", program)
        self.assertIn("Points(C : Bound := 1000)", program)
        self.assertNotIn('SetClassGroupBounds("GRH")', program)

    def test_top_nagao_candidates_use_the_frozen_shape(self) -> None:
        candidates = builder.load_nagao_cases(2)
        self.assertEqual([case.case_id for case in candidates], ["nagao-0001", "nagao-0002"])
        self.assertEqual([case.parameter for case in candidates], ["-5643/6760", "1452/7817"])
        self.assertTrue(all(len(case.generic_points) == 17 for case in candidates))
        self.assertTrue(all(not case.exceptional_points for case in candidates))

    def test_parser_classifies_exceptional_and_blind_recovery_spans(self) -> None:
        generic_rows = "\n".join(
            f"ELKIESR17REL2|stage=generic_class|index={i}|selmer_bits=[ {', '.join(['1' if j == i else '0' for j in range(1, 22)])} ]"
            for i in range(1, 18)
        )
        quotient_rows = "\n".join(
            f"ELKIESR17REL2|stage=quotient_basis|index={i}|selmer_bits=[ {', '.join(['1' if j == 17 + i else '0' for j in range(1, 22)])} ]|quotient_bits=[ {', '.join(['1' if j == i else '0' for j in range(1, 5)])} ]|alpha=a{i}"
            for i in range(1, 5)
        )
        covers = "\n".join(
            (
                f"ELKIESR17REL2|stage=blind_cover|index={i}|label=basis-{i}|quotient_bits=[ {', '.join(['1' if j == i else '0' for j in range(1, 5)])} ]|alpha=a{i}|quartic_f=x^4+{i}|quartic_h=0|construction_seconds=0.{i}|search_seconds=1.{i}|search_status="
                + (
                    f"point_found|cover_point=(1 : 1 : 1)|elliptic_point=({i} : {i} : 1)|recovered_quotient_bits=[ {', '.join(['1' if j == i else '0' for j in range(1, 5)])} ]"
                    if i <= 2
                    else "no_point_within_bound"
                )
            )
            for i in range(1, 5)
        )
        exceptional = "\n".join(
            f"ELKIESR17REL2|stage=exceptional_class|index={i}|selmer_bits=[ {', '.join(['1' if j == 17 + i else '0' for j in range(1, 22)])} ]|quotient_bits=[ {', '.join(['1' if j == i else '0' for j in range(1, 5)])} ]"
            for i in range(1, 5)
        )
        transcript = "\n".join(
            (
                "ELKIESR17REL2|version=1|stage=input|case=control-r21-t3_8|role=held-out-positive-control|parameter=3/8|generic_count=17|held_out_exceptional_count=4|search_bound=1000|enumerate_class_limit=1|provenance_sha256=x|magma=V2.29",
                "ELKIESR17REL2|stage=two_selmer|status=start|bound=-1|raw=true|hints=generic_only",
                "ELKIESR17REL2|stage=two_selmer|status=complete|seconds=12.5|total_dim=21|two_torsion_dim=0|generic_kummer_rank=17|residual_dim=4|factor_base_size=99",
                generic_rows,
                quotient_rows,
                "ELKIESR17REL2|stage=blind_plan|residual_dim=4|nonzero_class_count=15|enumerate_all=false|target_count=4",
                covers,
                "ELKIESR17REL2|stage=blind_end|seconds=9.0|target_count=4|recovered_class_count=2|recovered_quotient_rank=2",
                exceptional,
                "ELKIESR17REL2|stage=classification|status=complete|exceptional_count=4|exceptional_quotient_rank=4|known_realized_class_count=16|unexplained_dim=0|unrealized_class_count=0|blind_recovered_rank=2",
            )
        )
        expected = {
            "case_id": "control-r21-t3_8",
            "role": "held-out-positive-control",
            "parameter": "3/8",
            "global_minimal_model": ["model"],
            "certified_rank_lower_bound": 21,
            "held_out_exceptional_point_count": 4,
            "nagao_record": None,
        }
        parsed = result_parser.parse_case(transcript, expected)
        self.assertEqual(parsed["two_selmer"]["dimension"], 21)
        self.assertEqual(parsed["relative_quotient"]["dimension"], 4)
        self.assertEqual(parsed["held_out_exceptional_points"]["quotient_rank"], 4)
        self.assertEqual(parsed["unrealized_quotient"]["class_count"], 0)
        self.assertEqual(
            parsed["blind_cover_benchmark"]["recovered_known_exceptional_direction_rank"],
            2,
        )


if __name__ == "__main__":
    unittest.main()

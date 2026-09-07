from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_relative_2selmer_checkpointed as checkpointed  # noqa: E402
from build_elkies_2026_relative_2selmer_suite import load_record_pair_cases  # noqa: E402


class CheckpointedRelativeSelmerTests(unittest.TestCase):
    def test_mw29_backend_reverses_the_old_blind_order(self) -> None:
        source = (CAS / "run_mw29_relative_2selmer_from_bnf.sage").read_text()
        self.assertIn("post-discovery MW29-relative closure only", source)
        self.assertIn("prospective_mw17_candidate_gate_evidence", source)
        self.assertIn("ell2global_norm_envelope", source)
        self.assertIn("point_coordinates(", source)
        self.assertIn("known_norm_rows", source)
        self.assertIn("ell2allowed_at_place", source)
        runtime = source[source.index("def main()") :]
        self.assertLess(
            runtime.index("= point_coordinates("),
            runtime.index('pari("ell2allowed_at_place")'),
        )
        self.assertNotIn('pari("ell2selmer_basis_gen")', runtime)
        self.assertNotIn("QUOTIENT_COVER_WORKER", source)
        self.assertIn("if residual_upper == 0:\n            break", source)

    def test_mw29_backend_accepts_only_a_complete_unconditional_hecke_envelope(self) -> None:
        source = (CAS / "run_mw29_relative_2selmer_from_bnf.sage").read_text()
        provider = (
            CAS / "run_elkies_2026_curve356_hecke_s_squareclasses.jl"
        ).read_text()
        self.assertIn("--hecke-global-log", source)
        self.assertIn("parse_hecke_global_log", source)
        self.assertIn("stage=class_group|status=PASS|grh=false", source)
        self.assertIn("ell2allowed_at_place_nf", source)
        self.assertIn("class_group(O; GRH=false)", provider)
        self.assertIn("unit_group(O; GRH=false)", provider)
        self.assertIn("pselmer_group(2, S; algo=:raw)", provider)
        self.assertIn("norm_is_square", provider)

    def test_workers_preserve_blind_order(self) -> None:
        self.assertIn("prepared_bnf", checkpointed.BNF_WORKER)
        self.assertIn("bnfcertify", checkpointed.BNF_WORKER)
        self.assertIn("writebin", checkpointed.BNF_WORKER)
        self.assertIn("prepared_polredabs", checkpointed.BNF_WORKER)
        self.assertIn("prepared_simon", checkpointed.SELMER_WORKER)
        self.assertNotIn("pari.bnfcertify", checkpointed.SELMER_WORKER)
        self.assertIn("curve_theta_in_field", checkpointed.BNF_WORKER)
        self.assertIn("ell2selmer_basis_gen", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("elllocalimage_mapped", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("nfeltsign", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("Vec(localspaces)", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("Vec(localmeta)", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("descentprimes", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("matinverseimage(normspace,selmer)", checkpointed.SIMON_GP_FUNCTION)
        self.assertIn("local_allowed_subspaces", checkpointed.SELMER_WORKER)
        self.assertIn("selmer_in_normspace", checkpointed.SELMER_WORKER)
        self.assertIn("finite_local_condition_primes", checkpointed.SELMER_WORKER)
        self.assertIn("auxiliary_descent_primes", checkpointed.SELMER_WORKER)
        self.assertIn("elliptic_bad_local_data", checkpointed.SELMER_WORKER)
        self.assertIn("component_group_data", checkpointed.SELMER_WORKER)
        self.assertIn(
            "global_norm_square_subspace_basis_columns_in_s_squareclasses",
            checkpointed.SELMER_WORKER,
        )
        self.assertNotIn('payload["exceptional_points"]', checkpointed.SELMER_WORKER)
        self.assertNotIn('payload["generic_points"]', checkpointed.SELMER_WORKER)
        self.assertIn("cover_for", checkpointed.QUOTIENT_COVER_WORKER)
        self.assertNotIn('payload["exceptional_points"]', checkpointed.QUOTIENT_COVER_WORKER)

    def test_standard_extension(self) -> None:
        initial = [[1, 1, 0, 0], [0, 1, 1, 0]]
        extension = checkpointed.extend_standard(initial, 4)
        self.assertEqual(len(extension), 2)
        self.assertEqual(checkpointed.f2_rank([*initial, *extension]), 4)
        self.assertTrue(all(sum(row) == 1 for row in extension))

    def test_control_classification(self) -> None:
        class Case:
            exceptional_points = ((0, 0),) * 2

        selmer = {"two_selmer_dimension": 20}
        generic_rows = [
            [1 if column == row else 0 for column in range(20)]
            for row in range(17)
        ]
        exceptional_rows = [
            [1 if column == row else 0 for column in range(20)]
            for row in (17, 18)
        ]
        result = checkpointed.combine_results(
            Case(),
            selmer,
            {"point_selmer_rows": generic_rows},
            {"point_selmer_rows": exceptional_rows},
            {
                "quotient_basis": [
                    [1 if column == row else 0 for column in range(20)]
                    for row in (17, 18, 19)
                ],
                "enumeration_complete": True,
                "classes": [
                    {
                        "quotient_class_integer": mask,
                        "quotient_bits": [(mask >> index) & 1 for index in range(3)],
                        "blind_search": {"status": "no_point_within_bound"},
                    }
                    for mask in range(1, 8)
                ],
            },
        )
        self.assertEqual(result["quotient_dimension"], 3)
        self.assertEqual(result["exceptional_quotient_rank"], 2)
        self.assertEqual(result["classes_not_realized_by_known_exceptional_subgroup"], 4)
        labeled = {
            row["quotient_class_integer"]: row
            for row in result["quotient_cover_classification"]
        }
        self.assertTrue(labeled[3]["known_exceptional_subgroup_realizes_class"])
        self.assertFalse(labeled[4]["known_exceptional_subgroup_realizes_class"])

    def test_record_pair_inputs_and_rigid_quotient(self) -> None:
        cases = load_record_pair_cases()
        self.assertEqual([case.case_id for case in cases], ["record-r29-356", "record-r29-385"])
        self.assertEqual([len(case.exceptional_points) for case in cases], [12, 12])
        self.assertEqual(
            [case.rigid_complement_point_labels for case in cases],
            [
                ("P18", "P20", "P21", "P22", "P23", "P24", "P26", "P27", "P28", "P29"),
                ("P18", "P19", "P20", "P22", "P23", "P25", "P26", "P27", "P28", "P29"),
            ],
        )

        case = cases[0]
        generic_rows = [
            [1 if column == row else 0 for column in range(29)]
            for row in range(17)
        ]
        exceptional_rows = [
            [1 if column == row else 0 for column in range(29)]
            for row in range(17, 29)
        ]
        result = checkpointed.combine_results(
            case,
            {"two_selmer_dimension": 29},
            {"point_selmer_rows": generic_rows},
            {"point_selmer_rows": exceptional_rows},
        )
        rigid = result["rigid_character_quotient"]
        self.assertEqual(result["quotient_dimension"], 12)
        self.assertEqual(rigid["rigid_image_dimension"], 2)
        self.assertEqual(rigid["dimension_after_quotienting_rigid_plane"], 10)
        self.assertEqual(rigid["displayed_nonrigid_image_dimension"], 10)
        self.assertEqual(rigid["additional_dimension_beyond_all_twenty_nine_known_points"], 0)


if __name__ == "__main__":
    unittest.main()

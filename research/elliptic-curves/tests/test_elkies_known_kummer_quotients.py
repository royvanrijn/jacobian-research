from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from run_elkies_2026_relative_2selmer_open import f2_rank  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_kummer_quotients_suite_v1.json"
)
PROGRAM = CAS / "audit_elkies_2026_known_kummer_quotients.py"
COVER_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_exceptional_quotient_covers_v1.json"
)
BOTTLENECK_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_relative_2selmer_open_bottleneck_benchmarks_v2.json"
)
HECKE_MONITOR = CAS / "run_elkies_2026_s_class_hecke_monitor.jl"


class KnownKummerQuotientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())

    def test_suite_passes_all_pinned_cases(self) -> None:
        self.assertEqual(
            self.data["status"], "PASS_ALL_EXACT_KNOWN_KUMMER_QUOTIENTS"
        )
        self.assertEqual(self.data["run_count"], 15)
        self.assertEqual(self.data["passed_count"], 15)

    def test_control_ranks_and_quotients(self) -> None:
        expected = {
            "control-r21-t3_8": (21, 4),
            "control-r25": (25, 8),
            "control-r26": (26, 9),
            "control-r27": (27, 10),
            "control-r28": (28, 11),
        }
        records = {row["case_id"]: row for row in self.data["runs"]}
        for case_id, (known_rank, quotient_dimension) in expected.items():
            row = records[case_id]
            self.assertEqual(row["generic_kummer_rank"], 17)
            self.assertEqual(row["known_kummer_rank"], known_rank)
            self.assertEqual(
                row["known_exceptional_quotient_dimension"], quotient_dimension
            )
            self.assertEqual(
                row["known_realized_exceptional_quotient_class_count_including_zero"],
                2**quotient_dimension,
            )
            point_rows = [point["local_squareclass_row"] for point in row["points"]]
            self.assertEqual(f2_rank(point_rows), known_rank)

    def test_candidates_certify_only_the_generic_subgroup(self) -> None:
        candidates = [
            row
            for row in self.data["runs"]
            if row["role"] == "prospective-high-Nagao-candidate"
        ]
        self.assertEqual(len(candidates), 10)
        for row in candidates:
            self.assertEqual(row["generic_kummer_rank"], 17)
            self.assertEqual(row["known_kummer_rank"], 17)
            self.assertEqual(row["exceptional_point_count"], 0)
            self.assertEqual(row["known_exceptional_quotient_dimension"], 0)

    def test_method_is_class_group_free_and_claims_are_bounded(self) -> None:
        source = PROGRAM.read_text()
        self.assertNotIn("pari.nfinit", source)
        self.assertNotIn("bnfinit", source)
        boundaries = " ".join(self.data["claim_boundary"])
        self.assertIn("does not compute or upper-bound the full 2-Selmer group", boundaries)

    def test_all_known_exceptional_classes_have_explicit_covers(self) -> None:
        covers = json.loads(COVER_ARTIFACT.read_text())
        self.assertEqual(
            covers["status"],
            "PASS_ALL_CONTROL_KNOWN_EXCEPTIONAL_QUOTIENT_COVERS",
        )
        self.assertEqual(covers["total_nonzero_known_quotient_classes"], 3851)
        expected_dimensions = [4, 8, 9, 10, 11]
        self.assertEqual(
            [run["known_exceptional_quotient_dimension"] for run in covers["runs"]],
            expected_dimensions,
        )
        for run, dimension in zip(covers["runs"], expected_dimensions):
            self.assertEqual(
                run["nonzero_known_exceptional_quotient_class_count"],
                2**dimension - 1,
            )
            self.assertEqual(len(run["basis_class_records"]), dimension)
            self.assertTrue(
                all(
                    row["rational_cover_witness_verified"]
                    for row in run["basis_class_records"]
                )
            )
        full = covers["full_all_class_output"]
        full_path = Path(full["path"])
        self.assertTrue(full_path.is_file())
        self.assertEqual(sha256(full_path.read_bytes()).hexdigest(), full["sha256"])

    def test_specialized_s_class_failures_are_pinned_and_bounded(self) -> None:
        benchmark = json.loads(BOTTLENECK_ARTIFACT.read_text())
        specialized = benchmark["hecke_s_class_specialized_experiments"]
        self.assertEqual(
            sha256(HECKE_MONITOR.read_bytes()).hexdigest(),
            specialized["program_sha256"],
        )
        for key in (
            "all_S_ideals_augmented_factorbase",
            "direct_S_multiplier_collector",
        ):
            run = specialized[key]
            self.assertEqual(run["declared_S_columns_visible_in_factorbase"], 25)
            self.assertEqual(run["missing_S_ideals"], 0)
            self.assertEqual(run["factorbase_mod2_quotient_dimension"], 34)
            self.assertIn("zero-gain", run["outcome"])
            checkpoint = ROOT / run["checkpoint"]
            self.assertTrue(checkpoint.is_file())
            self.assertEqual(
                sha256(checkpoint.read_bytes()).hexdigest(),
                run["checkpoint_sha256"],
            )
        self.assertIn("not a class-group", specialized["warning"])


if __name__ == "__main__":
    unittest.main()

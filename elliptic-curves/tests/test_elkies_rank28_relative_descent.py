from __future__ import annotations

from pathlib import Path
import json
import inspect
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import build_elkies_2026_rank28_relative_descent_magma as builder  # noqa: E402
import build_elkies_2026_rank28_bad_place_ledger as bad_places  # noqa: E402
import parse_elkies_2026_rank28_relative_descent as parser  # noqa: E402
import run_elkies_2026_rank28_residual_selmer as supervised  # noqa: E402
from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402


class ElkiesRank28RelativeDescentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = builder.load_relative_input()
        cls.program = builder.build_magma(cls.source)

    def test_exact_generic_and_public_complement_replay(self) -> None:
        self.assertEqual(len(self.source.generic_points), 17)
        self.assertEqual(len(self.source.public_complement), 11)
        for point in self.source.generic_points + self.source.public_complement:
            self.assertTrue(is_on_weierstrass_curve(self.source.model, point))

    def test_complete_selmer_precedes_gate_and_covers(self) -> None:
        complete = self.program.index('stage=two_selmer|status=complete')
        classification = self.program.index("if unexplained_dim lt 4 then")
        covers = self.program.index("covers := TwoDescent(")
        self.assertLess(complete, classification)
        self.assertLess(classification, covers)
        self.assertIn("TwoSelmerGroup(E : Bound := -1)", self.program)
        self.assertIn(
            "RemoveGens := SequenceToSet(generic cat public_complement)",
            self.program,
        )
        self.assertIn("WithMaps := false", self.program)
        self.assertIn("assert residual_dim ge 11", self.program)
        self.assertIn("assert unexplained_dim ge 0", self.program)
        self.assertIn("assert #covers + 1 eq 2^unexplained_dim", self.program)

    def test_program_has_no_point_search_or_conditional_class_group(self) -> None:
        executable = "\n".join(self.program.splitlines()[3:])
        for forbidden in (
            "PointsQI",
            "PointSearch",
            "Ratpoints",
            "ratpoints",
            "FourDescent",
            'SetClassGroupBounds("GRH")',
        ):
            self.assertNotIn(forbidden, executable)

    def test_local_eclib_worker_is_strictly_selmer_only(self) -> None:
        worker = supervised.ECLIB_WORKER
        self.assertIn("selmer_only=True", worker)
        self.assertIn("first_limit=0", worker)
        self.assertIn("second_limit=0", worker)
        self.assertIn("backend.selmer_rank()", worker)
        self.assertNotIn("POINTS", worker)

    def test_factored_pari_worker_uses_only_proved_factor_hints(self) -> None:
        worker = supervised.PARI_FACTORED_WORKER
        self.assertIn("pari.allocatemem(PARI_STACK_BYTES)", worker)
        self.assertIn("pari.addprimes(factor_hint_primes)", worker)
        self.assertIn("ellrank(0, known)", worker)
        self.assertNotIn("ratpoints", worker.lower())
        self.assertEqual(len(bad_places.DISCRIMINANT_FACTORIZATION), 12)
        self.assertEqual(
            bad_places.factorization_product(),
            bad_places.load_input(bad_places.DEFAULT_SPECIALIZATION)["discriminant"],
        )

    def test_bad_place_ledger_is_complete_but_not_a_selmer_bound(self) -> None:
        ledger = json.loads(supervised.BAD_PLACE_LEDGER.read_text())
        self.assertEqual(
            ledger["status"],
            "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND",
        )
        self.assertTrue(ledger["factorization_product_verified"])
        self.assertTrue(ledger["factor_primality_proof_completed"])
        controls = json.loads(supervised.CONTROL_CERTIFICATE.read_text())
        self.assertEqual(
            int(ledger["descent_cubic_discriminant"]),
            256 * int(controls["fibres"][-1]["minimal_discriminant"]),
        )
        self.assertEqual(ledger["completed_local_block_count"], 13)
        self.assertEqual(ledger["combined_known_kummer_rank"], 15)
        certificate = supervised.validate_factor_hint_certificate(
            [str(value) for value in self.source.model]
        )
        self.assertEqual(certificate["factor_count"], 12)
        self.assertIn("not themselves a Selmer upper bound", certificate["claim_boundary"])

    def test_factored_pari_timeout_remains_fail_closed(self) -> None:
        directory = ROOT / "artifacts/generated-results/elliptic-curves"
        for name, timeout in (
            ("elkies_2026_rank28_residual_2selmer_pari_factored_8g_v1.json", 600),
            (
                "elkies_2026_rank28_residual_2selmer_pari_factored_8g_30min_v1.json",
                1800,
            ),
        ):
            with self.subTest(name=name):
                artifact = json.loads((directory / name).read_text())
                self.assertEqual(
                    artifact["status"],
                    "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN",
                )
                self.assertTrue(
                    artifact["descent_backend"]["factorization_supplied"]
                )
                self.assertIsNone(artifact["backend_result"])
                self.assertFalse(artifact["gate"]["expensive_search_authorized"])
                self.assertEqual(
                    artifact["supervisor"]["outcome"], "strict_wall_timeout"
                )
                self.assertEqual(artifact["supervisor"]["timeout_seconds"], timeout)
                self.assertEqual(
                    artifact["factor_hint_certificate"]["sha256"],
                    supervised.file_sha256(supervised.BAD_PLACE_LEDGER),
                )

    def test_factor_supplied_number_field_is_still_certified(self) -> None:
        source = inspect.getsource(bad_places._worker_setup)
        self.assertIn("pari.nfinit([pari(polynomial), ramified_primes])", source)
        self.assertIn("pari.nfcertify(nf)", source)

    def test_protocol_rejects_residual_fourteen(self) -> None:
        parsed = parser.parse_protocol(
            "\n".join(
                (
                    "ELKIESR28REL|version=1|stage=input|parameter=-9529/5471|generic=17|known_quotient_floor=11|known_rank_lower_bound=28|target_rank=32|required_beyond_known_rank28=4",
                    "ELKIESR28REL|stage=two_selmer|status=complete|total_selmer_dim=31|residual_dim=14|required_residual_dim=15|unexplained_dim=3|required_unexplained_dim=4",
                    "ELKIESR28REL|classification=REJECT_RANK32_BY_RESIDUAL_2_SELMER|total_selmer_dim=31|residual_dim=14|unexplained_dim=3|expensive_search_authorized=false",
                )
            )
        )
        self.assertFalse(parsed["gate"]["expensive_search_authorized"])

    def test_protocol_passes_residual_fifteen(self) -> None:
        parsed = parser.parse_protocol(
            "\n".join(
                (
                    "ELKIESR28REL|version=1|stage=input|parameter=-9529/5471|generic=17|known_quotient_floor=11|known_rank_lower_bound=28|target_rank=32|required_beyond_known_rank28=4",
                    "ELKIESR28REL|stage=two_selmer|status=complete|total_selmer_dim=32|residual_dim=15|required_residual_dim=15|unexplained_dim=4|required_unexplained_dim=4",
                    "ELKIESR28REL|classification=PASS_RANK32_RESIDUAL_2_SELMER_GATE|total_selmer_dim=32|residual_dim=15|unexplained_dim=4|expensive_search_authorized=true",
                )
            )
        )
        self.assertTrue(parsed["gate"]["expensive_search_authorized"])

    def test_old_low_complexity_raw_search_is_parked(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(
                    ROOT
                    / "elliptic-curves/scripts/search_q12o5867_low_complexity_x_ansatz.py"
                ),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("PARKED_BY_ELKIES_2026_DESCENT_FIRST_POLICY", completed.stderr)


if __name__ == "__main__":
    unittest.main()

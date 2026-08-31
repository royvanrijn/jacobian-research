from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_rank28_s_class_pari as s_class  # noqa: E402


class ElkiesRank28SClassPariTests(unittest.TestCase):
    def test_exact_rank28_inputs_and_factor_table_are_pinned(self) -> None:
        controls, ledger, coefficients, primes = s_class.validate_inputs()
        self.assertEqual(controls["fibres"][-1]["parameter"], "-9529/5471")
        self.assertEqual(coefficients[-1], 1)
        self.assertEqual(len(primes), 12)
        self.assertTrue(ledger["factor_primality_proof_completed"])

    def test_mod2_s_class_quotient_formula(self) -> None:
        cyclics = [2, 3, 4, 6]
        rows = [
            [1, 0, 1, 0],
            [0, 1, 1, 1],
        ]
        self.assertEqual(s_class.s_class_mod2_upper_bound(cyclics, rows), 1)
        with self.assertRaisesRegex(ValueError, "row length"):
            s_class.s_class_mod2_upper_bound(cyclics, [[1, 0]])

    def test_worker_uses_factor_supplied_nfinit_and_one_sided_certificate(self) -> None:
        _, _, coefficients, primes = s_class.validate_inputs()
        source = s_class.worker_source(
            coefficients=coefficients,
            primes=primes,
            stack_bytes=1_000_000_000,
            mode="class-quotient",
            tech=[0.01, 4.0, 20],
            debug=1,
        )
        self.assertIn("nfinit([polynomial, BAD_PRIMES])", source)
        self.assertIn("BNF_FLAG = 0", source)
        self.assertIn("CERTIFY_FLAG = 1", source)
        self.assertIn("bnfcertify(bnf, CERTIFY_FLAG)", source)
        self.assertIn("bnfisprincipal(bnf, prime_ideal, 0)", source)
        self.assertNotIn("ellsearch", source.lower())
        self.assertNotIn("ratpoints", source.lower())

    def test_protocol_preserves_last_completed_stage(self) -> None:
        result = {"s_class_group_mod2_dimension_upper_bound": 2}
        events, parsed = s_class.parse_protocol(
            "\n".join(
                (
                    "ELKIESR28SCLASS|stage=nfinit|status=start",
                    "ELKIESR28SCLASS|stage=nfinit|status=complete|seconds=0.1",
                    "ELKIESR28SCLASS|stage=bnfinit|status=start|flag=0",
                    "ELKIESR28SCLASS|result=" + json.dumps(result),
                )
            )
        )
        self.assertEqual(events[-1]["stage"], "bnfinit")
        self.assertEqual(parsed, result)

    def test_debug_progress_is_diagnostic_only(self) -> None:
        progress = s_class.parse_pari_progress(
            "\n".join(
                (
                    "LIMC = 1002, LIMC2 = 400879",
                    "Time factorbase (#subFB = 8) and ideal permutations: 125",
                    "#### Look for 158 relations in 243 ideals (small_norm)",
                    "#### Look for 153 relations in 137 ideals (small_norm)",
                )
            )
        )
        self.assertEqual(progress["minimum_candidate_ideal_count"], 137)
        self.assertIn("not certified", progress["interpretation"])

    def test_bnf_free_collector_has_exact_rank28_preset(self) -> None:
        source = (CAS / "run_fermigier_rank20_minkowski_specialq.py").read_text()
        self.assertIn('"--elkies-rank28"', source)
        self.assertIn("pari.addprimes(factor_hint_primes)", source)
        self.assertIn('"elkies-2026-rank28"', source)
        self.assertIn("exact_minkowski_ideal_relations_not_class_group_completion", source)
        for name in (
            "augment_bnf_free_canonical_principal_relations.py",
            "audit_bnf_free_s_class_quotient.py",
        ):
            with self.subTest(name=name):
                consumer = (CAS / name).read_text()
                self.assertIn("pari.addprimes(factor_hint_primes)", consumer)
                self.assertIn("validate_inputs()", consumer)

    def test_pinned_bnf_free_pilot_is_explicitly_uncertified(self) -> None:
        path = (
            ROOT
            / "artifacts/generated-results/elliptic-curves"
            / "elkies_2026_rank28_bnf_free_s_class_pilot_v1.json"
        )
        artifact = json.loads(path.read_text())
        self.assertEqual(artifact["classification"], "UNCERTIFIED_FACTOR_BASE")
        self.assertEqual(artifact["collector"]["curve_preset"], "elkies-2026-rank28")
        self.assertEqual(artifact["collector"]["special_ideal_mode"], "cycle-pairs")
        self.assertEqual(artifact["collector"]["sampled_generator_count"], 10460)
        self.assertEqual(artifact["collector"]["noncanonical_closed_relation_count"], 0)
        self.assertEqual(artifact["relation_rank"], 172)
        self.assertEqual(artifact["factor_base_quotient_dimension"], 141)
        self.assertIsNone(
            artifact["class_quotient_certification"][
                "remaining_dimension_upper_bound"
            ]
        )
        self.assertFalse(artifact["selmer_claim"]["expensive_search_authorized"])


if __name__ == "__main__":
    unittest.main()

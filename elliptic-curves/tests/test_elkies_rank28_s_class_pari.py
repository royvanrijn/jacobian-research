from __future__ import annotations

import json
from math import isqrt
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_rank28_s_class_pari as s_class  # noqa: E402
import merge_bnf_free_minkowski_relation_ledgers as merge_ledgers  # noqa: E402
import run_fermigier_rank20_minkowski_specialq as minkowski  # noqa: E402


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

    def test_reduced_field_model_preserves_an_explicit_generator_map(self) -> None:
        _, _, coefficients, primes = s_class.validate_inputs()
        source = s_class.worker_source(
            coefficients=coefficients,
            primes=primes,
            stack_bytes=1_000_000_000,
            mode="class-quotient",
            tech=[0.01, 4.0, 20],
            debug=1,
            field_model="polredabs",
        )
        self.assertIn("FIELD_MODEL = 'polredabs'", source)
        self.assertIn("pari.polredabs(original_polynomial, 1)", source)
        self.assertIn("original_generator_in_field_model", source)
        self.assertIn("pari.nfinit([polynomial, BAD_PRIMES])", source)
        self.assertIn('stage("polredabs",', source)

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

    def test_minkowski_collector_retains_mergeable_exact_partial_edges(self) -> None:
        source = (CAS / "run_fermigier_rank20_minkowski_specialq.py").read_text()
        self.assertIn("prime_range(2, args.trial_prime_bound + 1)", source)
        self.assertIn('"partial_relations": partial_relations', source)
        self.assertIn("primitive_projective_power_basis", source)
        self.assertIn("composite_above_bounded_factor_limit", source)
        self.assertIn("--batch-gcd-unresolved-cofactors", source)
        self.assertIn("--batch-gcd-engine", source)
        self.assertTrue((CAS / "merge_bnf_free_minkowski_relation_ledgers.py").is_file())
        self.assertTrue((CAS / "select_bnf_free_minkowski_feedback_specials.sage").is_file())

    def test_product_tree_batch_gcd_matches_pairwise_factor_splitting(self) -> None:
        values = [15, 21, 77]
        self.assertEqual(
            minkowski.product_tree_shared_divisors(values), [3, 21, 7]
        )
        product_parts, product_stats = minkowski.split_shared_cofactors(
            values, "product-tree"
        )
        pairwise_parts, _ = minkowski.split_shared_cofactors(values, "pairwise")
        self.assertEqual(product_parts, pairwise_parts)
        self.assertEqual(product_parts, [[3, 5], [3, 7], [7, 11]])
        self.assertEqual(product_stats["ambiguous_full_shared_cofactor_count"], 1)
        composite_aggregate_values = [3 * 5 * 11, 3 * 7, 5 * 13]

        def is_composite(value: int) -> bool:
            return any(value % divisor == 0 for divisor in range(2, isqrt(value) + 1))

        refined, refined_stats = minkowski.split_shared_cofactors(
            composite_aggregate_values,
            "product-tree",
            needs_refinement=is_composite,
        )
        pairwise_refined, _ = minkowski.split_shared_cofactors(
            composite_aggregate_values, "pairwise"
        )
        self.assertEqual(refined, pairwise_refined)
        self.assertEqual(refined[0], [3, 5, 11])
        self.assertEqual(refined_stats["composite_aggregate_fallback_count"], 1)
        with self.assertRaisesRegex(ValueError, "unknown batch-GCD engine"):
            minkowski.split_shared_cofactors(values, "unknown")

    def test_bounded_batch_gcd_selection_is_deterministic_and_cheap_first(self) -> None:
        records = [
            {"generator_index": 8, "cofactor": "15", "cofactor_bit_length": 4},
            {"generator_index": 3, "cofactor": "7", "cofactor_bit_length": 3},
            {"generator_index": 2, "cofactor": "5", "cofactor_bit_length": 3},
        ]
        selected = minkowski.select_batch_gcd_records(records, 2)
        self.assertEqual([index for index, _record in selected], [2, 1])
        self.assertEqual(
            minkowski.select_batch_gcd_records(records, 0), list(enumerate(records))
        )

    def test_generic_bnf_free_consumers_prove_declared_factor_hints(self) -> None:
        for name in (
            "augment_bnf_free_canonical_principal_relations.py",
            "audit_bnf_free_s_class_quotient.py",
        ):
            with self.subTest(name=name):
                source = (CAS / name).read_text()
                self.assertIn('ledger.get("selmer_rational_primes"', source)
                self.assertIn("value.is_prime()", source)
                self.assertIn("pari.addprimes(declared_primes)", source)

    @staticmethod
    def _synthetic_sparse_ledger(projective_key: list[str], row: str) -> dict:
        return {
            "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
            "curve_preset": None,
            "factor_hint_certificate": None,
            "defining_polynomial_ascending": ["1", "1", "0", "1"],
            "field_discriminant": "-31",
            "generator_coordinate_order": ["1", "theta", "theta^2"],
            "factor_base_bound": 3,
            "factor_base_completion": {
                "all_prime_ideals_above_rational_primes_through": 3,
                "materialized_complete_factor_base": True,
                "extra_declared_S_rational_primes": [],
            },
            "selmer_rational_primes": [],
            "factor_base": [
                {
                    "hnf": "P0",
                    "norm": 2,
                    "residue_degree": 1,
                    "rational_prime": 2,
                },
                {
                    "hnf": "P1",
                    "norm": 3,
                    "residue_degree": 1,
                    "rational_prime": 3,
                },
            ],
            "S_columns": [],
            "large_prime_merge_mode": "sparse-hypergraph",
            "generators": [
                {
                    "power_basis": projective_key,
                    "primitive_projective_power_basis": projective_key,
                }
            ],
            "partial_relations": [
                {
                    "generator_index": 0,
                    "fb_parity_mask_hex": row,
                    "large_prime_vertices": [[101, "Q101"], [103, "Q103"]],
                    "source": "synthetic",
                }
            ],
            "closed_relations": [],
            "unresolved_cofactors": [],
        }

    def test_cross_run_sparse_merge_finds_cycle_with_exact_provenance(self) -> None:
        left = self._synthetic_sparse_ledger(["1", "1", "0"], "0x1")
        right = self._synthetic_sparse_ledger(["1", "0", "1"], "0x2")
        merged = merge_ledgers.merge_loaded_ledgers([left, right])
        summary = merged["merged_relation_collection"]
        self.assertEqual(summary["new_cross_run_closed_relation_count"], 1)
        self.assertEqual(summary["cross_run_quotient_rank_gain"], 1)
        self.assertEqual(merged["large_prime_elimination"]["nullity"], 1)
        relation = merged["closed_relations"][0]
        self.assertEqual(relation["fb_parity_mask_hex"], "0x3")
        self.assertEqual(relation["generator_indices"], [0, 1])
        self.assertEqual(relation["source_ledger_indices"], [0, 1])

    def test_cross_run_sparse_merge_deduplicates_projective_generators(self) -> None:
        left = self._synthetic_sparse_ledger(["1", "1", "0"], "0x1")
        right = self._synthetic_sparse_ledger(["1", "1", "0"], "0x2")
        merged = merge_ledgers.merge_loaded_ledgers([left, right])
        summary = merged["merged_relation_collection"]
        self.assertEqual(summary["skipped_cross_run_projective_duplicate_count"], 1)
        self.assertEqual(summary["new_cross_run_closed_relation_count"], 0)
        self.assertEqual(merged["large_prime_elimination"]["edge_count"], 1)

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

    def test_pinned_reduced_field_model_is_exact_but_still_incomplete(self) -> None:
        path = (
            ROOT
            / "artifacts/generated-results/elliptic-curves"
            / "elkies_2026_rank28_s_class_pari_polredabs_v1.json"
        )
        artifact = json.loads(path.read_text())
        self.assertEqual(artifact["field_model"], "polredabs")
        self.assertEqual(
            artifact["status"], "INCOMPLETE_S_CLASS_COMPUTATION_SEARCH_FORBIDDEN"
        )
        events = artifact["supervisor"]["stage_events"]
        reduced = next(
            event
            for event in events
            if event.get("stage") == "polredabs" and event.get("status") == "complete"
        )
        nfinit = next(
            event
            for event in events
            if event.get("stage") == "nfinit" and event.get("status") == "complete"
        )
        self.assertEqual(
            reduced["polynomial"],
            "x^3 - 35676022072134269484503481261046298223875964999429256003*x "
            "- 81734190921553911625559669772737848345984148653181341176726216553622238508296306498",
        )
        self.assertTrue(reduced["original_generator"].startswith("Mod(-3*x + 1,"))
        ledger = json.loads(s_class.BAD_PLACE_LEDGER.read_text())
        self.assertEqual(
            int(ledger["descent_cubic_discriminant"]),
            729 * int(nfinit["polynomial_discriminant"]),
        )
        self.assertEqual(nfinit["defining_order_index"], "64023127168000")
        self.assertEqual(nfinit["signature"], "3:0")
        self.assertEqual(artifact["supervisor"]["outcome"], "strict_wall_timeout")
        self.assertEqual(artifact["supervisor"]["last_stage_event"]["stage"], "bnfinit")
        self.assertIsNone(artifact["backend_result"])
        self.assertEqual(
            artifact["supervisor"]["pari_progress"]["last_relation_request"],
            {
                "candidate_ideal_count": 153,
                "method": "rnd_rel",
                "requested_relations": 153,
            },
        )
        self.assertFalse(artifact["selmer_claim"]["expensive_search_authorized"])


if __name__ == "__main__":
    unittest.main()

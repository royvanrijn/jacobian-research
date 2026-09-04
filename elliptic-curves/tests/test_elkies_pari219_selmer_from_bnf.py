from __future__ import annotations

import argparse
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_pari219_selmer_from_bnf as worker  # noqa: E402
import run_elkies_2026_record_pari219_bnf as record_bnf  # noqa: E402


class Pari219SelmerFromBnfTests(unittest.TestCase):
    def test_record_bnf_requires_explicit_factor_certificate(self) -> None:
        self.assertEqual(record_bnf.parse_factor_primes("2,13,37"), (2, 13, 37))
        with self.assertRaises(argparse.ArgumentTypeError):
            record_bnf.parse_factor_primes("")
        certificate = (
            ROOT
            / "artifacts/generated-results"
            / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
        )
        factors = record_bnf.certified_factor_primes(certificate, 356)
        self.assertEqual(factors[:4], (2, 3, 5, 13))
        self.assertIn(28960331, factors)
        source = (CAS / "run_elkies_2026_record_pari219_bnf.py").read_text()
        self.assertIn("factor_certificate_sha256", source)
        self.assertIn("bnfcertify(b,{certify_flag})", source)
        self.assertIn("class-quotient-upper", source)
        self.assertIn("--relation-threads", source)
        self.assertIn("--random-seed", source)
        self.assertIn("setrand({args.random_seed})", source)
        self.assertIn("serial collector's early-abort strategies", source)
        self.assertIn("iferr(b=bnfinit", source)
        self.assertIn("reload_certified=1", source)
        self.assertIn('and "  ***" not in log_text', source)

    def test_record_bnf_one_sided_class_quotient_mode(self) -> None:
        self.assertEqual(
            record_bnf.mode_parameters("full-bnf"),
            (1, 0, "completed_certified_bnf"),
        )
        self.assertEqual(
            record_bnf.mode_parameters("class-quotient-upper"),
            (0, 1, "completed_certified_class_quotient_upper"),
        )
        with self.assertRaises(ValueError):
            record_bnf.mode_parameters("unsupported")
        parsed = record_bnf.parse_computed_class_group(
            "ELKIESR17RECORDPARI219BNF|stage=bnfinit|status=done|no=24|cyc=[12, 2]\n"
        )
        self.assertEqual(parsed["computed_class_group_order"], 24)
        self.assertEqual(parsed["computed_class_group_invariants"], [12, 2])
        self.assertEqual(parsed["computed_class_group_mod2_dimension"], 2)
        relation_search = record_bnf.parse_relation_search(
            "#### Look for 345 relations in 468 ideals (rnd_rel)\n"
            "#### Look for 345 relations in 340 ideals (rnd_rel)\n"
            "#### Look for 2456 relations in 2450 ideals (small_norm)\n"
        )["relation_search_by_strategy"]
        self.assertEqual(relation_search["rnd_rel"]["round_count"], 2)
        self.assertEqual(relation_search["rnd_rel"]["minimum_ideals_searched"], 340)
        self.assertEqual(
            relation_search["small_norm"]["latest"],
            {"relations_requested": 2456, "ideals_searched": 2450},
        )

    def test_complete_log_parser_retains_delete_one_ranks_and_basis(self) -> None:
        log = """\
ELKIESR17PARI219SELMER|stage=summary|selmer=29|sclass=44|norm=35|local_rank=6|bad=[2, 13]
ELKIESR17PARI219SELMER|stage=delete_one|place=-1|allowed=43|alone=34|omitted=30|rank=5
ELKIESR17PARI219SELMER|stage=delete_one|place=2|allowed=42|alone=33|omitted=29|rank=6
ELKIESR17PARI219SELMER|stage=basis|index=1|alpha=x^2 + 1
ELKIESR17PARI219SELMER|stage=norm_basis|basis=Mat([1, 0]~)
ELKIESR17PARI219SELMER|stage=local_basis|place=-1|basis=Mat([1, 0]~)
ELKIESR17PARI219SELMER|stage=local_basis|place=2|basis=Mat([0, 1]~)
"""
        summary, deleted, basis, norm_basis, local_bases = worker.parse_log(log)
        self.assertEqual(summary["two_selmer_dimension"], 29)
        self.assertEqual(summary["global_norm_square_subspace_dimension"], 35)
        self.assertEqual(
            [row["matrix_rank_after_deleting_this_place"] for row in deleted],
            [5, 6],
        )
        self.assertEqual(deleted[0]["place"], "infinity")
        self.assertEqual(basis[0]["field_squareclass_representative"], "x^2 + 1")
        self.assertEqual(norm_basis, "Mat([1, 0]~)")
        self.assertEqual([row["place"] for row in local_bases], ["infinity", "2"])

    def test_worker_reuses_the_audited_simon_definition(self) -> None:
        source = (CAS / "run_elkies_2026_pari219_selmer_from_bnf.py").read_text()
        self.assertIn("SIMON_GP_FUNCTION", source)
        self.assertIn("bnfcertify(b)", source)
        self.assertIn("stage=delete_one", source)


if __name__ == "__main__":
    unittest.main()

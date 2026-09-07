from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from overcomplete_quotient_bank import (  # noqa: E402
    build_overcomplete_quotient_bank,
    evaluate_candidate_batches_against_bank,
    evaluate_candidates_against_bank,
    load_overcomplete_quotient_bank,
    save_overcomplete_quotient_bank,
)


Q = Fraction


class OvercompleteQuotientBankTests(unittest.TestCase):
    def setUp(self) -> None:
        self.model = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        self.known = ((Q(0), Q(1)), (Q(1), Q(1)))
        self.candidate = (Q(3), Q(5))

    @staticmethod
    def fake_signature(model, points, reduction_prime, relation_prime):
        del model
        if relation_prime != 2 or reduction_prime not in (5, 7, 11):
            raise ValueError("outside the fixture")
        known_rows = {
            5: (1, 0),
            7: (0, 1),
            11: (1, 1),
        }
        candidate_coordinate = {5: 0, 7: 0, 11: 1}[reduction_prime]
        row = known_rows[reduction_prime]
        if len(points) == 3:
            row = (*row, candidate_coordinate)
        elif len(points) != 2:
            raise AssertionError("unexpected point count")
        return {
            "prime": reduction_prime,
            "group_order": {5: 7, 7: 9, 11: 13}[reduction_prime],
            "multiple_subgroup_order": 1,
            "quotient_dimension": 1,
            "rows": [list(row)],
        }

    @patch("overcomplete_quotient_bank.primes_up_to", return_value=(2, 5, 7, 11))
    @patch(
        "overcomplete_quotient_bank.finite_quotient_signature",
        side_effect=fake_signature.__func__,
    )
    def test_extra_rows_expose_a_candidate_column(self, _signature, _primes) -> None:
        bank = build_overcomplete_quotient_bank(
            self.model,
            self.known,
            moduli=(2,),
            prime_bound=11,
            row_target_per_modulus=3,
        )
        modulus_bank = bank.moduli[0]
        self.assertEqual(modulus_bank.baseline_rank, 2)
        self.assertEqual(modulus_bank.row_count, 3)
        self.assertEqual(len(modulus_bank.entries), 3)
        self.assertEqual(modulus_bank.no_torsion_witness, (5, 7))

        result = evaluate_candidates_against_bank(
            bank,
            self.known,
            (self.candidate,),
            candidate_labels=("new",),
        )
        profile = result["profile"]["profiles"][0]
        self.assertEqual(profile["baseline_rank"], 2)
        self.assertEqual(profile["combined_rank"], 3)
        self.assertEqual(profile["marginal_dimension"], 1)
        self.assertEqual(profile["individually_escaping_labels"], ["new"])

        batched = evaluate_candidate_batches_against_bank(
            bank,
            self.known,
            (self.candidate, self.candidate),
            candidate_labels=("first", "second"),
            batch_size=1,
        )
        self.assertEqual(batched["batch_count"], 2)
        self.assertEqual(batched["escaping_labels"], ["first", "second"])

    @patch("overcomplete_quotient_bank.primes_up_to", return_value=(2, 5, 7, 11))
    @patch(
        "overcomplete_quotient_bank.finite_quotient_signature",
        side_effect=fake_signature.__func__,
    )
    def test_bank_round_trip_is_bound_to_model_and_basis(self, _signature, _primes) -> None:
        bank = build_overcomplete_quotient_bank(
            self.model,
            self.known,
            moduli=(2,),
            prime_bound=11,
            row_target_per_modulus=3,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bank.json"
            save_overcomplete_quotient_bank(bank, path)
            replay = load_overcomplete_quotient_bank(
                path, model=self.model, known_points=self.known
            )
            self.assertEqual(replay, bank)
            with self.assertRaises(ValueError):
                load_overcomplete_quotient_bank(
                    path,
                    model=self.model,
                    known_points=(self.known[1], self.known[0]),
                )


if __name__ == "__main__":
    unittest.main()

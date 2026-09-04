from __future__ import annotations

from pathlib import Path
import random
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from r17_kummer_quotient_search import (  # noqa: E402
    BinaryRows,
    parse_strategies,
    projection_masks,
    select_companion_terms,
)


class R17KummerQuotientSearchTests(unittest.TestCase):
    def test_binary_rows_select_a_free_objective_column(self) -> None:
        rows = BinaryRows()
        self.assertTrue(rows.add(0b0011))
        self.assertFalse(rows.add(0b0011))
        self.assertEqual(rows.rank, 1)
        self.assertEqual(rows.free_column(3), 0)
        self.assertNotEqual(rows.reduce(1 << rows.free_column(3)), 0)
        self.assertEqual(rows.free_column(3, start=2), 2)

    def test_multi_column_target_stays_nonzero_in_the_quotient(self) -> None:
        rows = BinaryRows()
        rows.add(0b0011)
        selected = rows.free_column_combination(4, start=0, count=2)
        self.assertEqual(selected, (0, 2))
        mask = sum(1 << column for column in selected)
        self.assertNotEqual(rows.reduce(mask), 0)
        self.assertEqual(rows.free_column_combination(4, start=2, count=1), (2,))
        with self.assertRaises(ValueError):
            rows.free_column_combination(4, count=0)

    def test_strategy_cycle_realizes_single_pair_and_sparse_products(self) -> None:
        rng = random.Random(20260904)
        labels = tuple(f"P{index + 1}" for index in range(29))
        selected = [
            select_companion_terms(
                rng=rng,
                attempt=attempt,
                labels=labels,
                sparse_min=3,
                sparse_max=5,
                exponent_radius=2,
                signed_exponents=True,
            )
            for attempt in range(1, 7)
        ]
        self.assertEqual([strategy for strategy, _terms in selected], [
            "single", "pair", "sparse", "single", "pair", "sparse"
        ])
        self.assertEqual([len(terms) for _strategy, terms in selected[:2]], [1, 2])
        for strategy, terms in selected:
            self.assertEqual(len({term.index for term in terms}), len(terms))
            if strategy in ("single", "pair"):
                self.assertTrue(all(term.exponent == 1 for term in terms))
            else:
                self.assertTrue(all(0 < abs(term.exponent) <= 2 for term in terms))

    def test_exceptional_block_is_preferentially_sampled_but_not_exclusive(self) -> None:
        rng = random.Random(11)
        labels = tuple(f"P{index + 1}" for index in range(29))
        role_counts = {"generic_MW17": 0, "known_exceptional": 0}
        for attempt in range(1, 301):
            _strategy, terms = select_companion_terms(
                rng=rng,
                attempt=attempt,
                labels=labels,
                exceptional_weight=4,
            )
            for term in terms:
                role_counts[term.role] += 1
        self.assertGreater(role_counts["known_exceptional"], role_counts["generic_MW17"])
        self.assertGreater(role_counts["generic_MW17"], 0)

    def test_projection_retains_exceptional_classes_only_modulo_generic(self) -> None:
        generic, full_known = projection_masks(
            base_mask=0b101,
            exceptional_parity_mask=0b11,
            factor_base_width=3,
        )
        self.assertEqual(generic, 0b11101)
        self.assertEqual(full_known, 0b101)

    def test_strategy_parser_fails_closed(self) -> None:
        self.assertEqual(parse_strategies("single,pair"), ("single", "pair"))
        with self.assertRaises(ValueError):
            parse_strategies("single,norm-score")


if __name__ == "__main__":
    unittest.main()

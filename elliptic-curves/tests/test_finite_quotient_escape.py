#!/usr/bin/env python3
"""Exact tests for quotient-aware search promotion."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from finite_quotient_escape import (  # noqa: E402
    QuotientBlock,
    analyze_escape,
    analyze_multi_modulus_escape,
    rank_mod_prime,
)


class FiniteQuotientEscapeTests(unittest.TestCase):
    def test_exact_prime_field_rank_and_validation(self) -> None:
        rows = ((1, 2, 3), (2, 4, 1), (0, 1, 1))
        self.assertEqual(rank_mod_prime(rows, column_count=3, modulus=5), 2)
        with self.assertRaises(ValueError):
            rank_mod_prime(rows, column_count=3, modulus=6)
        with self.assertRaises(ValueError):
            rank_mod_prime(((1, 2),), column_count=3, modulus=5)

    def test_individual_and_collective_escape(self) -> None:
        block = QuotientBlock.build(
            modulus=3,
            rows=((1, 0, 1, 0), (0, 1, 1, 1)),
            column_count=4,
            source="toy-good-prime",
        )
        profile = analyze_escape(
            (block,),
            known_column_count=2,
            candidate_labels=("dependent", "escaping"),
        )
        self.assertEqual(profile.baseline_rank, 2)
        self.assertEqual(profile.combined_rank, 2)
        self.assertEqual(profile.marginal_dimension, 0)
        self.assertEqual(profile.individually_escaping_labels, ())

        escaping = QuotientBlock.build(
            modulus=3,
            rows=((1, 0, 1, 0), (0, 1, 1, 0), (0, 0, 0, 1)),
            column_count=4,
            source="second-good-prime",
        )
        profile = analyze_escape(
            (escaping,),
            known_column_count=2,
            candidate_labels=("dependent", "escaping"),
        )
        self.assertEqual(profile.marginal_dimension, 1)
        self.assertEqual(profile.individually_escaping_labels, ("escaping",))
        self.assertEqual(profile.independent_escape_basis_labels, ("escaping",))

    def test_real_R20_certificate_exposes_eight_prefix_escape_directions(self) -> None:
        artifact = json.loads(
            (
                ROOT
                / "artifacts/generated-results/"
                "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
            ).read_text()
        )
        certificate = artifact["imported_selected_twenty_basis"][
            "imported_ecsearch_cyclic_log_mod5_certificate"
        ]
        rows = tuple(tuple(record["logs"]) for record in certificate["rows"])
        labels = tuple(f"basis-column-{index}" for index in range(13, 21))
        profile = analyze_escape(
            (
                QuotientBlock.build(
                    modulus=5,
                    rows=rows,
                    column_count=20,
                    source="pinned-R20-mod5-certificate",
                ),
            ),
            known_column_count=12,
            candidate_labels=labels,
        )
        self.assertEqual(profile.baseline_rank, 12)
        self.assertEqual(profile.combined_rank, 20)
        self.assertEqual(profile.marginal_dimension, 8)
        self.assertEqual(profile.individually_escaping_labels, labels)
        self.assertEqual(profile.independent_escape_basis_labels, labels)

    def test_multi_modulus_priority_keeps_claim_boundary(self) -> None:
        labels = ("candidate",)
        profile = analyze_multi_modulus_escape(
            (
                QuotientBlock.build(
                    modulus=2,
                    rows=((1, 1),),
                    column_count=2,
                    source="mod2",
                ),
                QuotientBlock.build(
                    modulus=3,
                    rows=((1, 0), (0, 1)),
                    column_count=2,
                    source="mod3",
                ),
            ),
            known_column_count=1,
            candidate_labels=labels,
        )
        self.assertEqual(profile.moduli_with_escape, (3,))
        self.assertEqual(profile.aggregate_priority_key, (1, 1, 1))
        self.assertIn("not a calibrated rank probability", profile.to_record()["claim_boundary"])


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for the exact R20 Brumer--Kramer audit."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from audit_r20_brumer_kramer import (  # noqa: E402
    build_record,
    canonical_digest,
    monic_cubic_for_leading_scaled_root,
    weierstrass_discriminant,
)


class R20BrumerKramerTests(unittest.TestCase):
    def test_leading_scaled_monic_cubic(self) -> None:
        self.assertEqual(
            monic_cubic_for_leading_scaled_root((3, 5, 7, 4)),
            (48, 20, 7, 1),
        )
        with self.assertRaises(ValueError):
            monic_cubic_for_leading_scaled_root((1, 2, 3))

    def test_weierstrass_discriminant_replays_pinned_model(self) -> None:
        model = (
            1,
            1,
            1,
            -4437412060110743641525245114305,
            3586842216822165612930264910099076801587288127,
        )
        self.assertEqual(
            weierstrass_discriminant(model),
            int(
                "341580806733352802969707593326892076118994749446439889724290"
                "89096818399452131514624000000000"
            ),
        )

    def test_exact_cubic_field_and_brumer_kramer_terms(self) -> None:
        record = build_record()
        field = record["two_division_field"]
        brumer_kramer = record["brumer_kramer"]
        self.assertEqual(
            field["monic_defining_polynomial_coefficients_ascending"],
            [
                229557901876618599227536954246340915301586440144,
                -70998592961771898264403921828872,
                5,
                1,
            ],
        )
        self.assertEqual(field["irreducibility_witness"]["prime"], 37)
        self.assertEqual(field["irreducibility_witness"]["root_residues"], [])
        self.assertEqual(
            field["field_discriminant"],
            "17207612547621358265560224336784329653572551167050221201938192360",
        )
        self.assertEqual(field["power_order_index"], "712863540480000")
        self.assertEqual(field["signature"], {"real_places": 3, "complex_place_pairs": 0})
        self.assertEqual(
            brumer_kramer["phi_m_even_discriminant_valuation"],
            [3, 7, 13, 31, 79],
        )
        self.assertEqual(brumer_kramer["additive_primes"], [17])
        additive = brumer_kramer["additive_data"][0]
        self.assertEqual(additive["number_of_cubic_primes"], 1)
        self.assertEqual(
            additive["prime_decomposition"][0]["ramification_index"], 3
        )
        self.assertEqual(additive["prime_decomposition"][0]["residue_degree"], 1)
        self.assertEqual(brumer_kramer["u_term"], 2)
        self.assertEqual(brumer_kramer["n_term"], 5)
        self.assertEqual(brumer_kramer["class_group_two_rank_lower_bound"], 13)
        self.assertIn("unconditionally", brumer_kramer["conclusion"])

    def test_pinned_output_and_digest(self) -> None:
        record = build_record()
        declared = record.pop("result_sha256")
        self.assertEqual(canonical_digest(record), declared)
        record["result_sha256"] = declared
        pinned = ROOT / "artifacts/generated-results/elliptic_r20_brumer_kramer.json"
        if pinned.exists():
            self.assertEqual(json.loads(pinned.read_text(encoding="utf-8")), record)


if __name__ == "__main__":
    unittest.main()

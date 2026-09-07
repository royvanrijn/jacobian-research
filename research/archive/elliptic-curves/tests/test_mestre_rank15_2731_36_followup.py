#!/usr/bin/env python3
"""Focused checks for the exact T=2731/36 certificate and follow-ups."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
CERT_SCRIPT = CAS / "certify_mestre_rank15_2731_36.py"
FORMULA_SCRIPT = CAS / "explicit_formula_mestre_rank15_2731_36_delta22.py"
COVER_SCRIPT = CAS / "search_mestre_rank15_2731_36_alternate_covers.py"
CERT_ARTIFACT = GENERATED / "elliptic_mestre_rank15_2731_36.json"
FORMULA_ARTIFACT = GENERATED / "elliptic_mestre_rank15_2731_36_explicit_formula_delta22.json"
COVER_ARTIFACT = GENERATED / "elliptic_mestre_rank15_2731_36_alternate_covers.json"

EXPECTED_HASHES = {
    CERT_SCRIPT: "c872e93e8f9b951a40df0d6a4a95edbe359b5e2dab045b16d120b882b81e8cad",
    FORMULA_SCRIPT: "dff0f8f1ccc3346fd310ae997944cae3b0b8f54a8d4a6b3cc903a8f043ef2139",
    COVER_SCRIPT: "b620feba3a8c2347d3de2bc2731f3d04490b6e315609d2f590cb676dd2d119c4",
    CERT_ARTIFACT: "5f91987e9fd21887afbe0cd376e7b56844a37e0e70ade6fc713aaa3121e87c1a",
    FORMULA_ARTIFACT: "a0920e555b8a9ca32019d4d31f462969502863b41d67b00e795c138c00b56c10",
    COVER_ARTIFACT: "37dbc02c5468603edfc6f6608d7723daa68c5fb276878153cbb74c06b132bd97",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRank15T2731Over36FollowupTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = json.loads(CERT_ARTIFACT.read_text())
        cls.formula = json.loads(FORMULA_ARTIFACT.read_text())
        cls.covers = json.loads(COVER_ARTIFACT.read_text())

    def test_all_new_sources_and_artifacts_are_pinned(self) -> None:
        for path, expected in EXPECTED_HASHES.items():
            self.assertEqual(sha256(path), expected, path.name)
        self.assertEqual(
            self.certificate["provenance"]["script_sha256"],
            EXPECTED_HASHES[CERT_SCRIPT],
        )
        self.assertEqual(
            self.formula["provenance"]["script_sha256"],
            EXPECTED_HASHES[FORMULA_SCRIPT],
        )
        self.assertEqual(
            self.covers["provenance"]["script_sha256"],
            EXPECTED_HASHES[COVER_SCRIPT],
        )
        self.assertEqual(
            self.certificate["result_sha256"],
            "4d76dc4adf65376027135a1c6023c4cf9040e5d19aa2f8b84a1b79c7b8877252",
        )
        self.assertEqual(
            self.covers["result_sha256"],
            "be2aa57ccdc213062d3004568ffe3685b56c89140e5994bbba2885c8a879a4df",
        )

    def test_exact_curve_points_saturation_and_conductor(self) -> None:
        curve = self.certificate["curve"]
        self.assertEqual(curve["roots"], [0, 7, 93, 154, 161, 191])
        self.assertEqual(Fraction(curve["parameter"]), Fraction(2731, 36))
        self.assertEqual(Fraction(curve["sign_equivalent_parameter"]), Fraction(-2731, 36))
        self.assertEqual(
            curve["conductor"],
            "35180263184668233005022967592240992011410114878725318553939716560010",
        )
        self.assertEqual(curve["root_number"], -1)
        self.assertTrue(curve["exact_log_conductor_bound"]["strict_target_proved_exactly"])
        self.assertEqual(curve["exact_log_conductor_bound"]["deduced_log_conductor_upper_bound"], "3927/25")
        a1, a2, a3, coefficient_a, coefficient_b = map(
            Fraction, curve["weierstrass_coefficients"]
        )
        self.assertEqual((a1, a2, a3), (0, 0, 0))
        for collection in (self.certificate["input_points"], self.certificate["saturated_basis"]):
            self.assertEqual(len(collection), 15)
            for point in collection:
                x_value, y_value = Fraction(point["x"]), Fraction(point["y"])
                self.assertEqual(y_value**2, x_value**3 + coefficient_a * x_value + coefficient_b)

        source = self.certificate["point_source"]
        self.assertEqual(
            source["input_point_sha256"],
            "423d2e5c126501d0b0ed12f1bf8fadff6584ffffb20bc4e147e07d132ce8810e",
        )
        saturation = self.certificate["small_prime_saturation"]
        self.assertEqual(saturation["prime_bound_strict_upper_limit"], 20)
        self.assertEqual(saturation["returned_point_count"], 15)
        self.assertEqual(
            saturation["saturated_basis_sha256"],
            "6dd6cc43e9836b16af5874665dfa114e55762dbdb70284d7a67b0c17cebb903e",
        )
        self.assertEqual(Fraction(saturation["height_determinant_ratio"]), 2**28)

    def test_exact_finite_reduction_rank15_certificate(self) -> None:
        certificate = self.certificate["finite_reduction_certificate"]
        self.assertEqual(certificate["descent_modulus"], 3)
        self.assertEqual(certificate["point_count"], 15)
        self.assertEqual(certificate["combined_exact_rank_over_F3"], 15)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 15)
        self.assertEqual(
            certificate["certificate_primes"],
            [23, 31, 37, 47, 59, 73, 101, 109, 113, 131, 157, 163, 167, 179, 181],
        )
        self.assertEqual(certificate["rational_3_torsion_exclusion"]["prime"], 29)
        claim = self.certificate["claim"]
        self.assertFalse(claim["independence_uses_numerical_heights"])
        self.assertFalse(claim["independence_depends_on_ellsaturation_finite_index_hypothesis"])
        self.assertTrue(claim["does_not_claim_exact_mordell_weil_rank"])

    def test_delta22_narrowly_fails_to_close(self) -> None:
        self.assertEqual(
            self.formula["status"], "explicit_formula_does_not_close_below_17"
        )
        explicit = self.formula["explicit_formula"]
        self.assertEqual(explicit["delta"], "11/5")
        self.assertEqual(explicit["support_prime_limit"], 1_007_525)
        upper = float(explicit["conservative_explicit_formula_upper"])
        self.assertGreater(upper, 17)
        self.assertLess(upper, 17.3)
        self.assertFalse(explicit["strictly_below_17"])
        self.assertFalse(self.formula["conclusion"]["fixed_fiber_conditionally_closed"])
        self.assertEqual(
            self.formula["conclusion"]["next_action"],
            "run the declared fixed-fiber cover search",
        )

    def test_fixed_cover_budget_and_exact_negative_result(self) -> None:
        budget = self.covers["declared_budget"]
        self.assertEqual(budget["all_nonzero_certified_mod2_classes_identity_scored"], 32767)
        self.assertEqual(budget["low_weight_classes_built"], 120)
        self.assertEqual(budget["cover_classes_selected"], 20)
        self.assertEqual(budget["pilot_chart_count"], 60)
        self.assertEqual(budget["escalation_chart_count"], 8)
        self.assertEqual(budget["deep_chart_count"], 2)
        self.assertTrue(budget["one_pass_no_retry"])
        self.assertEqual(len(self.covers["cover_plans"]), 20)
        statuses = Counter((row["stage"], row["status"]) for row in self.covers["runs"])
        self.assertEqual(
            statuses,
            Counter(
                {
                    ("pilot", "completed"): 60,
                    ("escalation", "completed"): 8,
                    ("deep", "failed_or_timed_out"): 2,
                }
            ),
        )
        result = self.covers["results"]
        self.assertEqual(result["distinct_exact_affine_curve_points"], 74)
        self.assertEqual(result["nonbasis_candidate_points"], 59)
        self.assertEqual(result["exact_relations_in_certified_rank15_subgroup"], 59)
        self.assertEqual(result["unresolved_by_exact_relation_replay"], 0)
        self.assertEqual(result["combined_finite_reduction_rank"], 15)
        self.assertEqual(result["certified_new_directions"], 0)
        self.assertEqual(result["certified_rank_lower_bound_after_search"], 15)
        self.assertFalse(result["rank16_signal"])
        self.assertFalse(result["rank21_target_achieved"])


if __name__ == "__main__":
    unittest.main()

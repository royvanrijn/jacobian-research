from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json"
CAMPAIGN = ROOT / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-nagao-h1000-summary-v1.json"
SPECIALIZATIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json"
EXPECTED_HASHES = {
    CERTIFICATE: "7b94bf9b3dec377ce118c64e4455c6e2f089e88108cdabc9b3c5ca16c3a94b27",
    CAMPAIGN: "314fa5f91fd3ee86b5cbc0a25d60e81d3e1ddd3bf0eb3be0f29622b217bb9d45",
    SPECIALIZATIONS: "912d44bcbf9570d22a77300602875aaf4d1791e9e782c6f78f084f4a8f2fe562",
}
EXPECTED_FIBRES = {
    543: (12, 0, 0),
    544: (11, 0, 0),
    545: (11, 2, 2),
    531: (11, 3, 3),
    534: (11, 2, 2),
    535: (11, 0, 0),
    536: (11, 1, 1),
    537: (10, 0, 0),
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def evaluate_at_zero(record: dict) -> tuple[Fraction, Fraction]:
    numerator = [Fraction(value) for value in record["numerator_coefficients_low_to_high"]]
    denominator = [Fraction(value) for value in record["denominator_coefficients_low_to_high"]]
    value = numerator[0] / denominator[0]
    numerator_linear = numerator[1] if len(numerator) > 1 else Fraction(0)
    denominator_linear = denominator[1] if len(denominator) > 1 else Fraction(0)
    derivative = (
        numerator_linear * denominator[0] - numerator[0] * denominator_linear
    ) / denominator[0] ** 2
    return value, derivative


class R17ExtremeAnchoredMW18Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = json.loads(CERTIFICATE.read_text())
        cls.campaign = json.loads(CAMPAIGN.read_text())
        cls.specializations = json.loads(SPECIALIZATIONS.read_text())

    def test_frozen_bytes(self) -> None:
        for path, expected in EXPECTED_HASHES.items():
            self.assertEqual(digest(path), expected, str(path))

    def test_complete_refreshed_split_census(self) -> None:
        self.assertEqual(
            self.certificate["status"], "PASS_EXACT_EXTREME_ANCHORED_MW18_COVERS"
        )
        summary = self.certificate["summary"]
        self.assertEqual(summary["refreshed_fibre_count"], 8)
        self.assertEqual(summary["refreshed_split_cover_count"], 8)
        self.assertEqual(summary["refreshed_extreme_anchored_cover_count"], 8)
        self.assertEqual(summary["generic_rank_lower_bound_on_every_listed_cover"], 18)
        self.assertTrue(summary["generic_rank_is_not_proved_exact"])
        observed = {
            int(fibre["curve_id"]): (
                int(fibre["displayed_jump_over_MW17"]),
                int(fibre["split_count"]),
                int(fibre["anchored_exceptional_span_rank"]),
            )
            for chart in self.certificate["charts"]
            for fibre in chart["fibres"]
        }
        self.assertEqual(observed, EXPECTED_FIBRES)

    def test_every_refreshed_cover_is_exactly_anchored_and_normalized(self) -> None:
        for chart in self.certificate["charts"]:
            for fibre in chart["fibres"]:
                t0 = Fraction(fibre["native_parameter"])
                for cover in fibre["covers"]:
                    q = [Fraction(value) for value in cover["branch_quadratic_coefficients_low_to_high"]]
                    u0 = Fraction(cover["canonical_positive_square_root"])
                    self.assertEqual(sum(value * t0**i for i, value in enumerate(q)), u0**2)
                    self.assertTrue(cover["extreme_anchored"])
                    quotient = cover["exact_class_in_displayed_public_group_tensor_Q"]
                    self.assertTrue(quotient["nonzero_modulo_specialized_generic_MW17"])
                    self.assertTrue(any(Fraction(value) for value in quotient["exceptional_coordinates"]))
                    parameterization = cover["anchor_line_parameterization"]
                    self.assertEqual(parameterization["anchor_parameter"], "0")
                    self.assertEqual(parameterization["normalization"], "t(0)=t0 and dt/dr(0)=1")
                    value, derivative = evaluate_at_zero(parameterization["t_of_r"])
                    self.assertEqual(value, t0)
                    self.assertEqual(derivative, 1)

    def test_historical_rank28_cover_is_search_ready(self) -> None:
        historical = self.certificate["historical_rank28_anchor"]
        self.assertEqual(historical["label"], "orbit-15a68")
        self.assertTrue(historical["extreme_anchored"])
        self.assertEqual(historical["displayed_jump_over_MW17"], 11)
        value, derivative = evaluate_at_zero(
            historical["anchor_line_parameterization"]["t_of_r"]
        )
        self.assertEqual(value, Fraction(-9529, 5471))
        self.assertEqual(derivative, 1)

    def test_uniform_nagao_campaign(self) -> None:
        self.assertEqual(
            self.campaign["status"],
            "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_CAMPAIGN",
        )
        self.assertEqual(self.campaign["cover_count"], 9)
        self.assertEqual(self.campaign["total_stage_one_population_scored"], 10_950_912)
        self.assertEqual(self.campaign["total_final_survivor_count"], 178)
        search = self.campaign["uniform_search"]
        self.assertEqual((search["numerator_bound"], search["denominator_bound"]), (1000, 1000))
        priority = self.campaign["priority_order_by_top_nagao_score"]
        self.assertEqual(
            (priority[0]["anchor_id"], priority[0]["cover_label"]),
            ("curve-536", "08234-orbit-19188"),
        )
        self.assertEqual(len({(row["anchor_id"], row["cover_label"]) for row in priority}), 9)
        for row in priority:
            source = ROOT / row["source"]
            self.assertEqual(digest(source), row["source_sha256"])

    def test_every_finalist_has_exact_rank_at_least_18(self) -> None:
        audit = self.specializations
        self.assertEqual(
            audit["status"], "COMPLETE_EXACT_MW18_FINALIST_SPECIALIZATION_AUDIT"
        )
        self.assertEqual(audit["requested_finalist_count"], 178)
        self.assertEqual(audit["successful_specialization_count"], 178)
        self.assertEqual(audit["certified_rank_at_least_18_count"], 178)
        self.assertEqual(audit["structural_failure_count"], 0)
        self.assertEqual(audit["independence_unknown_count"], 0)
        self.assertTrue(any(row["r"] == "infinity" for row in audit["candidates"]))
        for row in audit["candidates"]:
            self.assertEqual(
                row["independence"]["status"],
                "CERTIFIED_INTEGRALLY_INDEPENDENT_RANK_AT_LEAST_18",
            )
            self.assertEqual(row["independence"]["selected_modulus"], 2)
            self.assertEqual(len(row["specialized_points"]["generic_R17"]), 17)
            self.assertTrue(
                row["specialized_points"]["all_section_identities_verified_exactly"]
            )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Pinned tests for the direct Nagao-root rank-16 research lane."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import certify_mestre_02557104116148_t62_35_rank16 as certificate_script  # noqa: E402
import search_mestre_02557104116148_direct_rational as direct  # noqa: E402
from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from search_mestre_root_tuple_scale import point_digest, point_on_short_curve  # noqa: E402


Q = Fraction
GENERATED = ROOT / "artifacts/generated-results"
DIRECT = GENERATED / "elliptic_mestre_02557104116148_direct_rational.json"
CERTIFICATE = GENERATED / "elliptic_mestre_02557104116148_t62_35_rank16_certificate.json"
EXPLICIT = GENERATED / "elliptic_mestre_02557104116148_t62_35_explicit_formula_delta22.json"
NEIGHBOR = GENERATED / "elliptic_mestre_02557104116148_t62_35_neighborhood.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class DirectNagaoRootProgramTests(unittest.TestCase):
    def test_even_family_formula(self) -> None:
        construction = SixRootMestreConstruction(tuple(Q(x) for x in direct.ROOTS))
        for parameter in (Q(1), Q(62, 35), Q(581, 58), Q(-17, 3)):
            self.assertEqual(
                construction.primitive_jacobian_coefficients(parameter),
                direct.family_coefficients(parameter),
            )
            self.assertEqual(
                direct.family_coefficients(parameter),
                direct.family_coefficients(-parameter),
            )

    def test_complete_direct_scan_artifact(self) -> None:
        self.assertEqual(
            sha256(DIRECT),
            "4874478c553c81ed69fffb49738b5975900a26a17d96f4dca9203a8244e75db6",
        )
        data = json.loads(DIRECT.read_text())
        self.assertEqual(
            data["result_sha256"],
            "c8d231506669c58fbb74b7e9d19b742a15564881f72b12c02d42fc9b1dadb687",
        )
        global_scan = next(
            row for row in data["modular_scan"]["strata"]
            if row["stratum"] == "global"
        )
        self.assertEqual(global_scan["primitive_population"], 18_244_819)
        self.assertEqual(global_scan["prior_specializations_excluded"], 668)
        self.assertEqual(global_scan["evaluated_population"], 18_244_151)
        self.assertEqual(
            data["prior_specialization_exclusion_audit"][
                "canonical_exclusion_lines_sha256"
            ],
            direct.EXPECTED_EXCLUSION_SHA256,
        )
        self.assertEqual(data["conductor_selection"]["selected_population"], 208)
        self.assertEqual(
            data["conductor_first_screen"],
            {
                "completed": 194,
                "errors": 0,
                "population_closed_before_any_point_or_rank_call": True,
                "selected_population": 208,
                "subtarget": 116,
                "timeouts": 14,
            },
        )
        self.assertEqual(
            data["point_search_protocol"]["maximum_stable_numerical_rank"], 16
        )
        rank16 = {
            row["parameter"]
            for row in data["selected_records"]
            if any(
                stage.get("stable_numerical_rank") == 16
                for stage in row.get("point_stages", {}).values()
            )
        }
        self.assertEqual(rank16, {"62/35", "581/58"})
        self.assertEqual(
            data["provenance"]["script_sha256"],
            "84fc344ed54f69b7a9d08e0635b08dac242014d87a1a5b5848c5b728f92ab05a",
        )
        self.assertEqual(
            data["provenance"]["scanner_sha256"],
            "3c3165b1a847106f6ef6ee0534e6fe9bf392885595aca4e7daa00aa8c7f25a54",
        )

    def test_exact_rank16_certificate(self) -> None:
        self.assertEqual(
            sha256(CERTIFICATE),
            "2c6d918546548227ac8f83287b3242e8d4261a98facd2665d506a8308f4c9fc7",
        )
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(
            data["theorem"]["certified_algebraic_rank_lower_bound"], 16
        )
        self.assertEqual(data["curve"]["parameter_T"], "62/35")
        self.assertEqual(data["curve"]["root_number"], 1)
        self.assertEqual(
            data["curve"]["conductor"], certificate_script.EXPECTED_CONDUCTOR
        )
        saturation = data["small_prime_saturation"]
        self.assertEqual(saturation["returned_point_count"], 16)
        self.assertEqual(
            saturation["saturated_basis_sha256"],
            certificate_script.EXPECTED_SATURATED_BASIS_SHA256,
        )
        points = tuple(
            (Q(row["jacobian_x"]), Q(row["jacobian_y"]))
            for row in saturation["saturated_basis"]
        )
        coefficients = direct.family_coefficients(Q(62, 35))
        self.assertEqual(point_digest(points), saturation["saturated_basis_sha256"])
        self.assertTrue(all(point_on_short_curve(coefficients, point) for point in points))
        finite = data["exact_finite_reduction_certificate"]
        self.assertEqual(finite["combined_exact_rank_over_F3"], 16)
        self.assertEqual(
            tuple(finite["certificate_primes"]),
            certificate_script.EXPECTED_CERTIFICATE_PRIMES,
        )
        self.assertEqual(finite["rational_3_torsion_exclusion"]["prime"], 23)
        self.assertEqual(
            data["provenance"]["script_sha256"],
            "d722996d7648bcc5bdb0bf0e6322d447bb296c64cb2f617f74cf5df477c7a06e",
        )

    def test_conditional_delta22_closure(self) -> None:
        self.assertEqual(
            sha256(EXPLICIT),
            "4352da34ecbc586a63042dc7c0c746be2f2d17d6ea751b8349cdac48956a56ad",
        )
        data = json.loads(EXPLICIT.read_text())
        explicit = data["explicit_formula"]
        self.assertEqual(explicit["delta"], "11/5")
        self.assertEqual(explicit["support_prime_limit"], 1_007_525)
        self.assertEqual(explicit["prime_count"], 79_057)
        self.assertEqual(explicit["prime_power_term_count"], 79_293)
        self.assertLess(
            Decimal(explicit["conservative_explicit_formula_upper"]), Decimal(18)
        )
        self.assertIn("no rank upper bound", data["conclusion"]["unconditional"])
        self.assertEqual(
            data["provenance"]["script_sha256"],
            "f6391b19c2f7070834bfa2dfa25ac41e4cd5cf76a7864c7590043a2a182d6b25",
        )

    def test_disjoint_neighborhood_artifact(self) -> None:
        self.assertEqual(
            sha256(NEIGHBOR),
            "7953c61b51b27e07840d2dbaaa6016b77d3f443c8b519fb746860f56b6b67af1",
        )
        data = json.loads(NEIGHBOR.read_text())
        audit = data["modular_scan"]["audit"]
        self.assertEqual(audit["primitive_population"], 2_918_494)
        self.assertEqual(audit["prior_excluded"], 6)
        self.assertEqual(audit["evaluated_population"], 2_918_488)
        self.assertEqual(audit["retained_union"], 8_192)
        self.assertEqual(data["conductor_selection"]["selected_population"], 160)
        self.assertEqual(data["conductor_first_screen"]["completed"], 120)
        self.assertEqual(data["conductor_first_screen"]["subtarget"], 0)
        self.assertEqual(
            data["point_search_protocol"]["maximum_stable_numerical_rank"], 12
        )
        self.assertEqual(data["target"]["hits"], [])
        self.assertEqual(
            data["provenance"]["script_sha256"],
            "30fc0d0318ed85212a7cd93e1f213db9fab4600637f07a607fac0678bf2819a5",
        )
        self.assertEqual(
            data["provenance"]["scanner_sha256"],
            "ebfd72712915799d4d2c4675cafd3fa835cb05f6c5c3b55a5e755d4ffdc8e05f",
        )


if __name__ == "__main__":
    unittest.main()

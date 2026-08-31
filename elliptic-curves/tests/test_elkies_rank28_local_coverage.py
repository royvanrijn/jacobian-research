from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
sys.path.insert(0, str(CAS))

import run_bnf_free_two_cover_local_supervisor as supervisor  # noqa: E402
import build_elkies_2026_rank28_public_selmer_controls as public_controls  # noqa: E402


class ElkiesRank28LocalCoverageTests(unittest.TestCase):
    def test_public_complement_is_an_exact_residual_selmer_positive_control(self) -> None:
        candidates = json.loads(
            (
                ARTIFACTS
                / "elkies_2026_rank28_public11_selmer_candidates_v1.json"
            ).read_text()
        )
        self.assertEqual(candidates, public_controls.build_controls())
        self.assertEqual(candidates["status"], "EXACT_MW_KUMMER_POSITIVE_CONTROLS")
        self.assertEqual(candidates["parameter"], "-9529/5471")
        self.assertEqual(candidates["candidate_count"], 11)
        self.assertEqual(candidates["residual_two_selmer_dimension_lower_bound"], 11)
        self.assertFalse(candidates["expensive_search_authorized"])
        self.assertEqual(
            candidates["independence_certificate"]["combined_rank"], 28
        )
        for candidate in candidates["candidates"]:
            with self.subTest(label=candidate["label"]):
                self.assertEqual(candidate["generator_coefficients"][1:], ["-1", "0"])
                self.assertEqual(
                    int(candidate["norm"]), int(candidate["norm_square_root"]) ** 2
                )
                self.assertEqual(
                    candidate["generator_coefficients"][0],
                    candidate["monic_cubic_point"][0],
                )
                self.assertEqual(candidate["rational_cover_witness"], ["1", "0", "0", "1"])

    def test_public_two_covers_replay_global_rational_witnesses(self) -> None:
        covers = json.loads(
            (
                ARTIFACTS
                / "elkies_2026_rank28_public11_two_cover_controls_v1.json"
            ).read_text()
        )
        audit = json.loads(
            (
                ARTIFACTS
                / "elkies_2026_rank28_public11_global_cover_witness_audit_v1.json"
            ).read_text()
        )
        self.assertEqual(len(covers["covers"]), 11)
        for cover in covers["covers"]:
            with self.subTest(label=cover["label"]):
                self.assertTrue(cover["rational_cover_witness_verified"])
                self.assertEqual(cover["rational_cover_witness"], ["1", "0", "0", "1"])
                self.assertEqual(
                    cover["rational_witness_affine_x"],
                    cover["alpha_coefficients"][0],
                )
        self.assertEqual(
            audit["status"], "GLOBAL_RATIONAL_POINT_POSITIVE_CONTROL_AUDIT"
        )
        self.assertEqual(audit["global_rational_witness_cover_count"], 11)
        self.assertEqual(len(audit["tested_rational_primes"]), 12)
        classifications = Counter(
            place["classification"]
            for cover in audit["covers"]
            for place in cover["finite_places"]
        )
        self.assertEqual(
            classifications, {"PROVED_QP_POINT_BY_GLOBAL_Q_WITNESS": 132}
        )
        self.assertTrue(
            any("no ambient class-group" in item for item in audit["claim_boundary"])
        )

    def test_rank28_signature_is_an_exact_positive_control_not_a_bound(self) -> None:
        signature = json.loads(
            (
                ARTIFACTS
                / "elkies_2026_rank28_generic17_local_signature_v1.json"
            ).read_text()
        )
        self.assertEqual(
            signature["schema"], "elliptic-curves.bnf-free-signature-map.v1"
        )
        self.assertEqual(
            signature["status"],
            "KNOWN_RANK28_BAD_PLACE_IMAGES_NOT_A_SELMER_BOUND",
        )
        self.assertEqual(signature["source"]["parameter"], "-9529/5471")
        self.assertEqual(len(signature["known_mw_images"]), 17)
        self.assertEqual(
            len(signature["public_positive_control_complement_images"]), 11
        )
        self.assertEqual(signature["local_dimension"], 53)
        self.assertEqual(signature["known_mw_local_rank"], 15)
        self.assertEqual(signature["combined_rank28_local_rank"], 15)
        self.assertEqual(signature["public_complement_bad_place_incremental_rank"], 0)
        self.assertEqual(len(signature["block_rank_comparison"]), 13)
        self.assertTrue(
            all(
                block["public_complement_incremental_rank"] == 0
                for block in signature["block_rank_comparison"]
            )
        )
        self.assertIsNone(
            signature["class_quotient_certification"][
                "remaining_dimension_upper_bound"
            ]
        )
        self.assertTrue(
            any("No class-group completeness" in item for item in signature["claim_boundary"])
        )

    def test_odd_and_real_coverage_is_partial_and_pinned(self) -> None:
        coverage = json.loads(
            (
                ARTIFACTS
                / "elkies_2026_rank28_generic17_local_coverage_v1.json"
            ).read_text()
        )
        self.assertEqual(coverage["protocol"], "BNFFREELOCALCOVERAGE-v1")
        self.assertEqual(coverage["signature_local_dimension"], 53)
        self.assertEqual(coverage["known_mw_image_count"], 17)
        self.assertEqual(
            coverage["signature_map"]["sha256"],
            "b86e7dd934def09faf9c307bd2167a223021c3b9984b99fd80880e84ab23cb19",
        )

        by_classification: dict[str, set[int]] = {}
        for place in coverage["odd_places"]:
            by_classification.setdefault(place["classification"], set()).add(
                place["rational_prime"]
            )
        self.assertEqual(
            by_classification["CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE"],
            {
                3,
                19,
                20650099,
                315574902691581877528345013999136728634663121,
            },
        )
        self.assertEqual(
            by_classification["UNRESOLVED_LOCAL_KUMMER_IMAGE_COVERAGE"],
            {
                5,
                7,
                11,
                13,
                17,
                48463,
                376018840263193489397987439236873583997122096511452343225772113000611087671413,
            },
        )
        self.assertEqual(
            coverage["real_place"]["classification"],
            "CERTIFIED_FULL_REAL_KUMMER_IMAGE_COVERAGE",
        )
        self.assertEqual(
            coverage["two_adic_place"]["classification"],
            "UNRESOLVED_TWO_ADIC_LOCAL_KUMMER_IMAGE_COVERAGE",
        )

    def test_norm_one_pilot_keeps_every_inconclusive_place_fail_closed(self) -> None:
        pilot = json.loads(
            (
                ARTIFACTS / "elkies_2026_rank28_norm_one_local_pilot12_v1.json"
            ).read_text()
        )
        self.assertEqual(pilot["protocol"], "BNFFREECOVERLOCAL-v1")
        self.assertEqual(
            pilot["status"], "SELECTED_FINITE_LOCAL_REDUCTION_AUDIT_ONLY"
        )
        self.assertEqual(pilot["covers_input"]["total_cover_count"], 49)
        self.assertEqual(pilot["covers_input"]["selected_cover_count"], 12)
        self.assertTrue(pilot["covers_input"]["selection_truncated"])
        self.assertEqual(pilot["tested_rational_primes"], [3, 5, 7, 11, 13, 17, 19])
        self.assertEqual(pilot["expected_place_count"], 84)
        self.assertEqual(pilot["completed_worker_place_count"], 79)
        self.assertEqual(pilot["incomplete_worker_place_count"], 5)
        self.assertEqual(pilot["certified_local_point_count"], 60)
        self.assertEqual(pilot["certified_local_obstruction_count"], 0)
        self.assertEqual(pilot["mathematically_inconclusive_place_count"], 24)
        classifications = Counter(
            place["classification"]
            for cover in pilot["covers"]
            for place in cover["finite_places"]
        )
        self.assertEqual(
            classifications,
            {
                "PROVED_QP_POINT_BY_SMOOTH_FP_LIFT": 60,
                "INCONCLUSIVE_SINGULAR_LIFT_STATE_CAP": 19,
                "INCONCLUSIVE_SUPERVISED_STRICT_WALL_TIMEOUT": 5,
            },
        )
        self.assertTrue(
            any("never authorizes point search" in item for item in pilot["claim_boundary"])
        )

    def test_supervisor_is_owned_resumable_and_parses_one_result(self) -> None:
        source = (CAS / "run_bnf_free_two_cover_local_supervisor.py").read_text()
        for required in (
            "start_new_session=True",
            "--retry-incomplete",
            "--overwrite-cache",
            "--timeout-per-place",
            "--rss-limit-bytes",
            "INCONCLUSIVE_SUPERVISED_",
        ):
            with self.subTest(required=required):
                self.assertIn(required, source)
        result = {"finite_places": [{"classification": "PROVED_QP_POINT"}]}
        stdout = supervisor.RESULT_PREFIX + json.dumps(result)
        self.assertEqual(supervisor.parse_worker(stdout), result)
        self.assertIsNone(supervisor.parse_worker(stdout + "\n" + stdout))


if __name__ == "__main__":
    unittest.main()

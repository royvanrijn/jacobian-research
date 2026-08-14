from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
)
from search_nagao_rank20_t5081_direction import (  # noqa: E402
    EXPECTED_INPUT_CERTIFICATE_SHA256,
    GENERIC_COMPANION_SECTIONS,
    INPUT_RANK,
    PARAMETER_T,
    CONSTRUCTION,
    exact_linear_combination,
    generic_companion_quartic_points,
    load_exact_basis,
    sha256_file,
    streaming_full_coset_frontier,
)
from search_nagao_u135_alternate_covers import (  # noqa: E402
    full_coset_identity_frontier,
)
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402


INPUT = GENERATED / "elliptic_nagao_rank20_t5081_rank20_certificate.json"
SCRIPT = CAS / "search_nagao_rank20_t5081_direction.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank20_t5081_direction.json"


class NagaoRank20T5081DirectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.coefficients, cls.basis, _ = load_exact_basis(INPUT)

    def test_pinned_input_and_companion_sections(self) -> None:
        self.assertEqual(sha256_file(INPUT), EXPECTED_INPUT_CERTIFICATE_SHA256)
        self.assertEqual(len(self.basis), INPUT_RANK)
        self.assertTrue(
            all(point_on_short_curve(self.coefficients, point) for point in self.basis)
        )
        quartic = primitive_quartic_coefficients(CONSTRUCTION, PARAMETER_T)
        companions = generic_companion_quartic_points(PARAMETER_T, quartic)
        self.assertEqual(len(companions), len(GENERIC_COMPANION_SECTIONS))
        for _, point in companions:
            self.assertEqual(point[1] ** 2, quartic_value(quartic, point[0]))
            image = quartic_point_to_short_jacobian(
                CONSTRUCTION, PARAMETER_T, point
            )
            self.assertTrue(point_on_short_curve(self.coefficients, image))

    def test_streaming_frontier_matches_reference_on_six_dimensions(self) -> None:
        retain = 20
        reference = full_coset_identity_frontier(
            self.coefficients, self.basis[:6], retain_count=retain
        )
        frontier, weights, metadata = streaming_full_coset_frontier(
            self.coefficients,
            self.basis[:6],
            retain_count=retain,
            progress_interval=0,
        )
        observed = tuple(
            (entry.score, entry.subset_indices) for entry in frontier
        )
        self.assertEqual(observed, reference)
        self.assertEqual(metadata["nonzero_classes_scored"], 63)
        self.assertEqual([entry.mask.bit_count() for entry in weights], list(range(1, 7)))

    def test_generated_artifact_replays_exactly_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the bounded direction artifact has not been generated")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["input"]["sha256"], EXPECTED_INPUT_CERTIFICATE_SHA256)
        scan = data["full_mod2_class_scan"]
        self.assertEqual(scan["nonzero_classes_scored"], 2**INPUT_RANK - 1)
        self.assertEqual(scan["weight_frontier_count"], INPUT_RANK)
        decontamination = data["generic_seed_decontamination"]
        self.assertEqual(decontamination["companion_section_count"], 6)
        self.assertTrue(
            decontamination[
                "all_specialized_images_replayed_in_certified_rank20_basis"
            ]
        )
        for record in decontamination["companions"]:
            self.assertTrue(record["fraction_group_law_replay"])
            point = (
                Q(record["jacobian_point"]["curve_x"]),
                Q(record["jacobian_point"]["curve_y"]),
            )
            self.assertEqual(
                exact_linear_combination(
                    Q(self.coefficients[3]), self.basis, record["basis_relation"]
                ),
                point,
            )
        results = data["results"]
        for record in results["candidate_points"]:
            if not record["exact_relation_in_certified_rank20_subgroup"]:
                continue
            point = Q(record["curve_x"]), Q(record["curve_y"])
            self.assertEqual(
                exact_linear_combination(
                    Q(self.coefficients[3]), self.basis, record["basis_relation"]
                ),
                point,
            )
        self.assertGreaterEqual(
            results["certified_rank_lower_bound_after_search"], INPUT_RANK
        )


if __name__ == "__main__":
    unittest.main()

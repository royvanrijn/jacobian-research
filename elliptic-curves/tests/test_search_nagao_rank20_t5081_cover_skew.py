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

from certify_nagao_rank20_t5081 import (  # noqa: E402
    CONSTRUCTION,
    PARAMETER_T,
    exact_curve_data,
)
from nagao_1994 import quartic_point_to_short_jacobian, quartic_value  # noqa: E402
from search_nagao_rank20_t5081_cover_skew import (  # noqa: E402
    EXPECTED_DIRECTION_ARTIFACT_SHA256,
    GENERIC_QUADRATIC_SECTIONS,
    SKEW_BOXES,
    generic_quadratic_quartic_points,
    load_direction_artifact,
    reconstruct_best_covers,
)
from search_nagao_rank20_t5081_direction import (  # noqa: E402
    INPUT_RANK,
    exact_linear_combination,
    load_exact_basis,
    sha256_file,
)
from triage_nagao_rank13_finalists import point_on_short_curve  # noqa: E402


CERTIFICATE = GENERATED / "elliptic_nagao_rank20_t5081_rank20_certificate.json"
DIRECTION = GENERATED / "elliptic_nagao_rank20_t5081_direction.json"
SCRIPT = CAS / "search_nagao_rank20_t5081_cover_skew.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank20_t5081_cover_skew.json"


class NagaoRank20T5081CoverSkewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.coefficients, cls.basis, _ = load_exact_basis(CERTIFICATE)
        cls.direction = load_direction_artifact(DIRECTION)

    def test_pinned_direction_and_cover_reconstruction(self) -> None:
        self.assertEqual(sha256_file(DIRECTION), EXPECTED_DIRECTION_ARTIFACT_SHA256)
        covers = reconstruct_best_covers(
            self.coefficients, self.basis, self.direction, 3
        )
        self.assertEqual(len(covers), 3)
        self.assertEqual(len(SKEW_BOXES), 2)
        self.assertEqual(SKEW_BOXES[0][3:], (1, 10))
        self.assertEqual(SKEW_BOXES[1][3:], (11, 100))

    def test_quadratic_generic_sections_specialize_exactly(self) -> None:
        quartic, _, _ = exact_curve_data()
        sections = generic_quadratic_quartic_points(PARAMETER_T, quartic)
        self.assertEqual(len(sections), len(GENERIC_QUADRATIC_SECTIONS))
        for _, point in sections:
            self.assertEqual(point[1] ** 2, quartic_value(quartic, point[0]))
            image = quartic_point_to_short_jacobian(
                CONSTRUCTION, PARAMETER_T, point
            )
            self.assertTrue(point_on_short_curve(self.coefficients, image))

    def test_generated_artifact_replays_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the cover-skew artifact has not been generated")
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(data["script_sha256"], sha256_file(SCRIPT))
        search = data["search"]
        self.assertEqual(
            len(search["runs"]),
            search["cover_count"] * search["boxes_per_cover"],
        )
        self.assertTrue(all(not record["retried"] for record in search["runs"]))
        generic = data["generic_quadratic_seed_decontamination"]
        self.assertEqual(generic["section_count"], 3)
        self.assertTrue(
            generic["all_specialized_images_replayed_in_certified_rank20_basis"]
        )
        for record in generic["sections"]:
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

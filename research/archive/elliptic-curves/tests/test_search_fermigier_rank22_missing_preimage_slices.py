from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from fermigier_mestre import FermigierMestreFamily  # noqa: E402
from search_fermigier_rank22_accidental_slices import (  # noqa: E402
    T0,
    build_slices,
    minimum_to_short,
    poly_evaluate,
    quartic_group_pullback,
)
from search_fermigier_rank22_missing_preimage_slices import (  # noqa: E402
    EXPECTED_EXTENDED78_SHA256,
    MISSING_LABELS,
    classify_hit,
    load_prior_state,
    missing_preimages,
)
from verify_fermigier_rank22_points import PUBLISHED_POINTS  # noqa: E402


class FermigierMissingPreimageSliceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.sources = missing_preimages()
        cls.slices = build_slices(cls.sources)
        cls.prior, cls.priority_sources, cls.manifest = load_prior_state(
            GENERATED
        )

    def test_exact_missing_preimages_and_group_pullbacks(self) -> None:
        self.assertEqual(tuple(label for label, _ in self.sources), MISSING_LABELS)
        self.assertEqual(len({point[0] for _, point in self.sources}), 5)
        for label, point in self.sources:
            self.assertEqual(
                point[1] ** 2,
                FermigierMestreFamily.quartic_value(T0, point[0]),
            )
            published_index = int(label[1:]) - 1
            self.assertEqual(
                quartic_group_pullback(T0, point),
                minimum_to_short(PUBLISHED_POINTS[published_index]),
            )

    def test_ten_source_slices_are_exact_nonsingular_quartics(self) -> None:
        self.assertEqual(len(self.slices), 10)
        self.assertEqual({item.slope for item in self.slices}, {-1, 1})
        for slice_data in self.slices:
            self.assertEqual(len(slice_data.coefficients) - 1, 4)
            self.assertEqual(
                poly_evaluate(slice_data.coefficients, T0),
                slice_data.source_point[1] ** 2,
            )
            self.assertEqual(slice_data.x_value(T0), slice_data.source_point[0])

    def test_pinned_prior78_and_four_survivor_fixture(self) -> None:
        self.assertEqual(len(self.prior), 78)
        self.assertEqual(len(self.priority_sources), 14)
        self.assertEqual(
            self.manifest["extended78_parameter_sha256"],
            EXPECTED_EXTENDED78_SHA256,
        )
        repeated = {
            Q("31331/104"),
            Q("22253/114"),
            Q("38633/138"),
        }
        survivors = {
            Q("100380/19589"),
            Q("184079/75"),
            Q("56441/810"),
            Q("56441/240"),
        }
        self.assertTrue(repeated <= self.prior)
        self.assertTrue(survivors.isdisjoint(self.prior))

    def test_exact_decontamination_classifications(self) -> None:
        generic_parameter = Q(2)
        generic_x = FermigierMestreFamily.known_quartic_points(
            generic_parameter
        )[0][0]
        classification, _ = classify_hit(
            generic_parameter,
            generic_x,
            prior_parameters=self.prior,
            priority_sources=self.priority_sources,
        )
        self.assertEqual(
            classification, "generic-or-T-sign-conjugate collision"
        )

        classification, _ = classify_hit(
            -Q("31331/104"),
            Q("565351/1976"),
            prior_parameters=self.prior,
            priority_sources=self.priority_sources,
        )
        self.assertEqual(classification, "previously searched Fermigier fiber")

        fixtures = (
            (Q("-100380/19589"), Q("11405580/19589")),
            (Q("-184079/75"), Q("-170677/75")),
            (Q("56441/810"), Q("1005913/4050")),
            (Q("56441/240"), Q("496637/1200")),
        )
        for raw_parameter, x_value in fixtures:
            classification, prior_slice_ids = classify_hit(
                raw_parameter,
                x_value,
                prior_parameters=self.prior,
                priority_sources=self.priority_sources,
            )
            self.assertEqual(classification, "genuinely new fiber incidence")
            self.assertEqual(prior_slice_ids, ())

    def test_pinned_bounded_artifact_frontier(self) -> None:
        path = (
            GENERATED
            / "elliptic_fermigier_rank22_missing_preimage_slices.json"
        )
        artifact = json.loads(path.read_text())
        summary = artifact["summary"]
        self.assertEqual(summary["missing_published_preimages"], 5)
        self.assertEqual(summary["genus_one_slices"], 10)
        self.assertEqual(summary["slice_searches_completed"], 10)
        self.assertEqual(summary["slice_search_timeouts"], 0)
        self.assertEqual(
            summary["new_parameters"],
            ["56441/810", "56441/240", "100380/19589", "184079/75"],
        )
        self.assertEqual(summary["genuinely_new_parameters"], 4)
        self.assertEqual(summary["conductor_calls_completed"], 2)
        self.assertEqual(summary["completed_conductors_below_strict_target"], 0)
        self.assertEqual(summary["specialized_height_50000_screens"], 0)
        self.assertEqual(summary["height_1000000_escalations"], 0)
        self.assertEqual(summary["target_hits"], [])
        self.assertEqual(summary["alternative_rank_hits"], [])
        script = CAS / "search_fermigier_rank22_missing_preimage_slices.py"
        self.assertEqual(
            hashlib.sha256(script.read_bytes()).hexdigest(),
            artifact["script_sha256"],
        )


if __name__ == "__main__":
    unittest.main()

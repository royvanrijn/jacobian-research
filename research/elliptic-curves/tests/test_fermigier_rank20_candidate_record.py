#!/usr/bin/env python3
"""Tests for the unified canonical Fermigier rank-20 candidate record."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import shutil
import sys
import unittest


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT / "cas"))
sys.path.insert(0, str(PROGRAM_ROOT))

from build_fermigier_rank20_candidate_record import (  # noqa: E402
    ANCHOR_U,
    FAMILY_ID,
    build_record,
)
from elliptic_candidate_record import (  # noqa: E402
    CANDIDATE_SCHEMA,
    WeierstrassChange,
    change_weierstrass_model,
    is_on_weierstrass_curve,
    model_from_record,
    point_from_record,
    point_sequence_sha256,
    sha256_file,
    source_point_to_target,
    stable_json_sha256,
    target_point_to_source,
    validate_candidate_identity,
    verify_finite_quotient_certificate,
)


Q = Fraction
ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)


class CommonCandidateInterfaceTests(unittest.TestCase):
    def test_exact_weierstrass_change_round_trip(self) -> None:
        source = (Q(1), Q(-3), Q(1), Q(5), Q(7))
        change = WeierstrassChange(Q(2, 3), Q(4, 5), Q(-1, 7), Q(9, 11))
        target = change_weierstrass_model(source, change)
        point = (Q(-1), Q(2))
        # The coordinate maps are inverse independently of curve membership.
        self.assertEqual(
            target_point_to_source(source_point_to_target(point, change), change),
            point,
        )
        self.assertEqual(len(target), 5)

    def test_canonical_identity_rejects_alias_as_key(self) -> None:
        record = {
            "schema": CANDIDATE_SCHEMA,
            "identity": {
                "family_id": FAMILY_ID,
                "canonical_parameter": {"name": "u", "value": "28917/20"},
                "sign_quotient": True,
                "candidate_key": f"{FAMILY_ID}:u=28917/20",
                "raw_parameter_strings_are_aliases_only": True,
            },
        }
        validate_candidate_identity(record)
        record["identity"]["candidate_key"] = f"{FAMILY_ID}:T=28917/10"
        with self.assertRaises(AssertionError):
            validate_candidate_identity(record)


@unittest.skipUnless(ARTIFACT.is_file(), "generated candidate artifact is absent")
class FermigierCandidateArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(ARTIFACT.read_text())

    def test_identity_result_digest_and_namespace_references(self) -> None:
        validate_candidate_identity(self.record)
        self.assertEqual(
            self.record["identity"]["candidate_key"],
            "fermigier-mestre-v1:u=28917/20",
        )
        self.assertEqual(Q(self.record["identity"]["canonical_parameter"]["value"]), ANCHOR_U)
        self.assertEqual(
            self.record["identity"]["aliases"],
            [
                {
                    "coordinate": "adapter_u",
                    "maps_to_canonical_u": "28917/20",
                    "role": "canonical",
                    "value": "28917/20",
                },
                {
                    "coordinate": "adapter_u",
                    "maps_to_canonical_u": "28917/20",
                    "role": "sign-quotient alias",
                    "value": "-28917/20",
                },
                {
                    "coordinate": "literal_symmetric_shift_s",
                    "map": "u=abs(s)/2",
                    "maps_to_canonical_u": "28917/20",
                    "role": "legacy literal-shift alias only",
                    "value": "28917/10",
                },
                {
                    "coordinate": "literal_symmetric_shift_s",
                    "map": "u=abs(s)/2",
                    "maps_to_canonical_u": "28917/20",
                    "role": "legacy sign alias only",
                    "value": "-28917/10",
                },
            ],
        )
        payload = deepcopy(self.record)
        claimed = payload.pop("result_sha256")
        payload.pop("generated_at_utc")
        self.assertEqual(stable_json_sha256(payload), claimed)
        bridge = self.record["artifact_namespace_bridge"]
        self.assertTrue(bridge["imported_ecsearch_namespace"])
        self.assertTrue(bridge["legacy_generated_results_namespace"])
        self.assertEqual(bridge["files_moved_or_rewritten"], [])
        self.assertTrue(
            {
                "elliptic_fermigier_bidegree21_all80.json",
                "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json",
                "elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json",
            }.issubset(
                {
                    Path(source["path"]).name
                    for source in bridge["legacy_generated_results_namespace"]
                }
            )
        )
        for namespace in (
            bridge["imported_ecsearch_namespace"],
            bridge["legacy_generated_results_namespace"],
        ):
            for source in namespace:
                self.assertEqual(
                    sha256_file(REPOSITORY_ROOT / source["path"]), source["sha256"]
                )

    def test_all_models_transforms_and_complete_pool(self) -> None:
        models = self.record["models"]
        canonical = model_from_record(
            models["canonical_generalized"]["coefficients_a1_a2_a3_a4_a6"]
        )
        legacy = model_from_record(
            models["legacy_normalized_short_jacobian"]["coefficients"]
        )
        raw = model_from_record(models["binary_quartic"]["raw_jacobian_coefficients"])
        minimum = model_from_record(models["global_minimal"]["coefficients"])
        changes = self.record["exact_transformations"]
        canonical_to_legacy = WeierstrassChange.from_values(
            changes["canonical_to_legacy_normalized_short"]["change_u_r_s_t"]
        )
        canonical_to_raw = WeierstrassChange.from_values(
            changes["canonical_to_raw_binary_quartic_jacobian"]["change_u_r_s_t"]
        )
        canonical_to_minimum = WeierstrassChange.from_values(
            changes["canonical_to_global_minimal"]["change_u_r_s_t"]
        )
        self.assertEqual(change_weierstrass_model(canonical, canonical_to_legacy), legacy)
        self.assertEqual(change_weierstrass_model(canonical, canonical_to_raw), raw)
        self.assertEqual(change_weierstrass_model(canonical, canonical_to_minimum), minimum)

        pool = self.record["complete_point_pool"]
        self.assertEqual(pool["abscissa_count"], len(pool["abscissas"]))
        self.assertEqual(pool["abscissa_count"], 58)
        self.assertEqual(pool["deduplicated_difference_count"], 115)
        canonical_points = tuple(
            point_from_record(entry["canonical_point"])
            for entry in pool["difference_pool"]
        )
        self.assertTrue(all(is_on_weierstrass_curve(canonical, point) for point in canonical_points))
        transported = {
            "canonical_generalized": canonical_points,
            "legacy_normalized_short": tuple(
                source_point_to_target(point, canonical_to_legacy)
                for point in canonical_points
            ),
            "raw_binary_quartic_jacobian": tuple(
                source_point_to_target(point, canonical_to_raw)
                for point in canonical_points
            ),
            "global_minimal": tuple(
                source_point_to_target(point, canonical_to_minimum)
                for point in canonical_points
            ),
        }
        model_lookup = {
            "canonical_generalized": canonical,
            "legacy_normalized_short": legacy,
            "raw_binary_quartic_jacobian": raw,
            "global_minimal": minimum,
        }
        for label, points in transported.items():
            replay = pool["transport_replay"][label]
            self.assertEqual(point_sequence_sha256(points), replay["point_sequence_sha256"])
            self.assertTrue(all(is_on_weierstrass_curve(model_lookup[label], point) for point in points))

    def test_selected_and_bounded_saturation_cross_certificates(self) -> None:
        legacy = model_from_record(
            self.record["models"]["legacy_normalized_short_jacobian"]["coefficients"]
        )
        selected = tuple(
            point_from_record(entry["points"]["legacy_normalized_short"])
            for entry in self.record["imported_selected_twenty_basis"]["basis"]
        )
        saturated = tuple(
            point_from_record(point)
            for point in self.record["bounded_saturation_status"][
                "returned_legacy_basis"
            ]
        )
        certificates = self.record["independent_cas_cross_certificates"]
        for key in ("mod_2", "mod_3", "mod_5"):
            verify_finite_quotient_certificate(
                legacy, selected, certificates["original_imported_basis"][key]
            )
            verify_finite_quotient_certificate(
                legacy,
                saturated,
                certificates["bounded_saturation_candidate_basis"][key],
            )
        original_ranks = {
            key: certificates["original_imported_basis"][key][
                "combined_rank_over_relation_field"
            ]
            for key in ("mod_2", "mod_3", "mod_5")
        }
        self.assertEqual(original_ranks, {"mod_2": 0, "mod_3": 19, "mod_5": 20})
        self.assertTrue(
            all(
                certificates["bounded_saturation_candidate_basis"][key][
                    "certified_independent"
                ]
                for key in ("mod_2", "mod_3", "mod_5")
            )
        )
        self.assertEqual(
            self.record["bounded_saturation_status"]["status"],
            "not-globally-proved-saturated",
        )

    def test_ledger_has_global_promotion_closure_and_both_exclusions(self) -> None:
        entries = {
            Path(entry["path"]).name: entry
            for entry in self.record["promotion_and_rejection_ledger"]["entries"]
        }
        self.assertEqual(
            entries["elliptic_fermigier_global.json"]["decision"],
            "enumerated-but-not-retained-in-discovery-union",
        )
        self.assertEqual(
            entries["fermigier_rank20_near_miss_v1.json"]["decision"],
            "promoted-to-exact-rank-at-least-20-anchor",
        )
        self.assertEqual(
            entries[
                "elliptic_fermigier_rank20_28917_20_explicit_formula_delta22.json"
            ]["decision"],
            "conditional-fixed-fiber-closure-not-identity-rejection",
        )
        self.assertEqual(
            entries[
                "elliptic_fermigier_rank20_adapter_neighborhood_audit.json"
            ]["decision"],
            "excluded-before-selection-as-prior-anchor",
        )
        self.assertEqual(
            entries["elliptic_fermigier_high_power_crt_gauss.json"]["decision"],
            "not-searched-as-a-fresh-CRT-Gauss-candidate",
        )
        exceptional = entries["elliptic_fermigier_exceptional_transport.json"]
        self.assertEqual(
            exceptional["decision"], "no-new-section-or-specialization-found"
        )
        self.assertEqual(exceptional["evidence"]["projective_mobius_pencils"], 88)
        self.assertEqual(exceptional["evidence"]["true_fiber_products"], 3160)
        quotient_ball = entries[
            "elliptic_fermigier_exceptional_quotient_ball.json"
        ]
        self.assertEqual(
            quotient_ball["decision"],
            "no-new-low-genus-transport-in-support-at-most-2-ball",
        )
        self.assertEqual(quotient_ball["evidence"]["maximum_support_weight"], 2)
        self.assertEqual(quotient_ball["evidence"]["signed_rank20_directions"], 128)
        self.assertEqual(quotient_ball["evidence"]["signed_E22_directions"], 200)
        pilot = entries[
            "elliptic_fermigier_bidegree21_p13_r20e1_pilot.json"
        ]
        self.assertEqual(
            pilot["decision"],
            "no-new-section-or-specialization-in-single-pair-pilot",
        )
        self.assertEqual(pilot["evidence"]["completed_pairs"], ["P13xR20E1"])
        self.assertEqual(pilot["evidence"]["completed_pair_count"], 1)
        self.assertEqual(pilot["evidence"]["possible_independent_pair_count"], 80)
        self.assertFalse(pilot["evidence"]["all_pairs_classified"])

        all80 = entries["elliptic_fermigier_bidegree21_all80.json"]
        self.assertEqual(
            all80["decision"],
            "no-new-section-or-specialization-in-all-80-finite-chart-classification",
        )
        self.assertTrue(
            all80["evidence"]["all_80_independent_pairs_classified"]
        )
        self.assertEqual(all80["evidence"]["pair_count"], 80)
        self.assertEqual(all80["evidence"]["irreducible_degree_32_pair_count"], 80)
        self.assertEqual(all80["evidence"]["unresolved_pair_count"], 0)

        simultaneous = entries[
            "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json"
        ]
        self.assertEqual(
            simultaneous["decision"],
            "no-third-parameter-in-genuine-simultaneous-square-H200000-box",
        )
        self.assertEqual(
            simultaneous["evidence"]["projective_height_bound"], 200000
        )
        self.assertEqual(simultaneous["evidence"]["direction_count"], 80)
        self.assertEqual(
            simultaneous["evidence"]["fiber_product_pair_count"], 3160
        )
        self.assertFalse(
            simultaneous["evidence"]["product_square_surrogate_used"]
        )
        self.assertEqual(simultaneous["evidence"]["new_third_parameter_count"], 0)
        self.assertFalse(simultaneous["evidence"]["outside_box_absence_proved"])

        nonlinear = entries[
            "elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json"
        ]
        self.assertEqual(
            nonlinear["decision"],
            "no-rational-point-on-single-degree32-component-in-H1024-box",
        )
        self.assertEqual(nonlinear["evidence"]["projective_height_bound"], 1024)
        self.assertTrue(nonlinear["evidence"]["affine_box_complete"])
        self.assertEqual(nonlinear["evidence"]["exact_hits"], 0)
        self.assertFalse(
            nonlinear["evidence"][
                "all_rational_points_on_degree32_component_classified"
            ]
        )
        self.assertFalse(
            nonlinear["evidence"]["other_exceptional_pairs_classified"]
        )

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is unavailable")
    def test_generator_replays_the_pinned_mathematical_record(self) -> None:
        replay = build_record()
        pinned = deepcopy(self.record)
        pinned.pop("generated_at_utc")
        self.assertEqual(replay, pinned)


if __name__ == "__main__":
    unittest.main()

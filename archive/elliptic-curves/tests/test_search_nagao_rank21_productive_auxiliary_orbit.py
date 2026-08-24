#!/usr/bin/env python3
"""Exact regression tests for the published-rank-21 auxiliary orbit."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from nagao_1994 import RANK21_PUBLISHED_POINTS  # noqa: E402
from search_nagao_rank21_productive_auxiliary_orbit import (  # noqa: E402
    PRODUCTIVE_INTERCEPT,
    PRODUCTIVE_SLOPE,
    PRODUCTIVE_SOURCE_X,
    Q,
    Slice,
    PointedQuartic,
    build_slices,
    minimum_to_short,
    select_reconstruction_convention,
    short_to_minimum,
    slice_polynomial,
    validate_bivariate_quartic,
)
from search_nagao_rank21_accidental_slices import normalize_slice  # noqa: E402


SCRIPT = CAS / "search_nagao_rank21_productive_auxiliary_orbit.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_rank21_productive_auxiliary_orbit.json"
)


class ProductiveAuxiliaryOrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_exact_model_and_bivariate_replays(self) -> None:
        validate_bivariate_quartic()
        self.assertEqual(
            tuple(short_to_minimum(minimum_to_short(point)) for point in RANK21_PUBLISHED_POINTS),
            RANK21_PUBLISHED_POINTS,
        )
        polynomial = slice_polynomial(PRODUCTIVE_SLOPE, PRODUCTIVE_INTERCEPT)
        self.assertEqual(len(polynomial) - 1, 4)
        normalized = normalize_slice(polynomial)
        self.assertEqual(normalized.genus, 1)
        self.assertEqual(
            normalized.normalized_coefficients,
            (
                4239687950187868058329,
                -230514202948783942032,
                984650398674293128,
                79837193272671744,
                -63470653224944,
            ),
        )

    def test_unique_published_preimage_convention(self) -> None:
        reconstruction, trials = select_reconstruction_convention()
        visible_matches = max(record["generic_abscissa_matches"] for record in trials)
        self.assertEqual(visible_matches, 10)
        accidental = tuple(
            point
            for point in reconstruction
            if point[0]
            not in {
                Fraction(root) + sign * Fraction(14721, 188)
                for root in (0, 4, 47, 352, 380, 399)
                for sign in (-1, 1)
            }
        )
        self.assertEqual(len(accidental), 11)
        self.assertIn(PRODUCTIVE_SOURCE_X, {point[0] for point in accidental})

    def test_pointed_quartic_round_trip(self) -> None:
        source_y = Q(1483060942564315, 170720356)
        normalized = normalize_slice(slice_polynomial(-1, PRODUCTIVE_INTERCEPT))
        item = Slice(5, (PRODUCTIVE_SOURCE_X, source_y), -1, PRODUCTIVE_INTERCEPT, normalized)
        auxiliary = PointedQuartic.from_slice(item)
        normalized_point = Q(950, 139), Q(52302717723)
        image = auxiliary.forward(normalized_point)
        self.assertEqual(auxiliary.inverse(image), normalized_point)

    def test_artifact_frontier_and_hash(self) -> None:
        data = self.data
        self.assertFalse(data["target_hit"])
        self.assertEqual(data["published_point_reconstruction"]["accidental_preimage_count"], 11)
        self.assertEqual(data["slice_search"]["slice_count"], 22)
        self.assertEqual(data["slice_search"]["completed_slice_count"], 22)
        self.assertEqual(data["slice_search"]["distinct_forced_nonvisible_parameter_count"], 5)
        self.assertEqual(
            data["productive_auxiliary_curve"]["height_selection"]["stable_numerical_rank"],
            8,
        )
        orbit = data["orbit"]
        self.assertEqual(orbit["enumerated_nonzero_coefficient_vector_count"], 116448)
        self.assertEqual(orbit["distinct_nonzero_canonical_parameter_count"], 105938)
        self.assertEqual(orbit["nonvisible_proxy_below_gate_parameters"], ["12122/417"])
        self.assertEqual(
            data["reproduction"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )

    def test_every_conductor_misses_strict_target(self) -> None:
        conductors = self.data["exact_conductors"]
        self.assertEqual(
            set(conductors),
            {"4745/248", "12122/417", "36449/16393", "50221/2085", "112537/2180"},
        )
        self.assertTrue(all(record["status"] == "completed" for record in conductors.values()))
        self.assertTrue(
            all(not record["below_strict_log_conductor_target"] for record in conductors.values())
        )
        self.assertEqual(
            conductors["12122/417"]["log_conductor"],
            "184.919703675014813309434588107747671757372855900775692755238",
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import unittest


CAS = Path(__file__).resolve().parents[1] / "cas"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CAS))

from search_fermigier_published_pair_fiber_products import (  # noqa: E402
    exact_slices,
    published_accidentals,
)
from search_fermigier_rank22_auxiliary_orbits import (  # noqa: E402
    T0,
    generic_intersections,
    make_auxiliary_slice,
    orbit_vectors,
    prior_parameter_manifest,
    radical_upper_proxy,
)


SCRIPT = CAS / "search_fermigier_rank22_auxiliary_orbits.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_rank22_auxiliary_orbits.json"
)
L1_SEVEN_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_rank22_auxiliary_orbits_l1_7.json"
)
PRIMARY = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_fermigier_rank22_accidental_slices.json"
)
EXPECTED_SCRIPT_SHA256 = (
    "d8b2e69b409cb0e5097324c804b1c065f3a0ec28f1cf2d3e7abd4e4d9ef4c4e5"
)
EXPECTED_ARTIFACT_SHA256 = (
    "ef840dce82dcc7b70bd0aa72995aa9c16a562a8b84cafad68cad59a55589d4b1"
)
EXPECTED_L1_SEVEN_ARTIFACT_SHA256 = (
    "38e128b3b08e5a4929f536a6bf3fa7b72b5990dd3896460f59a91b2a2bea6112"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FermigierAuxiliaryOrbitTests(unittest.TestCase):
    def test_sign_quotient_vector_counts(self) -> None:
        self.assertEqual(
            len(
                orbit_vectors(
                    6,
                    maximum_absolute_coefficient=3,
                    maximum_l1_norm=5,
                )
            ),
            1754,
        )
        self.assertEqual(
            len(
                orbit_vectors(
                    7,
                    maximum_absolute_coefficient=3,
                    maximum_l1_norm=5,
                )
            ),
            3493,
        )
        self.assertEqual(
            len(
                orbit_vectors(
                    6,
                    maximum_absolute_coefficient=3,
                    maximum_l1_norm=7,
                )
            ),
            8088,
        )
        self.assertEqual(
            len(
                orbit_vectors(
                    7,
                    maximum_absolute_coefficient=3,
                    maximum_l1_norm=7,
                )
            ),
            20987,
        )

    def test_pointed_quartic_map_round_trip(self) -> None:
        primary = json.loads(PRIMARY.read_text())
        source_slice = exact_slices(published_accidentals(primary))[0]
        auxiliary = make_auxiliary_slice(source_slice)
        intersections = generic_intersections(auxiliary)
        self.assertEqual(len(intersections), 7)
        parameter = next(value for value, _ in intersections if value != T0)
        ordinate_square = auxiliary.quartic_value(parameter)
        numerator = isqrt(ordinate_square.numerator)
        denominator = isqrt(ordinate_square.denominator)
        ordinate = Fraction(numerator, denominator)
        self.assertEqual(ordinate**2, ordinate_square)
        image = auxiliary.forward((parameter, ordinate))
        self.assertEqual(auxiliary.inverse(image), (parameter, ordinate))

    def test_prior_manifest_and_record_proxy(self) -> None:
        prior, record = prior_parameter_manifest(
            ROOT / "artifacts" / "generated-results"
        )
        self.assertEqual(len(prior), 590)
        self.assertEqual(record["parameter_sha256"], (
            "64c09a13b427938a44251a91f74a116f7f9e685aed07c6159550e7ec3ea51291"
        ))
        proxy = radical_upper_proxy(Fraction(39508, 39), prime_bound=97)
        self.assertAlmostEqual(
            proxy["log_radical_upper_proxy"], 179.8916976065812, places=12
        )

    def test_generated_artifact_is_pinned(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        artifact = json.loads(ARTIFACT.read_text())
        self.assertEqual(artifact["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertFalse(artifact["target"]["target_hit"])
        self.assertEqual(artifact["outcome"]["unique_orbit_parameters"], 71439)
        self.assertEqual(artifact["outcome"]["multiple_source_parameters"], 0)
        self.assertEqual(artifact["outcome"]["subtarget_conductors"], 0)
        self.assertEqual(
            artifact["proxy_selection"]["summary"][
                "minimum_log_radical_upper_proxy"
            ],
            193.88691890362722,
        )
        conductor_rows = artifact["conductor_and_rank_tranche"]
        self.assertEqual(len(conductor_rows), 1)
        self.assertEqual(conductor_rows[0]["parameter_t"], "330449/190")
        self.assertEqual(
            conductor_rows[0]["conductor_probe"]["log_conductor"],
            "189.195571021398057056294598360395604148857131780223475108370",
        )

    def test_l1_seven_extension_is_pinned(self) -> None:
        self.assertEqual(
            sha256(L1_SEVEN_ARTIFACT), EXPECTED_L1_SEVEN_ARTIFACT_SHA256
        )
        artifact = json.loads(L1_SEVEN_ARTIFACT.read_text())
        self.assertEqual(artifact["script_sha256"], EXPECTED_SCRIPT_SHA256)
        self.assertEqual(artifact["parameters"]["maximum_l1_norm"], 7)
        self.assertFalse(artifact["target"]["target_hit"])
        self.assertEqual(artifact["outcome"]["unique_orbit_parameters"], 419023)
        self.assertEqual(artifact["outcome"]["multiple_source_parameters"], 0)
        self.assertEqual(artifact["outcome"]["subtarget_conductors"], 0)
        self.assertEqual(
            artifact["proxy_selection"]["summary"][
                "minimum_log_radical_upper_proxy"
            ],
            193.88691890362722,
        )
        self.assertIn("--max-l1-norm 7", artifact["reproducing_command"])


if __name__ == "__main__":
    unittest.main()

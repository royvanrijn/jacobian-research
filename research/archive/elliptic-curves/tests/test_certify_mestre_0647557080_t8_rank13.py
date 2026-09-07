#!/usr/bin/env python3
"""Focused exact checks for the Mestre T=8 rank-13 certificate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    mod2_reduction_signature,
)
from search_mestre_root_tuple_scale import (  # noqa: E402
    point_digest,
    point_on_short_curve,
    quartic_point_to_jacobian,
    quartic_value,
)


SCRIPT = CAS / "certify_mestre_0647557080_t8_rank13.py"
ARTIFACT = GENERATED / "elliptic_mestre_0647557080_t8_rank13_certificate.json"
INPUT = GENERATED / "elliptic_mestre_root_tuple_scale_max100.json"
EXPECTED_SCRIPT_SHA256 = (
    "b4eb1e0c911310bf7d97e5835532c0c83d5b7ff363f526488a331fa8c4734747"
)
EXPECTED_ARTIFACT_SHA256 = (
    "cef8f6f0279b48a114459be7c36f9c9c96cd2416d45be840a8d74d58518e0250"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Mestre0647557080T8Rank13CertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text())
        cls.construction = SixRootMestreConstruction(
            tuple(Fraction(value) for value in (0, 6, 47, 55, 70, 80))
        )
        cls.parameter = Fraction(8)
        cls.coefficients = cls.construction.primitive_jacobian_coefficients(
            cls.parameter
        )
        cls.basis = tuple(
            (Fraction(point["jacobian_x"]), Fraction(point["jacobian_y"]))
            for point in cls.data["small_prime_saturation"]["saturated_basis"]
        )

    def test_pinned_files_and_input(self) -> None:
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            self.data["provenance"]["script_sha256"], EXPECTED_SCRIPT_SHA256
        )
        self.assertEqual(
            self.data["input"]["max100_artifact_sha256"], sha256(INPUT)
        )
        self.assertEqual(self.data["input"]["roots"], [0, 6, 47, 55, 70, 80])
        self.assertEqual(self.data["input"]["parameter"], 8)
        self.assertEqual(
            self.data["input"]["conductor_phase"]["log_conductor"],
            "82.3515440580102287639528785405710145305597823305329922776592",
        )

    def test_visible_and_accidental_rank_split_is_explicit(self) -> None:
        reconstruction = self.data["H5000_reconstruction"]
        self.assertEqual(reconstruction["displayed_point_count"], 12)
        self.assertEqual(reconstruction["signless_quartic_point_count"], 61)
        self.assertEqual(reconstruction["pool_point_count_modulo_inverse"], 61)
        self.assertEqual(
            reconstruction["pool_point_sha256"],
            "5cbb40eac95c97f6da9fcd6e6c7be57da78fd9b5c2f9b7daed29924f05596681",
        )
        self.assertEqual(
            [
                record["stable_numerical_rank"]
                for record in reconstruction["prefix_rank_replay"]
            ],
            [10, 11, 12, 13, 13],
        )
        self.assertEqual(reconstruction["displayed_points_span_stable_numerical_rank"], 10)
        expected_x = [Fraction(75, 2), Fraction(175, 37), Fraction(243, 4)]
        accidental = reconstruction["first_three_accidental_directions"]
        self.assertEqual(
            [Fraction(record["quartic_point"]["x"]) for record in accidental],
            expected_x,
        )
        quartic_coefficients = self.construction.primitive_quartic_coefficients(
            self.parameter
        )
        for record in accidental:
            quartic_point = (
                Fraction(record["quartic_point"]["x"]),
                Fraction(record["quartic_point"]["y"]),
            )
            self.assertEqual(
                quartic_point[1] ** 2,
                quartic_value(quartic_coefficients, quartic_point[0]),
            )
            jacobian = (
                Fraction(record["jacobian_point"]["x"]),
                Fraction(record["jacobian_point"]["y"]),
            )
            self.assertEqual(
                quartic_point_to_jacobian(
                    self.construction, self.parameter, quartic_point
                ),
                jacobian,
            )

    def test_saturated_basis_has_exact_finite_reduction_rank_13(self) -> None:
        self.assertEqual(len(self.basis), 13)
        self.assertTrue(
            all(point_on_short_curve(self.coefficients, point) for point in self.basis)
        )
        self.assertEqual(
            point_digest(self.basis),
            "1a038326d2caff9bac0310cc484f21921270666c20f9bfdbb2d48cf8abd7975f",
        )
        certificate = self.data["exact_finite_reduction_certificate"]
        primes = [5, 13, 19, 23, 29, 31, 41, 47, 67, 71, 73]
        self.assertEqual(certificate["certificate_primes"], primes)
        signatures = tuple(
            mod2_reduction_signature(self.coefficients, self.basis, prime)
            for prime in primes
        )
        self.assertEqual(combined_mod2_rank(signatures, len(self.basis)), 13)
        self.assertEqual(
            [list(map(list, signature.rows)) for signature in signatures],
            [record["rows"] for record in certificate["signatures"]],
        )
        self.assertEqual(certificate["two_torsion_certificate_prime"], 11)
        self.assertEqual(certificate["certified_algebraic_rank_lower_bound"], 13)
        self.assertFalse(
            self.data["conclusion"][
                "depends_on_PARIs_finite_index_saturation_hypothesis"
            ]
        )

    def test_T_sign_quotient_replays_exactly(self) -> None:
        symmetry = self.data["symmetry"]
        self.assertEqual(symmetry["exact_quotient"], "T is identified with -T")
        for value in symmetry[
            "quartic_and_jacobian_coefficients_replayed_equal_on_panel"
        ]:
            parameter = Fraction(value)
            self.assertEqual(
                self.construction.primitive_quartic_coefficients(parameter),
                self.construction.primitive_quartic_coefficients(-parameter),
            )
            self.assertEqual(
                self.construction.primitive_jacobian_coefficients(parameter),
                self.construction.primitive_jacobian_coefficients(-parameter),
            )
        self.assertEqual(self.data["provenance"]["owned_processes_remaining"], 0)


if __name__ == "__main__":
    unittest.main()

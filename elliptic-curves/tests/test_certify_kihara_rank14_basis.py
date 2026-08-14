from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS_DIRECTORY = ROOT / "elliptic-curves" / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from alternate_quartic_covers import point_on_short_curve, short_add  # noqa: E402
from certify_kihara_rank14_basis import (  # noqa: E402
    SATURATED_SUBSET_INDICES,
    covariant_images,
    signature_from_record,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)


ARTIFACT = ROOT / "artifacts" / "generated-results" / "elliptic_kihara_rank14_basis.json"


def point(record: list[str]) -> tuple[Q, Q]:
    return Q(record[0]), Q(record[1])


class KiharaRank14BasisCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        cls.coefficients = tuple(Q(value) for value in cls.data["short_jacobian_coefficients"])

    def test_covariant_images_and_differences_recompute_exactly(self) -> None:
        coefficients, images, differences = covariant_images()
        self.assertEqual(coefficients, self.coefficients)
        self.assertEqual(tuple(point(row) for row in self.data["covariant_images"]), images)
        self.assertEqual(tuple(point(row) for row in self.data["difference_images"]), differences)
        self.assertEqual(len(set(images)), 15)

    def test_every_stored_point_is_exact_and_saturation_relations_hold(self) -> None:
        differences = tuple(point(row) for row in self.data["difference_images"])
        halves = tuple(point(row) for row in self.data["first_halves"])
        basis = tuple(point(row) for row in self.data["saturated_basis"])
        for stored_point in (*differences, *halves, *basis):
            self.assertTrue(point_on_short_curve(self.coefficients, stored_point))
        for difference, half in zip(differences, halves):
            self.assertEqual(short_add(self.coefficients, half, half), difference)
        subset_sum = None
        for index in SATURATED_SUBSET_INDICES:
            subset_sum = short_add(self.coefficients, subset_sum, halves[index])
        self.assertEqual(short_add(self.coefficients, basis[0], basis[0]), subset_sum)
        self.assertEqual(basis[1:], halves[1:])

    def test_finite_reductions_certify_rank_fourteen(self) -> None:
        signatures = tuple(
            signature_from_record(record)
            for record in self.data["mod2_reduction_signatures"]
        )
        self.assertEqual(combined_mod2_rank(signatures, 14), 14)
        self.assertEqual(self.data["combined_mod2_rank"], 14)
        prime = int(self.data["no_rational_2_torsion_prime"])
        self.assertTrue(
            short_curve_has_no_rational_2_torsion_modular_certificate(
                self.coefficients, prime
            )
        )
        self.assertEqual(self.data["exact_specialized_rank_lower_bound"], 14)


if __name__ == "__main__":
    unittest.main()


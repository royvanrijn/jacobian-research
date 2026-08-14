from __future__ import annotations

from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from derive_kihara_rank14_identities import derive  # noqa: E402


class KiharaRank14IdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = derive()

    def test_center_obstruction_and_all_extra_squares(self) -> None:
        self.assertTrue(self.result["all_identities_exact"])
        center = self.result["center_geometry"]
        self.assertEqual(center["mean_zero_elementary_relation"], "2*e5=e2*e3")
        self.assertTrue(center["kihara_relation_exact"])
        self.assertTrue(
            all(
                record["identity_exact"]
                for record in self.result["extra_sections"].values()
            )
        )

    def test_large_specialized_square_factors_are_pinned(self) -> None:
        expected = {
            "P14": (
                58,
                "b6b5bb73a584cacddcbe1e45d0fb7e839f1489ffaffa3f2db14e606d0784326a",
            ),
            "P15": (
                44,
                "bedb18dc27516b36336fc7690319a05e2b32b1052ca4068c5cd0412860fcd585",
            ),
        }
        for name, (degree, digest) in expected.items():
            factors = self.result["extra_sections"][name]["factorization"][
                "numerator_factors"
            ]
            high = [factor for factor in factors if factor["degree"] > 8]
            self.assertEqual(len(high), 1)
            self.assertEqual(high[0]["degree"], degree)
            self.assertEqual(high[0]["primitive_coefficients_sha256"], digest)


if __name__ == "__main__":
    unittest.main()

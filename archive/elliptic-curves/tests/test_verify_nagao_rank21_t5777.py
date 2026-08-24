from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from ek_k3 import valuation  # noqa: E402
from nagao_1994 import RANK21_CONSTRUCTION  # noqa: E402
from verify_nagao_rank21_t5777 import (  # noqa: E402
    EXPECTED_HEIGHT_SUBSET,
    EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS,
    EXPECTED_RESIDUES,
    PARAMETER_T,
)


Q = Fraction


class NagaoRank21T5777VerifierTests(unittest.TestCase):
    def test_parameter_and_corrected_residues(self) -> None:
        self.assertEqual(PARAMETER_T, Q(5777, 32))
        residues = {
            prime: PARAMETER_T.numerator
            * pow(PARAMETER_T.denominator, -1, prime)
            % prime
            for prime in EXPECTED_RESIDUES
        }
        self.assertEqual(residues, EXPECTED_RESIDUES)
        self.assertEqual(residues[13], 3)
        self.assertEqual((-residues[13]) % 13, 10)

    def test_primitive_discriminant_local_profile(self) -> None:
        discriminant = RANK21_CONSTRUCTION.primitive_discriminant_value(PARAMETER_T)
        self.assertEqual(
            {
                prime: valuation(discriminant, prime)
                for prime in EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS
            },
            EXPECTED_PRIMITIVE_DISCRIMINANT_VALUATIONS,
        )

    def test_pinned_height_subset_has_seventeen_indices(self) -> None:
        self.assertEqual(len(EXPECTED_HEIGHT_SUBSET), 17)
        self.assertEqual(EXPECTED_HEIGHT_SUBSET[:11], tuple(range(1, 12)))


if __name__ == "__main__":
    unittest.main()

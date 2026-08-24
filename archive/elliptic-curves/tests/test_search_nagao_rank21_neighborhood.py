from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from crt_lattice import gauss_reduce  # noqa: E402
from nagao_1994 import (  # noqa: E402
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_PARAMETER,
)
from search_nagao_rank21_neighborhood import (  # noqa: E402
    DISCRIMINANT_POLYNOMIAL,
    LOCAL_CONDITIONS,
    classify_condition,
    discriminant_valuation_at_rational,
    enumerate_neighborhood,
    fixed_discriminant_valuation,
    forced_valuation,
    profile_crt,
    search_conditions,
)


class NagaoRank21NeighborhoodTests(unittest.TestCase):
    def test_discriminant_and_all_claimed_local_powers_are_exact(self) -> None:
        self.assertEqual(len(DISCRIMINANT_POLYNOMIAL) - 1, 20)
        self.assertTrue(
            all(
                DISCRIMINANT_POLYNOMIAL[index] == 0
                for index in range(1, 21, 2)
            )
        )
        self.assertEqual(fixed_discriminant_valuation(23), 2)
        self.assertEqual(
            {condition.label: forced_valuation(condition) for condition in LOCAL_CONDITIONS},
            {
                "p5-broad": 3,
                "p5-record-lift": 4,
                "p7-record": 4,
                "p13-record": 4,
                "p17-record": 3,
                "p23-record": 3,
            },
        )

    def test_each_search_ball_is_uniformly_split_multiplicative(self) -> None:
        for condition in search_conditions():
            classification = classify_condition(condition)
            self.assertEqual(classification["reduction"], "split multiplicative")
            self.assertEqual(classification["conductor_exponent"], 1)
            self.assertTrue(classification["split_multiplicative"])
            self.assertEqual(classification["tangent_legendre_symbol"], 1)

    def test_record_profile_crt_and_gauss_basis_are_pinned(self) -> None:
        residue, modulus = profile_crt()
        self.assertEqual((residue, modulus), (354_942, 889_525))
        self.assertEqual(
            gauss_reduce((modulus, 0), (residue, 1)),
            ((-763, -614), (-938, 411)),
        )
        self.assertEqual(
            RANK21_CONSTRUCTOR_PARAMETER, 2 * RANK21_PUBLISHED_PARAMETER
        )
        self.assertEqual(RANK21_CONSTRUCTOR_PARAMETER, Q(14_721, 188))
        self.assertEqual(
            (
                RANK21_CONSTRUCTOR_PARAMETER.numerator
                - residue * RANK21_CONSTRUCTOR_PARAMETER.denominator
            )
            % modulus,
            0,
        )

    def test_published_constructor_parameter_replays_record_valuations(self) -> None:
        observed = {
            condition.prime: discriminant_valuation_at_rational(
                RANK21_CONSTRUCTOR_PARAMETER.numerator,
                RANK21_CONSTRUCTOR_PARAMETER.denominator,
                condition.prime,
            )
            for condition in search_conditions()
        }
        self.assertEqual(observed, {5: 4, 7: 4, 13: 4, 17: 3, 23: 3})

    def test_lattice_box_is_exhaustive_and_contains_the_record(self) -> None:
        radius_ten = enumerate_neighborhood(10)
        radius_twelve = enumerate_neighborhood(12)
        self.assertEqual(len(radius_ten), 76)
        self.assertEqual(len(radius_twelve), 110)
        self.assertEqual(radius_twelve[0].parameter, Q(763, 614))
        self.assertEqual(radius_twelve[0].height, 763)
        self.assertTrue(
            any(
                candidate.parameter == RANK21_CONSTRUCTOR_PARAMETER
                for candidate in radius_ten
            )
        )
        residue, modulus = profile_crt()
        self.assertEqual(
            len({(candidate.numerator, candidate.denominator) for candidate in radius_twelve}),
            len(radius_twelve),
        )
        for candidate in radius_twelve:
            self.assertGreater(candidate.denominator, 0)
            self.assertEqual(
                (candidate.numerator - residue * candidate.denominator) % modulus,
                0,
            )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from nagao_1994 import RANK21_CONSTRUCTION, short_jacobian_coefficients  # noqa: E402
from search_nagao_rank21_unbiased import (  # noqa: E402
    CALIBRATION_PARAMETERS,
    INVARIANT_I,
    INVARIANT_J,
    build_residue_tables,
    polynomial_value,
    prefilter_population,
    primitive_population_count,
    projective_index,
    residue_score,
)


class NagaoRank21UnbiasedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = build_residue_tables()

    def test_invariant_polynomials_replay_the_exact_constructor(self) -> None:
        for parameter in (Q(1), Q(3, 31), Q(531, 2), Q(956, 9)):
            coefficients = short_jacobian_coefficients(
                RANK21_CONSTRUCTION, parameter
            )
            self.assertEqual(-27 * polynomial_value(INVARIANT_I, parameter), coefficients[3])
            self.assertEqual(-27 * polynomial_value(INVARIANT_J, parameter), coefficients[4])

    def test_projective_infinity_and_calibration_scores_are_pinned(self) -> None:
        self.assertEqual(projective_index(1, 5, 5), 5)
        expected = {
            Q(531, 2): (26.10514955621011, 40, 4),
            Q(956, 9): (25.3590985503976, 40, 4),
            Q(1637, 12): (20.83904798709409, 37, 7),
            Q(5777, 32): (21.789166497095458, 37, 7),
            Q(1, 5): (19.54446207713391, 38, 6),
        }
        for parameter, (score, good, bad) in expected.items():
            observed = residue_score(
                parameter.numerator, parameter.denominator, self.tables
            )
            self.assertAlmostEqual(observed[0], score, places=12)
            self.assertEqual(observed[1:], (good, bad))

    def test_full_population_count_and_calibration_exclusion_are_exact(self) -> None:
        self.assertEqual(primitive_population_count(10_000, 100), 608_337)
        retained, audit = prefilter_population(
            a_max=100,
            b_max=20,
            keep_count=25,
            tables=self.tables,
            excluded_parameters=CALIBRATION_PARAMETERS,
        )
        self.assertEqual(audit["complete_primitive_population_count"], 1_249)
        # None of the four large calibration parameters lies in this small box.
        self.assertEqual(audit["eligible_population_count"], 1_249)
        self.assertEqual(len(retained), 25)
        self.assertEqual(len({candidate.parameter for candidate in retained}), 25)
        self.assertEqual(retained[0].parameter, Q(5, 16))


if __name__ == "__main__":
    unittest.main()

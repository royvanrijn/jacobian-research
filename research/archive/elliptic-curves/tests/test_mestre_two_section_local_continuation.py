import unittest

from probe_mestre_two_section_local_continuation import small_pade_models


class MestreTwoSectionLocalContinuationTest(unittest.TestCase):
    def test_pade_requires_a_holdout_coefficient(self):
        # With 17 coefficients, a (8,8) fit has exactly enough unknowns to
        # interpolate but no independent coefficient on which to verify it.
        # It must not be reported as a recognized rational model.
        series = {"c1": [str((index + 1) ** 5 + 3 * index) for index in range(17)]}
        models = small_pade_models(series, max_degree=8)
        self.assertNotIn(
            {"coordinate": "c1", "numerator_degree": 8, "denominator_degree": 8},
            models,
        )


if __name__ == "__main__":
    unittest.main()

"""Do not convert an absent marking witness into a rational-point exclusion."""

import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py'
SPEC = importlib.util.spec_from_file_location('arithmetic_foundry', SCRIPT)
foundry = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(foundry)


class ArithmeticMarkingBoundary(unittest.TestCase):
    def row(self, classification, tests=()):
        return {'classification': classification, 'arithmetic_tests': list(tests),
                'full_discriminant_marking_curve': {'status': 'UNKNOWN'},
                'easy_quotient_maps': []}

    def test_absence_is_unknown_for_unscreened_and_unknown_rows(self):
        for row in (None, self.row('UNKNOWN')):
            with self.subTest(row=row):
                self.assertIsNone(foundry.full_curve_summary(row)['rational_non_CM_point'])

    def test_exclusion_and_positive_witness_are_distinct_decisions(self):
        excluded = foundry.full_curve_summary(self.row('ARITHMETICALLY_EXCLUDED'))
        possible = foundry.full_curve_summary(self.row('ARITHMETICALLY_POSSIBLE',
            [{'status': 'PASS_EXACT_FULL_RATIONAL_MARKING_WITNESS'}]))
        self.assertIs(excluded['rational_non_CM_point'], False)
        self.assertIs(possible['rational_non_CM_point'], True)


if __name__ == '__main__':
    unittest.main()

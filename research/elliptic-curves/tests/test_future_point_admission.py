import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from future_point_admission import FinitePointAdmission
from audit_recorded_point_mod2_rank_v3 import signature
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore


class AdmissionTests(unittest.TestCase):
    def test_exact_columns_against_reference_rows(self):
        model = (0, 0, 0, 0, -2)
        p = (F(3), F(5))
        admission = FinitePointAdmission(model, [p], prime_bound=97)
        cache = ReductionCache(MemoryFactStore())
        col, offset = 0, 0
        for prime in admission.primes:
            sig = signature(cache, model, [p], prime)
            for row in sig.rows:
                col |= int(row[0]) << offset
                offset += 1
        self.assertEqual(admission.columns[p], col)
        self.assertEqual(admission.consider((3, -5))['status'], 'KNOWN_POINT_UP_TO_SIGN')
        # 2P is exact, yet a mod-2 zero is not labelled rational dependence.
        twice = (F(129, 100), F(-383, 1000))
        self.assertEqual(admission.consider(twice)['status'], 'UNKNOWN_FINITE_COLUMN_IN_SPAN')
        with self.assertRaises(ValueError):
            admission.consider((3, 6))


if __name__ == '__main__':
    unittest.main()

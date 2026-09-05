from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState
from search_nagao_u42_skew_height import short_multiply


class QuotientOnlyTests(unittest.TestCase):
    def test_same_signatures_with_eviction_and_nonintegral_points(self):
        model = (0, 0, 0, -7, 10)
        points = [short_multiply(F(-7), (F(1), F(2)), n) for n in range(1, 8)]
        old = ReductionCache(MemoryFactStore())
        new = QuotientOnlyReductionCache(MemoryFactStore(), point_cache_limit=3)
        for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
            try:
                expected = old.signature(model, points, p)
            except ValueError:
                continue
            self.assertEqual(new.signature(model, points, p), expected)
            self.assertLessEqual(len(new._points), 3)
        self.assertTrue(any(q.denominator > 1 for point in points for q in point))
        self.assertEqual({f['record']['key']['namespace'] for f in new.store.snapshot()['facts']},
                         {'finite-field/mod2-quotient'})
        # The fifth multiple has x denominator 25 and reduces to O at5.
        q = points[4]
        self.assertEqual(q[0].denominator % 5, 0)
        self.assertEqual(new.point_signature(model, q, 5)[0], 0)

    def test_general_model_and_portable_state_replay(self):
        # y^2=x^3-7x+10 via y=Y+x gives [2,-1,0,-7,10].
        model = (2, -1, 0, -7, 10)
        points = [(1, 1), (2, 0)]
        old = ReductionCache(MemoryFactStore())
        new = QuotientOnlyReductionCache(MemoryFactStore(), point_cache_limit=0)
        a = raw_state(model, points, cache=old, prime_bound=31)
        b = raw_state(model, points, cache=new, prime_bound=31)
        self.assertEqual(a.record(), b.record())
        fresh = MemoryFactStore()
        fresh.import_snapshot(new.store.snapshot())
        self.assertEqual(MWState.from_record(b.record(), cache=ReductionCache(fresh)), b)
        self.assertEqual(len(new._points), 0)
        with self.assertRaises(ValueError):
            new.point_signature(model, (1, 3), 5)


if __name__ == '__main__':
    unittest.main()

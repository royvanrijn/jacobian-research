from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.cached_observation_state import CachedObservationMWState, _contains
from research_runtime.mw_state import MWState, PointObservation
from research_runtime.arithmetic import CurveModel
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.search_state import raw_state
from search_nagao_u42_skew_height import short_multiply


class CachedObservationTests(unittest.TestCase):
    def test_records_survive_known_ambiguous_and_new_admissions(self):
        model = (0, 0, 0, -7, 10)
        old_cache = ReductionCache(MemoryFactStore())
        new_cache = ReductionCache(MemoryFactStore())
        old = raw_state(model, [(1, 2)], cache=old_cache, prime_bound=31)
        new_cache.store.import_snapshot(old_cache.store.snapshot())
        new = CachedObservationMWState.from_record(old.record(), cache=new_cache)
        twice = short_multiply(F(-7), (F(1), F(2)), 2)
        for P in [(1, 2), (1, -2), twice, (2, 2), twice, (2, -2)]:
            old = old.adjoin(P, cache=old_cache, extra_primes=(5, 7, 11, 13, 17, 19, 23, 29, 31))
            new = new.adjoin(P, cache=new_cache, extra_primes=(5, 7, 11, 13, 17, 19, 23, 29, 31))
            self.assertEqual(new.record(), old.record())
            self.assertIsInstance(new, CachedObservationMWState)
        self.assertEqual(new.rank, 2)
        self.assertEqual(MWState.from_record(new.record(), cache=new_cache).record(), old.record())

    def test_cache_never_accepts_off_curve_observation_or_changed_model(self):
        cache = ReductionCache(MemoryFactStore())
        old = raw_state((0, 0, 0, -7, 10), [(1, 2)], cache=cache, prime_bound=31)
        state = CachedObservationMWState.from_record(old.record(), cache=cache)
        good = PointObservation(('1', '2'), 'KNOWN_POINT_UP_TO_SIGN', None)
        bad = PointObservation(('1', '3'), 'KNOWN_POINT_UP_TO_SIGN', None)
        replace(state, observations=(good,))
        for _ in range(2):
            with self.assertRaises(ValueError):
                replace(state, observations=(good, bad))
        self.assertTrue(_contains(CurveModel((0, 0, 0, -7, 10)), ('1', '2')))
        self.assertFalse(_contains(CurveModel((0, 0, 0, -7, 9)), ('1', '2')))


if __name__ == '__main__':
    unittest.main()

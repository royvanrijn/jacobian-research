import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
from research_runtime.preloaded_prime_state import preload
from research_runtime.rotated_observation_state import rotate
from alternate_quartic_covers import short_add

class RotationTests(unittest.TestCase):
    def test_independent_dependent_and_archived_history(self):
        model=(0,0,0,0,-2);P=(3,5);points=(P,short_add(model,P,P),(3,-5));histories=[]
        for rotating in (False,True):
            cache=QuotientOnlyReductionCache(MemoryFactStore());s,_=preload(raw_state(model,(),cache=cache,prime_bound=43),cache,43);history=[]
            for p in points:
                if rotating:
                    old=s; s,archive=rotate(s)
                    self.assertEqual(archive,old.record());self.assertEqual(s.basis,old.basis);self.assertEqual(s.parent_state,archive['key']);self.assertFalse(s.observations)
                s=s.adjoin(p,cache=cache)
                history.append((s.rank,s.basis,s.observations[-1].point,s.observations[-1].status,s.observations[-1].finite_relation_mask))
            histories.append(history)
        self.assertEqual(*histories)
if __name__=='__main__':unittest.main()

import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
from research_runtime.preloaded_prime_state import preload
from alternate_quartic_covers import short_add
from mod2_reduction_independence import _is_prime

class PrimeBankTests(unittest.TestCase):
    def test_independent_then_dependent_and_sign(self):
        model=(0,0,0,0,-2);P=(3,5);twice=short_add(model,P,P);histories=[]
        for full in (False,True):
            cache=QuotientOnlyReductionCache(MemoryFactStore());s=raw_state(model,(),cache=cache,prime_bound=43)
            if full:s,record=preload(s,cache,43);self.assertTrue(record['basis_unchanged'])
            history=[]
            for p in (P,twice,(3,-5)):
                s=s.adjoin(p,cache=cache,extra_primes=() if full else tuple(p for p in range(3,44) if _is_prime(p)))
                history.append((s.rank,s.basis,s.observations[-1].status,s.observations[-1].finite_relation_mask))
            histories.append(history)
        self.assertEqual(histories[0],histories[1]);self.assertEqual(histories[0][0][0],1)
    def test_rejects_invalid_bound(self):
        with self.assertRaises(ValueError):preload(None,None,1001)
if __name__=='__main__':unittest.main()

import sys,unittest,json
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from memory_rank_certificate import checked_rank
import certify_compact_r17_candidates as original
class MemoryRankTests(unittest.TestCase):
    def test_exact_backend_agreement_and_false_claim_rejection(self):
        E=tuple(map(F,(0,0,0,0,-2)));P=(F(3),F(5))
        self.assertEqual(json.dumps(checked_rank(E,[P],[5,7,11],7),sort_keys=True),json.dumps(original.checked_rank(E,[P],[5,7,11],7),sort_keys=True))
        with self.assertRaises(ArithmeticError):checked_rank(E,[P,P],[5,7,11],7)
        with self.assertRaises(ArithmeticError):checked_rank(E,[(F(3),F(4))],[5,7,11],7)
if __name__=='__main__':unittest.main()

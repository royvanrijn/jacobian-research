import sys,random,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from audit_retained_cloud_modl import insert,finite
class StreamedRows(unittest.TestCase):
    def test_agrees_with_independent_dense_elimination(self):
        rng=random.Random(17)
        for ell in (3,5):
            for _ in range(80):
                rows=[[rng.randrange(ell) for _ in range(9)] for _ in range(rng.randrange(1,12))];basis={}
                for row in rows:insert(basis,row,ell)
                self.assertEqual(sorted(basis),finite.pivots(rows,ell))
if __name__=='__main__':unittest.main()

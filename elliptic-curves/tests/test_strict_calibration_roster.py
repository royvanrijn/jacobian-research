import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from replay_backend_calibrations_strict import roster
class RosterTests(unittest.TestCase):
    def test_rejects_empty_duplicate_and_reordered(self):
        roster([{'chart':0},{'chart':1}],2)
        for rows in ([],[{'chart':0},{'chart':0}],[{'chart':1},{'chart':0}]):
            with self.assertRaises(ArithmeticError):roster(rows,2)
if __name__=='__main__':unittest.main()

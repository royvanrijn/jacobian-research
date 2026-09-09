"""Factor-free shared-cofactor hint regressions; hints are not primality proofs."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from run_inventory_conductor_long_v3 import shared_splits


class SharedSplits(unittest.TestCase):
    def test_exact_splits(self):
        values={'a':15,'b':21,'c':77,'d':13}
        result=shared_splits(values)
        self.assertEqual(result,{'a':{3},'b':{3,7},'c':{7},'d':set()})
        for identifier,splits in result.items():
            for split in splits:self.assertEqual(values[identifier]%split,0)

    def test_equal_and_dividing_cofactors(self):
        self.assertEqual(shared_splits({'a':15,'b':15,'c':3}),{'a':{3},'b':{3},'c':set()})


if __name__=='__main__':unittest.main()

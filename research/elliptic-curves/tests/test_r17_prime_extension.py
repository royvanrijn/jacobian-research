import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from extend_retained_r17_prime_scores import choose,sums
class PrimeExtensionTests(unittest.TestCase):
    def test_validation_does_not_select(self):
        rows=[dict(retained_index=i,combined_selection_units=i,combined_good=10,denominator=1,numerator=i,validation_units=-i) for i in range(8)]
        a=choose(rows)
        for r in rows:r['validation_units']=10**80*(8-r['retained_index'])
        self.assertEqual(a,choose(rows));self.assertEqual(a['extended_top_two'],[7,6])
        base=sums([[32749,1],[32771,2]])
        change=sums([[32749,1],[32771,-200]])
        self.assertEqual(base['extension_selection_units'],change['extension_selection_units']);self.assertNotEqual(base['validation_units'],change['validation_units'])
if __name__=='__main__':unittest.main()

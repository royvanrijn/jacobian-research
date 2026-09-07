"""Regression: columns remain the same generic sections across specializations."""
import json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import mestre_parent_adapter as adapter
from search_mestre_root_tuple_scale_max200 import mod_l_reduction_signature,gf_l_rank_and_pivots

class CoherentParentAdapterTests(unittest.TestCase):
    def test_all_historical_probes_match_independently_audited_columns(self):
        data=json.loads(adapter.SOURCE.read_text());rows=[]
        for r in data['probes']:
            model,points=adapter.specialize(r['outer_u'],r['T'])
            selected=points[:11]+points[12:]
            signature=mod_l_reduction_signature(model,selected,r['prime'],modulus=3)
            self.assertEqual([list(row) for row in signature.rows],r['coherent_rows'])
            rows.extend(signature.rows)
        self.assertEqual(gf_l_rank_and_pivots(rows,13,3)[0],11)

    def test_new_parent_inputs_have_fourteen_exact_labelled_points(self):
        for u in (11,13,17,19,23,29):
            model,points=adapter.specialize(u,1)
            self.assertEqual(len(points),14)
            self.assertTrue(all(y*y==x*x*x+model[3]*x+model[4] for x,y in points))

    def test_parent_and_fibre_poles_are_not_admitted(self):
        for u,t in ((0,1),(-2,1),(11,0)):
            with self.assertRaises((ValueError,ZeroDivisionError)):
                adapter.specialize(u,t)

if __name__=='__main__':unittest.main()

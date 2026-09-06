import copy
import importlib.util
import unittest
import retrospective as r
import local_collision as lc
import scalar_cup as scalar
import strict_cup_obstruction as nonscalar
import verify_scalar_cup_control as independent


class ScalarCup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(scalar.OUTPUT)
        cls.control=r.read(scalar.CONTROL)

    def test_explicit_obstruction_planes_and_radicals(self):
        pairs={0:[(1,2),(6,11),(30,132),(37,94)],
               4:[(1,2),(7,8),(45,139)],5:[(1,2),(5,17),(14,34)]}
        for row in self.data["production_cases"]:
            M=row["scalar_cup_matrix"];planes=pairs[row["case_index"]]
            vectors=[x for pair in planes for x in pair]
            radical=row["annihilator_masks"];n=len(M)
            self.assertEqual(r.rank(vectors+radical),n)
            for i,v in enumerate(vectors):
                for j,w in enumerate(vectors):
                    self.assertEqual(lc.pairing(v,w,M),int((i^1)==j))
            for v in radical:
                self.assertTrue(all(lc.pairing(v,1<<j,M)==0 for j in range(n)))
            self.assertEqual(row["detected_scalar_cup_rank"],len(vectors))

    def test_low_gain_control_also_has_full_obstruction(self):
        low=next(x for x in self.data["production_cases"] if x["case_index"]==5)
        self.assertEqual(low["detected_scalar_cup_rank"],low["strict_dimension"])
        self.assertEqual(low["retained_space_necessary_twist_solubility_dimension_upper_bound"],0)
        self.assertIn("No whole-curve rank upper bound",self.data["boundary"])

    def test_scalar_minus_one_is_not_parameter_minus_one(self):
        a=self.data["fixed_cubic_strict_case"]["scalar_cup_matrix"]
        b=next(x for x in r.read(nonscalar.OUTPUT)["cases"] if x["u"]==-1)["strict_CT_matrix"]
        self.assertEqual(a[0][1],1)
        self.assertEqual(b[0][1],0)
        self.assertNotEqual(a,b)

    @unittest.skipUnless(importlib.util.find_spec("sage"),"independent arithmetic replay uses Sage")
    def test_independent_nonzero_control(self):
        self.assertEqual(independent.verify(self.control),[[0,1],[1,0]])

    @unittest.skipUnless(importlib.util.find_spec("sage"),"independent arithmetic replay uses Sage")
    def test_corrupted_norm_witness_is_rejected(self):
        wrong=copy.deepcopy(self.control)
        wrong["norm_witnesses"][0]["norm_a"][0]="0"
        with self.assertRaises(AssertionError):
            independent.verify(wrong)

    @unittest.skipUnless(importlib.util.find_spec("sage"),"independent arithmetic replay uses Sage")
    def test_fabricated_zero_cup_is_rejected(self):
        wrong=copy.deepcopy(self.control)
        wrong["independent_norm_cup_matrix"]=[[0,0],[0,0]]
        with self.assertRaises(AssertionError):
            independent.verify(wrong)


if __name__=="__main__":
    unittest.main()

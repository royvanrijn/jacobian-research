import unittest
from fractions import Fraction
from math import isqrt
import retrospective as r
import component_kummer_gate as ex
from cubic_bridge import Cubic


class ComponentKummerGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(ex.OUTPUT)
        cls.inputs=r.read(ex.INPUT)["cases"]

    def test_S3_hypotheses_on_both_anchors(self):
        for row in self.inputs:
            A,B=map(Fraction,row["model"][3:]);disc=-4*A**3-27*B*B
            self.assertNotEqual(B*disc,0)
            square=disc>0 and isqrt(disc.numerator)**2==disc.numerator and isqrt(disc.denominator)**2==disc.denominator
            self.assertFalse(square)
            Bsquare=B>0 and isqrt(B.numerator)**2==B.numerator and isqrt(B.denominator)**2==B.denominator
            self.assertFalse(Bsquare)
            witnesses=[p for p in r.primes(97) if r.roots_at(str(A),str(B),p)==()]
            self.assertTrue(witnesses)

    def test_obstruction_has_an_exact_nonsquare_witness(self):
        for row,input_row in zip(self.data["cases"],self.inputs):
            A,B=map(Fraction,input_row["model"][3:])
            xs=[Fraction(P[0]) for P in input_row["generic_points"]]
            for rec in row["mixed_characters"]:
                i,j=rec["obstructing_indices"];w=rec["witness"];p=w["p"];root=w["root"]
                self.assertEqual((root**3+r.mod(A,p)*root+r.mod(B,p))%p,0)
                value=(r.mod(xs[i],p)-root)*(r.mod(xs[j],p)-root)%p
                self.assertEqual(pow(value,(p-1)//2,p),p-1)
                self.assertEqual(rec["arithmetic_generic_rank"],0)

    def test_same_class_control_preserves_free_independence(self):
        for row,input_row in zip(self.data["cases"],self.inputs):
            K=Cubic(*map(Fraction,input_row["model"][3:]))
            P0=tuple(map(Fraction,input_row["generic_points"][0]))
            beta0=K.sub(K.scalar(P0[0]),K.theta)
            control=row["same_class_control"]
            self.assertEqual(control["index_in_original_three_point_lattice"],4)
            for c in control["derived_points"]:
                x,y=map(Fraction,c["point"]);root=tuple(map(Fraction,c["square_root_of_class_ratio"]))
                self.assertEqual(K.mul(beta0,K.square(root)),K.sub(K.scalar(x),K.theta))
            self.assertEqual(control["finite_additive_component_gate"],"PASS")
            self.assertEqual(control["mixed_triple_rational_section_exists"],"UNKNOWN")
            infinity=next(x for x in row["pairwise_squareclass_witnesses"] if x["indices"]==[0,3])
            self.assertTrue(infinity["unequal_squareclasses"])

    def test_component_invariants_leave_only_one_bit(self):
        c=self.data["component_module"]
        self.assertEqual(c["even_subset_invariants"],[0])
        self.assertEqual(c["I2_component_flag_invariants"],[0,7])
        self.assertEqual(c["kappa_valuation_rows"],[[0,1,1],[1,0,1],[1,1,0]])

    def test_exact_arithmetic_rank_does_not_bound_production_fibres(self):
        for row in self.data["cases"]:
            self.assertEqual(row["pair_base_arithmetic_generic_rank"],2)
            self.assertEqual(row["triple_base_arithmetic_generic_rank"],3)
            self.assertEqual(row["production_curve_rank"],"UNKNOWN")
        replay=r.read(r.OUT/"rank_jump_component_kummer_gate_verification_v1.json")
        self.assertEqual(replay["status"],"PASS")
        self.assertEqual(replay["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(replay["verifier_sha256"],r.digest((ex.HERE/"verify_component_kummer_gate.py").read_bytes()))
        self.assertEqual(self.data["bindings"],ex.bindings())


if __name__=="__main__":
    unittest.main()

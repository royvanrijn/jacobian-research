import base64
from fractions import Fraction
import unittest
import retrospective as r
import triple_character as ex


class TripleCharacter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(ex.OUTPUT)
        cls.geometry=r.read(r.OUT/"rank_jump_triple_character_geometry_v1.json")

    def test_projection_uses_only_generic_points(self):
        panel={x["id"]:x for x in r.read(r.INPUT)["rows"]}
        for row in r.read(ex.INPUT)["cases"]:
            src=panel[row["id"]]
            model,points=r.short(src["model"],src["generic_points"][:3])
            self.assertEqual(row["generic_point_indices"],[0,1,2])
            self.assertEqual((row["model"],row["generic_points"]),(model,points))
            self.assertNotIn("points",row)
        self.assertEqual(r.read(ex.PROTOCOL)["limits"]["new_rational_point_searches"],0)

    def test_zero_ranks_have_two_place_witnesses(self):
        for case,row in enumerate(self.data["bounds"]):
            chars={x["mask"]:x for x in row["characters"]}
            for mask in (5,6):
                self.assertEqual(chars[mask]["geometric_rank_interval"],[0,0])
                primes=chars[mask]["zero_rank_witness_primes"]
                self.assertEqual(len(primes),2)
                witnesses=[x for x in self.data["reductions"]
                           if x["case"]==case and x["mask"]==mask and x["p"] in primes]
                self.assertTrue(all(x["reduction_geometric_Picard_rank"]==18 for x in witnesses))
                self.assertEqual(len({x["NS_discriminant_squareclass"] for x in witnesses}),2)
            for mask in (3,7):
                self.assertEqual(chars[mask]["geometric_rank_interval"],[0,1])
                self.assertIsNone(chars[mask]["zero_rank_witness_primes"])
            self.assertEqual(row["full_base_arithmetic_rank_interval"],[3,5])
            self.assertEqual(row["production_curve_rank"],"UNKNOWN")

    def test_smooth_infinity_changes_the_trace(self):
        row=next(x for x in self.data["reductions"] if (x["case"],x["mask"],x["p"])==(0,7,17))
        finite=[x["finite_sum"] for x in row["fields"]]
        self.assertEqual(finite,[28,-186,6355])
        self.assertEqual([x["infinity_sum"] for x in row["fields"]],[3,25,-126])
        self.assertEqual(row["traces"],[31,-161,6229])
        # Omitting infinity fails the reciprocal five-dimensional reconstruction.
        solutions=0
        for e in (-1,1):
            p=17;s=finite[0]-e*p;c=Fraction(s*s-(finite[1]-p*p),2)
            solutions+=finite[2]-e*p**3==s**3-3*s*c+3*s*p*p
        self.assertEqual(solutions,0)

    def test_independent_trace_chain_is_bound(self):
        replay=r.read(r.OUT/"rank_jump_triple_character_verification_v1.json")
        self.assertEqual(replay["status"],"PASS")
        self.assertEqual(replay["counts_sha256"],r.digest(ex.RAW.read_bytes()))
        self.assertEqual(replay["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(replay["new_directly_counted_base_parameters"],2926)
        self.assertEqual(replay["previous_independently_verified_base_parameters"],20056)
        self.assertEqual(replay["verifier_sha256"],r.digest((ex.HERE/"verify_triple_character.py").read_bytes()))
        self.assertEqual(self.geometry["verifier_sha256"],
                         r.digest((ex.HERE/"verify_triple_character_geometry.py").read_bytes()))

    def test_auxiliary_nontorsion_witnesses_are_actual_points(self):
        self.assertEqual(self.geometry["base_genus"],1)
        self.assertEqual(self.geometry["norm_map_degree"],4)
        for row in self.geometry["auxiliary_bases"]:
            self.assertTrue(row["point_is_nontorsion"])
            a1,a2,a3,a4,a6=map(Fraction,row["model"])
            self.assertEqual((a1,a3),(0,0))
            for point in (row["point"],row["bound_multiple"]):
                x,y=map(Fraction,point)
                self.assertEqual(y*y,x**3+a2*x*x+a4*x+a6)


if __name__=="__main__":
    unittest.main()

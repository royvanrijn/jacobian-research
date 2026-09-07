import unittest
from fractions import Fraction as F
from math import comb
from pathlib import Path
import retrospective as r
import nonscalar_cup_control as ex
import verify_nonscalar_cup_control as verify
import nonscalar_unramified_lifts as lifts
import nonscalar_cup_orbit as orbit


class NonscalarCupControl(unittest.TestCase):
    def test_nonscalar_and_opposite_obstruction_maps(self):
        data=r.read(ex.OUTPUT)
        self.assertEqual(data["bindings"],ex.bindings())
        self.assertEqual(data["bilinearity_status"],"PASS")
        self.assertEqual([x["matrix"] for x in data["records"]],[[[0,1],[1,0]],[[0,0],[0,0]]])
        f=lambda t:sum(c*t**i for i,c in enumerate(verify.POLY))
        for left,right in [(-2,-1),(-1,0),(12,13)]:self.assertLess(f(left)*f(right),0)
        self.assertEqual([verify.norm(list(map(F,x["gamma"]))) for x in data["records"]],[-1,1])

    def test_replay_has_complete_support_and_bound_sources(self):
        result=r.read(verify.OUTPUT)
        self.assertEqual(result["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(result["verifier_sha256"],r.digest(Path(verify.__file__).read_bytes()))
        support={p for x in result["records"] for w in x["witnesses"] for p in w["support_primes"]}
        self.assertEqual(support,{5,13,37})
        self.assertEqual({p["extension"] for x in result["records"] for w in x["witnesses"] for p in w["places"]},{"split","inert"})

    def test_principal_corrections_construct_both_lifts(self):
        lifts.build(check=True)
        data=r.read(lifts.OUTPUT)
        self.assertEqual(data["norm_image_dimension"],2)
        self.assertEqual(data["elliptic_rational_solubility"],"UNKNOWN")
        for row in data["lifts"]:
            self.assertTrue(all(v%2==0 for p in row["corrected_valuation_checks"] for v in p["corrected_valuations"]))

    def test_unit_orbit_exhausts_eight_squareclasses(self):
        orbit.build(check=True)
        data=r.read(orbit.OUTPUT)
        self.assertEqual(len(data["records"]),8)
        self.assertEqual([r["unramified_norm_image_dimension_on_U"] for r in data["records"]],[2,2,2,2,0,0,0,0])
        theta_image=list(map(F,data["automorphism_theta"]))
        self.assertEqual(verify.mul(theta_image,[F(1),F(1),F(0)]),[-1,0,0])

    def test_elliptic_quotients_are_self_and_scalar_twist(self):
        # Exact expansion f(z-1), followed by reciprocal reversal.
        shifted=[sum(verify.POLY[i]*comb(i,j)*(-1)**(i-j) for i in range(j,4)) for j in range(4)]
        self.assertEqual(shifted,[1,11,-14,1])
        reciprocal=shifted[::-1]
        self.assertEqual(reciprocal,[1,-14,11,1])
        minus=[reciprocal[i]*(-1)**(3-i) for i in range(4)]
        self.assertEqual(minus,verify.POLY)


if __name__=="__main__":
    unittest.main()

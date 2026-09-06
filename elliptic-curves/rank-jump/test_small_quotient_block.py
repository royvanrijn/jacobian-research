import unittest
from fractions import Fraction as F
import retrospective as r
import small_quotient_block as ex
import verify_small_quotient_block as verifier
import small_quotient_covers as covers


class SmallQuotientBlock(unittest.TestCase):
    def test_exact_rank_and_full_selmer_accounting(self):
        data=r.read(ex.OUTPUT)
        self.assertEqual(data["bindings"],ex.bindings())
        self.assertEqual(data["exact_MW_ranks"],[1,3])
        self.assertEqual(data["full_2_Selmer_dimensions"],[3,3])
        self.assertEqual(data["Sha_2_dimensions"],[2,0])
        self.assertEqual(data["strict_rational_dimensions"],[0,2])
        for rank,sha,selmer in zip(data["exact_MW_ranks"],data["Sha_2_dimensions"],data["full_2_Selmer_dimensions"]):
            self.assertEqual(rank+sha,selmer)

    def test_rank_proof_does_not_use_new_descent_upper_bounds(self):
        verifier.verify(check=True)
        result=r.read(verifier.OUTPUT)
        self.assertFalse(result["new_descent_upper_bounds_used_as_proof"])
        self.assertEqual(result["original_strict_CT_rank"],2)
        self.assertEqual(result["norm_kernel_dimension"]-result["real_condition_codimension"],3)

    def test_every_cover_witness_checks_the_actual_quadrics(self):
        covers.build(check=True)
        for row in r.read(covers.OUTPUT)["records"]:
            z=list(map(F,row["twist_rational_point"]));h=z.pop()
            self.assertEqual(covers.evaluate(list(map(F,row["Q2"])),z),0)
            self.assertEqual(covers.evaluate(list(map(F,row["Q1"])),z),h*h)
            self.assertNotEqual(covers.evaluate(list(map(F,row["Q1"])),z)+h*h,0)
            self.assertEqual(row["original_rational_solubility"],"NO: nonzero Sha[2] class, exact rank/CT certificate")

    def test_cover_coefficients_expand_in_the_cubic_algebra(self):
        for row in r.read(covers.OUTPUT)["records"]:
            beta=list(map(F,row["beta"]))
            for z in ([F(1),F(2),F(3)],[F(-2),F(1,3),F(4)]):
                expanded=covers.a.mul(beta,covers.a.mul(z,z))
                expected=[covers.evaluate(list(map(F,row[key])),z) for key in ("Q0","Q1","Q2")]
                self.assertEqual(expanded,expected)


if __name__=="__main__":
    unittest.main()

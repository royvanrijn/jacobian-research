import unittest
from fractions import Fraction
from itertools import permutations
from math import isqrt
from pathlib import Path
import retrospective as r
import equal_class_high_followup as ex


class HighFollowup(unittest.TestCase):
    def test_coefficient_only_prime_rule(self):
        protocol=r.read(ex.PROTOCOL)
        row=r.read(ex.ex.INPUT)["cases"][0]
        decisions=[{"p":p,"good":ex.ex.previous.good(row,p,7)} for p in r.primes(29) if p>17]
        self.assertEqual(decisions,protocol["eligibility"])
        self.assertEqual([v["p"] for v in decisions if v["good"]],[29])

    def test_match_leaves_high_unknown_and_low_proof_intact(self):
        report=r.read(ex.OUTPUT);old=r.read(ex.OLD_OUTPUT)
        self.assertEqual(report["bounds"],old["bounds"])
        self.assertEqual(report["reductions"][:-1],old["reductions"])
        new=report["reductions"][-1]
        self.assertEqual(new["traces"],[-25,335,7865])
        self.assertEqual(new["NS_discriminant_squareclass"],-1)
        self.assertEqual(new["reduction_geometric_Picard_rank"],18)
        self.assertEqual(report["bounds"][0]["arithmetic_mixed_rank_interval"],[0,1])

    def test_independent_counts_and_provenance(self):
        result=r.read(ex.VERIFICATION)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["bindings"],ex.bindings())
        self.assertEqual(result["counts_sha256"],r.digest(ex.RAW.read_bytes()))
        self.assertEqual(result["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(sum(v["base_parameters"] for v in result["new_record"]["fields"]),25259)
        self.assertEqual(sum(v["direct_orbit_character_sums"] for v in result["new_record"]["fields"]),8613)

    def test_trivial_lattice_and_contact_square_equation(self):
        D4=[[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]
        determinant=0
        for perm in permutations(range(4)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
            term=sign
            for i in range(4):term*=D4[i][perm[i]]
            determinant+=term
        trivial=-determinant**3*(-2)**3
        self.assertEqual(trivial,512)
        allowed=[]
        for m in range(64):
            height=4+2*m-3*Fraction(1,2)
            self.assertEqual(Fraction(trivial)*height/(5+4*m),256)
            if isqrt(5+4*m)**2==5+4*m:allowed.append(m)
        self.assertEqual(allowed,[1,5,11,19,29,41,55])
        self.assertEqual(4+2*allowed[0]-Fraction(3,2),Fraction(9,2))
        for m in allowed:
            height=Fraction(5+4*m,2)
            twice_contact=3+4*m
            self.assertEqual(4*height,4+2*twice_contact)
            self.assertEqual(isqrt(twice_contact+2)**2,twice_contact+2)


if __name__=="__main__":
    unittest.main()

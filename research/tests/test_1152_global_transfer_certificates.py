"""Arithmetic tamper controls and explicit CAS-assurance boundaries."""
from copy import deepcopy
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,ROOT/"elkies-k3/scripts"/filename)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


F=load("fibration_distinction","verify_1152_fibration_distinction.py")
J=load("paired_jacobians","verify_1152_pair_jacobians.py")
A=load("first_pair_arithmetic","verify_1152_first_pair_arithmetic.py")


class GlobalTransferTamperTests(unittest.TestCase):
    def certificate(self,module):
        return json.loads(module.DEFAULT.read_text())

    def test_critical_polynomial_constant_tamper_is_rejected(self):
        d=self.certificate(F)
        d["models"][0]["critical_value_characteristic"][0]^=1
        with self.assertRaisesRegex(ValueError,"complete critical-value"):
            F.verify(d)

    def test_surface_identity_overclaim_is_rejected(self):
        d=self.certificate(F);d["surface_identity"]="ISOMORPHIC"
        with self.assertRaisesRegex(ValueError,"surface identity"):
            F.verify(d)

    def test_jacobian_twist_sign_is_rejected(self):
        d=self.certificate(J)
        row=d["pairs"][0]
        row["jacobian_a_invariants"][1]=str(-Q(row["jacobian_a_invariants"][1]))
        with self.assertRaisesRegex(ValueError,"Jacobian coefficients"):
            J.verify(d)

    def test_wrong_isogeny_quotient_is_rejected(self):
        d=self.certificate(J)
        row=d["pairs"][0];row["quotient_a_invariants"][3]=str(Q(row["quotient_a_invariants"][3])+1)
        with self.assertRaisesRegex(ValueError,"quotient coefficients"):
            J.verify(d)

    def test_missing_bad_prime_is_rejected(self):
        d=self.certificate(A);d["local_witnesses"].pop()
        with self.assertRaisesRegex(ValueError,"all bad primes"):
            A.verify(d)

    def test_false_2adic_witness_is_rejected(self):
        d=self.certificate(A)
        next(z for z in d["local_witnesses"] if z["prime"]=="2")["integer_parameter"]=0
        with self.assertRaisesRegex(ValueError,"local square"):
            A.verify(d)

    def test_false_large_prime_certificate_is_rejected(self):
        d=self.certificate(A)
        key=max(d["prime_certificates"],key=int)
        d["prime_certificates"][key]["full_order_witness"]=1
        with self.assertRaisesRegex(ValueError,"Lucas full-order"):
            A.verify(d)

    def test_false_quartic_map_is_rejected(self):
        d=self.certificate(A)
        z=d["selmer_basis_covers"][0]["map_x_y_quadratic_algebra"][0][0]
        z["num"][0]=str(Q(z["num"][0])+1)
        with self.assertRaisesRegex(ValueError,"quartic-to-Jacobian"):
            A.verify(d)

    def test_arithmetic_replay_does_not_promote_global_completeness(self):
        r=A.verify(self.certificate(A))
        self.assertTrue(r["all_places_locally_soluble"])
        self.assertFalse(r["independent_selmer_completeness"])
        self.assertFalse(r["independent_exact_rank"])
        self.assertEqual(r["anti_invariant_heights_on_paired_base"],["24","24"])


if __name__=="__main__":unittest.main()

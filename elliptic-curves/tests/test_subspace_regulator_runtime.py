"""Cheap regressions for the two theorem-directed pre-search gates."""

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"elliptic-curves/cas"))
from research_runtime.arithmetic import ArithmeticContext,CurveModel,TwoTorsionContext
from research_runtime.subspace import GlobalSquareclasses,SubspaceDescent,local_intersection,restricted_radical
from research_runtime.regulator import Surface,VerifiedReduction,frobenius_invariants,pre_search_gate


class SubspaceTests(unittest.TestCase):
    def test_locally_bad_generators_can_have_good_sum(self):
        self.assertEqual(local_intersection(3, [[[1],[1],[0]],[[0],[0],[1]]]),(3,))

    def test_actual_fixed_field_matrix(self):
        source=ROOT/"artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
        ct=source.with_name("fixed_cubic_u_minus1_cassels_tate_v1.json")
        if not(source.exists() and ct.exists()):self.skipTest("retained fixed-field witnesses absent")
        data=json.loads(source.read_text())
        for row in data["runs"]:
            maps=[r["known_span_quotient_rows"] for r in row["finite_local_conditions"]]
            maps.append(row["real_local_condition"]["known_span_quotient_rows"])
            masks=local_intersection(20,maps)
            self.assertEqual(list(masks),[r["mask"] for r in row["W_u_basis"]])
            if row["parameter_u"]=="-1":
                result=restricted_radical(masks,json.loads(ct.read_text())["arithmetic"]["matrix"],global_dimension=20)
                self.assertEqual(result["pairing_rank"],16)
                self.assertEqual(result["radical_dimension"],2)
                self.assertEqual(result["obstructed_class_count"],262140)
                self.assertEqual(result["nonzero_compatible_class_count"],3)
                self.assertIsNone(result["full_curve_rank_upper"])

    def test_bad_matrices_fail_closed(self):
        for pairing in ([[1,0],[0,0]],[[0,1],[0,0]],[[0]]):
            with self.assertRaises(ValueError):restricted_radical([1,2],pairing,global_dimension=2)
        with self.assertRaises(ValueError):local_intersection(2,[[[1],[0,1]]])

    def test_witness_replay_never_calls_discovery(self):
        model=CurveModel((0,0,0,-1,1));algebra=TwoTorsionContext(model.two_division_polynomial)
        context=ArithmeticContext(model,model,(1,0,0,0),((2,4),(23,1)),algebra,(0,1,0))
        classes=GlobalSquareclasses(algebra.key,(("1","1","0"),("2","1","0")),"toy-independent-input")
        class Backend:
            allow_discovery=True
            def verify_global(self,*args):return True
            def required_places(self,*args):return (2,23,"infinity")
            def local_map(self,*args):
                assert self.allow_discovery
                return {"quotient_rows":[[],[]]}
            def verify_local(self,*args):return True
            def cover(self,*args):
                assert self.allow_discovery
                return {"mask":args[-1]}
            def verify_cover(self,*args):return args[-1]["mask"]==args[-2]
            def ct_pairing(self,*args):
                assert self.allow_discovery
                return {"matrix":[[0,1],[1,0]]}
            def verify_ct(self,*args):return True
        backend=Backend();pipeline=SubspaceDescent(context,classes,backend)
        witness=pipeline.run();backend.allow_discovery=False
        self.assertEqual(pipeline.run(retained=witness),witness)
        tampered=deepcopy(witness);tampered["radical"]["obstructed_class_count"]=0
        with self.assertRaises(ArithmeticError):pipeline.run(retained=tampered)
        backend.verify_local=lambda *args:False
        with self.assertRaises(ArithmeticError):pipeline.run(retained=witness)


class RegulatorTests(unittest.TestCase):
    def setUp(self):self.surface=Surface((1,1),(1,0,1),(1,0,1))

    def row(self,prime,rank,value,geometric=4,isometry=True):
        return VerifiedReduction(self.surface.key,prime,rank,geometric,str(value),str(prime),isometry)

    def test_rank_one_and_higher_rank_obstructions(self):
        for rank in [1,2,3]:
            gate=pre_search_gate(self.surface,[self.row(131,rank,18),self.row(137,rank,75)],candidate_rank=rank)
            self.assertEqual(gate["arithmetic_rank_upper"],rank-1)
            self.assertTrue(gate["rank_target_excluded"])
            self.assertEqual(gate["geometric_rank_upper"],4)

    def test_compatible_unknown_and_higher_reduction_ineligible(self):
        gate=pre_search_gate(self.surface,[self.row(131,1,18),self.row(137,1,50)],candidate_rank=1)
        self.assertEqual(gate["status"],"UNKNOWN")
        self.assertTrue(gate["section_search_eligible"])
        gate=pre_search_gate(self.surface,[self.row(131,2,18),self.row(137,1,75)],candidate_rank=1)
        self.assertEqual(gate["two_prime_regulator_test"],"NOT_APPLICABLE")
        self.assertEqual(gate["arithmetic_rank_upper"],1)

    def test_regulator_sharpens_a_bound_below_the_original_target(self):
        gate=pre_search_gate(self.surface,[self.row(131,1,18),self.row(137,1,75)],
                             candidate_rank=3,candidate_regulator=7)
        self.assertEqual(gate['arithmetic_rank_upper'],0)
        self.assertEqual(gate['height_tests'],[])  # A rank-three Gram cannot be compared with rank one.

    def test_distinct_primes_and_height_isometry_required(self):
        with self.assertRaises(ValueError):
            pre_search_gate(self.surface,[self.row(131,1,18),self.row(131,1,75)],candidate_rank=1)
        gate=pre_search_gate(self.surface,[self.row(131,1,18),self.row(137,1,75,isometry=False)],candidate_rank=1)
        self.assertEqual(gate["status"],"UNKNOWN")
        gate=pre_search_gate(self.surface,[self.row(131,1,18)],candidate_rank=1,candidate_regulator=10)
        self.assertTrue(gate["candidate_regulator_excluded"])
        self.assertFalse(gate["rank_target_excluded"])

    def test_actual_rank_two_frobenius_and_moment_tamper(self):
        path=ROOT/"artifacts/generated-results/elkies-k3-r17-product-alternate-orbit-0f82c--alternate-orbit-025be-p131-toric-frobenius-v1.json"
        if not path.exists():self.skipTest("retained Frobenius witness absent")
        data=json.loads(path.read_text());row=data["elliptic_L"]
        facts=frobenius_invariants(row["frobenius_characteristic_coefficients_low_to_high"],131,
                                  expected_degree=28,moments=row["power_sums_n1_n2"])
        self.assertEqual((facts["arithmetic_rank_upper"],facts["geometric_rank_upper"]),(2,2))
        with self.assertRaises(ValueError):
            frobenius_invariants(row["frobenius_characteristic_coefficients_low_to_high"],131,expected_degree=28,moments=[0,0])

    def test_no_cyclotomic_cutoff(self):
        # F(Z)=(Z-p)(Z+p)(Z²+pZ+p²) has arithmetic rank 1 and geometric rank 4.
        p=5
        facts=frobenius_invariants([-p**4,-p**3,0,p,1],p,expected_degree=4)
        self.assertEqual(facts["arithmetic_rank_upper"],1)
        self.assertEqual(facts["geometric_rank_upper"],4)


if __name__ == "__main__":unittest.main()

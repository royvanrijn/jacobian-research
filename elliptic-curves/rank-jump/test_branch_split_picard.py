import base64
import struct
import unittest
from pathlib import Path
import retrospective as r
import branch_split_picard as gate
from cubic_bridge import Cubic


class BranchSplitPicard(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(gate.OUTPUT)
        cls.raw=r.read(gate.INPUT)
        cls.geometry=r.read(r.OUT/"rank_jump_branch_split_geometry_v1.json")
        cls.replay=r.read(r.OUT/"rank_jump_branch_split_picard_verification_v1.json")

    def test_independent_replay_binds_every_trace(self):
        self.assertEqual(self.raw["bindings"],gate.bindings())
        self.assertEqual(self.data["input_sha256"],r.digest(gate.INPUT.read_bytes()))
        for record in self.replay["verification"]:
            self.assertEqual(record["status"],"PASS")
            self.assertEqual(record["input_sha256"],r.digest(gate.INPUT.read_bytes()))
            self.assertEqual(record["analysis_sha256"],r.digest(gate.OUTPUT.read_bytes()))
            self.assertEqual(record["verifier_sha256"],
                             r.digest((gate.HERE/"verify_branch_split_picard.py").read_bytes()))
            fields=next(x["fields"] for x in self.raw["primes"] if x["p"]==record["prime"])
            for row,verified in zip(fields,record["fields"]):
                traces=gate.decode(row)
                packed=struct.pack("<"+"h"*len(traces),*traces)
                self.assertEqual(r.digest(packed),verified["trace_sha256"])
                self.assertEqual(row["residual_H2_traces"],verified["residual_H2_traces"])
                self.assertEqual(len(traces),verified["base_parameters"])

    def test_corrupt_trace_length_is_rejected(self):
        row=dict(self.raw["primes"][0]["fields"][0])
        row["trace_values_i16_le_base64"]=base64.b64encode(b"\0\0").decode()
        with self.assertRaises(AssertionError):
            gate.decode(row)

    def test_picard_endpoint_requires_distinct_squareclasses(self):
        witnesses=self.data["two_place_witness"]
        self.assertEqual({x["p"] for x in witnesses},{23,59})
        self.assertEqual({x["NS_discriminant_squareclass"] for x in witnesses},{-19,-23})
        for row in witnesses:
            p=row["p"];t=-row["quadratic_factor_ascending"][1]
            self.assertEqual(row["t2"],t*t-p*p)
            self.assertEqual(row["t1"],row["algebraic_eigenvalue"]+t)
            self.assertNotIn(t,(-2*p,-p,0,p,2*p))
            self.assertEqual(gate.squarefree(t*t-4*p*p),row["NS_discriminant_squareclass"])
        unknown=[x for x in self.data["reductions"] if x["status"]=="UNKNOWN"]
        self.assertEqual(len(unknown),2)
        for row in unknown:
            self.assertNotIn(row,witnesses)
            self.assertNotIn("geometric_picard_rank",row)

    def test_full_split_control_uses_no_exceptional_point(self):
        anchor=r.read(gate.lc.INPUT)["anchor"]
        A,B=map(r.F,anchor["short_model_ainvariants"][3:])
        K=Cubic(A,B);c=self.geometry["split_control"]
        u=r.F(c["parameter"]);z=tuple(map(r.F,c["square_root_coordinates"]))
        self.assertNotEqual(u,0)
        self.assertEqual(K.square(z),K.sub(K.one,K.scale(K.theta,u)))
        self.assertEqual(K.norm(z)**2,1+A*u*u+B*u**3)
        self.assertNotEqual(K.norm(z),0)
        self.assertEqual(c["specialized_curve_rank"],"UNKNOWN")
        self.assertEqual(self.geometry["verifier_sha256"],
                         r.digest((gate.HERE/"verify_branch_split_geometry.py").read_bytes()))
        self.assertEqual(self.geometry["picard_certificate_sha256"],r.digest(gate.OUTPUT.read_bytes()))

    def test_character_rank_sum_does_not_bound_fibres(self):
        self.assertEqual(self.geometry["character_ranks"],[1,0,0,0,0,0,0,0])
        for row in self.geometry["character_geometry"]:
            euler=sum(8 if x=="I2*" else 2 for x in row["finite"])
            euler+=6 if row["infinity"]=="I0*" else 0
            self.assertEqual(euler,12*row["chi"])
        self.assertEqual(self.data["full_branch_cover_geometric_MW_rank"],1)
        self.assertIn("No specialized curve rank",self.data["boundary"])
        limits=r.read(gate.PROTOCOL)["limits"]
        self.assertEqual(limits["new_parameter_searches"],0)
        self.assertEqual(limits["new_rational_point_searches"],0)


if __name__=="__main__":
    unittest.main()

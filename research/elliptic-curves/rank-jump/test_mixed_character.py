import base64
import struct
import unittest
from pathlib import Path
import retrospective as r
import mixed_character as ex
from verify_mixed_character import Field


class MixedCharacter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.counts=r.read(ex.COUNTS)
        cls.data=r.read(ex.OUTPUT)
        cls.geom=r.read(r.OUT/"rank_jump_mixed_character_geometry_v1.json")
        cls.verified=r.read(r.OUT/"rank_jump_mixed_character_verification_v1.json")["records"]

    def test_no_exceptional_point_in_projection(self):
        panel=r.read(r.INPUT);indexed={x["id"]:x for x in panel["rows"]}
        projection=r.read(ex.INPUT)
        self.assertEqual(projection["source_sha256"],r.digest(r.INPUT.read_bytes()))
        for row in projection["cases"]:
            src=indexed[row["id"]]
            model,points=r.short(src["model"],src["generic_points"][:2])
            self.assertEqual(row["generic_point_indices"],[0,1])
            self.assertEqual(row["model"],model)
            self.assertEqual(row["generic_points"],points)
            self.assertNotIn("points",row)
        self.assertEqual(r.read(ex.PROTOCOL)["limits"]["point_searches"],0)

    def test_all_arrays_have_independent_replay(self):
        self.assertEqual(len(self.verified),4)
        total=0
        for verified in self.verified:
            self.assertEqual(verified["status"],"PASS")
            self.assertEqual(verified["counts_sha256"],r.digest(ex.COUNTS.read_bytes()))
            self.assertEqual(verified["verifier_sha256"],
                             r.digest((ex.HERE/"verify_mixed_character.py").read_bytes()))
            raw=next(x for x in self.counts["reductions"]
                     if (x["case"],x["p"])==(verified["case"],verified["p"]))
            for a,b in zip(raw["fields"],verified["fields"]):
                packed=base64.b64decode(a["trace_values_i16_le_base64"])
                self.assertEqual(len(packed),2*a["q"])
                self.assertEqual(r.digest(packed),b["trace_sha256"])
                self.assertEqual(a["residual_H2_trace"],b["residual_H2_trace"])
                total+=a["q"]
        self.assertEqual(total,20056)

    def test_matching_discriminants_do_not_close_rank(self):
        usable=[x for x in self.data["reductions"] if x["status"]=="RHO_18_REDUCTION"]
        self.assertEqual(len(usable),3)
        self.assertEqual({x["NS_discriminant_squareclass"] for x in usable},{-1})
        for row in self.data["cases"]:
            self.assertIsNone(row["two_place_witness"])
            self.assertEqual(row["geometric_Picard_rank"],"UNKNOWN")
            self.assertEqual(row["mixed_character_geometric_MW_rank"],"UNKNOWN")
        extra=next(x for x in self.data["reductions"] if x["case"]==1 and x["p"]==17)
        self.assertEqual(extra["reduction_geometric_Picard_rank"],20)
        self.assertEqual(extra["status"],"UNKNOWN")

    def test_bounds_are_on_the_new_function_field(self):
        self.assertEqual(self.geom["base_genus"],0)
        for row in self.geom["bounds"]:
            self.assertEqual(row["mixed_geometric_rank_interval"],[0,1])
            self.assertEqual(row["full_pair_base_geometric_rank_interval"],[3,4])
            self.assertEqual(row["full_pair_base_arithmetic_rank_interval"],[2,3])
            self.assertEqual(row["production_curve_rank"],"UNKNOWN")
        self.assertEqual(self.geom["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(self.geom["verifier_sha256"],
                         r.digest((ex.HERE/"verify_mixed_character_geometry.py").read_bytes()))

    def test_nonsplit_node_regression(self):
        # This prime caused the first count attempt to stop. The cubic has only
        # one rational root; its two nonrational roots still coalesce rationally.
        row=r.read(ex.INPUT)["cases"][1];p=19
        A,B=map(lambda x:r.mod(x,p),row["model"][3:])
        roots=[t for t in range(p) if (t**3+A*t+B)%p==0]
        self.assertEqual(len(roots),1)
        theta=roots[0];u=pow(theta,-1,p)
        double=(-A*u-theta)%p;single=2*theta%p
        self.assertNotEqual(double,single)
        F=Field(p,1,[0,1])
        vals=[]
        for x in range(p):
            f=(x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3)%p
            self.assertEqual(f,((x-double)**2*(x-single))%p)
            vals.append(0 if f==0 else (1 if pow(f,(p-1)//2,p)==1 else -1))
        raw=next(x for x in self.counts["reductions"] if x["case"]==1 and x["p"]==19)["fields"][0]
        values=struct.unpack("<"+"h"*p,base64.b64decode(raw["trace_values_i16_le_base64"]))
        self.assertEqual(sum(vals),values[u])


if __name__=="__main__":
    unittest.main()

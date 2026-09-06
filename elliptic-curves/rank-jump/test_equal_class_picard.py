import unittest
import base64
from pathlib import Path
import retrospective as r
import equal_class_picard as ex


class EqualClassPicard(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=r.read(ex.OUTPUT)
        cls.replay=r.read(r.OUT/"rank_jump_equal_class_verification_v1.json")

    def test_low_mismatch_is_the_zero_rank_witness(self):
        row=self.report["bounds"][1]
        self.assertEqual(row["zero_rank_witness_primes"],[17,31])
        reductions=[v for v in self.report["reductions"] if v["case"]==1]
        self.assertEqual([v["reduction_geometric_Picard_rank"] for v in reductions],[18,18])
        self.assertEqual([v["NS_discriminant_squareclass"] for v in reductions],[-1,-113])
        self.assertEqual(row["geometric_mixed_rank_interval"],[0,0])
        self.assertEqual(row["full_base_arithmetic_rank_interval"],[3,3])

    def test_matching_discriminants_leave_high_unknown(self):
        row=self.report["bounds"][0]
        self.assertIsNone(row["zero_rank_witness_primes"])
        self.assertEqual(row["geometric_mixed_rank_interval"],[0,1])
        self.assertEqual(row["arithmetic_mixed_rank_interval"],[0,1])
        self.assertEqual(row["full_base_arithmetic_rank_interval"],[3,4])
        for row in self.report["bounds"]:
            self.assertEqual(row["production_curve_rank"],"UNKNOWN")

    def test_counts_sources_and_independent_coverage(self):
        replay=self.replay
        self.assertEqual(replay["status"],"PASS")
        self.assertEqual(self.report["bindings"],ex.bindings())
        self.assertEqual(replay["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(replay["counts_sha256"],r.digest(ex.RAW.read_bytes()))
        self.assertEqual(replay["verifier_sha256"],r.digest((ex.HERE/"verify_equal_class_picard.py").read_bytes()))
        new=reused=direct=0
        for record,raw in zip(replay["records"],r.read(ex.RAW)["records"]):
            self.assertEqual((record["case"],record["p"]),(raw["case"],raw["p"]))
            for checked,saved in zip(record["fields"],raw["fields"]):
                self.assertEqual(checked["trace_sha256"],r.digest(base64.b64decode(saved["trace_values_i16_le_base64"])))
                if checked["reused"]:reused+=checked["base_parameters"]
                else:new+=checked["base_parameters"]
                direct+=checked["direct_orbit_character_sums"]
        self.assertEqual((new,reused,direct),(30783,12817,10478))

    def test_smooth_infinity_is_included(self):
        for record in self.report["reductions"]:
            self.assertEqual(record["traces"],[v["finite_sum"]+v["infinity_sum"] for v in record["fields"]])
            self.assertTrue(any(v["infinity_sum"] for v in record["fields"]))


if __name__=="__main__":
    unittest.main()

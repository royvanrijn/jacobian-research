import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))
import wide_arithmetic_profile_core as core


class ProfileCoreTests(unittest.TestCase):
    def test_extract_ainvs_shapes(self):
        self.assertEqual(core.extract_ainvs({"ainvs":[1,2,3,4,6]}), [1,2,3,4,6])
        self.assertEqual(core.extract_ainvs({"curve":{"A":-7,"B":11}}), [0,0,0,-7,11])
        self.assertEqual(core.extract_ainvs({"model":{"a1":1,"a2":0,"a3":1,"a4":-2,"a6":3}}), [1,0,1,-2,3])
        self.assertIsNone(core.extract_ainvs({"equation":"y^2=x^3-x"}))

    def test_normalize_nested_row(self):
        row = core.normalize_row({
            "live_case_id":"b-test",
            "parameter":{"numerator":3,"denominator":5},
            "result":{"curve":{"A":-7,"B":11}, "rank_lower_bound":"≥ 23"},
            "metadata":{"initial_stage_rank_lower_bound":21,"generic_rank":17},
        }, "x.json", 4)
        self.assertEqual(row["id"], "b-test")
        self.assertEqual(row["t"], "3/5")
        self.assertEqual(row["ainvs"], [0,0,0,-7,11])
        self.assertEqual(row["declared_initial_stage_rank_lower_bound"], 21)
        self.assertEqual(row["generic_rank_lower_bound"], 17)
        self.assertEqual(row["final_rank_lower_bound"], 23)

    def test_aggregate_continuations(self):
        base = {"ainvs":[0,0,0,-7,11],"t":"3/5","source":"a","source_ordinal":0,"id":"b"}
        rows = [
            dict(base, declared_initial_stage_rank_lower_bound=20, generic_rank_lower_bound=17, stage="initial", final_rank_lower_bound=22),
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage="continuation-1", final_rank_lower_bound=23, source="c"),
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage="continuation-2", final_rank_lower_bound=24, source="d"),
        ]
        got = core.aggregate_rows(rows)
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0]["initial_rank_lower_bound"], 20)
        self.assertEqual(got[0]["final_rank_lower_bound"], 24)
        self.assertTrue(got[0]["improved_since_initial"])
        self.assertEqual(got[0]["observed_rank_lower_bounds"], [22,23,24])

    def test_aggregate_earliest_observed_label(self):
        base = {"ainvs":[0,0,0,-9,13],"t":None,"source":"a","source_ordinal":0,"id":"x"}
        got = core.aggregate_rows([
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage=None, final_rank_lower_bound=19),
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage=None, final_rank_lower_bound=21, source="b"),
        ])[0]
        self.assertEqual(got["initial_rank_lower_bound"], 19)
        self.assertEqual(got["initial_rank_source"], "earliest_observed_rank")

    def test_brumer_kramer_term(self):
        self.assertEqual(core.brumer_kramer_local_term(1, 4, [3]), {"u":2,"n":6,"local_term":8})
        self.assertEqual(core.brumer_kramer_local_term(-1, 1, [1,2]), {"u":1,"n":2,"local_term":3})
        with self.assertRaises(ValueError):
            core.brumer_kramer_local_term(1, 0, [4])

    def test_rank_buckets(self):
        self.assertEqual(core.rank_bucket(17), "17-18")
        self.assertEqual(core.rank_bucket(22), "21-22")
        self.assertEqual(core.rank_bucket(23), "23")
        self.assertEqual(core.rank_bucket(24), "24")
        self.assertEqual(core.rank_bucket(25), "25+")

    def test_spearman_ties(self):
        self.assertAlmostEqual(core.spearman([1,2,3],[10,20,30]), 1.0)
        self.assertAlmostEqual(core.spearman([1,2,3],[30,20,10]), -1.0)
        self.assertIsNone(core.spearman([1,1,1],[1,2,3]))

    def test_deterministic_controls(self):
        rows=[]
        for i, r in enumerate([17,17,18,19,20,21,22,22,23,24,25]):
            rows.append({"curve_key":f"k{i}","final_rank_lower_bound":r})
        a=core.deterministic_controls(rows, high_threshold=23, per_bucket=1, salt="s")
        b=core.deterministic_controls(list(reversed(rows)), high_threshold=23, per_bucket=1, salt="s")
        self.assertEqual(a,b)
        self.assertEqual(len(a["high"]),3)
        self.assertEqual(len(a["controls"]),3)

    def test_flatten_forced_g_boundary(self):
        row={"curve_key":"x","id":"x","t":None,"initial_rank_lower_bound":22,"final_rank_lower_bound":25,"improved_since_initial":True}
        base={"status":"PASS","delta_sign":1,"log2_abs_minimal_discriminant":100,"log2_abs_polynomial_discriminant":90,"rational_2torsion_rank":0,"cubic_irreducible":True}
        local={"status":"PASS","bk_local_term":9,"bk_u":2,"bk_n":7,"phi_m_count":4,"phi_a_count":1,"field_signature":[3,0],"field_discriminant":"123","log2_abs_field_discriminant":7,"field_ramified_prime_count":4,"root_number":-1,"conductor":"999","log2_conductor":10}
        got=core.flatten_profile(row,base,local,None)
        self.assertEqual(got["forced_class_2rank_lower_from_known_rank"],16)
        self.assertEqual(got["g_upper_needed_to_close_known_rank"],16)

    def test_summary_never_calls_lower_bound_exact(self):
        flat=[
            {"curve_key":"a","id":"a","final_rank_lower_bound":23,"initial_rank_lower_bound":22,"improved_since_initial":True,"rank_bucket":"23","bk_local_term":5},
            {"curve_key":"b","id":"b","final_rank_lower_bound":17,"initial_rank_lower_bound":17,"improved_since_initial":False,"rank_bucket":"17-18","bk_local_term":2},
            {"curve_key":"c","id":"c","final_rank_lower_bound":18,"initial_rank_lower_bound":17,"improved_since_initial":True,"rank_bucket":"17-18","bk_local_term":3},
        ]
        s=core.summarize_profiles(flat)
        self.assertIn("lower bounds", s["interpretation_boundary"])
        md=core.markdown_summary(s,flat)
        self.assertIn("certified lower bounds, not exact ranks",md)
        self.assertNotIn("exact rank distribution",md.lower())

    def test_read_json_nested_and_dedupe(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"results.json"
            p.write_text(json.dumps({"results":[
                {"case_id":"a","curve":{"ainvs":[0,0,0,-7,11]},"rank_lower_bound":22},
                {"case_id":"a","curve":{"ainvs":[0,0,0,-7,11]},"rank_lower_bound":23},
            ]}))
            rows=core.read_rows_from_file(p)
            self.assertGreaterEqual(len(rows),2)
            agg=core.aggregate_rows(rows)
            self.assertEqual(len(agg),1)
            self.assertEqual(agg[0]["final_rank_lower_bound"],23)


if __name__ == "__main__":
    unittest.main()

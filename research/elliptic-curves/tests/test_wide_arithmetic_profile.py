import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import subprocess

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))
import wide_arithmetic_profile_core as core
import run_wide_arithmetic_profile as controller


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

    def test_aggregate_unlabelled_observations_do_not_bind_initial(self):
        base = {"ainvs":[0,0,0,-9,13],"t":None,"source":"a","source_ordinal":0,"id":"x"}
        got = core.aggregate_rows([
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage=None, final_rank_lower_bound=19),
            dict(base, declared_initial_stage_rank_lower_bound=None, generic_rank_lower_bound=17, stage=None, final_rank_lower_bound=21, source="b"),
        ])[0]
        self.assertIsNone(got["initial_rank_lower_bound"])
        self.assertIsNone(got["initial_rank_source"])
        self.assertIsNone(got["improved_since_initial"])
        self.assertEqual(got["observed_rank_lower_bounds"], [19,21])

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


class ControllerTests(unittest.TestCase):
    def campaign_fixture(self, root):
        def put(path,obj):
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(core.stable_json(obj))
            return core.sha256_file(path)
        model=[0,0,1,-1,0];item={'id':'b-test','model':model,'parameter':'1/2','backend':'native'}
        qh=put(root/'queue.json',{'rows':[item]})
        ph=put(root/'plan.json',{});mh=put(root/'manifest.json',{})
        final=None
        for n,rank in [(0,18),(1,19)]:
            folder=root/'runtime/research/broad-cases/b-test'/f'batch-{n:03d}'
            packet=folder/'packet.json';pkh=put(packet,{'curve':model,'rank_lower_bound':rank})
            put(folder/'packet-verified.json',{'status':'PASS_TWO_FINITE_IMPLEMENTATIONS','packet_sha256':pkh})
            sh=put(folder/'broad-state.json',{'id':'b-test','round':n,'status':'CERTIFIED','rank':rank,
                'packet':str(packet.relative_to(root/'runtime/research')),'packet_sha256':pkh})
            final={'id':'b-test','final_round':n,'state_sha256':sh,'packet_sha256':pkh}
        put(root/'COMPLETE.json',{'status':'COMPLETE_BOUNDED_CAMPAIGN','queue_sha256':qh})
        put(root/'COMPLETION_REVIEW.json',{'status':'COMPLETE_BOUNDED_CAMPAIGN','queue_sha256':qh,
            'plan_sha256':ph,'manifest_sha256':mh,'unknown_or_censored':0,'still_eligible':0,
            'final_fibres_checked':1,'final_certificate_bindings':[final],
            'initial_lower_bound_histogram':{'18':1},'final_lower_bound_histogram':{'19':1},
            'fibres_improved_in_continuations':1,'sum_of_certified_lower_bound_increases':1})

    def test_completed_campaign_binds_actual_initial_endpoint(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.campaign_fixture(root)
            _,rows,_=core.read_completed_campaign(root)
            row=core.aggregate_rows(rows)[0]
            self.assertEqual(row['initial_rank_lower_bound'],18)
            self.assertEqual(row['generic_rank_lower_bound'],17)
            self.assertTrue(row['improved_since_initial'])

    def test_completed_campaign_rejects_packet_tampering(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.campaign_fixture(root)
            p=root/'runtime/research/broad-cases/b-test/batch-000/packet.json'
            obj=json.loads(p.read_text());obj['rank_lower_bound']=17;p.write_text(core.stable_json(obj))
            with self.assertRaisesRegex(RuntimeError,'hash mismatch'):core.read_completed_campaign(root)

    def test_completed_campaign_rejects_replay_failure(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.campaign_fixture(root)
            p=root/'runtime/research/broad-cases/b-test/batch-000/packet-verified.json'
            obj=json.loads(p.read_text());obj['status']='UNKNOWN';p.write_text(core.stable_json(obj))
            with self.assertRaisesRegex(RuntimeError,'replay not bound'):core.read_completed_campaign(root)

    def fixture(self, out):
        for name in ['inputs','base','local','logs']:(out/name).mkdir()
        row={'curve_key':'abc','id':'test','ainvs':[0,0,1,-1,0]}
        (out/'inputs/abc.json').write_text(core.stable_json(row))
        return row

    def test_worker_command_memory_and_resume(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);row=self.fixture(out)
            class Child:
                returncode=0
                def communicate(self, timeout):
                    (out/'base/abc.tmp.json').write_text(core.stable_json({'curve_key':'abc','status':'PASS'}))
                    return 'worker output',''
            with patch.object(controller.subprocess,'Popen',return_value=Child()) as p:
                self.assertEqual(controller.run_one(out,row,'base',60,2.5,'sage'),'PASS')
                self.assertEqual(p.call_args.args[0][-1],'2.5')
                self.assertTrue(p.call_args.kwargs['start_new_session'])
                self.assertEqual(controller.run_one(out,row,'base',60,2.5,'sage'),'PASS')
                self.assertEqual(p.call_count,1)
            self.assertEqual((out/'logs/base-abc.txt').read_text(),'worker output')

    def test_timeout_checkpoints_unknown_and_kills_process_group(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);row=self.fixture(out)
            class Child:
                pid=123456789
                def communicate(self, timeout=None):
                    if timeout is not None:raise subprocess.TimeoutExpired('sage',timeout)
                    return '', ''
            with patch.object(controller.subprocess,'Popen',return_value=Child()),patch.object(controller.os,'killpg') as kill:
                self.assertEqual(controller.run_one(out,row,'base',1,0,'sage'),'UNKNOWN_TIMEOUT')
                kill.assert_called_once()
            result=json.loads((out/'base/abc.json').read_text())
            self.assertEqual(result['input_sha256'],core.sha256_file(out/'inputs/abc.json'))

    def test_wrong_curve_output_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td);row=self.fixture(out)
            class Child:
                returncode=0
                def communicate(self,timeout):
                    (out/'base/abc.tmp.json').write_text(core.stable_json({'curve_key':'wrong','status':'PASS'}))
                    return '', ''
            with patch.object(controller.subprocess,'Popen',return_value=Child()):
                self.assertEqual(controller.run_one(out,row,'base',1,0,'sage'),'UNKNOWN_CONTROLLER_FAILURE')

    def test_missing_completion_receipt_is_not_inferred(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(FileNotFoundError):core.read_completed_campaign(Path(td))


if __name__ == "__main__":
    unittest.main()

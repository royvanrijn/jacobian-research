"""Deterministic contracts/controller tests. Mocked workers prove NO mathematics."""
import copy
from dataclasses import dataclass
from fractions import Fraction as F
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import det1092_v3_contract as c
import run_det1092_v3_trial as runner
import v3_warm_support as s


def population():
    rows = []
    for band, power in ((10,214),(11,237)):
        for k in range(24):
            i = len(rows)
            model = ['0','0','0',str(-(2**power+2*k+1)),str(2*k+3)]
            j = c.j_invariant(model)
            rows.append({'id':f'scale-{i:07d}','family':'det1092-reduced','parameter':str(F(i+1,101)),
              'model':model,'height_bin':band,'stratum':'strong' if k < 16 else 'moderate' if k < 20 else 'lower_fixed',
              'score_units':10000-100*k-band,'j_numerator_bits':abs(j.numerator).bit_length(),
              'j_denominator_bits':j.denominator.bit_length()})
    return {'status':'PASS','selected':rows}


def small_terminal(folder, count=2):
    folder.mkdir(parents=True,exist_ok=True)
    wd = folder/'replay-M17/epoch-00'; wd.mkdir(parents=True)
    s.atomic(folder/'protocol.json',{'initial_rank':17})
    centres = [{'point':[str(i),'2']} for i in range(count)]
    sel = {'rank':17,'basis':[['1','2']]*17,'centres':centres,'full_cosets_scored':32}
    s.atomic(wd/'selection.json',sel)
    for i,centre in enumerate(centres):
        s.atomic(wd/f'chart-{i:03d}.json',{'index':i,'centre':centre,'search':{'status':'bounded_search_complete'}})
    audit = wd/f'mod2-{count-1:03d}.json'; s.atomic(audit,{'fixture':'not an arithmetic certificate'})
    stage = {'epoch':0,'before':17,'after':17,'charts':count,'stale_charts_cancelled':0,
             'audit':audit.name,'audit_sha256':s.sha(audit),'full_cosets_scored':32,
             'stop_reason':'COMPLETE_FINITE_NO_GAIN','censored_charts':0}
    s.atomic(wd/'stage.json',stage)
    terminal = {'status':'TERMINAL_BOUNDED_TRANSFER','initial_rank':17,'final_rank_lower_bound':17,
                'protocol_sha256':s.sha(folder/'protocol.json'),'stages':[stage],'charts':count,
                'stop_reason':'COMPLETE_FINITE_NO_GAIN'}
    s.atomic(folder/'replay-M17/terminal.json',terminal)
    return terminal


class SelectionTests(unittest.TestCase):
    def setUp(self): self.data = population()

    def test_sizes_roles_and_disjoint_reserve(self):
        pilot,reserve = c.choose(self.data)
        self.assertEqual((len(pilot),len(reserve)),(8,40))
        self.assertEqual([p['pilot_role'] for p in pilot],['strong']*4+['middle']*2+['sha-control']*2)
        self.assertEqual(len({p['id'] for p in pilot+reserve}),48)

    def test_row_order_invariance(self):
        old = c.choose(self.data); self.data['selected'].reverse()
        self.assertEqual(c.choose(self.data),old)

    def test_outcome_annotations_cannot_select(self):
        old = c.choose(self.data)
        for i,row in enumerate(self.data['selected']):
            row.update(known_rank=100-i,points=[[i,99]],prior_gains=i%2,exceptional_label='ignored')
        self.assertEqual(c.choose(self.data),old)

    def test_four_strongest_frozen_scores(self):
        pilot,_ = c.choose(self.data)
        self.assertEqual([r['id'] for r in pilot[:4]], [r['id'] for r in sorted(self.data['selected'],key=c.score_key)[:4]])

    def test_middle_and_sha_cover_both_bands(self):
        pilot,_ = c.choose(self.data)
        self.assertEqual([r['height_bin'] for r in pilot[4:]], [10,11,10,11])

    def test_missing_case_no_refill(self):
        self.data['selected'].pop()
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_unfinished_intake(self):
        self.data['status'] = 'RUNNING'
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_duplicate_id(self):
        self.data['selected'][1]['id'] = self.data['selected'][0]['id']
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_duplicate_parameter(self):
        self.data['selected'][1]['parameter'] = self.data['selected'][0]['parameter']
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_duplicate_curve(self):
        for k in ('model','j_numerator_bits','j_denominator_bits'):
            self.data['selected'][1][k] = self.data['selected'][0][k]
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_path_injection(self):
        self.data['selected'][0]['id'] = '../../other'
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_floating_score_rejected(self):
        self.data['selected'][0]['score_units'] = 1.0
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_j_height_rechecked(self):
        self.data['selected'][0]['j_numerator_bits'] += 1
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_strata_not_refilled(self):
        self.data['selected'][0]['stratum'] = 'moderate'
        with self.assertRaises(ValueError): c.choose(self.data)

    def test_policy_numerical_fields_unchanged(self):
        old = dict(anchors_per_shell=16,canonical_per_shell=25,exact_cvp_node_limit=900,
                   height=125000,seconds_per_chart=10,gp_sha256='abc',max_epochs=16,max_charts=4096,
                   selection='frozen',target_rank=31,sources={'numerical':'fixed'})
        new = c.trial_policy(old)
        for k,v in old.items():
            if k != 'target_rank': self.assertEqual(new[k],v)
        self.assertEqual(new['target_rank'],32)
        self.assertEqual(old['target_rank'],31)

    def test_missing_policy_field(self):
        with self.assertRaises(ValueError): c.trial_policy({})

    def test_no_gain_not_expansion(self):
        rows = [{'gain':0,'stop_reason':'COMPLETE_FINITE_NO_GAIN'}]*8
        r = c.result_summary(rows)
        self.assertTrue(r['all_eight_uncensored_terminal']); self.assertFalse(r['expansion_eligible'])

    def test_censored_not_clean_even_with_gain(self):
        rows = [{'gain':0,'stop_reason':'COMPLETE_FINITE_NO_GAIN'}]*7 + [{'gain':1,'stop_reason':'CENSORED_SEARCH'}]
        self.assertFalse(c.result_summary(rows)['expansion_eligible'])

    def test_gain_eligible_but_not_automatic(self):
        rows = [{'gain':1,'stop_reason':'COMPLETE_FINITE_NO_GAIN'}]*8
        r = c.result_summary(rows)
        self.assertTrue(r['expansion_eligible']); self.assertFalse(r['automatic_expansion'])

    def test_partial_roster_not_closed(self):
        self.assertFalse(c.result_summary([{'gain':1,'stop_reason':'COMPLETE_FINITE_NO_GAIN'}])['all_eight_uncensored_terminal'])


    def test_packed_selection_fallback_preserves_raw_hash(self):
        import gzip,hashlib
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); raw=json.dumps(self.data).encode(); packed=root/'selection.json.gz'
            packed.write_bytes(gzip.compress(raw))
            with patch.object(c,'SELECTION',root/'missing.json'),patch.object(c,'PACKED_SELECTION',packed),patch.object(c,'SELECTION_SHA256',hashlib.sha256(raw).hexdigest()):
                actual,path=c.selection_input()
                self.assertEqual(actual,self.data); self.assertEqual(path,packed)

    def test_bad_local_input_not_silently_replaced_by_pack(self):
        import gzip,hashlib
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); raw=json.dumps(self.data).encode(); packed=root/'selection.json.gz'
            packed.write_bytes(gzip.compress(raw)); local=root/'selection.json'; local.write_text('{}')
            with patch.object(c,'SELECTION',local),patch.object(c,'PACKED_SELECTION',packed),patch.object(c,'SELECTION_SHA256',hashlib.sha256(raw).hexdigest()):
                with self.assertRaises(ValueError): c.selection_input()

    def test_exact_git_blob_hash(self):
        import hashlib
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'blob'; p.write_bytes(b'abc')
            self.assertEqual(c.git_blob(p),hashlib.sha1(b'blob 3\0abc').hexdigest())


class RationalAndCheckpointTests(unittest.TestCase):
    def test_exact_homogeneous_specialization(self):
        rat = lambda a:{'numerator':a,'denominator':['1']}
        parent = {'a_invariants':[rat([]),rat([]),rat([]),rat([-2]),rat([1])],
                  'basis_weierstrass_coordinates':[[rat([1]),rat([0])]]*17}
        row = {'parameter':'2/3','model':[0,0,0,-2*3**8,3**12]}
        model,points = c.specialize(parent,row)
        self.assertEqual(points[0],(F(3**4),F(0)))
        self.assertTrue(c.on_curve(model,points[0]))
        # These deliberately dependent fixture points are membership tests only.

    def test_wrong_equation_rejected(self):
        rat = lambda a:{'numerator':a,'denominator':['1']}
        parent = {'a_invariants':[rat([])]*3+[rat([-2]),rat([1])],
                  'basis_weierstrass_coordinates':[[rat([1]),rat([0])]]*17}
        with self.assertRaises(ValueError): c.specialize(parent,{'parameter':'2/3','model':[0,0,0,-2,1]})

    def test_pole_rejected(self):
        with self.assertRaises(ValueError): c.evaluate({'numerator':[1],'denominator':[-1,1]},F(1))

    def test_exact_value_normalization_preserves_sign_order(self):
        self.assertEqual(s.point_tuple([['2/4','-3']]),((F(1,2),F(-3)),))
        self.assertNotEqual(s.point_tuple([['1','2'],['3','4']]),s.point_tuple([['3','4'],['1','2']]))

    def test_1508_numeric_order(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            for i in range(1508): (p/f'chart-{i:03d}.json').write_text('{}')
            paths = s.indexed_paths(p,expected=1508)
            self.assertEqual(paths[101].name,'chart-101.json'); self.assertEqual(paths[1000].name,'chart-1000.json')

    def test_alias_index_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); (p/'chart-0.json').write_text('{}'); (p/'chart-000.json').write_text('{}')
            with self.assertRaises(ValueError): s.indexed_paths(p)

    def test_mutated_centre_rejected(self):
        with self.assertRaises(ValueError):
            s.check_chart({'index':0,'centre':{'a':1}},{'centres':[{'a':2}]},0)

    def test_immutable_writes(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'record.json'; s.atomic(p,{'a':1},immutable=True); s.atomic(p,{'a':1},immutable=True)
            with self.assertRaises(ValueError): s.atomic(p,{'a':2},immutable=True)

    def test_sealed_terminal_does_not_mean_certificate(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)/'case'; small_terminal(folder)
            self.assertIsNotNone(s.terminal_structure(folder))
            self.assertFalse((folder/'trial-verified.json').exists())

    def test_partial_schedule_not_complete_negative(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)/'case'; small_terminal(folder)
            p = folder/'replay-M17/epoch-00/selection.json'; data = s.read(p); data['centres'].append({'point':['99','2']}); s.atomic(p,data)
            with self.assertRaises(ValueError): s.terminal_structure(folder)


    def test_preflight_skips_zero_orbit_and_wrong_shell(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'orbits.tsv'
            p.write_text('orbit_mask\tminimum_norm\tparent_MW17_w\n0\t0\t'+'0 '*17+'\n1\t4\t'+'1 '*17+'\n2\t8\t'+'1 '+'0 '*16+'\n')
            self.assertEqual(c.preflight_word(p,17),[1]+[0]*16)

    def test_preflight_pads_extra_directions(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'orbits.tsv'
            p.write_text('orbit_mask\tminimum_norm\tparent_MW17_w\n2\t10\t'+'1 '+'0 '*16+'\n')
            self.assertEqual(c.preflight_word(p,27),[1]+[0]*26)

    def test_preflight_rejects_table_without_allowed_shell(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'orbits.tsv'; p.write_text('orbit_mask\tminimum_norm\tparent_MW17_w\n0\t0\t'+'0 '*17+'\n')
            with self.assertRaises(ValueError): c.preflight_word(p,17)


class ControllerTests(unittest.TestCase):
    def test_every_numeric_action_uses_sage(self):
        for action in ('prepare','preflight','search','replay'):
            argv = runner.worker_command('/fake/sage',action,None if action == 'prepare' else 'scale-0000000',Path('/session.json'))
            self.assertEqual(argv[:3],['/fake/sage','-python','-u'])
            self.assertIn(action,argv)

    def test_invalid_worker_action(self):
        with self.assertRaises(ValueError): runner.worker_command('/sage','expand',None,'x')

    def test_status_standard_python_no_sage_import(self):
        script = 'import sys; import run_det1092_v3_trial; print("SAGE_FREE", not any(k == "sage" or k.startswith("sage.") for k in sys.modules))'
        result = subprocess.run([sys.executable,'-c',script],cwd=CAS,capture_output=True,text=True,timeout=10)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('SAGE_FREE True',result.stdout)

    def test_budget_validation_before_launch(self):
        for hours in (0,-1,float('inf'),float('nan'),25):
            with self.assertRaises(ValueError): runner.launch(hours)

    def test_zombie_and_reused_pid(self):
        with patch.object(s,'process_info',return_value={'pid':99,'state':'Z','start_token':'1'}):
            self.assertFalse(runner.alive({'pid':99,'start_token':'1'}))
        with patch.object(s,'process_info',return_value={'pid':99,'state':'S','start_token':'2'}):
            self.assertFalse(runner.alive({'pid':99,'start_token':'1'}))

    def test_live_lifetime_lock(self):
        with tempfile.TemporaryDirectory() as d, patch.object(c,'AUTO',Path(d)),patch.object(runner,'LOCK',Path(d)/'lock'):
            fd = runner.take_lock()
            try:
                with self.assertRaises(RuntimeError): runner.take_lock()
            finally: os.close(fd)
            fd = runner.take_lock(); os.close(fd)

    def test_supervised_backend_error_is_not_no_gain(self):
        with tempfile.TemporaryDirectory() as d:
            module = types.ModuleType('research_runtime.supervisor')
            module.Limits = lambda *args: args
            module.run = lambda *args,**kwargs:{'outcome':'backend_failure'}
            report = types.SimpleNamespace(update=lambda *a,**k:None)
            root = Path(d); session = root/'session.json'
            with patch.dict(sys.modules,{'research_runtime.supervisor':module}),patch.object(c,'ROOT',root),patch.object(c,'D',root/'data'):
                with self.assertRaisesRegex(RuntimeError,'backend_failure'):
                    runner.run_step('search','case',session,{'sage':'/sage'},runner.time.monotonic()+60,report)

    def test_resource_stop_not_negative(self):
        with tempfile.TemporaryDirectory() as d:
            module = types.ModuleType('research_runtime.supervisor'); module.Limits = lambda *a:a
            module.run = lambda *a,**k:{'outcome':'strict_rss_limit'}
            root=Path(d); report=types.SimpleNamespace(update=lambda *a,**k:None)
            with patch.dict(sys.modules,{'research_runtime.supervisor':module}),patch.object(c,'ROOT',root),patch.object(c,'D',root/'data'):
                with self.assertRaises(runner.ResourceStop):
                    runner.run_step('search','case',root/'session.json',{'sage':'/sage'},runner.time.monotonic()+60,report)

    def test_sealed_search_routes_to_replay_after_abnormal_exit(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); directory=root/'data'; small_terminal(directory/'case')
            module=types.ModuleType('research_runtime.supervisor'); module.Limits=lambda *a:a
            module.run=lambda *a,**k:{'outcome':'backend_failure'}
            events=[]; report=types.SimpleNamespace(update=lambda *a,**k:events.append(a[0]))
            with patch.dict(sys.modules,{'research_runtime.supervisor':module}),patch.object(c,'ROOT',root),patch.object(c,'D',directory):
                runner.run_step('search','case',root/'session.json',{'sage':'/sage'},runner.time.monotonic()+60,report)
            self.assertIn('SEALED_SEARCH_REUSED',events)
            self.assertFalse((directory/'case/trial-verified.json').exists())


    def test_controller_completes_eight_not_forty_reserves(self):
        # Mocked arithmetic: verifies dispatch/order only, not any rank statement.
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); directory=root/'data'; directory.mkdir()
            s.atomic(directory/'roster.json',{'fixture':True})
            session=root/'session.json'; s.atomic(session,{'sources':{},'inputs':{},'hours':1,'sage':'/sage'})
            lock=root/'lock'; fd=os.open(lock,os.O_CREAT|os.O_RDWR,0o600)
            ids=[f'scale-{i:07d}' for i in range(8)]; results={}; calls=[]; events=[]
            class Report:
                def __init__(self,*args): pass
                def __enter__(self): return self
                def __exit__(self,*args): pass
                def update(self,status=None,**values): events.append((status,values))
            def run_step(action,case,*args):
                calls.append((action,case))
                if action=='replay': results[case]={'case':case,'gain':0,'rank_lower_bound':17,'stop_reason':'COMPLETE_FINITE_NO_GAIN'}
            with patch.object(c,'ROOT',root),patch.object(c,'D',directory),patch.object(c,'own_sources',return_value={}),patch.object(c,'validate_roster',return_value={'pilot':[{'id':x} for x in ids],'reserve':[{}]*40}),patch.object(runner,'LOCK',lock),patch.object(runner,'bindings'),patch.object(runner,'Reporter',Report),patch.object(runner,'verified',side_effect=lambda case:results.get(case)),patch.object(runner,'run_step',side_effect=run_step),patch.object(runner.signal,'signal'):
                runner.worker(session,fd)
            self.assertEqual(len(calls),24)
            self.assertEqual([case for action,case in calls if action=='search'],ids)
            self.assertEqual(events[-1][0],'COMPLETE_EIGHT_PILOT')
            self.assertEqual(len(results),8)


if __name__ == '__main__': unittest.main()

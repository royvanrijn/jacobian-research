"""Small deterministic regressions; no Sage installation or point-search campaign.

Mock arithmetic tests exercise orchestration ONLY, not mathematical correctness.
Native Sage preflight in the controller separately verifies each actual seed/map.
"""
from __future__ import annotations

import contextlib
from fractions import Fraction as F
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace, ModuleType
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import v3_warm_support as support
import v3_warm_engine as engine
import v3_warm_replay as replay
import run_v3_warm_start_overnight as controller


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def fake_points(rank=27):
    return tuple((F(i+1), F(2*i+1)) for i in range(rank))


MODEL = (F(0), F(0), F(0), F(-1), F(1))


class FakeState:
    """State-shaped fixture, NOT an elliptic-curve rank certificate."""
    def __init__(self, points):
        self.basis = tuple(tuple(map(str, p)) for p in points)
        self.rank = len(points)
        self.model = SimpleNamespace(coefficients=MODEL)
        self.key = 'mock-' + str(self.rank)
    def record(self):
        return {'key': self.key, 'state': {'reductions': {'points': [list(p) for p in self.basis]}}}


def terminal_fixture(folder, count=1508):
    policy = {'initial_rank': 27}
    put(folder/'protocol.json', policy)
    wd = folder/'replay-M17/epoch-00'
    centres = [{'point': [str(i), '1'], 'representative': [i]} for i in range(count)]
    selection = {'rank': 27, 'basis': [list(map(str, p)) for p in fake_points()],
                 'centres': centres, 'full_cosets_scored': 32768}
    put(wd/'selection.json', selection)
    for i, centre in enumerate(centres):
        put(wd/f'chart-{i:03d}.json', {'index': i, 'centre': centre,
                                     'search': {'status': 'bounded_search_complete'}})
    audit = wd/f'mod2-{count-1:03d}.json'
    put(audit, {'fixture': 'structural only'})
    stage = {'epoch': 0, 'before': 27, 'after': 27, 'charts': count,
             'stale_charts_cancelled': 0, 'audit': audit.name, 'audit_sha256': support.sha(audit),
             'full_cosets_scored': 32768, 'stop_reason': 'COMPLETE_FINITE_NO_GAIN', 'censored_charts': 0}
    put(wd/'stage.json', stage)
    terminal = {'status': 'TERMINAL_BOUNDED_TRANSFER', 'initial_rank': 27,
                'final_rank_lower_bound': 27, 'charts': count, 'stages': [stage],
                'protocol_sha256': support.sha(folder/'protocol.json'), 'stop_reason': 'COMPLETE_FINITE_NO_GAIN'}
    put(folder/'replay-M17/terminal.json', terminal)
    return wd, selection


class IndexAndBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
    def tearDown(self):
        self.temp.cleanup()

    def test_original_1508_chart_sort_failure_is_reproduced(self):
        wd, selection = terminal_fixture(self.root)
        lexical = sorted(wd.glob('chart-*.json'))
        mismatch = next(i for i, p in enumerate(lexical) if support.read(p)['index'] != i)
        self.assertEqual(mismatch, 101)
        self.assertEqual(lexical[mismatch].name, 'chart-1000.json')
        with self.assertRaisesRegex(ValueError, 'index mismatch'):
            support.check_chart(support.read(lexical[mismatch]), selection, mismatch)
        numeric = support.indexed_paths(wd, expected=1508)
        for i, path in enumerate(numeric):
            support.check_chart(support.read(path), selection, i)
        self.assertEqual(support.terminal_structure(self.root)['charts'], 1508)

    def test_numeric_boundaries_82_1000_1001_1508_4096(self):
        for count in (82, 1000, 1001, 1508, 4096):
            folder = self.root/str(count)
            folder.mkdir()
            for i in range(count):
                (folder/f'chart-{i:03d}.json').write_text('{}')
            paths = support.indexed_paths(folder, expected=count)
            self.assertEqual([int(p.stem.split('-')[1]) for p in paths], list(range(count)))

    def test_missing_index_rejected(self):
        put(self.root/'chart-000.json', {})
        put(self.root/'chart-002.json', {})
        with self.assertRaisesRegex(ValueError, 'noncontiguous'):
            support.indexed_paths(self.root)

    def test_alias_index_rejected(self):
        put(self.root/'chart-000.json', {})
        put(self.root/'chart-0000.json', {})
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            support.indexed_paths(self.root)

    def test_payload_index_rejected(self):
        with self.assertRaisesRegex(ValueError, 'index mismatch'):
            support.check_chart({'index': 1, 'centre': {}}, {'centres': [{}]}, 0)

    def test_wrong_centre_still_rejected(self):
        with self.assertRaisesRegex(ValueError, 'centre differs'):
            support.check_chart({'index': 0, 'centre': {'x': 1}}, {'centres': [{'x': 2}]}, 0)

    def test_boolean_index_rejected(self):
        with self.assertRaises(ValueError):
            support.check_chart({'index': False, 'centre': {}}, {'centres': [{}]}, 0)

    def test_terminal_audit_tamper_rejected(self):
        wd, _ = terminal_fixture(self.root, 3)
        put(wd/'mod2-002.json', {'tampered': True})
        with self.assertRaisesRegex(ValueError, 'audit changed'):
            support.terminal_structure(self.root)

    def test_terminal_missing_last_chart_rejected(self):
        wd, _ = terminal_fixture(self.root, 3)
        (wd/'chart-002.json').unlink()
        with self.assertRaisesRegex(ValueError, 'count'):
            support.terminal_structure(self.root)

    def test_censored_chart_cannot_be_clean_no_gain(self):
        wd, _ = terminal_fixture(self.root, 3)
        chart = support.read(wd/'chart-002.json')
        chart['search']['status'] = 'bounded_search_timeout'
        put(wd/'chart-002.json', chart)
        with self.assertRaisesRegex(ValueError, 'censored'):
            support.terminal_structure(self.root)

    def test_rational_strings_are_not_fraction_tuples(self):
        points = ((F(3, 7), F(-5, 11)),)
        state = FakeState(points)
        self.assertNotEqual(tuple(state.basis), points)  # original failing assertion
        support.assert_basis(state, MODEL, points, 1)

    def test_real_point_change_not_normalized_away(self):
        state = FakeState(((F(3, 7), F(-5, 11)),))
        with self.assertRaisesRegex(ValueError, 'coordinates differ'):
            support.assert_basis(state, MODEL, ((F(3, 7), F(5, 11)),), 1)

    def test_rank_and_curve_mismatch_rejected(self):
        state = FakeState(fake_points())
        with self.assertRaisesRegex(ValueError, 'rank'):
            support.assert_basis(state, MODEL, fake_points(), 28)
        with self.assertRaisesRegex(ValueError, 'model'):
            support.assert_basis(state, (0,0,0,-2,1), fake_points(), 27)

    def test_atomic_immutable_no_clobber(self):
        p = self.root/'record.json'
        support.atomic(p, {'v': 1}, immutable=True)
        digest = support.sha(p)
        support.atomic(p, {'v': 1}, immutable=True)
        with self.assertRaises(ValueError):
            support.atomic(p, {'v': 2}, immutable=True)
        self.assertEqual(support.sha(p), digest)

    def test_path_escape_rejected(self):
        with self.assertRaises(ValueError):
            support.within(self.root, '../escape')

    def test_bounded_log_tail(self):
        p = self.root/'log'
        p.write_bytes(b'a'*100000+b'last')
        result = support.tail(p, 100)
        self.assertEqual(len(result), 100)
        self.assertTrue(result.endswith('last'))


class SearchRestartTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name)
        self.executed = []
        self.calls = []
        self.policy = {'initial_rank': 27, 'max_charts': 10, 'max_epochs': 3,
                       'height': 125000, 'seconds_per_chart': 10, 'gp_sha256': 'mock', 'target_rank': 32}
        put(self.folder/'protocol.json', self.policy)
        self.gain_at = None
        self.mapper = SimpleNamespace(pari=SimpleNamespace(allocatemem=lambda *a, **k: None),
            mapping=lambda model,basis,c: {'centre': c, 'coordinate_policy': {'kind': 'raw'}})
        self.numeric = SimpleNamespace(v1=SimpleNamespace(READS=set()), load=lambda *a: self.mapper,
                                       landscape=self.landscape)
        self.ctx = engine.Context('warm-11952-41', self.folder, self.folder, self.policy,
                                  self.numeric, MODEL, FakeState(fake_points()), {})
        backend = ModuleType('pari_pointed_backend')
        backend.execute, backend.replay = self.execute, lambda *a: ()
        pointed = ModuleType('pointed_quartic_search')
        pointed.PointedQuarticSearch = lambda **kw: SimpleNamespace(**kw)
        self.stack = contextlib.ExitStack()
        self.stack.enter_context(patch.dict(sys.modules, {'pari_pointed_backend': backend, 'pointed_quartic_search': pointed}))
        self.stack.enter_context(patch.object(engine, '_audit', side_effect=self.audit))
        self.stack.enter_context(patch.object(engine, '_odd', side_effect=self.odd))
        self.stack.enter_context(patch.object(engine, 'certified_state', side_effect=lambda m,p,proof:FakeState(p)))
        self.stack.enter_context(patch.object(engine, 'restore_state', side_effect=lambda r,m,p:FakeState(p)))
        self.stack.enter_context(patch.object(replay, 'verify_landscape', return_value={}))
        self.stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
    def tearDown(self):
        self.stack.close()
        self.temp.cleanup()
    def landscape(self, model, basis, tested, wd, policy):
        centres = [{'representative': [j]+[0]*(len(basis)-1), 'point': [str(j+len(basis)*10),'1']} for j in range(3)]
        put(wd/'selection.json', {'rank':len(basis), 'basis':[list(map(str,p)) for p in basis],
                                  'full_cosets_scored':32, 'centres':centres})
        return centres
    def execute(self, search, mapping, height, seconds, gp):
        self.executed.append(mapping['centre']['point'])
        return {'status':'bounded_search_complete','finite_curve_points':[]}, ()
    def audit(self, snapshot, path):
        snap=support.read(snapshot)
        basis=point_tuple_local(snap['final_state']['state']['reductions']['points'])
        points = basis
        if self.gain_at and len(basis) == 27 and len(snap['charts']) == self.gain_at:
            points = basis+((F(999),F(1000)),)
        result = {'independent_points':[list(map(str,p)) for p in points], 'rank_lower_bound':len(points),
                  'rank_certificate':{}, 'points':[list(map(str,p)) for p in points]}
        put(path,result)
        return result
    def odd(self, audit, path):
        rank=support.read(audit)['rank_lower_bound']
        put(path,{'fixture_rank':rank})
        return {'3':rank,'5':rank}

    def test_complete_then_resume_does_not_execute_again(self):
        engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),3)
        first=support.sha(self.folder/'replay-M17/terminal.json')
        engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),3)
        self.assertEqual(support.sha(self.folder/'replay-M17/terminal.json'),first)

    def test_crash_after_chart_before_cloud_resumes_without_search_repeat(self):
        real=engine.atomic
        stopped=[False]
        def fail(path, value, **kw):
            real(path,value,**kw)
            if path.name=='chart-000.json' and not stopped[0]:
                stopped[0]=True
                raise RuntimeError('injected crash after chart')
        with patch.object(engine,'atomic',side_effect=fail):
            with self.assertRaisesRegex(RuntimeError,'injected'):
                engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),1)
        first=support.sha(self.folder/'replay-M17/epoch-00/chart-000.json')
        engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),3)
        self.assertEqual(support.sha(self.folder/'replay-M17/epoch-00/chart-000.json'),first)

    def test_crash_after_stage_before_terminal_resumes_without_search_repeat(self):
        real=engine.atomic
        def fail(path,value,**kw):
            real(path,value,**kw)
            if path.name=='stage.json':raise RuntimeError('injected stage crash')
        with patch.object(engine,'atomic',side_effect=fail):
            with self.assertRaises(RuntimeError):engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),3)
        engine.run_search(self.ctx)
        self.assertEqual(len(self.executed),3)
        self.assertIsNotNone(support.terminal_structure(self.folder))

    def test_rank_gain_rebuilds_and_preserves_prefix(self):
        self.gain_at=2
        engine.run_search(self.ctx)
        terminal=support.terminal_structure(self.folder)
        self.assertEqual([(s['before'],s['after']) for s in terminal['stages']],[(27,28),(28,28)])
        self.assertEqual(terminal['stages'][0]['stale_charts_cancelled'],1)
        self.assertEqual(len(self.executed),5)

    def test_chart_cap_is_cumulative_on_resume(self):
        self.policy['max_charts']=2
        put(self.folder/'protocol.json',self.policy)
        engine.run_search(self.ctx)
        t=support.terminal_structure(self.folder)
        self.assertEqual(t['charts'],2)
        self.assertEqual(t['stop_reason'],'CHART_BUDGET_EXHAUSTED')

    def test_odd_rank_disagreement_is_not_no_gain(self):
        with patch.object(engine,'_odd',return_value={'3':28,'5':27}):
            with self.assertRaisesRegex(ValueError,'odd-prime rank'):
                engine.run_search(self.ctx)
        self.assertFalse((self.folder/'replay-M17/terminal.json').exists())


def point_tuple_local(points):
    return tuple(tuple(F(x) for x in p) for p in points)


class ControllerTests(unittest.TestCase):
    def test_worker_commands_use_sage_for_every_action(self):
        for action in ('preflight','search','replay'):
            cmd=controller.worker_command('/chosen/sage',action,'warm-11952-41',Path('/session.json'))
            self.assertEqual(cmd[:3],['/chosen/sage','-python','-u'])
            self.assertIn('--session',cmd)
            self.assertNotEqual(cmd[0],sys.executable)

    def test_orphan_worker_blocks_concurrent_restart(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            put(root/'warm-11952-41/warm-run.supervisor.json',{'pid':12,'start_token':'token'})
            with patch.object(controller,'D',root),patch.object(controller,'AUTO',root/'auto'),patch.object(controller,'same_process',return_value=True):
                with self.assertRaisesRegex(RuntimeError,'still alive'):
                    controller.ensure_no_live_workers()

    def test_invalid_case_rejected(self):
        with self.assertRaises(ValueError):
            controller.worker_command('/sage','search','curve302',Path('/session'))

    def test_import_does_not_require_sage(self):
        code = '''import importlib.abc,sys,runpy
class NoSage(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname == "sage" or fullname.startswith("sage."):
   raise RuntimeError("unexpected Sage import")
sys.meta_path.insert(0,NoSage())
sys.path.insert(0,sys.argv[1])
runpy.run_path(sys.argv[1]+"/run_v3_warm_start_overnight.py",run_name="import_test")
'''
        p=subprocess.run([sys.executable,'-c',code,str(CAS)],capture_output=True,text=True,timeout=10)
        self.assertEqual(p.returncode,0,p.stderr)

    def test_lock_is_inherited_until_child_exits(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td)
            with patch.object(controller,'AUTO',folder),patch.object(controller,'LOCK',folder/'lock'):
                fd=controller.take_lock()
                child=subprocess.Popen([sys.executable,'-c','import sys; print("ready",flush=True); sys.stdin.readline()'],
                    pass_fds=(fd,),stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
                try:
                    self.assertEqual(child.stdout.readline().strip(),'ready')
                    os.close(fd)
                    with self.assertRaisesRegex(RuntimeError,'owns the lock'):controller.take_lock()
                finally:
                    child.communicate('\n',timeout=5)
                fd=controller.take_lock();os.close(fd)

    def test_process_token_and_zombie_parse(self):
        with tempfile.TemporaryDirectory() as td:
            folder=Path(td);(folder/'12').mkdir()
            fields=['Z','1']+['0']*17+['123456']
            (folder/'12/stat').write_text('12 (worker (name)) '+' '.join(fields))
            info=support.process_info(12,folder)
            self.assertEqual(info['state'],'Z')
            self.assertEqual(info['start_token'],'123456')

    def test_reused_pid_token_is_not_alive(self):
        row={'pid':12,'state':'S','start_token':'new','ppid':1}
        with patch.object(support,'process_info',return_value=row):
            self.assertFalse(support.same_process(12,'old'))
            self.assertTrue(support.same_process(12,'new'))

    def test_dead_status_does_not_claim_running(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);p=root/'state.json'
            put(p,{'status':'SEARCHING','pid':987654321,'start_token':'x'})
            with patch.object(controller,'STATE',p),patch.object(controller,'D',root),contextlib.redirect_stdout(io.StringIO()) as out:
                controller.status()
            self.assertIn('STOPPED_REVIEW_REQUIRED',out.getvalue())
            self.assertEqual(support.read(p)['status'],'SEARCHING')  # read-only inspection

    def test_worker_failure_records_terminal_status(self):
        self.run_fake_controller(fail_replay=True)

    def test_sealed_case_skips_search_then_advances(self):
        self.run_fake_controller(fail_replay=False)

    def run_fake_controller(self,fail_replay):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);auto=root/'auto';auto.mkdir();d=root/'d';session=auto/'session.json'
            dummy=root/'source.py';dummy.write_text('test')
            hashes={'source.py':support.sha(dummy)}
            put(session,{'sources':hashes,'inputs':hashes,'hours':1})
            terminal_fixture(d/'warm-11952-41',3)
            calls=[];finished={}
            def step(action,case,*args):
                calls.append((action,case))
                if action=='replay':
                    if fail_replay:raise RuntimeError('injected replay mismatch')
                    finished[case]={'rank_lower_bound':27,'gain':0}
            with contextlib.ExitStack() as stack:
                for name,value in [('ROOT',root),('AUTO',auto),('STATE',auto/'state.json'),('LOCK',auto/'lock'),('D',d)]:
                    stack.enter_context(patch.object(controller,name,value))
                stack.enter_context(patch.object(controller,'sources',return_value=hashes))
                stack.enter_context(patch.object(controller,'run_step',side_effect=step))
                stack.enter_context(patch.object(controller,'verify_result',side_effect=lambda case:finished.get(case)))
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                fd=controller.take_lock()
                before=support.sha(d/'warm-11952-41/replay-M17/epoch-00/chart-002.json')
                if fail_replay:
                    with self.assertRaises(RuntimeError):controller.worker(session,fd)
                    self.assertEqual(support.read(auto/'state.json')['status'],'STOPPED_REVIEW_REQUIRED')
                    self.assertNotIn(('preflight','warm-11952-72'),calls)
                else:
                    controller.worker(session,fd)
                    self.assertEqual(support.read(auto/'state.json')['status'],'COMPLETE_WARM_ROSTER')
                    self.assertIn(('search','warm-11952-72'),calls)
                    self.assertIn(('search','warm-11952-186'),calls)
                self.assertNotIn(('search','warm-11952-41'),calls)
                self.assertEqual(support.sha(d/'warm-11952-41/replay-M17/epoch-00/chart-002.json'),before)


if __name__ == '__main__':
    unittest.main()

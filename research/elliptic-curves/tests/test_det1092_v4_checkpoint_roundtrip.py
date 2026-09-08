"""Regression for V4 tuple/list checkpoints and zero-chart source recovery.

CVP tests use the REAL small-dimensional solver. Search/replay integration uses
mock arithmetic (not rank certificates) with the real worker/control flow,
JSON publication, checkpoint files and schedule comparisons. No Sage needed
except that the real CVP test requires the already-used NumPy dependency.
"""
from contextlib import ExitStack, redirect_stdout
import ast
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import det1092_v4_bootstrap_contract as c
import det1092_v4_bootstrap_worker as w
import det1092_v4_bootstrap_replay as replay
import run_det1092_v4_bootstrap as controller
from v3_warm_support import atomic, read, sha, point_tuple


def raw_selection():
    # This matches the actual solver's nested tuple representation.
    proof = {'norm': 4, 'minima': [(-1, 0), (1, 0)], 'nodes': 4,
             'certificate': 'complete exact rational LDL ellipsoid enumeration'}
    return {'rank': 17, 'exact_cvp_classes': 2, 'cvp_censored_masks': [], 'centres': [
        {'representative': [1, 0]+[0]*15, 'point': ['1', '2'], 'cvp': copy.deepcopy(proof)},
        {'representative': [0, 1]+[0]*15, 'point': ['3', '4'], 'cvp': copy.deepcopy(proof)}]}


class SelectionTests(unittest.TestCase):
    def test_actual_cvp_exhibits_original_failure(self):
        try:
            from visibility_lattice_v2 import ExactParity
        except ModuleNotFoundError as exc:
            if exc.name == 'numpy':
                self.skipTest('run under sage -python for the real CVP regression')
            raise
        proof = ExactParity([[4, 1], [1, 4]]).solve([1, 0], [1, 0])
        raw = {'centres': [{'cvp': proof}]}
        saved = json.loads(json.dumps(raw))
        self.assertIsInstance(proof['minima'][0], tuple)
        self.assertNotEqual(saved, raw)
        self.assertEqual(saved, c.selection_record(raw))

    def test_normalization_is_idempotent_and_does_not_mutate_solver_output(self):
        raw = raw_selection(); original = copy.deepcopy(raw)
        normalized = c.selection_record(raw)
        self.assertEqual(normalized, c.selection_record(normalized))
        self.assertEqual(raw, original)
        self.assertIsInstance(raw['centres'][0]['cvp']['minima'][0], tuple)

    def test_first_publication_readback_and_second_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'selection.json'
            generated = w.finish_selection(path, raw_selection(), True)
            self.assertEqual(generated, read(path))
            before = path.read_bytes()
            self.assertEqual(generated, w.finish_selection(path, raw_selection(), True))
            self.assertEqual(before, path.read_bytes())

    def test_pre_fix_saved_selection_bytes_are_compatible(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'selection.json'
            atomic(path, raw_selection(), immutable=True)
            before = path.read_bytes()
            w.finish_selection(path, raw_selection(), True)
            self.assertEqual(before, path.read_bytes())

    def test_large_integers_and_coordinate_strings_are_exact(self):
        original = {'n': 2**1000+1, 'point': ('-1/3', '9007199254740993/7')}
        result = c.selection_record(original)
        self.assertEqual(result, json.loads(json.dumps(original)))
        self.assertEqual(result['n'], 2**1000+1)

    def test_invalid_types_are_not_silently_coerced(self):
        for value in (1.0, float('nan'), float('inf'), F(1, 2), object(), {1}):
            with self.subTest(value=repr(value)), self.assertRaises(TypeError):
                c.selection_record({'value': value})
        with self.assertRaisesRegex(ValueError, 'non-string selection key'):
            c.selection_record({1: 'wrong'})

    def test_changed_cvp_coordinate_has_exact_diff_path(self):
        saved = c.selection_record(raw_selection())
        generated = copy.deepcopy(saved)
        generated['centres'][0]['cvp']['minima'][0][0] = -3
        with self.assertRaisesRegex(ValueError, r'\$\.centres\[0\]\.cvp\.minima\[0\]\[0\]'):
            c.require_same_selection(saved, generated)

    def test_changed_norm_is_not_a_serialization_fix(self):
        saved = c.selection_record(raw_selection()); generated = copy.deepcopy(saved)
        generated['centres'][0]['cvp']['norm'] += 1
        with self.assertRaisesRegex(ValueError, r'cvp\.norm'):
            c.require_same_selection(saved, generated)

    def test_reordered_centres_are_rejected(self):
        saved = c.selection_record(raw_selection()); other = copy.deepcopy(saved)
        other['centres'].reverse()
        with self.assertRaises(ValueError):
            c.require_same_selection(saved, other)

    def test_reordered_minima_are_rejected(self):
        saved = c.selection_record(raw_selection()); other = copy.deepcopy(saved)
        other['centres'][0]['cvp']['minima'].reverse()
        with self.assertRaises(ValueError):
            c.require_same_selection(saved, other)

    def test_bool_and_int_are_not_interchangeable(self):
        with self.assertRaisesRegex(ValueError, 'types bool/int'):
            c.require_same_selection({'n': True}, {'n': 1})

    def test_missing_and_extra_fields_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'keys differ'):
            c.require_same_selection({'a': 1}, {'a': 1, 'b': 0})

    def test_failed_immutable_write_keeps_original_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'selection.json'
            w.finish_selection(path, raw_selection(), True); before = path.read_bytes()
            changed = raw_selection(); changed['rank'] = 18
            with self.assertRaisesRegex(ValueError, r'\$\.rank'):
                w.finish_selection(path, changed, True)
            self.assertEqual(before, path.read_bytes())

    def test_actual_builder_uses_shared_roundtrip_boundary(self):
        tree = ast.parse((CAS/'det1092_v4_bootstrap_worker.py').read_text())
        func = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'bootstrap_selection')
        self.assertIsInstance(func.body[-1], ast.Return)
        self.assertEqual(func.body[-1].value.func.id, 'finish_selection')


class FakeState:
    def __init__(self, points):
        self.basis = point_tuple(points); self.rank = len(self.basis)
    def record(self):
        return {'state': {'reductions': {'points': [list(map(str, p)) for p in self.basis]}}}


class WorkerCheckpointTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack(); self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        self.folder = self.root/'case'; self.folder.mkdir()
        state = FakeState([(i+1, i+101) for i in range(17)])
        self.mapper = types.SimpleNamespace(pari=types.SimpleNamespace(allocatemem=lambda *a, **k: None))
        self.mapper.mapping = lambda model, basis, centre: {'centre': centre, 'coordinate_policy': {'kind': 'raw'}}
        self.ctx = types.SimpleNamespace(root=self.root, folder=self.folder, case='scale-test', state=state,
            model=(F(0), F(0), F(0), F(-1), F(1)), policy={'height':125000, 'seconds_per_chart':10, 'gp_sha256':'test'},
            engine=types.SimpleNamespace(load=lambda name: self.mapper), sources={'mock-arithmetic': 'not-a-proof'})
        self.execute_count = 0; self.interrupt_at = None
        backend = types.ModuleType('pari_pointed_backend')
        def execute(search, mapping, height, seconds, digest):
            self.execute_count += 1
            if self.execute_count == self.interrupt_at:
                raise RuntimeError('injected worker interruption')
            return {'status':'bounded_search_complete', 'finite_curve_points':[],
                    'height_bound':height, 'timeout_seconds':seconds, 'gp_binary_sha256':digest}, ()
        backend.execute = execute; backend.replay = lambda *a: (); backend.validate_map = lambda *a: None
        pointed = types.ModuleType('pointed_quartic_search')
        pointed.PointedQuarticSearch = lambda **kw: types.SimpleNamespace(**kw)
        eng = types.ModuleType('v3_warm_engine')
        eng.certified_state = lambda model, points, proof: FakeState(points)
        def audit(snapshot, output):
            record = read(snapshot)
            points = record['final_state']['state']['reductions']['points']
            result = {'status':'COMPLETE_DECLARED_FINITE_AUDIT','rank_lower_bound':17,
                'input_sha256':sha(snapshot),'input_path':str(snapshot.relative_to(self.root)),
                'curve':record['curve'],'points':points,'independent_points':points,
                'rank_certificate':{'MOCK':True}}
            atomic(output, result, immutable=True); return result
        eng._audit = audit
        def odd(audit_path, output):
            points = read(audit_path)['points']
            atomic(output, {'input_sha256':sha(audit_path), 'points':points,
                'audits':[{'modulus':3,'finite_column_rank':17},{'modulus':5,'finite_column_rank':17}]}, immutable=True)
            return {'3':17, '5':17}
        eng._odd = odd
        mod2 = types.ModuleType('audit_recorded_point_mod2_rank_v3'); mod2.check = lambda p: None
        modl = types.ModuleType('audit_retained_cloud_modl'); modl.check = lambda p: None
        self.stack.enter_context(patch.dict(sys.modules, {m.__name__:m for m in (backend,pointed,eng,mod2,modl)}))
        self.raw = raw_selection()
        def selection(ctx, publish=False):
            return w.finish_selection(ctx.folder/'bootstrap/selection.json', self.raw, publish)
        self.stack.enter_context(patch.object(w, 'bootstrap_selection', side_effect=selection))
        self.stack.enter_context(patch.object(c, 'FRESH_CHARTS', 2))
        self.stack.enter_context(redirect_stdout(io.StringIO()))
        for name, content in (('protocol.json',self.ctx.policy),('seed-input.json',{'parameter':'1'}),('seed-proof.json',{'MOCK':True})):
            atomic(self.folder/name, content)

    def test_preflight_exercises_the_save_read_compare_path(self):
        w.preflight(self.ctx)
        self.assertEqual(read(self.folder/'bootstrap/selection.json'), c.selection_record(self.raw))
        self.assertEqual(self.execute_count, 0)

    def test_preflight_search_and_independent_replay_roundtrip(self):
        w.preflight(self.ctx); w.run_search(self.ctx)
        result = replay.replay_case(self.ctx)
        self.assertEqual(self.execute_count, 2)
        self.assertEqual(result['bootstrap_charts'], 2)
        self.assertEqual(result['stop_reason'], 'COMPLETE_FRESH_BOOTSTRAP_NO_GAIN')

    def test_interruption_after_one_chart_reuses_it(self):
        self.interrupt_at = 2
        with self.assertRaisesRegex(RuntimeError, 'injected'):
            w.run_search(self.ctx)
        self.assertTrue((self.folder/'bootstrap/chart-000.json').exists())
        self.execute_count = 0; self.interrupt_at = None
        w.run_search(self.ctx)
        self.assertEqual(self.execute_count, 1)
        replay.replay_case(self.ctx)

    def test_sealed_search_is_not_run_again(self):
        w.run_search(self.ctx); self.execute_count = 0
        w.run_search(self.ctx)
        self.assertEqual(self.execute_count, 0)

    def test_real_selection_change_blocks_search(self):
        w.preflight(self.ctx)
        self.raw['centres'][0]['cvp']['norm'] = 6
        with self.assertRaisesRegex(ValueError, 'selection differs'):
            w.run_search(self.ctx)
        self.assertEqual(self.execute_count, 0)

    def test_real_selection_change_blocks_replay(self):
        w.run_search(self.ctx)
        self.raw['centres'].reverse()
        with self.assertRaisesRegex(ValueError, 'independent V4 atlas selection differs'):
            replay.replay_case(self.ctx)


class RepairTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack(); self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        self.directory = self.root/'campaign'; self.directory.mkdir()
        self.auto = self.root/'controller'; self.auto.mkdir()
        self.stack.enter_context(patch.object(c, 'D', self.directory))
        self.stack.enter_context(patch.object(controller, 'ROOT', self.root))
        self.stack.enter_context(patch.object(controller, 'AUTO', self.auto))
        self.stack.enter_context(patch.object(controller, 'controller_alive', return_value=False))
        self.stack.enter_context(redirect_stdout(io.StringIO()))
        self.cases = [f'scale-{i:07d}' for i in range(8)]
        rows = [{'id':case} for case in self.cases]
        self.stack.enter_context(patch.object(c, 'validate_v3_panel', return_value=({},rows)))
        self.old = {'status':'STOPPED_REVIEW_REQUIRED','error':'saved V4 selection differs'}
        numeric = self.root/'numeric.py'; numeric.write_text('# frozen\n')
        jobs = {}
        for case in self.cases:
            folder = self.directory/case; folder.mkdir()
            seed = folder/'seed-input.json'; atomic(seed, {'points':[]})
            p = folder/'protocol.json'
            atomic(p, {'inputs':{str(seed.relative_to(self.root)):sha(seed)}, 'sources':{'numeric.py':sha(numeric)}})
            jobs[case] = sha(p)
        atomic(self.directory/'roster.json', {'status':'READY_V4_EIGHT','cases':rows,'jobs':jobs})
        self.selection = self.directory/self.cases[0]/'bootstrap/selection.json'
        atomic(self.selection, raw_selection())

    def test_repair_preserves_all_bytes_and_carries_selection_without_rehashing(self):
        before = {str(p.relative_to(self.directory)):p.read_bytes() for p in self.directory.rglob('*') if p.is_file()}
        archive = controller.archive_unsearched_preparation(self.old)
        for name, data in before.items():
            self.assertEqual((archive/'campaign'/name).read_bytes(), data)
        self.assertEqual(self.selection.read_bytes(), before[str(self.selection.relative_to(self.directory))])
        self.assertFalse((self.directory/'roster.json').exists())
        self.assertEqual(read(archive/'manifest.json')['previous_controller_state'], self.old)

    def test_repair_refuses_even_one_chart(self):
        chart = self.selection.parent/'chart-000.json'; atomic(chart, {})
        with self.assertRaisesRegex(ValueError, 'repair refuses search evidence'):
            controller.archive_unsearched_preparation(self.old)
        self.assertTrue(chart.exists()); self.assertTrue((self.directory/'roster.json').exists())

    def test_repair_refuses_completed_result(self):
        atomic(self.directory/self.cases[0]/'v4-verified.json', {})
        with self.assertRaisesRegex(ValueError, 'repair refuses search evidence'):
            controller.archive_unsearched_preparation(self.old)

    def test_repair_refuses_changed_bound_input(self):
        atomic(self.directory/self.cases[0]/'seed-input.json', {'changed':True})
        with self.assertRaisesRegex(ValueError, 'changed bound file'):
            controller.archive_unsearched_preparation(self.old)
        self.assertTrue((self.directory/'roster.json').exists())

    def test_repair_refuses_live_controller(self):
        with patch.object(controller,'controller_alive',return_value=True), self.assertRaisesRegex(ValueError,'still alive'):
            controller.archive_unsearched_preparation(self.old)

    def test_repair_refuses_live_supervised_worker(self):
        atomic(self.auto/'sessions/s/attempts/c/search-1/supervisor.json', {'pid':123,'start_token':'5'})
        with patch.object(controller,'process_info',return_value={'state':'S'}), patch.object(controller,'same_process',return_value=True):
            with self.assertRaisesRegex(ValueError,'worker may still be alive'):
                controller.archive_unsearched_preparation(self.old)

    def test_repair_ignores_dead_supervised_worker(self):
        atomic(self.auto/'sessions/s/attempts/c/search-1/supervisor.json', {'pid':123,'start_token':'5'})
        with patch.object(controller,'process_info',return_value=None):
            controller.archive_unsearched_preparation(self.old)

    def test_repair_refuses_symlinks(self):
        (self.directory/'linked-input').symlink_to(self.root/'numeric.py')
        with self.assertRaisesRegex(ValueError,'symlink'):
            controller.archive_unsearched_preparation(self.old)

    def test_legacy_directory_names_and_roster_policy_stay_unchanged(self):
        text = (CAS/'det1092_v4_bootstrap_contract.py').read_text()
        self.assertIn("D = LOCAL/'det1092-v4-wide-bootstrap-v2'", text)
        self.assertIn("DOMAIN = 'det1092-v4-wide-bootstrap-v2'", text)
        self.assertEqual(c.FRESH_CHARTS, 512)


if __name__ == '__main__':
    unittest.main()

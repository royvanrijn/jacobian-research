"""Retained curve40 regression: timeout, exact fallback, and forged-map rejection."""
import json
from fractions import Fraction as F
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import bounded_map_receipts as maps
from v3_warm_engine import certified_state

LOCAL = CAS.parents[1]/'artifacts/local/elliptic-curves'


class BoundedMapReceiptsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        seed_path = LOCAL/'curve40-short-v3-discovery-v1/seed.json'
        input_path = LOCAL/'curve40-bounded-map-diagnostic-v1/factor_free-input.json'
        if not seed_path.exists() or not input_path.exists() or not shutil.which('sage'):
            raise unittest.SkipTest('requires retained curve40 inputs and Sage')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.epoch = Path(cls.tmp.name)
        seed = json.loads(seed_path.read_text())
        cls.model = tuple(map(F, seed['curve']))
        cls.points = tuple(tuple(map(F, p)) for p in seed['points'])
        cls.centre = json.loads(input_path.read_text())['centre']
        cls.state = certified_state(cls.model, cls.points, seed['proof'])
        cls.protocol = dict(map_python=shutil.which('sage'), seconds_per_map=5,
                            map_rss_bytes=1024**3)
        cls.minimized = cls.obtain('quartic_minimized')
        cls.factor_free = cls.obtain('factor_free')

    @classmethod
    def obtain(cls, policy, replay=False):
        return maps.obtain(cls.epoch, 0, policy, cls.model, cls.points, cls.centre,
                           cls.state, cls.protocol, replay=replay)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_timeout_replay_never_restarts_worker(self):
        self.assertIsNone(self.minimized[0])
        self.assertEqual(self.minimized[2], 'strict_wall_timeout')
        with patch.object(maps, 'run', side_effect=AssertionError('worker started during replay')):
            self.assertEqual(self.obtain('quartic_minimized', replay=True), self.minimized)

    def test_factor_free_replay_is_exact_and_worker_free(self):
        self.assertIsNone(self.factor_free[2])
        golden = json.loads((LOCAL/'curve40-bounded-map-diagnostic-v1/factor_free-result.json').read_text())
        self.assertEqual(self.factor_free[0], golden['mapping'])
        with patch.object(maps, 'run', side_effect=AssertionError('worker started during replay')):
            self.assertEqual(self.obtain('factor_free', replay=True), self.factor_free)

    def test_corrupt_map_fails_even_with_rebound_container_hashes(self):
        base = self.epoch/'map-0000-factor_free'
        result_path, sup_path = base/'result.json', base/'supervisor.json'
        receipt_path = base.with_suffix('.json')
        original = {p:p.read_bytes() for p in (result_path, sup_path, receipt_path)}
        try:
            result = maps.read(result_path)
            result['mapping']['reduced_P'][0] = str(F(result['mapping']['reduced_P'][0])+1)
            result_path.write_text(json.dumps(result))
            supervisor = maps.read(sup_path)
            supervisor['worker_result'] = result
            supervisor['worker_result_sha256'] = maps.sha(result_path)
            sup_path.write_text(json.dumps(supervisor))
            receipt = maps.read(receipt_path)
            receipt.update(result_sha256=maps.sha(result_path), supervisor_sha256=maps.sha(sup_path))
            receipt_path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ArithmeticError, 'reduced quartic identity failed'):
                self.obtain('factor_free', replay=True)
        finally:
            for path, raw in original.items():
                path.write_bytes(raw)


if __name__ == '__main__':
    unittest.main()

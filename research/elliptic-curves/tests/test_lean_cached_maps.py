"""Inherited-map routing preserves original receipts and never starts workers."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import run_lean_cached_seed_v3 as cached
import lean_preconditioned_map_receipts as maps


class CachedMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = CAS.parents[1]
        cls.source = cls.root/'artifacts/local/elliptic-curves/curve48-lean-preconditioned-v3-discovery-v1'
        if not (cls.source/'verified.json').exists():
            raise unittest.SkipTest('requires retained verified curve48 pass')
        read = lambda path: json.loads(path.read_text())
        seed = read(cls.source/'seed.json')
        cls.model = tuple(map(F, seed['curve']))
        cls.points = tuple(tuple(map(F, p)) for p in seed['points'])
        cls.state = cached.parent.certified_state(cls.model, cls.points, seed['proof'])
        cls.centre = read(cls.source/'epoch-00/landscape/selection.json')['centres'][0]
        cls.protocol = read(cls.source/'protocol.json')
        cls.key = 'epoch-00/map-0000-factor_free'
        cls.protocol['inherited_map_locations'] = {cls.key: str((cls.source/'epoch-00').relative_to(cls.root))}

    def obtain(self, protocol=None, centre=None):
        return cached.obtain_map(Path('/unused-continuation/epoch-00'), 0, 'factor_free',
            self.model, self.points, centre or self.centre, self.state, protocol or self.protocol)

    def test_original_map_reused_without_worker(self):
        with patch.object(maps, 'run', side_effect=AssertionError('worker launched')):
            mapping, seal, limited = self.obtain()
        original = json.loads((self.source/'epoch-00/map-0000-factor_free/result.json').read_text())
        self.assertEqual(mapping, original['mapping'])
        self.assertEqual(seal, maps.sha(self.source/'epoch-00/map-0000-factor_free.json'))
        self.assertIsNone(limited)

    def test_missing_inherited_receipt_does_not_start_worker(self):
        protocol = {**self.protocol, 'inherited_map_locations': {self.key: 'missing-map-regression/epoch-00'}}
        with patch.object(maps, 'run', side_effect=AssertionError('worker launched')):
            with self.assertRaisesRegex(ArithmeticError, 'missing map receipt during replay'):
                self.obtain(protocol=protocol)

    def test_changed_centre_rejected(self):
        centre = {**self.centre, 'unexpected_change': True}
        with patch.object(maps, 'run', side_effect=AssertionError('worker launched')):
            with self.assertRaisesRegex(ArithmeticError, 'map input differs'):
                self.obtain(centre=centre)


if __name__ == '__main__':
    unittest.main()

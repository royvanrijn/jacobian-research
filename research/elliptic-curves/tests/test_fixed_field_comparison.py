"""Comparison replay, transport and fail-closed gates; run with sage -python."""
import copy
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
from run_fixed_field_comparison import ROOT, read, setup, verify, search_queue
from research_runtime.store import FactStore
from research_runtime.sage_subspace import SageSubspaceBackend


class ComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = read(ROOT/'artifacts/generated-results/elliptic-curves/fixed_field_comparison_u1_v1.json.gz')

    def replay(self, evidence):
        with TemporaryDirectory() as directory:
            return verify(evidence, FactStore(directory))

    def test_empty_cache_replay_without_discovery(self):
        with patch.object(SageSubspaceBackend, 'cover', side_effect=AssertionError('conic search')), \
             patch('research_runtime.sage_subspace.quartic_local_witness', side_effect=AssertionError('local search')):
            result = self.replay(self.evidence)
        self.assertEqual(result['profile'], [13, 12, 1, 0])
        self.assertIsNone(result['full_curve_rank_upper'])
        self.assertEqual(result['radical_point_or_sha_status'], 'UNKNOWN')

    def test_missing_pair_cannot_be_filled_with_zero(self):
        bad = copy.deepcopy(self.evidence)
        bad['ct']['pairs'].pop()
        with self.assertRaisesRegex(ArithmeticError, 'incomplete CT'):
            self.replay(bad)

    def test_wrong_u_fails_binding(self):
        bad = copy.deepcopy(self.evidence)
        bad['parameter_u'] = '2'
        with self.assertRaises((ArithmeticError, FileNotFoundError)):
            self.replay(bad)

    def test_normalized_cover_map_is_checked(self):
        bad = copy.deepcopy(self.evidence)
        bad['covers'][0]['d_over_quartic_y'] = '0'
        with self.assertRaisesRegex(ArithmeticError, 'two-quadric'):
            self.replay(bad)

    def test_outside_radical_search_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        bad['searches'][0]['mask'] = bad['admissible_masks'][0]
        with self.assertRaisesRegex(ArithmeticError, 'schedule|radical'):
            self.replay(bad)

    def test_fabricated_realization_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        bad['searches'][0]['points'] = [{'quartic_point': ['0','1','1'], 'raw_point': ['0','0']}]
        # Update only the alleged count to exercise the actual point gate.
        from run_fixed_field_comparison import summarize
        bad['summary'] = summarize(bad['admissible_masks'],bad['ct'],bad['searches'],'1')
        with self.assertRaisesRegex(ArithmeticError, 'invalid quartic point'):
            self.replay(bad)

    def test_queue_is_bounded_without_exponential_pool(self):
        result = list(search_queue({'radical_global_masks': [1<<i for i in range(20)]}, 3))
        self.assertEqual(result, [1,2,4])

    def test_control_all_190_entries_vanish(self):
        evidence = read(ROOT/'artifacts/generated-results/elliptic-curves/fixed_field_comparison_u0_v1.json.gz')
        self.assertEqual(len(evidence['ct']['pairs']), 190)
        self.assertEqual(self.replay(evidence)['profile'], [20,0,20,20])


if __name__ == '__main__':
    unittest.main()

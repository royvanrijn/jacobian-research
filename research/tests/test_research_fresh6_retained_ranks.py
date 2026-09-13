"""Small witness controls; never restart a cohort search or its full replay."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import verify_fresh6_retained_ranks as replay


class PrimeSelectionTests(unittest.TestCase):
    def test_keep_sufficient_prime_blocks_and_drop_dependent_ones(self):
        blocks = [{'prime': p, 'rows': [row]} for p, row in
                  [(3, [1, 0]), (5, [1, 0]), (7, [1, 1]), (11, [0, 1])]]
        original = copy.deepcopy(blocks)
        self.assertEqual([s['prime'] for s in replay.select_signatures(blocks, 2)], [3, 7])
        self.assertEqual(blocks, original)

    def test_deficient_saved_rows_do_not_prove_independence(self):
        with self.assertRaisesRegex(ValueError, 'do not span'):
            replay.select_signatures([{'prime': 3, 'rows': [[1, 1]]}], 2)

    def test_duplicate_prime_is_rejected_even_after_full_span(self):
        with self.assertRaisesRegex(ValueError, 'duplicate prime'):
            replay.select_signatures([{'prime': 3, 'rows': [[1]]}]*2, 1)

    def test_binary_rows_require_exact_width_and_integer_bits(self):
        for row in ([1], [1, 0, 0], [1, 2], [1, True]):
            with self.subTest(row=row), self.assertRaisesRegex(ValueError, 'malformed binary'):
                replay.select_signatures([{'prime': 3, 'rows': [row]}], 2)

    def test_expanded_modulus_budget_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'small-prime protocol'):
            replay.select_signatures([{'prime': 1009, 'rows': [[1]]}], 1)


class RetainedEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(replay.WITNESS.read_text())
        cls.atlas = json.loads(replay.ATLAS.read_text())

    def setUp(self):
        # One endpoint suffices for these boundary controls; not all 37 replays.
        self.row = copy.deepcopy(self.data['rows'][0])
        self.family = next(f for f in self.atlas['families'] if f['family'] == self.row['family'])

    def test_pinned_input_bytes(self):
        self.assertEqual(hashlib.sha256(replay.WITNESS.read_bytes()).hexdigest(), replay.WITNESS_SHA256)
        self.assertEqual(hashlib.sha256(replay.ATLAS.read_bytes()).hexdigest(), self.data['atlas_sha256'])

    def test_one_literal_endpoint(self):
        self.assertEqual(replay.replay_row(self.row, self.family)[1], 18)

    def test_changed_extra_point_fails_membership(self):
        self.row['packet']['points'][-1][1] = str(replay.Fraction(self.row['packet']['points'][-1][1])+1)
        with self.assertRaisesRegex(ArithmeticError, 'point membership'):
            replay.replay_row(self.row, self.family)

    def test_changed_native_parameter_fails(self):
        self.row['parameter'] = '1/2'
        with self.assertRaisesRegex(ValueError, 'native family or section prefix'):
            replay.replay_row(self.row, self.family)

    def test_rank_header_is_not_trusted(self):
        self.row['packet']['rank_lower_bound'] += 1
        with self.assertRaisesRegex(ValueError, 'rank header'):
            replay.replay_row(self.row, self.family)

    def test_retained_finite_group_order_is_recomputed(self):
        self.row['packet']['proof']['signatures'][0]['group_order'] += 1
        with self.assertRaisesRegex(ValueError, 'recomputed finite witness'):
            replay.replay_row(self.row, self.family)

    def test_selected_composite_modulus_fails(self):
        self.row['packet']['proof']['signatures'][0]['prime'] = 9
        with self.assertRaisesRegex(ValueError, 'certificate primes'):
            replay.replay_row(self.row, self.family)

    def test_missing_torsion_witness_fails(self):
        del self.row['packet']['proof']['no_rational_2_torsion_prime']
        with self.assertRaises(KeyError):
            replay.replay_row(self.row, self.family)

    def test_incomplete_endpoint_roster_fails_before_arithmetic(self):
        with self.assertRaisesRegex(ValueError, 'endpoint roster'):
            replay.verify({**self.data, 'rows': self.data['rows'][:-1]}, self.atlas)


class ProvenanceTests(unittest.TestCase):
    def test_missing_sources_are_not_reconstructed(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(replay, 'ROOT', Path(folder)):
            with self.assertRaises(FileNotFoundError):
                replay.check_sources({'source_inputs': {'missing.json': '0'*64}, 'rows': []})
            self.assertEqual(list(Path(folder).iterdir()), [])

    def test_changed_source_bytes_fail(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(replay, 'ROOT', Path(folder)):
            Path(folder, 'source.json').write_text('{}\n')
            with self.assertRaisesRegex(ValueError, 'original input changed'):
                replay.check_sources({'source_inputs': {'source.json': '0'*64}, 'rows': []})


if __name__ == '__main__':
    unittest.main()

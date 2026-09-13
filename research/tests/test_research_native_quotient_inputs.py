"""Validate replay proposals and immutable inputs without importing Sage."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'elkies-k3/scripts'))
import replay_r17_norm12_native_icarm_quotient_audit as replay


class RetainedQuotientInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(replay.CERTIFICATE.read_text())

    def test_generation_files_are_preserved(self):
        self.assertEqual(hashlib.sha256(replay.PRODUCER.read_bytes()).hexdigest(), replay.PRODUCER_SHA256)
        self.assertEqual(hashlib.sha256(replay.CERTIFICATE.read_bytes()).hexdigest(), replay.CERTIFICATE_SHA256)

    def test_every_retained_matrix_has_exact_integer_dimensions(self):
        for record in self.certificate['fibres']:
            with self.subTest(curve=record['curve_id']):
                before = deepcopy(record)
                transport = replay.coordinate_rows(record, with_covers=False)
                complete = replay.coordinate_rows(record, with_covers=True)
                self.assertEqual(len(transport), record['displayed_point_count'])
                self.assertEqual(len(transport[0]), 17)
                self.assertEqual(len(complete[0]), 17 + len(record['alternate_q80_cover_audit']['splits']))
                self.assertEqual([row[:17] for row in complete], transport)
                complete[0][0] += 1
                self.assertEqual(record, before)

    def test_missing_rows_or_columns_are_rejected(self):
        for mutate in (lambda rows: rows.pop(), lambda rows: rows[0].pop()):
            rows = [[1, 0], [0, 1]]
            mutate(rows)
            with self.assertRaisesRegex(ArithmeticError, 'exact dimensions'):
                replay.integer_matrix(rows, 2, 2)

    def test_integer_coercion_cannot_hide_corruption(self):
        for value in (True, 1.0, 1.5, '1', None):
            with self.subTest(value=value), self.assertRaises(ArithmeticError):
                replay.integer_matrix([[value]], 1, 1)

    def test_duplicate_split_proposals_are_rejected(self):
        record = deepcopy(self.certificate['fibres'][1])
        record['alternate_q80_cover_audit']['splits'][1]['label'] = record['alternate_q80_cover_audit']['splits'][0]['label']
        with self.assertRaisesRegex(ArithmeticError, 'duplicate retained split'):
            replay.coordinate_rows(record, with_covers=True)

    def test_selected_subset_is_distinct_and_keeps_certificate_order(self):
        self.assertEqual(replay.selected_curves(None), replay.CURVE_IDS)
        self.assertEqual(replay.selected_curves([378, 12]), (12, 378))
        for selection in ([], [12, 12], [999]):
            with self.assertRaises(ValueError):
                replay.selected_curves(selection)

    def test_changed_input_bytes_do_not_enter_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.json'
            path.write_text('{}')
            with self.assertRaisesRegex(ArithmeticError, 'pinned input changed'):
                replay.load_pinned(path, replay.CERTIFICATE_SHA256)

    def test_interrupted_checkpoint_keeps_the_previous_verified_subset(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'receipt.json'
            first = {'status': 'INCOMPLETE_CHECKPOINT', 'fibres': [{'curve_id': 12}]}
            replay.write_checkpoint(path, first)
            with patch.object(replay.os, 'replace', side_effect=OSError('interrupted publication')):
                with self.assertRaises(OSError):
                    replay.write_checkpoint(path, {'status': 'PASS_SELECTED_NATIVE_COMPONENTS'})
            self.assertEqual(json.loads(path.read_text()), first)
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_replay_requires_an_explicit_component_choice(self):
        with patch.object(sys, 'argv', [str(replay.__file__)]), patch('sys.stderr'):
            with self.assertRaises(SystemExit) as failure:
                replay.main()
        self.assertEqual(failure.exception.code, 2)


if __name__ == '__main__':
    unittest.main()

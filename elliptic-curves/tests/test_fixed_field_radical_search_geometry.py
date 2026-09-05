"""Reject corrupted bounded-exclusion certificates (sage -python)."""
import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
import audit_fixed_field_radical_search_geometry as audit


class SearchGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = audit.inputs()
        cls.certificate = json.loads(audit.OUTPUT.read_text())

    def test_exact_replay(self):
        self.assertEqual(audit.verify(self.records, self.certificate), 37944)

    def test_missing_real_band_is_rejected(self):
        bad = copy.deepcopy(self.certificate['models'][0])
        bad['root_bands'].pop()
        with self.assertRaises(AssertionError):
            audit.verify_row(self.records[0], bad)

    def test_shifted_band_is_rejected(self):
        bad = copy.deepcopy(self.certificate['models'][0])
        bad['root_bands'][0] = ['-2', '-3/2']
        with self.assertRaises(AssertionError):
            audit.verify_row(self.records[0], bad)

    def test_false_uniform_lower_bound_is_rejected(self):
        bad = copy.deepcopy(self.certificate['models'][1])
        bad['uniform_f_lower_bound'] = str(10**80)
        with self.assertRaises(AssertionError):
            audit.verify_row(self.records[1], bad)

    def test_primitive_lattice_cannot_change_silently(self):
        bad = copy.deepcopy(self.records[0])
        bad['quadric_model']['variable_transform'][3] = '2'
        with self.assertRaises(AssertionError):
            audit.geometry(bad)

    def test_bounded_exclusion_cannot_become_sha(self):
        bad = copy.deepcopy(self.certificate)
        bad['point_or_sha']['596921'] = 'SHA'
        with self.assertRaises(AssertionError):
            audit.verify(self.records, bad)


if __name__ == '__main__':
    unittest.main()

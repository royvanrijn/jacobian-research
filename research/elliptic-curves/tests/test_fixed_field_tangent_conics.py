"""Fail-closed checks for the genuine-lift construction gate (sage -python)."""
import copy
import gzip
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
import run_fixed_field_tangent_conics as gate


class TangentConics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = json.loads(gzip.decompress(gate.EVIDENCE.read_bytes()))

    def test_exact_equivalences_and_claim_boundary(self):
        result = gate.audit(self.evidence)
        self.assertEqual(result['exact_reduced_conic_maps'], 8)
        self.assertEqual(result['genuine_higher_covers'], 0)
        self.assertEqual(set(result['point_or_sha'].values()), {'UNKNOWN'})

    def test_wrong_pencil_root_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        bad['conics'][0]['lambda'][0] = '1'
        with self.assertRaises(AssertionError):
            gate.audit(bad)

    def test_wrong_reduction_map_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        bad['local_runs'][0]['checkpoint']['map_to_initial_norm_conic'][0] = ['1', '0', '0']
        with self.assertRaises(AssertionError):
            gate.audit(bad)

    def test_fabricated_auxiliary_point_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        bad['local_runs'][0]['result']['point_on_auxiliary_conic'] = [['1', '0', '0'], ['0']*3, ['0']*3]
        with self.assertRaises(AssertionError):
            gate.audit(bad)

    def test_done_after_timeout_is_not_a_witness(self):
        raw = '<calculator><headers><warning>time limit</warning></headers><results><line>CONIC_OK true</line><line>DONE_CONIC</line></results></calculator>'
        self.assertEqual(gate.remote_status(raw), 'TIME_LIMIT_NO_WITNESS')

    def test_unexamined_cas_completion_is_rejected(self):
        bad = copy.deepcopy(self.evidence)
        raw = '<calculator><results><line>DONE_CONIC</line></results></calculator>'
        bad['remote_attempts'][0].update(xml=raw, status=gate.remote_status(raw))
        with self.assertRaisesRegex(AssertionError, 'completed CAS job'):
            gate.audit(bad)

    def test_new_indefinite_probe_witness_requires_replay(self):
        bad = copy.deepcopy(self.evidence)
        bad['indefinite_probe']['result']['point_on_auxiliary_conic'] = [['1'], ['0'], ['0']]
        with self.assertRaisesRegex(AssertionError, 'probe witness requires'):
            gate.audit(bad)


if __name__ == '__main__':
    unittest.main()

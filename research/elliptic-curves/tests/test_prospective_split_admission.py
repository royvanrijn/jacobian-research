"""Admission/replay safety regressions; synthetic curves, no research searches."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from sage.all import QQ, EllipticCurve
from split_seed_descent import build_frame, point_record
from prospective_split_admission import ProspectiveAdmission, frame_from_packet
from replay_split_admission import ReplayFrame


class ProspectiveAdmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.E = EllipticCurve(QQ, [-4, 1])
        cls.P, cls.Q = cls.E([0, 1]), cls.E([2, 1])
        cls.frame = build_frame(cls.E, [cls.P], 101)
        cls.engine = ProspectiveAdmission(cls.frame)
        cls.replay = ReplayFrame(cls.frame)

    def decision(self, P, steps=8):
        result = self.engine.consider(point_record(P), max_steps=steps)
        self.assertEqual(self.replay.replay(result)['status'], result['status'])
        return result

    def test_in_span_is_not_dependence(self):
        result = self.decision(self.Q)
        self.assertEqual(result['status'], 'NEW_INDEPENDENT_DIRECTION')
        self.assertEqual(result['reason'], 'GLOBAL_NONHALVING_ESCAPE')

    def test_independence_hidden_by_two_halves(self):
        result = self.decision(4*self.Q)
        self.assertEqual(result['status'], 'NEW_INDEPENDENT_DIRECTION')
        self.assertEqual(result['steps'], 2)

    def test_odd_index_cycle(self):
        frame = build_frame(self.E, [3*self.P], 101)
        result = ProspectiveAdmission(frame).consider(point_record(self.P))
        self.assertEqual(result['relation_multiplier'], '3')
        self.assertEqual(ReplayFrame(frame).replay(result)['status'], 'INHERITED_RATIONAL_SPAN')

    def test_integral_relation_and_zero(self):
        for P in [-3*self.P, self.E(0)]:
            self.assertEqual(self.decision(P)['status'], 'INHERITED_RATIONAL_SPAN')

    def test_unknown_retains_pending_unique_word(self):
        result = self.decision(self.P, 0)
        self.assertEqual(result['status'], 'UNKNOWN')
        self.assertEqual(result['pending']['parity_word'], [1])
        self.assertIsNone(result['pending']['target'])

    def test_finite_escape(self):
        full = build_frame(self.E, [self.P, self.Q], 101)
        frame = dict(full, basis=full['basis'][:1], inherited_rank=1,
                     finite_rows=[row[:1] for row in full['finite_rows']],
                     records=[dict(row, codes=row['codes'][:1]) for row in full['records']])
        result = ProspectiveAdmission(frame).consider(point_record(self.Q))
        self.assertEqual(result['reason'], 'FINITE_FOOTPRINT_ESCAPE')
        self.assertEqual(ReplayFrame(frame).replay(result)['status'], 'NEW_INDEPENDENT_DIRECTION')

    def test_deficient_reference_rejected(self):
        frame = copy.deepcopy(self.frame)
        frame['basis'] *= 2
        frame['finite_rows'] = [row*2 for row in frame['finite_rows']]
        for row in frame['records']:
            row['codes'] *= 2
        with self.assertRaisesRegex(ArithmeticError, 'inherited mod2 independence'):
            ProspectiveAdmission(frame)

    def test_packet_adapter_preserves_sealed_places(self):
        signatures = []
        for record in self.frame['records']:
            dim = record['dimension']
            signatures.append(dict(prime=record['prime'], group_order=record['order'],
                doubled_subgroup_order=record['order']//2**dim, quotient_dimension=dim,
                rows=[[(code >> j) & 1 for code in record['codes']] for j in range(dim)]))
        packet = dict(curve=self.frame['curve'], points=self.frame['basis'], proof=dict(
            signatures=signatures, no_rational_2_torsion_prime=self.frame['no_two_torsion_prime']))
        converted = frame_from_packet(packet)
        self.assertEqual(converted['finite_rows'], self.frame['finite_rows'])
        self.assertEqual(converted['records'], self.frame['records'])
        ProspectiveAdmission(converted)
        packet['proof']['signatures'].append(signatures[0])
        with self.assertRaisesRegex(ArithmeticError, 'duplicate footprint place'):
            frame_from_packet(packet)

    def test_relation_and_halving_tampering_rejected(self):
        result = self.decision(-3*self.P)
        bad = copy.deepcopy(result)
        bad['relation_word'][0] = '0'
        with self.assertRaises(ArithmeticError):
            self.replay.replay(bad)
        bad = copy.deepcopy(result)
        bad['history'][0]['next'] = point_record(self.P)
        with self.assertRaises(ArithmeticError):
            self.replay.replay(bad)

    def test_no_artifact_reads_during_admission_or_replay(self):
        with patch.object(Path, 'open', side_effect=AssertionError('unexpected artifact read')):
            result = self.engine.consider(point_record(self.Q))
            self.replay.replay(result)


if __name__ == '__main__':
    unittest.main()

"""Finite controls for the period-group proof correction; no CAS or search."""
import importlib.util
import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('ns0031_group', ROOT/'elkies-k3/scripts/verify_ns0031_period_group.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class PeriodGroupTests(unittest.TestCase):
    def setUp(self):
        self.witness = json.loads(replay.WITNESS.read_text())
        self.original = json.loads(replay.ORIGINAL.read_text())

    def test_literal_counterwitness(self):
        self.assertIn('existence remains UNKNOWN', replay.verify(self.witness, self.original))

    def test_changed_reflection_fails(self):
        self.witness['reflection'][2][0] = 3
        with self.assertRaisesRegex(ValueError, 'reflection formula'):
            replay.verify(self.witness, self.original)

    def test_nontrivial_discriminant_action_fails(self):
        self.witness['discriminant_kernel_integral_matrix'][2][2] = 0
        with self.assertRaisesRegex(ValueError, 'discriminant-kernel'):
            replay.verify(self.witness, self.original)

    def test_orientation_witness_must_be_fixed(self):
        self.witness['fixed_positive_plane_basis'][0][2] = 0
        with self.assertRaisesRegex(ValueError, 'plane is not fixed'):
            replay.verify(self.witness, self.original)

    def test_wrong_projective_representative_fails(self):
        self.witness['projective_spin_matrix'] = [[1,0],[0,37]]
        with self.assertRaisesRegex(ValueError, 'adjoint action'):
            replay.verify(self.witness, self.original)

    def test_wrong_order_transport_fails(self):
        self.witness['order_conjugation_coordinates'][0][0] = 2
        with self.assertRaisesRegex(ValueError, 'normalizer identity'):
            replay.verify(self.witness, self.original)

    def test_wrong_clifford_identification_fails(self):
        self.witness['trace_zero_clifford_images'][1][0][1] = 3
        with self.assertRaisesRegex(ValueError, 'identification mismatch'):
            replay.verify(self.witness, self.original)


class CorrectionPropagationTests(unittest.TestCase):
    surface = 'K3-d1b1381f87d69f1c'

    def load(self, filename):
        return json.loads((ROOT/'artifacts/generated-results'/filename).read_text())

    def test_old_certificate_label_cannot_restore_exclusion(self):
        data = self.load('elkies-k3-rank19-arithmetic-marking-classifier-v1.json')
        row = next(r for r in data['candidates'] if r['surface_id'] == self.surface)
        self.assertEqual(row['classification'], 'UNKNOWN')
        self.assertFalse(row['equation_agent_eligible'])
        self.assertFalse(row['different_ns_foundry_equation_eligible'])
        self.assertEqual(row['easy_quotient_maps'], [])
        curve = row['full_discriminant_marking_curve']
        self.assertIsNone(curve.get('genus'))
        self.assertEqual(curve['coarse_norm_one_curve']['genus'], 23)
        self.assertTrue(any('counterwitness' in c['path'] for c in row['certificate_replay']))

    def test_unknown_returns_to_research_but_not_equation_queue(self):
        data = self.load('elkies-k3-rank7-determinant-aware-ranking-v1.json')
        self.assertTrue(any(r['surface_id'] == self.surface for r in data['candidates']))
        self.assertFalse(any(r['surface_id'] == self.surface for r in data['arithmetic_marking_rejections']))
        self.assertNotIn(self.surface, json.dumps(data['expensive_equation_scoring_queue']))

    def test_research_queue_preserves_known_arithmetic_and_exact_remaining_gate(self):
        data = self.load('elkies-k3-arithmetic-first-marked-t-foundry-v1.json')
        row = next(r for r in data['curve_identification_queue'] if r['surface_id'] == self.surface)
        self.assertEqual(row['arithmetic_classification'], 'UNKNOWN')
        self.assertIsNone(row['full_marking_curve']['rational_non_CM_point'])
        self.assertIsNone(row['full_marking_curve']['genus'])
        self.assertEqual(row['coarse_curve_diagnostic']['genus'], 23)
        self.assertEqual(row['phase_one_priority_key'][2], 23)
        self.assertIn('determinant-minus-one', row['next_gate'])
        self.assertEqual(data['new_positive_NS_rootless_handoff'], [])

    def test_new_witness_does_not_claim_rational_k3_nonexistence(self):
        entries = {e['id']: e for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
        self.assertEqual(entries['EC-K3-NS0031-QQ-MARKING-OBSTRUCTION']['state'], 'partial')
        counter = entries['EC-K3-NS0031-PERIOD-GROUP-COUNTERWITNESS']
        self.assertEqual(counter['state'], 'proved')
        self.assertIn('existence, which remains UNKNOWN', counter['scope'])
        self.assertNotIn('EC-K3-NS0031-QQ-MARKING-OBSTRUCTION', counter['dependencies'])

    def test_pre_correction_evidence_is_preserved_byte_for_byte(self):
        receipt = json.loads((ROOT/'archive/repository-cleanup-2026-09-12/ns0031-period-group-review/REVIEW.json').read_text())
        self.assertIn('MATH_STATUS.json', receipt['before_files'])
        self.assertIn('elkies-k3/NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md', receipt['before_files'])
        for source, record in receipt['before_files'].items():
            with self.subTest(source=source):
                self.assertEqual(hashlib.sha256((ROOT/record['preserved_path']).read_bytes()).hexdigest(), record['sha256'])


if __name__ == '__main__':
    unittest.main()

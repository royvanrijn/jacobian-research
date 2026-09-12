"""Rebinding container hashes must not authorize changed research protocols."""

import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'rescue_commitment', ROOT / 'elliptic-curves/cas/audit_mw17_rescue_commitment.py',
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class RescueProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads((ROOT / 'artifacts/generated-results/elkies-k3-mw17-jump-v2-zero-gain-rescue-arm-v1.json').read_text())

    def rebound(self, change):
        changed = copy.deepcopy(self.original)
        change(changed)
        changed['protocol_definition_sha256'] = audit.canonical_hash(
            {k: v for k, v in changed.items() if k != 'protocol_definition_sha256'},
        )
        return changed

    def test_rehashed_budget_change_is_rejected(self):
        changed = self.rebound(lambda d: d['rescue_detector'].update(additional_budget_chart_count=302))
        with self.assertRaisesRegex(ValueError, 'budget differs'):
            audit.compare_definition(self.original, changed, (), lambda _: b'')

    def test_rehashed_assignment_change_is_rejected(self):
        changed = self.rebound(lambda d: d['assignments'][0].update(
            assigned_to_rescue_arm=not d['assignments'][0]['assigned_to_rescue_arm'],
        ))
        with self.assertRaisesRegex(ValueError, 'assignment'):
            audit.compare_definition(self.original, changed, (), lambda _: b'')

    def test_unapproved_source_change_is_rejected(self):
        source = next(iter(self.original['implementation_hashes']))
        changed = self.rebound(lambda d: d['implementation_hashes'].update({source: '0' * 64}))
        with self.assertRaisesRegex(ValueError, 'source change outside'):
            audit.compare_definition(self.original, changed, (), lambda _: b'')


if __name__ == '__main__':
    unittest.main()

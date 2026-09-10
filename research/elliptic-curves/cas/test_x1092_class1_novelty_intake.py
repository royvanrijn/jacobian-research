"""Closed-gate and frozen-input integrity regressions; no search exposure."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import freeze_x1092_class1_arithmetic as freeze
import x1092_class1_novelty_intake as gate

class FrozenIntakeTests(unittest.TestCase):
    def test_no_checker_cannot_release_any_stage(self):
        r=gate.decision()
        self.assertFalse(r['V3_release'])
        self.assertFalse(r['parameter_panel_enabled'])
        self.assertEqual(r['allowed_class_indices'],[1])
        self.assertEqual(r['gates']['cover_construction'],'BLOCKED')

    def test_forged_pass_or_independence_does_not_bypass_checker(self):
        for candidate in ({'status':'PASS','rank':32},
                          {'constructed_strict_class':'CERTIFIED_NONZERO',
                           'independent_quotient_direction':'CERTIFIED_INDEPENDENT'},
                          {'class_index':19,'V3_release':True}):
            r=gate.decision(candidate)
            self.assertFalse(r['V3_release'])
            self.assertEqual(r['gates']['prospective_class_novelty'],'UNAVAILABLE')

    def test_write_once_rejects_changed_existing_parent(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'parent.json';freeze.write_once(p,b'original')
            freeze.write_once(p,b'original')
            with self.assertRaises(ValueError):freeze.write_once(p,b'mutated')
            self.assertEqual(p.read_bytes(),b'original')

    def test_source_parent_hash_pin_is_enforced(self):
        with patch.object(freeze,'PINNED_PARENT','0'*64):
            with self.assertRaisesRegex(ValueError,'parent changed'):freeze.frozen_manifest()

    def test_frozen_copy_tamper_is_rejected(self):
        manifest=freeze.frozen_manifest()
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)
            for name in manifest['files']:(out/name).write_bytes((freeze.DEST/name).read_bytes())
            (out/'freeze.json').write_text(json.dumps(manifest))
            (out/'parent.json').write_text('{}')
            with patch.object(freeze,'DEST',out):
                with self.assertRaisesRegex(ValueError,'frozen input changed'):freeze.check()

if __name__=='__main__':unittest.main()

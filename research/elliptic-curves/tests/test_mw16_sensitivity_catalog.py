from hashlib import sha256
from importlib.machinery import SourceFileLoader
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/'scripts/audit_artifact_catalog.py'
audit=SourceFileLoader('sensitivity_catalog_audit',str(path)).load_module()


class EmbeddedEvidenceTests(unittest.TestCase):
    name='artifacts/local/mw16-catalogue-test-absent.json'

    def bundle(self):
        text='{"status":"COMPLETE"}\n'
        return {'schema':'elliptic-curves.mw16-sensitivity-evidence.v1',
            'files':{self.name:{'text':text,'sha256':sha256(text.encode()).hexdigest()}}}

    def test_embedded_checkpoint_does_not_require_ignored_file(self):
        self.assertFalse((audit.ROOT/self.name).exists())
        self.assertEqual(audit.audit_references({'bundle':self.bundle(),'summary':{'input':self.name}}),1)

    def test_missing_or_corrupted_embedded_input_fails_closed(self):
        with self.assertRaises(AssertionError): audit.audit_references({'summary':{'input':self.name}})
        bundle=self.bundle();bundle['files'][self.name]['text']='corrupted'
        with self.assertRaises(AssertionError): audit.audit_references({'bundle':bundle})


if __name__=='__main__': unittest.main()

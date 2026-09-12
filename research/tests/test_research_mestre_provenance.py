"""Historical source preservation must not hide a changed mathematical result."""

import copy
from hashlib import sha256
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'elliptic-curves/cas'))
from mestre_replay_provenance import check_retained_replay


class RetainedMestreTests(unittest.TestCase):
    def test_changed_payload_is_rejected(self):
        old = {'rank': 11, 'sources': {}}
        with self.assertRaisesRegex(ValueError, 'mathematical payload'):
            check_retained_replay({'rank': 13, 'sources': {}}, old, Path('.'))

    def test_unlisted_source_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'helper.py').write_bytes(b'changed')
            old = {'rank': 11, 'sources': {'helper.py': sha256(b'original').hexdigest()}}
            current = copy.deepcopy(old)
            current['sources']['helper.py'] = sha256(b'changed').hexdigest()
            with self.assertRaisesRegex(ValueError, 'undeclared'):
                check_retained_replay(current, old, root)

    def test_unmodified_replay_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'helper.py').write_bytes(b'original')
            value = {'rank': 11, 'sources': {'helper.py': sha256(b'original').hexdigest()}}
            check_retained_replay(value, value, root)

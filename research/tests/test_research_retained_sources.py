"""A historical source view must never mask a changed data input."""

import hashlib
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "elkies-k3/scripts"))
from audit_retained_r17_singleton_po0 import source_digest_view


class RetainedSourceTests(unittest.TestCase):
    def test_source_relocation_preserves_data_drift_and_restores_checker(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            active, archived, data = (root / name for name in ("producer.py", "archived.py", "input.json"))
            active.write_text("new producer")
            archived.write_text("generation-time producer")
            data.write_text("original witness")
            digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            module = SimpleNamespace(digest=digest)
            original_data_hash = digest(data)
            with self.assertRaisesRegex(RuntimeError, "simulated audit failure"):
                with source_digest_view(module, {active.resolve(): archived}):
                    self.assertEqual(module.digest(active), digest(archived))
                    self.assertEqual(module.digest(data), original_data_hash)
                    data.write_text("changed witness")
                    self.assertNotEqual(module.digest(data), original_data_hash)
                    raise RuntimeError("simulated audit failure")
            self.assertIs(module.digest, digest)
            self.assertEqual(module.digest(active), digest(active))

    def test_change_during_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            data = Path(directory) / "input.json"
            data.write_text("original")

            def unstable_digest(path):
                value = hashlib.sha256(path.read_bytes()).hexdigest()
                path.write_text("changed during read")
                return value

            module = SimpleNamespace(digest=unstable_digest)
            with source_digest_view(module, {}):
                with self.assertRaisesRegex(ArithmeticError, "changed while hashing"):
                    module.digest(data)


if __name__ == "__main__":
    unittest.main()

"""Current input declarations must not silently drift with edited proof files."""

from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research


class CurrentInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.root = self.repository / "research"
        self.directory = self.root / "result"
        self.directory.mkdir(parents=True)
        for name, value in (("ROOT", self.root), ("REPOSITORY", self.repository)):
            changed = patch.object(research, name, value)
            changed.start()
            self.addCleanup(changed.stop)
        self.proof = self.directory / "PROOF.md"
        self.proof.write_text("The retained proof.\n")
        self.digest = hashlib.sha256(self.proof.read_bytes()).hexdigest()

    def record(self, name, data):
        path = self.directory / name
        path.write_text(json.dumps(data))
        return path

    def test_edited_proof_invalidates_context_and_local_source_pins(self):
        records = [
            self.record("context.json", {"retained_sha256": {"result/PROOF.md": self.digest}}),
            self.record("sources.json", {"local_sources": [
                {"path": "result/PROOF.md", "sha256": self.digest}]}),
        ]
        self.assertEqual(research.check_current_inputs(records), 2)
        self.proof.write_text("An edited proof.\n")
        for record in records:
            with self.subTest(record=record.name):
                with self.assertRaisesRegex(AssertionError, "changed current input"):
                    research.check_current_inputs([record])

    def test_historical_receipt_is_separate_from_current_replay_inputs(self):
        record = self.record("proof_verification.json", {
            "record_scope": "historical",
            "files": {"unbundled-old-input.json": "0" * 64},
            "current_replay": {"files": {"PROOF.md": self.digest}},
        })
        self.assertEqual(research.check_current_inputs([record]), 1)
        self.proof.write_text("An edited proof.\n")
        with self.assertRaisesRegex(AssertionError, "changed current input"):
            research.check_current_inputs([record])

    def test_non_markdown_primary_mirror_is_still_checksum_checked(self):
        primary = self.directory / "primary.pdf"
        primary.write_bytes(b"%PDF-1.4\noriginal primary source\n")
        record = self.record("source_ledger.json", {"cached_primary_sources": [{
            "file": "primary.pdf", "url": "https://example.org/primary.pdf",
            "sha256": hashlib.sha256(primary.read_bytes()).hexdigest(),
        }]})
        self.assertEqual(research.source_markdown_mirrors([record]), set())
        primary.write_bytes(b"%PDF-1.4\naltered primary source\n")
        with self.assertRaisesRegex(AssertionError, "changed pinned upstream"):
            research.source_markdown_mirrors([record])


if __name__ == "__main__":
    unittest.main()

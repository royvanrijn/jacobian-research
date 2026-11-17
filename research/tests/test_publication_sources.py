"""Extracted published sources must match their already checksum-validated archive."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research


class ExtractedPublicationSourceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.directory = self.repository / "results" / "paper"
        self.directory.mkdir(parents=True)
        self.contents = {
            "main.tex": b"\\documentclass{article}\r\n\\begin{document}Exact source\\end{document}\r\n",
            "main.bbl": b"\\begin{thebibliography}{1}\r\n\\end{thebibliography}\r\n",
        }
        self.archive_name = "official-source.zip"
        self.members = {name: "arxiv-upload/" + name for name in self.contents}
        with ZipFile(self.directory / self.archive_name, "w") as archive:
            for name, contents in self.contents.items():
                archive.writestr(self.members[name], contents)
                (self.directory / name).write_bytes(contents)
        self.meta = {
            "source": "main.tex",
            "source_bundle": self.archive_name,
            "source_bundle_members": dict(self.members),
        }
        repository_patch = patch.object(research, "REPOSITORY", self.repository)
        repository_patch.start()
        self.addCleanup(repository_patch.stop)

    def test_exact_primary_source_and_bibliography_match_the_sealed_archive(self):
        self.assertIsNone(research.check_extracted_sources(
            self.meta, self.directory, {self.archive_name}))

    def test_altered_extracted_primary_or_bibliography_is_rejected(self):
        for name, contents in self.contents.items():
            with self.subTest(source=name):
                (self.directory / name).write_bytes(contents.replace(b"\r\n", b"\n"))
                with self.assertRaisesRegex(AssertionError, "extracted publication source changed"):
                    research.check_extracted_sources(
                        self.meta, self.directory, {self.archive_name})
                (self.directory / name).write_bytes(contents)

    def test_existing_but_unsealed_archive_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "published source lacks a sealed archive"):
            research.check_extracted_sources(self.meta, self.directory, {"paper.pdf"})

    def test_primary_source_must_be_declared_as_an_archive_member(self):
        self.meta["source_bundle_members"] = {"main.bbl": self.members["main.bbl"]}
        with self.assertRaisesRegex(AssertionError, "published source is not bound to an archive member"):
            research.check_extracted_sources(
                self.meta, self.directory, {self.archive_name})


if __name__ == "__main__":
    unittest.main()

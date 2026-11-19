"""Deposited editions stay sealed while current revisions can change independently."""

import hashlib
import json
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



class CurrentAndDepositedEditionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name)
        self.current = self.repository / "results" / "paper"
        self.deposited = self.current / "published"
        self.draft = self.repository / "unpublished" / "other-paper"
        for directory in (self.current, self.deposited, self.draft):
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "README.md").write_text("Paper edition\n")
        self.current_meta = self.make_revision(self.current)
        self.current_meta["published_edition"] = "published/"
        self.write_json(self.current / "metadata.json", self.current_meta)
        self.write_json(self.draft / "metadata.json", self.make_revision(self.draft))
        self.deposited_contents = {
            "main.tex": b"Exact deposited source\r\n",
            "paper.pdf": b"%PDF-deposited-edition\n",
        }
        for name, contents in self.deposited_contents.items():
            (self.deposited / name).write_bytes(contents)
        self.deposited_meta = {
            "title": "Deposited paper title",
            "summary": "The precise deposited result.",
            "author": {"name": "An Author", "orcid": "0009-0000-2808-0645"},
            "status": "published-preprint",
            "doi": "10.5281/zenodo.12345",
            "zenodo": "https://zenodo.org/records/12345",
            "publication_date": "2026-07-23",
            "source": "main.tex",
            "pdf": "paper.pdf",
            "current_edition": "../",
            "verification": {
                "scope": "Written proof.",
                "limits": "No independent review claimed.",
                "evidence": ["main.tex"],
            },
        }
        self.write_json(self.deposited / "metadata.json", self.deposited_meta)
        self.seal = {
            "record": self.deposited_meta["zenodo"],
            "files": [
                {
                    "file": name,
                    "sha256": hashlib.sha256(contents).hexdigest(),
                    "md5": hashlib.md5(contents).hexdigest(),
                }
                for name, contents in self.deposited_contents.items()
            ],
        }
        self.write_json(self.deposited / "deposited-files.json", self.seal)
        repository_patch = patch.object(research, "REPOSITORY", self.repository)
        repository_patch.start()
        self.addCleanup(repository_patch.stop)

    @staticmethod
    def write_json(path, data):
        path.write_text(json.dumps(data))

    def make_revision(self, directory):
        (directory / "main.tex").write_bytes(b"Current revised source\n")
        (directory / "paper.pdf").write_bytes(b"%PDF-current-revision\n")
        return {
            "title": "A paper",
            "summary": "A precise result.",
            "author": {"name": "An Author", "orcid": "0009-0000-2808-0645"},
            "status": "revision",
            "published": False,
            "source": "main.tex",
            "pdf": "paper.pdf",
            "verification": "Source reviewed; external review remains open.",
            "build": "latexmk -pdf main.tex",
            "next_steps": ["Independent mathematical review."],
        }

    def test_current_revision_can_differ_from_its_separately_sealed_deposit(self):
        self.assertNotEqual(
            (self.current / "main.tex").read_bytes(),
            (self.deposited / "main.tex").read_bytes(),
        )
        self.assertEqual(research.check_portfolio(), (1, 1))
        (self.current / "main.tex").write_bytes(b"A further local revision\n")
        (self.current / "paper.pdf").write_bytes(b"%PDF-further-revision\n")
        self.assertEqual(research.check_portfolio(), (1, 1))

    def test_tampered_deposited_source_or_pdf_is_rejected(self):
        for name, contents in self.deposited_contents.items():
            with self.subTest(file=name):
                (self.deposited / name).write_bytes(contents + b"changed")
                with self.assertRaisesRegex(AssertionError, "deposit SHA mismatch"):
                    research.check_portfolio()
                (self.deposited / name).write_bytes(contents)

    def test_pending_submission_can_live_in_results_without_a_public_archive(self):
        submitted = self.current.parent / "pending"
        submitted.mkdir()
        (submitted / "README.md").write_text("Pending submission\n")
        meta = self.make_revision(submitted)
        meta.update(status="submission-pending", submission={
            "venue": "Hexagon", "reserved_identifier": "hexagon:2611.00001v1",
        })
        self.write_json(submitted / "metadata.json", meta)
        self.assertEqual(research.check_portfolio(), (2, 1))

        meta["submission"] = {}
        self.write_json(submitted / "metadata.json", meta)
        with self.assertRaisesRegex(AssertionError, "missing submission venue"):
            research.check_portfolio()

    def test_pending_submission_does_not_hide_an_existing_archive(self):
        self.current_meta.pop("published_edition")
        self.current_meta.update(status="submission-pending", submission={"venue": "Hexagon"})
        self.write_json(self.current / "metadata.json", self.current_meta)
        with self.assertRaisesRegex(AssertionError, "missing published archive declaration"):
            research.check_portfolio()

    def test_missing_published_archive_declaration_is_rejected(self):
        self.current_meta.pop("published_edition")
        self.write_json(self.current / "metadata.json", self.current_meta)
        with self.assertRaisesRegex(AssertionError, "missing published archive"):
            research.check_portfolio()

    def test_missing_published_archive_directory_is_rejected(self):
        self.deposited.rename(self.current / "misplaced-archive")
        with self.assertRaisesRegex(AssertionError, "missing reference: published/"):
            research.check_portfolio()

    def test_current_revision_cannot_claim_to_be_a_public_deposit(self):
        self.current_meta["published"] = True
        self.write_json(self.current / "metadata.json", self.current_meta)
        with self.assertRaisesRegex(AssertionError, "current revision mislabeled as deposited"):
            research.check_portfolio()

    def test_current_revision_cannot_use_published_status_with_false_flag(self):
        self.current_meta["status"] = "published-preprint"
        self.write_json(self.current / "metadata.json", self.current_meta)
        with self.assertRaises(AssertionError):
            research.check_portfolio()

    def test_archived_edition_must_keep_its_published_status(self):
        self.deposited_meta["status"] = "revision"
        self.write_json(self.deposited / "metadata.json", self.deposited_meta)
        with self.assertRaisesRegex(AssertionError, "archive is not a deposited edition"):
            research.check_portfolio()

    def test_archived_title_must_be_nonempty(self):
        for title in ("", "   "):
            with self.subTest(title=title):
                self.deposited_meta["title"] = title
                self.write_json(self.deposited / "metadata.json", self.deposited_meta)
                with self.assertRaises(AssertionError):
                    research.check_portfolio()

    def test_archived_author_must_have_name_and_valid_orcid(self):
        for author in (
            {"name": "", "orcid": "0009-0000-2808-0645"},
            {"name": "An Author", "orcid": "invalid-orcid"},
        ):
            with self.subTest(author=author):
                self.deposited_meta["author"] = author
                self.write_json(self.deposited / "metadata.json", self.deposited_meta)
                with self.assertRaises(AssertionError):
                    research.check_portfolio()

    def test_seal_must_identify_the_same_public_record(self):
        self.seal["record"] = "https://zenodo.org/records/99999"
        self.write_json(self.deposited / "deposited-files.json", self.seal)
        with self.assertRaisesRegex(AssertionError, "deposit record mismatch"):
            research.check_portfolio()

    def make_hexagon_deposit(self):
        self.deposited_meta.pop("doi")
        self.deposited_meta.pop("zenodo")
        self.deposited_meta["hexagon"] = "https://hexagonmath.org/2611.00596v1"
        self.seal["record"] = self.deposited_meta["hexagon"]
        self.write_json(self.deposited / "metadata.json", self.deposited_meta)
        self.write_json(self.deposited / "deposited-files.json", self.seal)

    def test_hexagon_deposit_without_a_zenodo_doi_is_accepted(self):
        self.make_hexagon_deposit()
        self.assertEqual(research.check_portfolio(), (1, 1))

    def test_hexagon_deposit_requires_an_exact_versioned_record(self):
        self.make_hexagon_deposit()
        for record in (
            "https://hexagonmath.org/2611.00596",
            "https://hexagonmath.org/2611.00596v0",
            "https://example.com/2611.00596v1",
            "https://hexagonmath.org/2611.00596v1/other",
        ):
            with self.subTest(record=record):
                self.deposited_meta["hexagon"] = record
                self.write_json(self.deposited / "metadata.json", self.deposited_meta)
                with self.assertRaisesRegex(AssertionError, "invalid versioned Hexagon"):
                    research.check_portfolio()

    def test_hexagon_seal_must_identify_the_same_version(self):
        self.make_hexagon_deposit()
        self.seal["record"] = "https://hexagonmath.org/2611.00596v2"
        self.write_json(self.deposited / "deposited-files.json", self.seal)
        with self.assertRaisesRegex(AssertionError, "deposit record mismatch"):
            research.check_portfolio()

    def test_hexagon_deposited_files_remain_sealed(self):
        self.make_hexagon_deposit()
        (self.deposited / "paper.pdf").write_bytes(b"%PDF-changed")
        with self.assertRaisesRegex(AssertionError, "deposit SHA mismatch"):
            research.check_portfolio()

    def test_hexagon_link_does_not_bypass_zenodo_doi_validation(self):
        self.deposited_meta.update(
            doi="invalid-doi", hexagon="https://hexagonmath.org/2611.00596v1")
        self.write_json(self.deposited / "metadata.json", self.deposited_meta)
        with self.assertRaises(AssertionError):
            research.check_portfolio()

    def test_deposited_pdf_must_be_sealed(self):
        self.seal["files"] = [item for item in self.seal["files"] if item["file"] != "paper.pdf"]
        self.write_json(self.deposited / "deposited-files.json", self.seal)
        with self.assertRaisesRegex(AssertionError, "published PDF is not the sealed deposited file"):
            research.check_portfolio()


if __name__ == "__main__":
    unittest.main()

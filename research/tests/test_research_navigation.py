"""Retrieval must preserve authority, coverage and non-execution boundaries."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research
from research_documents import document_paths
from render_status import CORE_ORDER, render


def claim(identifier, state="proved", source="elliptic-curves/result.md", **updates):
    entry = dict(id=identifier, state=state, kind="example", title="Pointed quartic result",
                 scope="Bounded miss; the exact rank remains UNKNOWN.", canonical_source=source,
                 dependencies=[], replaced_by=[], checker="never-execute.py")
    entry.update(updates)
    return entry


class RetrievalTests(unittest.TestCase):
    def test_exact_id_retains_falsification_and_replacement(self):
        old = claim("OLD-Q8", state="falsified", replaced_by=["NEW-Q8"])
        new = claim("NEW-Q8")
        with patch.object(research, "source_documents", return_value=[]):
            hits = research.search("OLD-Q8", [new, old], [])
        self.assertEqual(hits[0]["id"], "OLD-Q8")
        self.assertEqual(hits[0]["state"], "falsified")
        self.assertEqual(hits[0]["replaced_by"], ["NEW-Q8"])
        self.assertIn("UNKNOWN", research.show("OLD-Q8", [old, new], [])["entry"]["scope"])

    def test_search_reads_edits_and_labels_unregistered_and_archived_notes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "elliptic-curves/new.md"
            history = root / "archive/elliptic-curves/old.md"
            source.parent.mkdir(parents=True)
            history.parent.mkdir(parents=True)
            source.write_text("# New investigation\nPointed quartic timeout remains UNKNOWN.\n")
            history.write_text("# Old investigation\nPointed quartic history.\n")
            with patch.object(research, "ROOT", root), patch.object(research, "source_documents", return_value=[source, history]):
                hits = research.search("pointed quartic", [], [])
                self.assertEqual({r["kind"] for r in hits}, {"reference", "historical"})
                source.write_text("# New investigation\nAnnihilator optimization.\n")
                hits = research.search("annihilator", [], [])
                self.assertEqual([r["source"] for r in hits], ["elliptic-curves/new.md"])

    def test_historical_phrase_never_hides_current_closure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / "archive/elliptic-curves/old.md"
            old.parent.mkdir(parents=True)
            old.write_text("# Curve302 parent remains UNKNOWN\n")
            current = claim("RECOVERED", title="Recovered parent of Curve302", scope="Explicit MW17 parent is complete.")
            with patch.object(research, "ROOT", root), patch.object(research, "source_documents", return_value=[old]):
                hits = research.search("curve302 parent", [current], [])
            self.assertEqual(hits[0]["kind"], "claim")
            self.assertEqual(hits[0]["id"], "RECOVERED")
            self.assertEqual(hits[1]["kind"], "historical")

    def test_programme_archive_keeps_open_state_but_leaves_default_search(self):
        archived = claim("PLANE-OPEN", state="open", kind="open_problem",
                         source="archive/non-elliptic/plane-jc/result.md", programme_status="archived")
        current = claim("EC-OPEN", state="open", kind="open_problem")
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(research, "ROOT", Path(directory)), \
                patch.object(research, "source_documents", return_value=[]):
            hits = research.search("pointed quartic", [archived, current], [])
            history = research.search("pointed quartic", [archived, current], [], history=True)
        self.assertEqual([r['id'] for r in hits], ['EC-OPEN'])
        self.assertEqual({r['id'] for r in history}, {'EC-OPEN', 'PLANE-OPEN'})
        self.assertEqual(research.show('PLANE-OPEN', [archived, current], [])['entry']['state'], 'open')
        self.assertEqual(research.area(archived['canonical_source']), 'plane-jc')
        self.assertNotIn('PLANE-OPEN', render({'entries': [archived, current]}))

    def test_active_claim_cannot_use_archived_programme_proof(self):
        from research_programmes import programme_path_allowed
        path = 'archive/non-elliptic/plane-jc/result.md'
        self.assertFalse(programme_path_allowed(path, claim('C')))
        self.assertTrue(programme_path_allowed(path, claim('C', programme_status='archived')))
        self.assertFalse(programme_path_allowed('archive/retired-proof.md', claim('C', programme_status='archived')))

    def test_default_work_excludes_archived_programmes(self):
        import research_work
        entries = [claim('OLD', state='partial', programme_status='archived'), claim('EC', state='partial')]
        self.assertEqual([r['id'] for r in research_work.records({}, {}, entries, 'partial')], ['EC'])
        self.assertEqual(len(research_work.records({}, {}, entries, 'partial', history=True)), 2)

    def test_full_show_does_not_import_or_execute_checker(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "checker.py"
            p.write_text("raise RuntimeError('checker must never execute during retrieval')\n")
            item = claim("C", checker=str(p), scope="A " + "full scope " * 1000 + "UNKNOWN")
            result = research.show("C", [item], [])
            self.assertEqual(result["entry"], item)
            self.assertTrue(result["entry"]["scope"].endswith("UNKNOWN"))
            with self.assertRaises(ValueError):
                research.show("MISSING", [item], [])

    def test_every_registered_claim_is_in_exactly_one_area_catalogue(self):
        entries = [claim("A", source="verified/core.md"), claim("B", source="HC4_TEST.md"),
                   claim("C", state="partial"), claim("D", state="archived"),
                   claim("E", state="parked", kind="open_problem", source="plane-jc/front.md")]
        with patch.object(research, "source_documents", return_value=[]):
            views = research.render_outputs(entries, [])
        catalogues = [text for p, text in views.items() if p.name in {a + ".md" for a in research.AREAS}]
        for entry in entries:
            self.assertEqual(sum(text.count(f"| `{entry['id']}` |") for text in catalogues), 1)
        self.assertIn("UNKNOWN", research.show("C", entries, [])["entry"]["scope"])

    def test_failed_routes_keep_scope_and_witness(self):
        proof = claim("BOUNDED", scope="Only one fixed chart is excluded.")
        problem = claim("OPEN", kind="open_problem", state="open",
                        forbidden_attack_classes=[dict(attack="Repeat fixed chart", reason="Only that chart is excluded.", witnesses=["BOUNDED"])])
        with patch.object(research, "source_documents", return_value=[]):
            views = research.render_outputs([proof, problem], [])
        text = views[research.ROOT / "knowledge/FAILED_ROUTES.md"]
        self.assertIn("Only that chart is excluded.", text)
        self.assertIn("[BOUNDED]", text)

    def test_lesson_with_unknown_evidence_is_rejected(self):
        item = {key: "text" for key in ("id", "title", "when", "use", "avoid", "boundary", "reopen_when")}
        item.update(area="core", claims=["MISSING"], sources=["missing.md"], implementation=[], keywords=[], reviewed_claims={})
        with self.assertRaisesRegex(AssertionError, "unknown claim"):
            research.validate_lessons([item], [])

    def test_changed_scope_marks_lesson_for_review(self):
        old = claim("C", scope="One chart is open.")
        item = dict(id="L", claims=["C"], reviewed_claims={"C": research.claim_fingerprint(old)})
        self.assertEqual(research.stale_claims(item, [old]), [])
        corrected = {**old, "scope": "The chart is closed by an exact obstruction."}
        self.assertEqual(research.stale_claims(item, [corrected]), ["C"])

    def test_fresh_render_does_not_approve_stale_lesson(self):
        old = claim("C", scope="Proposed rank upper bound.")
        new = {**old, "state": "partial", "scope": "Upper bound withdrawn; UNKNOWN."}
        item = {key: "text" for key in ("id", "title", "when", "use", "avoid", "boundary", "reopen_when")}
        item.update(area="core", claims=["C"], sources=["canonical.md"], implementation=[], keywords=[],
                    reviewed_claims={"C": research.claim_fingerprint(old)})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "canonical.md").write_text("Updated canonical proof boundary.\n")
            with patch.object(research, "ROOT", root), patch.object(research, "source_documents", return_value=[]):
                research.validate_lessons([item], [new])
                generated = research.render_outputs([new], [item])
                self.assertIn("Review needed", generated[root / "knowledge/ALGORITHMS.md"])
                with self.assertRaisesRegex(AssertionError, "review changed claim scopes"):
                    research.validate_lessons([item], [new], check_reviews=True)

    def test_compact_status_keeps_queue_and_falsification_anchors(self):
        entries = [claim(i) for i in CORE_ORDER]
        entries += [claim("OPEN", state="open", kind="open_problem"),
                    claim("CLOSED", state="parked", kind="open_problem", replaced_by=["F1"]),
                    claim("FALSE", state="falsified")]
        text = render({"entries": entries})
        for phrase in ("## Active open problems", "## Parked problems", "## Falsified claims",
                       "`OPEN`", "`CLOSED`", "replaced by F1", "`FALSE`"):
            self.assertIn(phrase, text)
        self.assertLess(len(text), 15000)


class DocumentDiscoveryTests(unittest.TestCase):
    def test_keeps_pinned_ignored_docs_and_skips_untracked_run_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / ".gitignore").write_text("cache/\n")
            (root / "cache").mkdir()
            (root / "cache/pinned.md").write_text("pinned proof\n")
            (root / "cache/worker-output.md").write_text("transient\n")
            (root / "new note.md").write_text("new source\n")
            subprocess.run(["git", "add", "-f", "cache/pinned.md"], cwd=root, check=True)
            paths = {p.relative_to(root).as_posix() for p in document_paths(root)}
            self.assertEqual(paths, {"cache/pinned.md", "new note.md"})


if __name__ == "__main__":
    unittest.main()
